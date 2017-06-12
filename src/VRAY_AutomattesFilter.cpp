/*
 */
#include <UT/UT_DSOVersion.h>
#include <VRAY/VRAY_SpecialChannel.h>
#include <UT/UT_Args.h>
#include <UT/UT_StackBuffer.h>
#include <SYS/SYS_Floor.h>
#include <SYS/SYS_Math.h>
#include <OpenEXR/half.h>
#include <IMG/IMG_DeepShadow.h>
#include <IMG/IMG_File.h>
#include <PXL/PXL_Raster.h>
#include <PXL/PXL_DeepSampleList.h>

#include "VRAY_AutomattesFilter.hpp"

#include <iostream>
#include <map>


using namespace HA_HDK;

VRAY_PixelFilter *
allocPixelFilter(const char *name)
{
    // NOTE: We could use name to distinguish between multiple pixel filters,
    //       in the same library, but we only have one.
    return new VRAY_AutomatteFilter();
}


VRAY_AutomatteFilter::VRAY_AutomatteFilter()
    : mySamplesPerPixelX(1) 
    , mySamplesPerPixelY(1) 
    , mySortByPz(1)
    , myOffset(0)
    , myFilterWidth(2)
    , myGaussianAlpha(1)
    , myGaussianExp(0)
    , myRank(0)
    , myHashChannel("m3hash")
{

}


VRAY_AutomatteFilter::~VRAY_AutomatteFilter()
{

}

VRAY_PixelFilter *
VRAY_AutomatteFilter::clone() const
{
    // In this case, all of our members can be default-copy-constructed,
    // so we don't need to write a copy constructor implementation.
    VRAY_AutomatteFilter *pf = new VRAY_AutomatteFilter(*this);
    return pf;
}

void
VRAY_AutomatteFilter::setArgs(int argc, const char *const argv[])
{
    UT_Args args;
    args.initialize(argc, argv);
    args.stripOptions("w:o:");

    if (args.found('w')) { myFilterWidth = args.fargp('w'); }
    if (args.found('o')) { myOffset      = args.fargp('o'); }
    // if (args.found('h')) { myHashChannel = args.argp('h'); }
}

void
VRAY_AutomatteFilter::getFilterWidth(float &x, float &y) const
{
    float filterwidth = myFilterWidth;
    x = filterwidth;
    y = filterwidth;
}

void
VRAY_AutomatteFilter::addNeededSpecialChannels(VRAY_Imager &imager)
{
    // if (myUseOpID)
    //     addSpecialChannel(imager, VRAY_SPECIAL_OPID);
    addSpecialChannel(imager, VRAY_SPECIAL_PZ);

}

namespace {
    float VRAYcomputeSumX2(int samplesperpixel, float width, int &halfsamplewidth)
    {
      float sumx2 = 0;
      if (samplesperpixel & 1 ) {
          halfsamplewidth = (int)SYSfloor(float(samplesperpixel)*0.5f*width);
          for (int i = -halfsamplewidth; i <= halfsamplewidth; ++i) {
              float x = float(i)/float(samplesperpixel);
              sumx2 += x*x;
          }
      } else {
          halfsamplewidth = (int)SYSfloor(float(samplesperpixel)*0.5f*width + 0.5f);
          for (int i = -halfsamplewidth; i < halfsamplewidth; ++i) {
              float x = (float(i)+0.5f)/float(samplesperpixel);
              sumx2 += x*x;
          }
      }
      return sumx2;
    }

}

void
VRAY_AutomatteFilter::prepFilter(int samplesperpixelx, int samplesperpixely)
{
    mySamplesPerPixelX = samplesperpixelx;
    mySamplesPerPixelY = samplesperpixely;

    myOpacitySumX2 = VRAYcomputeSumX2(mySamplesPerPixelX, myFilterWidth, myOpacitySamplesHalfX);
    myOpacitySumY2 = VRAYcomputeSumX2(mySamplesPerPixelY, myFilterWidth, myOpacitySamplesHalfY);
    myGaussianExp  = SYSexp(-myGaussianAlpha * myFilterWidth * myFilterWidth);

}

void
VRAY_AutomatteFilter::filter(
    float *destination,
    int vectorsize,
    const VRAY_SampleBuffer &source,
    int channel,
    int sourcewidth,
    int sourceheight,
    int destwidth,
    int destheight,
    int destxoffsetinsource,
    int destyoffsetinsource,
    const VRAY_Imager &imager) const
{

    // It's not technically necessery, but some convention needs to be taken.
    UT_ASSERT(vectorsize == 4);

    const float *const colordata = getSampleData(source, channel);
    // const float *const pzdata    = getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_PZ));

    for (int desty = 0; desty < destheight; ++desty) 
    {
        for (int destx = 0; destx < destwidth; ++destx)
        {
            // First, compute the sample bounds of the pixel
            const int sourcefirstx = destxoffsetinsource + destx*mySamplesPerPixelX;
            const int sourcefirsty = destyoffsetinsource + desty*mySamplesPerPixelY;
            const int sourcelastx = sourcefirstx + mySamplesPerPixelX-1;
            const int sourcelasty = sourcefirsty + mySamplesPerPixelY-1;
            // Find the first sample to read for opacity and Pz
            const int sourcefirstox = sourcefirstx + (mySamplesPerPixelX>>1) - myOpacitySamplesHalfX;
            const int sourcefirstoy = sourcefirsty + (mySamplesPerPixelY>>1) - myOpacitySamplesHalfY;
            // // Find the last sample to read for colour and z gradients
            const int sourcelastox = sourcefirstx + ((mySamplesPerPixelX-1)>>1) + myOpacitySamplesHalfX;
            const int sourcelastoy = sourcefirsty + ((mySamplesPerPixelY-1)>>1) + myOpacitySamplesHalfY;

            int sourcefirstrx = sourcefirstox;
            int sourcefirstry = sourcefirstoy;
            int sourcelastrx = sourcelastox;
            int sourcelastry = sourcelastoy;
          
            
            UT_StackBuffer<float> sample(vectorsize);
            for (int i = 0; i < vectorsize; ++i)
                sample[i] = 0.f;

            HashMap hash_map;
            float gaussianNorm = 0;

            for (int sourcey = sourcefirstry; sourcey <= sourcelastry; ++sourcey)
            {
                for (int sourcex = sourcefirstrx; sourcex <= sourcelastrx; ++sourcex)
                {
                    const int sourceidx = sourcex + sourcewidth*sourcey;
                    if(sourcex >= sourcefirstox && sourcex <= sourcelastox &&\
                      sourcey >= sourcefirstoy && sourcey <= sourcelastoy) 
                    {

                        // Find (x,y) of sample relative to *middle* of pixel
                        const float x = (float(sourcex)-0.5f*float(sourcelastx + sourcefirstx))\
                            / float(mySamplesPerPixelX);
                        const float y = (float(sourcey)-0.5f*float(sourcelasty + sourcefirsty))\
                            / float(mySamplesPerPixelY);

                        // TODO: remove magic number
                        const float gaussianWeight = gaussianFilter(x*1.66667, y*1.66667, myGaussianExp, \
                            myGaussianAlpha);
                        gaussianNorm += gaussianWeight;

                        // for (int i = 0; i < vectorsize; ++i) {
                        //     const float c = colordata[vectorsize*sourceidx+i];
                        //     uint seed  = static_cast<const uint>(c);
                        //     sample[i] += gaussianWeight* SYSfastRandom(seed);
                        // }
                        // For now this is our coverage sample (Af), which means 
                        //  no transparency support with precomposed pixel samples.
                        const float alpha = colordata[vectorsize*sourceidx+3] * gaussianWeight; 
                        // R channel of ID export. Not sure atm how to manage three types of IDs /
                        // We can export from a shader hashes for objects and materiale names,
                        // groupid doesn't work nor would it have much sense anyway.
                        const float object_id   = colordata[vectorsize*sourceidx];    // R -> object_id
                        // const float material_id = colordata[vectorsize*sourceidx+1];  // G -> object_id

                        uint seed  = static_cast<const uint>(object_id);
                        sample[1] += gaussianWeight* SYSfastRandom(seed);

                        if (hash_map.find(object_id) == hash_map.end()) {
                            hash_map.insert(std::pair<float, float>(object_id, alpha)); 
                        }
                        else {
                            hash_map[object_id] += alpha;
                        } 
                    }
                }
            }


            
            HashMap coverage_map;
            // sort by coverage (alpha here)
            HashMap::const_iterator it(hash_map.begin());
            for(; it != hash_map.end(); ++it) {
                const float alpha     = it->second;// / gaussianNorm;
                const float object_id = it->first;// ? alpha != 0.0f: 0.f;
                coverage_map.insert(std::pair<float, float>(alpha, object_id));
            }

            HashMap::const_reverse_iterator rit(coverage_map.rbegin());

            if (myOffset == 0) {
                sample[2] = rit->second * gaussianNorm;
                for (int i = 0; i< 4; ++i, ++destination) {
                    *destination  = sample[i] / gaussianNorm; 
                }
                    
            } else {
        
                const size_t id_offset = (static_cast<size_t>(myOffset) - 1) * 2; 
                std::advance(rit, id_offset);

                for (int i = id_offset; i < 2; ++rit,  ++i) {
                    if (rit == coverage_map.rend()) {
                        destination += i*2;
                        break;
                    }
                    destination[0] = rit->second; // object_id
                    destination[1] = rit->first / gaussianNorm; // alpha
                    destination += 2;
                }   
                       
            }
        }
    }
}



