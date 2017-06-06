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
    , myUseOpID(true)
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
    args.stripOptions("i:w:h:");

    if (args.found('i')) { myUseOpID  = true; }
    if (args.found('w')) { myFilterWidth = args.fargp('w'); }
    if (args.found('h')) { myHashChannel = args.argp('h'); }
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
    if (myUseOpID)
        addSpecialChannel(imager, VRAY_SPECIAL_OPID);
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

    const float *const colordata = getSampleData(source, channel);
    const float *const pzdata    = getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_PZ));

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

            // HashMap hash_map;
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
                        const float x = (float(sourcex) - 0.5f*float(sourcelastx + sourcefirstx))/float(mySamplesPerPixelX);
                        const float y = (float(sourcey) - 0.5f*float(sourcelasty + sourcefirsty))/float(mySamplesPerPixelY);

                        // TODO: remove magic number
                        const float gaussianWeight = gaussianFilter(x*1.66667, y*1.66667, myGaussianExp, myGaussianAlpha);
                        gaussianNorm += gaussianWeight;

                        for (int i = 0; i < vectorsize; ++i) {
                            sample[i] += gaussianWeight*colordata[vectorsize*sourceidx+i];
                        }

                        // const float alpha = colordata[vectorsize*sourceidx+3] * gaussianWeight; // TODO: move to opacitySamples?
                        // const float hash = m3hashdata[sourceidx];

                        // if (hash_map.find(hash) == hash_map.end()) {
                        //     hash_map.insert(std::pair<float, float>(hash, alpha)); 
                        // }
                        // else {
                        //     hash_map[hash] += alpha;
                        // } 
                    }
                }
            }

            // const int pixelIndex = (destxoffsetinsource + destx*mySamplesPerPixelX) + \
            // sourcewidth*(destyoffsetinsource + desty*mySamplesPerPixelY);
            // uint px, py; px = py = 0;

            // if (pixeldata) {
            //     px = static_cast<uint>(pixeldata[3*pixelIndex]);
            //     py = static_cast<uint>(pixeldata[3*pixelIndex+1]);
            //     px = SYSmin(px, myXRes-1);
            //     py = SYSmin(py, myYRes-1);
            // }


           // { //
           //      std::map<float, float>   hashOrderedByCoverage;
    
           //      float combinedHash = 0;
           //      // two ids per raster == 2*3+first technical raster (all RGBA)
           //      HashSamples::const_iterator it(hashMap.begin());
           //      for (uint i =1; i < myRasters.size()*2; ++i, ++it) {
           //          if (it!=hashMap.end()) {   
           //              const float alpha = SYSmax(it->second/gaussianNorm, 0.f);
           //              const float hash = it->first ? alpha != 0.0f: 0;
           //              hashOrderedByCoverage.insert(std::pair<float, float>(alpha, hash));
           //              combinedHash += hash; // TODO: should I add floats or ints?
           //          } 
           //      }

           //      // First 'technical' raster RGBA
           //      float vals[4];
           //      vals[0] = vals[1] = vals[2] = vals[3] = 0.f;
           //      vals[0] = combinedHash; // TODO: ?
           //      uint32_t s = static_cast<uint32_t>(combinedHash);
           //      vals[1] = ((float) ((s << 8)) /  (float) UINT32_MAX);
           //      vals[2] = ((float) ((s << 16)) / (float) UINT32_MAX);
           //      vals[3] = 0.0f;
           //      myRasters(0)->setPixelValue(px, py, vals);

                
           //      // Then 3 rasters with two mattes each:
           //      std::map<float, uint32_t>::const_reverse_iterator jt(hashOrderedByCoverage.rbegin());
           //      for (uint i=1; i<myRasters.size(); ++i, ++jt) {
           //          vals[0] = vals[1] = vals[2] = vals[3] = 0.f;
           //          for (uint j=0; j<2; ++j, ++jt) {
           //              if (jt==hashOrderedByCoverage.rend())
           //                  break;
           //              float val  = jt->first;
           //              float hash = hash_to_float(jt->second);
           //              vals[2*j]   = hash;
           //              vals[2*j+1] = val;
                        
           //          }
           //          myRasters(i)->setPixelValue(px, py, vals);
           //      }
           //  } // hide symbols

            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = sample[i] / gaussianNorm;
           
        }
    }
}
