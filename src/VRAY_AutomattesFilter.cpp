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


inline float gaussian(float d, float expv, float alpha)
{
    return SYSmax(0.f, float(SYSexp(-alpha*d*d) - expv));
}

inline float gaussianFilter(float x, float y, float expv, float alpha)
{
    return gaussian(x, expv, alpha) * gaussian(y, expv, alpha);
}

inline void packFloats(const float a, const float b, float &store)
{
    const half first = half(a); const half second = half(b);
    int16_t sh1 = *reinterpret_cast<int16_t*>((void*) &first);
    int16_t sh2 = *reinterpret_cast<int16_t*>((void*) &second);
    int32_t tmp = ( sh2 << 16) | sh1;
          store = *reinterpret_cast<float*>((void*)&(tmp)); 
}


inline void unpackFloats(const float store, float &a, float &b)
{
    int16_t unpack16a = *reinterpret_cast<int16_t*>((void*)&store);
    int16_t unpack16b = *reinterpret_cast<int32_t*>((void*)&store) >> 16;
    a = static_cast<float>(*reinterpret_cast<half*>((void*)&unpack16a));
    b = static_cast<float>(*reinterpret_cast<half*>((void*)&unpack16b));
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

    const float *const zdata = mySortByPz
        ? getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_PZ))
        : NULL;

    const float *const opiddata = myUseOpID
        ? getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_OPID))
        : NULL;

    const float *const colourdata = \
    getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_CFAF));

    const float * m3hashdata = !myUseOpID
        ? getChannelIdxByName(imager, myHashChannel)
        : NULL;

    
    UT_ASSERT(mySortByPz == (zdata != NULL));
    UT_ASSERT(myUseOpID == (opiddata != NULL));
    UT_ASSERT(myUseOpID != (m3hashdata != NULL));


    // IMG_DeepPixelWriter writer(*myDsm);

    for (int desty = 0; desty < destheight; ++desty)
    {
        for (int destx = 0; destx < destwidth; ++destx)
        {

            // Thanks SESI for HDK example...
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
                sample[i] = 0;

            HashSamples hashMap;
            float gaussianNorm = 0;
            for (int sourcey = sourcefirstry; sourcey <= sourcelastry; ++sourcey)
            {
                for (int sourcex = sourcefirstrx; sourcex <= sourcelastrx; ++sourcex)
                {
                    const int sourcei = sourcex + sourcewidth*sourcey;
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
                            sample[i] += gaussianWeight*colourdata[vectorsize*sourcei+i];
                        }

                        const float   alpha = colourdata[vectorsize*sourcei+3] * gaussianWeight; // TODO: move to opacitySamples?
                        const uint32_t hash = m3hashdata[sourcei];

                        if (hashMap.find(hash) == hashMap.end()) {
                            hashMap.insert(std::pair<uint32_t, float>(hash, alpha)); 
                        }
                        else {
                            hashMap[hash] += alpha;
                        } 
                    }
                }
            }

            const int pixelIndex = (destxoffsetinsource + destx*mySamplesPerPixelX) + \
            sourcewidth*(destyoffsetinsource + desty*mySamplesPerPixelY);
            uint px, py; px = py = 0;

            if (pixeldata) {
                px = static_cast<uint>(pixeldata[3*pixelIndex]);
                py = static_cast<uint>(pixeldata[3*pixelIndex+1]);
                px = SYSmin(px, myXRes-1);
                py = SYSmin(py, myYRes-1);
            }


           { 
                
                std::map<float, uint32_t>   hashOrderedByCoverage;
    
                uint32_t combinedHash = 0;
                // two ids per raster == 2*3+first technical raster (all RGBA)
                HashSamples::const_iterator it(hashMap.begin());
                for (uint i =1; i < myRasters.size()*2; ++i, ++it) {
                    if (it!=hashMap.end()) {   
                        const float   alpha = SYSmax(it->second/gaussianNorm, 0.f);
                        const uint32_t hash = it->first ? alpha != 0.0f: 0;
                        hashOrderedByCoverage.insert(std::pair<float, uint32_t>(alpha, hash));
                        combinedHash += hash; // TODO: should I add floats or ints?
                    } 
                    // else {
                    //     // FIXME: This is bug.
                    //     hashOrderedByCoverage.insert(std::pair<float, uint32_t>(0.f, 0));
                    // }
                }

                // First 'technical' raster RGBA
                float vals[4];
                vals[0] = vals[1] = vals[2] = vals[3] = 0.f;
                vals[0] = hash_to_float(combinedHash); // TODO: ?
                vals[1] = ((float) ((combinedHash << 8)) /  (float) UINT32_MAX);
                vals[2] = ((float) ((combinedHash << 16)) / (float) UINT32_MAX);
                vals[3] = 0.0f;
                myRasters(0)->setPixelValue(px, py, vals);

                
                // Then 3 rasters with two mattes each:
                std::map<float, uint32_t>::const_reverse_iterator jt(hashOrderedByCoverage.rbegin());
                for (uint i=1; i<myRasters.size(); ++i, ++jt) {
                    vals[0] = vals[1] = vals[2] = vals[3] = 0.f;
                    for (uint j=0; j<2; ++j, ++jt) {
                        if (jt==hashOrderedByCoverage.rend())
                            break;
                        float val  = jt->first;
                        float hash = hash_to_float(jt->second);
                        vals[2*j]   = hash;
                        vals[2*j+1] = val;
                        
                    }
                    myRasters(i)->setPixelValue(px, py, vals);
                }
            }


            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = sample[i] / gaussianNorm;
           
        }
    }
}
