/*
 */

#include <UT/UT_DSOVersion.h>

#include "VRAY_MagicMaskFilter.hpp"
#include <VRAY/VRAY_SpecialChannel.h>
#include <UT/UT_Args.h>
#include <UT/UT_StackBuffer.h>
#include <SYS/SYS_Floor.h>
#include <SYS/SYS_Math.h>
#include <iostream>
#include <map>
#include <OpenEXR/half.h>

using namespace HA_MMask;

VRAY_PixelFilter *
allocPixelFilter(const char *name)
{
    // NOTE: We could use name to distinguish between multiple pixel filters,
    //       in the same library, but we only have one.
    return new VRAY_MagicMaskFilter();
}


VRAY_MagicMaskFilter::VRAY_MagicMaskFilter()
    : mySamplesPerPixelX(1) 
    , mySamplesPerPixelY(1) 
    , mySortByOpacity(false)
    , mySortByPz(true)
    , myUseOpID(true)
    , myMaskNumber(4)
    , myFilterWidth(2)
    , myGaussianAlpha(1)
    , myGaussianExp(0)
{
}

VRAY_MagicMaskFilter::~VRAY_MagicMaskFilter()
{
}

VRAY_PixelFilter *
VRAY_MagicMaskFilter::clone() const
{
    // In this case, all of our members can be default-copy-constructed,
    // so we don't need to write a copy constructor implementation.
    VRAY_MagicMaskFilter *pf = new VRAY_MagicMaskFilter(*this);
    return pf;
}

void
VRAY_MagicMaskFilter::setArgs(int argc, const char *const argv[])
{
    UT_Args args;
    args.initialize(argc, argv);
    args.stripOptions("i:z:o:n:w:");

    if (args.found('i')) {
        myUseOpID = true;
    }
    if (args.found('z')) {
        mySortByPz = true;
    }
    if (args.found('o')) {
        mySortByOpacity = true;
    }
    if (args.found('n')) {
       myMaskNumber =  args.iargp('n');
    }
    if (args.found('w')) {
        myFilterWidth = args.fargp('w');
    }


}

void
VRAY_MagicMaskFilter::getFilterWidth(float &x, float &y) const
{
    float filterwidth = myFilterWidth;
    x = filterwidth;
    y = filterwidth;
}

void
VRAY_MagicMaskFilter::addNeededSpecialChannels(VRAY_Imager &imager)
{
    if (myUseOpID)
        addSpecialChannel(imager, VRAY_SPECIAL_OPID);
    if (mySortByPz)
        addSpecialChannel(imager, VRAY_SPECIAL_PZ);

    addSpecialChannel(imager, VRAY_SPECIAL_OPACITYSAMPLES);
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
VRAY_MagicMaskFilter::prepFilter(int samplesperpixelx, int samplesperpixely)
{
    mySamplesPerPixelX = samplesperpixelx;
    mySamplesPerPixelY = samplesperpixely;

    myOpacitySumX2 = VRAYcomputeSumX2(mySamplesPerPixelX, myFilterWidth, myOpacitySamplesHalfX);
    myOpacitySumY2 = VRAYcomputeSumX2(mySamplesPerPixelY, myFilterWidth, myOpacitySamplesHalfY);
    myGaussianExp  = SYSexp(-myGaussianAlpha * myFilterWidth * myFilterWidth);

    // std::cout << "myOpacitySamplesHalf: " << myOpacitySamplesHalfX  <<", " << myOpacitySamplesHalfY << std::endl; 
    // std::cout << "myOpacitySum*2: " << myOpacitySumX2  <<", " << myOpacitySumY2 << std::endl; 
}

void
VRAY_MagicMaskFilter::filter(
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

    const float *const opacitydata = \
    getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_OPACITYSAMPLES));
    const float *const zdata = mySortByPz
        ? getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_PZ))
        : NULL;

    const float *const opiddata = myUseOpID
        ? getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_OPID))
        : NULL;

    const float *const colourdata = \
    getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_CFAF));

    UT_ASSERT(opacitydata != NULL);
    UT_ASSERT(mySortByPz == (zdata != NULL));
    UT_ASSERT(myUseOpID == (opiddata != NULL));

    // const float normalizer = 1.0f/mySamplesPerPixelX*mySamplesPerPixelY;
    for (int desty = 0; desty < destheight; ++desty)
    {
        for (int destx = 0; destx < destwidth; ++destx)
        {

            // I am giving up my own attempts and copy/paste from HDK example...
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
            // sourcefirstrx = SYSmin(sourcefirstrx, sourcefirstox);
            // sourcefirstry = SYSmin(sourcefirstry, sourcefirstoy);
            // sourcelastrx  = SYSmax(sourcelastrx, sourcelastox);
            // sourcelastry  = SYSmax(sourcelastry, sourcelastoy);
            
           
            UT_StackBuffer<float> sample(vectorsize);
            for (int i = 0; i < vectorsize; ++i)
                sample[i] = 0;

            std::map<int, float> alphaMap;
            std::map<int, float> counterMap;

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
                        // Why is that magic number?
                        const float gaussianWeight = gaussianFilter(x*1.66667, y*1.66667, myGaussianExp, myGaussianAlpha);
                        gaussianNorm += gaussianWeight;
                        for (int i = 0; i < vectorsize; ++i) {
                            sample[i] += gaussianWeight*colourdata[vectorsize*sourcei+i];
                        }
                        const int id = opiddata[sourcei];
                        if (alphaMap.find(id) == alphaMap.end() && id != -1) {
                            alphaMap.insert(std::pair<int, float>(id, 0.f));
                            counterMap.insert(std::pair<int, int>(id, 0));
                        } 
                        // const float alpha = sample[vectorsize-1];  
                        if (id != -1) {
                            alphaMap[id]   +=  gaussianWeight*colourdata[vectorsize*sourcei+3];
                            counterMap[id] += 1;
                        }
                    }
                }
            }

            for(std::map<int, float>::iterator it(alphaMap.begin()); \
                it != alphaMap.end(); ++it) {
                it->second /= gaussianNorm;
            }
            
            /* make mask indexing */

            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = sample[i];
           
        }
    }
}
