/*
 */

#include <UT/UT_DSOVersion.h>

#include "VRAY_AutomattesFilter.hpp"
#include <VRAY/VRAY_SpecialChannel.h>
#include <UT/UT_Args.h>
#include <UT/UT_StackBuffer.h>
#include <SYS/SYS_Floor.h>
#include <SYS/SYS_Math.h>
#include <iostream>
#include <map>
#include <OpenEXR/half.h>
#include <IMG/IMG_DeepShadow.h>
#include <PXL/PXL_DeepSampleList.h>

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
    , mySortByOpacity(false)
    , mySortByPz(true)
    , myUseOpID(true)
    , myMaskNumber(4)
    , myFilterWidth(2)
    , myGaussianAlpha(1)
    , myGaussianExp(0)
    , myXRes(1280)
    , myYRes(720)
    , myDeepImagePath("./test.rat")
{
}


VRAY_AutomatteFilter::~VRAY_AutomatteFilter()
{
    IMG_DeepPixelWriter writer(*myDsm);
    for (int y=0; y<myYRes; ++y) {
        for(int x=0; x<myXRes; ++x) {
            const IdSamples mask = mySamples->get(x, y);
            if (mask.size() == 0)
                continue; 
            writer.open(x, y);
            // DEBUG_PRINT("pixel %i, %i. Samples:", x, y);
            IdSamples::const_iterator it(mask.begin());
            for (; it != mask.end(); ++it) {
                const float z = it->first;
                float v[3];
                v[0] = v[1] = v[2] = it->second;
                // DEBUG_PRINT("%i at %f, ", it->first, v[0]);
                writer.write(z, v, 3, PXL_DeepSampleList::MATTE_SURFACE, -1, 0);
            }
            // std::cout << "\n";
            writer.close();
        }
    }
    
    myDsm->close();
    delete myDsm;
    delete mySamples;
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
    args.stripOptions("i:z:o:n:w:x:y:p:");

    if (args.found('i')) { myUseOpID  = true; }
    if (args.found('z')) { mySortByPz = true; }
    if (args.found('o')) { mySortByOpacity = true; }
    if (args.found('n')) { myMaskNumber  = args.iargp('n'); }
    if (args.found('w')) { myFilterWidth = args.fargp('w'); }
    if (args.found('p')) { myDeepImagePath = args.argp('p');}
    if (args.found('x')) { myXRes = args.iargp('x'); }
    if (args.found('y')) { myYRes = args.iargp('y'); }
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
VRAY_AutomatteFilter::prepFilter(int samplesperpixelx, int samplesperpixely)
{
    mySamplesPerPixelX = samplesperpixelx;
    mySamplesPerPixelY = samplesperpixely;

    myOpacitySumX2 = VRAYcomputeSumX2(mySamplesPerPixelX, myFilterWidth, myOpacitySamplesHalfX);
    myOpacitySumY2 = VRAYcomputeSumX2(mySamplesPerPixelY, myFilterWidth, myOpacitySamplesHalfY);
    myGaussianExp  = SYSexp(-myGaussianAlpha * myFilterWidth * myFilterWidth);
    myDsm          = new IMG_DeepShadow();
    myDsm->setOption("deepcompression", "1");
    myDsm->setOption("zbias", "0.05");
    myDsm->setOption("depth_planes", "Pz,Zback");
    myDsm->setOption("compositing", 1);
    
    myDsm->create(myDeepImagePath, myXRes, myYRes, 1  /*mySamplesPerPixelX*/,1  /*mySamplesPerPixelY*/); // !!!
    // DEBUG_PRINT("%s", "Starting new filter.");

    mySamples = new AutomatteSamples();
    mySamples->init(myXRes, myYRes);   

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

    const float * pixeldata = NULL;
    const int pixelnumindex = getChannelIdxByName(imager, "pixelnum");
    if (pixelnumindex != -1)
        pixeldata = getSampleData(source, pixelnumindex );

    

    UT_ASSERT(opacitydata != NULL);
    UT_ASSERT(mySortByPz == (zdata != NULL));
    UT_ASSERT(myUseOpID == (opiddata != NULL));
    UT_ASSERT(pixeldata != NULL);


    IMG_DeepPixelWriter writer(*myDsm);

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

            IdSamples sampleMap;
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

                        const float alpha   = colourdata[vectorsize*sourcei+3]; // TODO: move to opacitySamples?
                        const int   idMatte = opiddata[sourcei];

                        if (sampleMap.find(idMatte) == sampleMap.end())
                            sampleMap.insert(std::pair<int, float>(idMatte, alpha));
                        else 
                            sampleMap[idMatte] += (alpha * gaussianWeight);
                    }
                }
            }

            
            const int pixelIndex = (destxoffsetinsource + destx*mySamplesPerPixelX) + \
            sourcewidth*(destyoffsetinsource + desty*mySamplesPerPixelY);
            int px, py; px = py = 0;

            if (pixeldata) {
                px = static_cast<int>(pixeldata[3*pixelIndex]);
                py = static_cast<int>(pixeldata[3*pixelIndex+1]);
            }

            // IMG_DeepPixelWriter can't handle this..? 
            if (sampleMap.size())
                mySamples->write(px, py, gaussianNorm, sampleMap);

            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = sample[i] / gaussianNorm;
           
        }
    }
}
