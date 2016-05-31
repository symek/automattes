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
    , myUseOpID(false)
    , myMaskNumber(4)
    , myFilterWidth(2)
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

}

void
VRAY_MagicMaskFilter::prepFilter(int samplesperpixelx, int samplesperpixely)
{
    mySamplesPerPixelX = samplesperpixelx;
    mySamplesPerPixelY = samplesperpixely;

    myOpacitySumX2 = VRAYcomputeSumX2(mySamplesPerPixelX, myFilterWidth, myOpacitySamplesHalfX);
    myOpacitySumY2 = VRAYcomputeSumX2(mySamplesPerPixelY, myFilterWidth, myOpacitySamplesHalfY);
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

     const float *const colourdata = getSampleData(source, VRAY_SPECIAL_CFAF);

    UT_ASSERT(opacitydata != NULL);
    UT_ASSERT(mySortByPz == (zdata != NULL));
    UT_ASSERT(myUseOpID == (opiddata != NULL));

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
            int sourcefirstrx = sourcefirstx;
            int sourcefirstry = sourcefirsty;
            int sourcelastrx = sourcelastx;
            int sourcelastry = sourcelasty;
            // sourcefirstrx = SYSmin(sourcefirstrx, sourcefirstox);
            // sourcefirstry = SYSmin(sourcefirstry, sourcefirstoy);
            // sourcelastrx  = SYSmax(sourcelastrx, sourcelastox);
            // sourcelastry  = SYSmax(sourcelastry, sourcelastoy);
            
           
            UT_StackBuffer<float> opacitySample(vectorsize);
            for (int i = 0; i < vectorsize; ++i)
                opacitySample[i] = 0;

            int counter=1;
            float value = 0;
            for (int sourcey = sourcefirstry; sourcey <= sourcelastry; ++sourcey)
            {
                for (int sourcex = sourcefirstrx; sourcex <= sourcelastrx; ++sourcex)
                {
                    if(sourcex >= sourcefirstox && sourcex <= sourcelastox &&\
                      sourcey >= sourcefirstoy && sourcey <= sourcelastoy) 
                    {

                        // Find (x,y) of sample relative to *middle* of pixel
                        // float x = (float(sourcex) - 0.5f*float(sourcelastx + sourcefirstx))/float(mySamplesPerPixelX);
                        // float y = (float(sourcey) - 0.5f*float(sourcelasty + sourcefirsty))/float(mySamplesPerPixelY);
                        const int sourcei = sourcex + sourcewidth*sourcey;
                        for (int i = 0; i < vectorsize; ++i) {
                            opacitySample[i] = colourdata[vectorsize*sourcei+i];
                            value += colourdata[vectorsize*sourcei+3];
                        }
                        counter++;
                    }
                }
            }

            for (int i = 0; i < vectorsize; ++i)
                opacitySample[i] /= counter;

            value /= counter;

            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = value;
        }
    }
}
