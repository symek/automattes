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
    : mySamplesPerPixelX(1) // Initialized just in case; value shouldn't be used
    , mySamplesPerPixelY(1) // Initialized just in case; value shouldn't be used
    , mySortByOpacity(false)
    , mySortByPz(true)
    , myUseOpID(false)
    , myMaskNumber(4)
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
    args.stripOptions("i:z:o:n:");

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
}

void
VRAY_MagicMaskFilter::getFilterWidth(float &x, float &y) const
{
    float filterwidth = 1.0f;
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



void
VRAY_MagicMaskFilter::prepFilter(int samplesperpixelx, int samplesperpixely)
{
    mySamplesPerPixelX = samplesperpixelx;
    mySamplesPerPixelY = samplesperpixely;
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

    UT_ASSERT(opacitydata != NULL);
    UT_ASSERT(mySortByPz == (zdata != NULL));
    UT_ASSERT(myUseOpID == (opiddata != NULL));

   
    const int lwidth  = mySamplesPerPixelY * sourcewidth * vectorsize;
    const int swidth  = mySamplesPerPixelX * vectorsize;
    const int reoffx  = destxoffsetinsource / mySamplesPerPixelX;
    const int reoffy  = destyoffsetinsource / mySamplesPerPixelY;

    for (int desty = 0; desty < destheight; ++desty)
    {
        for (int destx = 0; destx < destwidth; ++destx)
        {

            const int dxf = destx + reoffx;
            const int dyf = desty + reoffy;
                  int sx  = dyf * lwidth + dxf * swidth;

          
            float value = 0;
            for (int y=0; y<mySamplesPerPixelY; ++y) {
                sx += y*sourcewidth*vectorsize;
                for( int x=0; x<mySamplesPerPixelX; ++x) {
                    value += zdata[sx+x];
                }
            }

            value /= (mySamplesPerPixelX*mySamplesPerPixelY);

            if (value > 1000)
                value = 0.0;

            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = value; //SYSmin(value, 100.0f);
        }
    }
}
