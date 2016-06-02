/*
 Taken stright from HDK example.
 */

#pragma once

#ifndef __VRAY_MagicMaskFilter__
#define __VRAY_MagicMaskFilter__

#include <VRAY/VRAY_PixelFilter.h>
#include <VRAY/VRAY_Procedural.h>

class VRAY_Imager;
class VRAY_SampleBuffer;
class IMG_DeepShadow;

namespace HA_MMask {

class VRAY_MagicMaskFilter : public VRAY_PixelFilter {
public:
    VRAY_MagicMaskFilter();
    virtual ~VRAY_MagicMaskFilter();

    virtual VRAY_PixelFilter *clone() const;

    /// setArgs is called with the options specified after the pixel filter
    /// name in the Pixel Filter parameter on the Mantra ROP.
    /// -i 1 use OpId istead of 'mask' raster
    /// -z 1 sort by Pz
    /// -o 1 sort by opacity
    /// -n 8 size of stored sampels per pixel. 
 
    virtual void setArgs(int argc, const char *const argv[]);

    /// getFilterWidth is called after setArgs when Mantra needs to know
    /// how far to expand the render region.
    virtual void getFilterWidth(float &x, float &y) const;

    /// addNeededSpecialChannels is called after setArgs so that this filter
    /// can indicate that it depends on having special channels like z-depths
    /// or Op IDs.
    virtual void addNeededSpecialChannels(VRAY_Imager &imager);

    /// prepFilter is called after setArgs so that this filter can
    /// precompute data structures or values for use in filtering that
    /// depend on the number of samples per pixel in the x or y directions.
    virtual void prepFilter(int samplesperpixelx, int samplesperpixely);

    /// filter is called for each destination tile region with a source
    /// that is at least as large as is needed by this filter, based on
    /// the filter widths returned by getFilterWidth.
    virtual void filter(
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
        const VRAY_Imager &imager) const;

private:
    /// These must be saved in prepFilter.
    /// Each pixel has mySamplesPerPixelX*mySamplesPerPixelY samples.
    /// @{
    int mySamplesPerPixelX;
    int mySamplesPerPixelY;
    /// @}

    /// true if generate mask using Operator ID
    bool myUseOpID;

    /// true if sorting by Pz is enabled
    bool mySortByPz;

    /// true if sorting by Opacity is enabled
    bool mySortByOpacity;

    /// How many object we will store.
    int myMaskNumber;

    /// Filter width (default 2)
    float myFilterWidth;

    float myOpacitySumX2;
    float myOpacitySumY2;

    int myOpacitySamplesHalfX;
    int myOpacitySamplesHalfY;

    /// Gaussians
    float myGaussianExp;
    float myGaussianAlpha;

    IMG_DeepShadow *myDsm;

    int myXRes;
    int myYRes;

    const char *myDeepImagePath;
};


} // End HA_MMask namespace

#endif
