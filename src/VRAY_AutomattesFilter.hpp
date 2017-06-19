/*
 Partialy taken from HDK examples:
 http://www.sidefx.com/docs/hdk15.5/_v_r_a_y_2_v_r_a_y__demo_edge_detect_filter_8h-example.html

 Contains snippets from PsyOp
 [1] - Jonah Friedman, Andrew C. Jones, Fully automatic ID mattes with support for motion blur and transparency.
 */

#pragma once

#ifndef __VRAY_AutomatteFilter__
#define __VRAY_AutomatteFilter__

#include <VRAY/VRAY_PixelFilter.h>
#include <VRAY/VRAY_Procedural.h>

#define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif


typedef  std::map<float, float>  HashMap;

namespace HA_HDK {

// static const char* plane_names[] = {"CryptoObject",   "CryptoObject00", 
                                    // "CryptoObject01", "CryptoObject02"};

enum Automatte_HashType {
    MANTRA,
    CRYPTO
};


enum Automatte_IdType {
    ASSET, // not supported yet
    OBJECT,
    MATERIAL,
    GROUP // not supported yet
};

inline float gaussian(float d, float expv, float alpha) {
    return SYSmax(0.f, float(SYSexp(-alpha*d*d) - expv));
}

inline float gaussianFilter(float x, float y, float expv, float alpha) {
    return gaussian(x, expv, alpha) * gaussian(y, expv, alpha);
}

inline void packFloats(const float a, const float b, float &store) {
    const half first = half(a); const half second = half(b);
    int16_t sh1 = *reinterpret_cast<int16_t*>((void*) &first);
    int16_t sh2 = *reinterpret_cast<int16_t*>((void*) &second);
    int32_t tmp = ( sh2 << 16) | sh1;
          store = *reinterpret_cast<float*>((void*)&(tmp)); 
}

inline void unpackFloats(const float store, float &a, float &b) {
    int16_t unpack16a = *reinterpret_cast<int16_t*>((void*)&store);
    int16_t unpack16b = *reinterpret_cast<int32_t*>((void*)&store) >> 16;
    a = static_cast<float>(*reinterpret_cast<half*>((void*)&unpack16a));
    b = static_cast<float>(*reinterpret_cast<half*>((void*)&unpack16b));
}

// From Cryptomatte specification[1]
float hash_to_float(uint32_t hash)
{
    uint32_t mantissa = hash & (( 1 << 23) - 1);
    uint32_t exponent = (hash >> 23) & ((1 << 8) - 1);
    exponent = std::max(exponent, (uint32_t) 1);
    exponent = std::min(exponent, (uint32_t) 254);
    exponent = exponent << 23;
    uint32_t sign = (hash >> 31);
    sign = sign << 31;
    uint32_t float_bits = sign | exponent | mantissa;
    float f;
    std::memcpy(&f, &float_bits, 4);
    return f;
}

class VRAY_AutomatteFilter : public VRAY_PixelFilter {
public:
    VRAY_AutomatteFilter();
    virtual ~VRAY_AutomatteFilter();

    virtual VRAY_PixelFilter *clone() const;

    /// setArgs is called with the options specified after the pixel filter
    /// name in the Pixel Filter parameter on the Mantra ROP.
    /// -i 1 use OpId istead of 'mask' raster

    virtual void setArgs(int argc, const char *const argv[]);

   
    virtual void getFilterWidth(float &x, float &y) const;

    virtual void addNeededSpecialChannels(VRAY_Imager &imager);

    virtual void prepFilter(int samplesperpixelx, int samplesperpixely);

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
   
    int mySamplesPerPixelX;
    int mySamplesPerPixelY;

    float myOpacitySumX2;
    float myOpacitySumY2;

    int myOpacitySamplesHalfX;
    int myOpacitySamplesHalfY;

    // 0: filtered pseudo color, 
    // 1: for 0-1 ranks, 2: 2-3 and so on...
    int myRank; 

    const char* myIdTypeName;    // user string flag 
    Automatte_IdType myIdType;   // corresponding enumerator.

    const char* myHashTypeName;  // user string flag 
    Automatte_HashType myHashType; // corresponding enumerator.

    int mySortByPz;

    // Filter width (default 2)
    float myFilterWidth;
    //  Gaussians parms
    float myGaussianExp;
    float myGaussianAlpha;


};


} // End HA_HDK namespace

#endif
