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
    , myImagePath("./test.exr")
{
}


VRAY_AutomatteFilter::~VRAY_AutomatteFilter()
{
    myDsm->close();
    myImage->writeImages(myRasters, true); // free myRasters;
    myImage->close();
    delete myImage;
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
    if (args.found('d')) { myDeepImagePath = args.argp('d');}
    if (args.found('p')) { myImagePath     = args.argp('p'); }
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

    IMG_Stat stat = IMG_Stat(myXRes, myYRes); 
    for (uint i=0; i<4; ++i) {
        IMG_Plane *plane   = stat.addPlane(plane_names[i], IMG_FLOAT32, IMG_RGBA);
        PXL_Raster *raster = new PXL_Raster(PACK_RGBA, PXL_FLOAT32, myXRes, myYRes); // note: will free at writeImages();
        myRasters.append(raster);
    }
    myImage = IMG_File::create(myImagePath, stat);
    
    DEBUG_PRINT("%s", "Created exr image.");

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

    const float * m3hashdata = NULL;
    const int m3hashindex = getChannelIdxByName(imager, "m3hash");
    if (m3hashindex != -1)
        m3hashdata = getSampleData(source, m3hashindex);
    else
        // TODO: make option to use either opid or m3hash.
        m3hashdata = opiddata;

    

    UT_ASSERT(opacitydata != NULL);
    UT_ASSERT(mySortByPz == (zdata != NULL));
    UT_ASSERT(myUseOpID == (opiddata != NULL));
    UT_ASSERT(pixeldata != NULL);


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

            IdSamples sampleMap;
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

                        const float alpha = colourdata[vectorsize*sourcei+3] * gaussianWeight; // TODO: move to opacitySamples?
                        // const uint32_t  idMatte   = static_cast<uint32_t>(opiddata[sourcei]);
                        const float     hashMatte = m3hashdata[sourcei];

                        if (hashMap.find(hashMatte) == hashMap.end()) {
                            // sampleMap.insert(std::pair<uint32_t, float>(idMatte, alpha));
                            hashMap.insert(std::pair<float, float>(hashMatte, alpha)); 
                        }
                        else {
                            // sampleMap[idMatte] += (alpha * gaussianWeight);
                            hashMap[hashMatte] += alpha;
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

            #if 0

            IMG_DeepPixelWriter writer(*myDsm);
            writer.open(px, py);

            IdSamples::const_iterator it(sampleMap.begin());
            for (; it != sampleMap.end(); ++it) {
                const float z = it->first;
                float v[3];
                v[0] = v[1] = v[2] = SYSmax(it->second/gaussianNorm, 0.f);;
                writer.write(z, v, 3, PXL_DeepSampleList::MATTE_SURFACE, -1, 0);
            }

            writer.close();

            #else

           { 
                // uint32_t combinedId = 0;
                float  combinedHash = 0;
                // std::map<float, uint32_t> idsOrderedByCoverage;
                std::map<float, float>    hashOrderedByCoverage;
                // IdSamples::const_iterator   it(sampleMap.begin());
                HashSamples::const_iterator kt(hashMap.begin());
    
                // two ids per raster == 2*3+first technical raster (all RGBA)
                for (uint i =1; i < myRasters.size()*2; ++i, /*++it,*/ ++kt) {
                    if (kt!=hashMap.end()) {   
                        const float alpha = SYSmax(kt->second/gaussianNorm, 0.f);
                        const float hash = kt->first ? alpha != 0.0f: 0.0f;
                        // idsOrderedByCoverage.insert(std::pair<float, uint32_t>(alpha, it->first));
                        hashOrderedByCoverage.insert(std::pair<float, float>(alpha,   hash));
                        // combinedId   += it->first;
                        combinedHash += hash; // TODO: can I add hash as  floats?
                    } else {
                        // idsOrderedByCoverage.insert(std::pair<float, uint32_t>(0.0f, 0));
                        hashOrderedByCoverage.insert(std::pair<float, float>(0.f, 0.f));
                    }
                }

                // TODO: implement hashing generation most probably in VEX, where
                // we have an access to objects' names;
                float vals[4];
                float m3hash = combinedHash; //tmp
                vals[0] = m3hash;
                // vals[1] = ((float) ((m3hash << 8)) /  (float) UINT32_MAX);
                // vals[2] = ((float) ((m3hash << 16)) / (float) UINT32_MAX);
                uint32_t tmp;
                std::memcpy(&tmp, &m3hash, 4); // TODO: is it enough for back-to-uint?
                vals[1] = ((float) ((tmp << 8)) /  (float) UINT32_MAX);
                vals[2] = ((float) ((tmp << 16)) / (float) UINT32_MAX);
                vals[3] = 0.0f;
                myRasters(0)->setPixelValue(px, py, vals);
                 
                // std::map<float, uint32_t>::const_iterator jt(idsOrderedByCoverage.end());
                std::map<float, float>::const_iterator tt(hashOrderedByCoverage.end());
                for (uint i=1; i<myRasters.size(); ++i, /*--jt,*/ --tt) {
                    // for (uint j=0; j<2; ++j, --jt, --tt) {
                    float val  = tt->first;
                    float hash = tt->second;
                    vals[2*0]   = hash;
                    vals[2*0+1] = val;
                    --tt;
                    val  = tt->first;
                    hash = tt->second;
                    vals[2*1]   = hash;
                    vals[2*1+1] = val;

                    myRasters(i)->setPixelValue(px, py, vals);
                    // }
            
                }
            }
            #endif

            // if (sampleMap.size())
            //     mySamples->write(px, py, gaussianNorm, sampleMap);

            for (int i = 0; i < vectorsize; ++i, ++destination)
                *destination = sample[i] / gaussianNorm;
           
        }
    }
}
