/*
 Partialy taken from HDK examples:
 http://www.sidefx.com/docs/hdk15.5/_v_r_a_y_2_v_r_a_y__demo_edge_detect_filter_8h-example.html

 Contains snippets from PsyOp
 [1] - Jonah Friedman, Andrew C. Jones, Fully automatic ID mattes with support for motion blur and transparency.
 */
//HDK
#include <UT/UT_DSOVersion.h>
#include <VRAY/VRAY_SpecialChannel.h>
#include <UT/UT_Args.h>
#include <UT/UT_StackBuffer.h>
#include <SYS/SYS_Floor.h>
#include <SYS/SYS_Math.h>
#include <UT/UT_Thread.h>
#include <UT/UT_PointGrid.h>
#include <GU/GU_Detail.h>

//STD
#include <iostream>
#include <map>
#include <memory>
#include <limits>
#include <cmath>
#include <queue>

//OWN
#include "VRAY_AutomattesFilter.hpp"
#include "AutomattesHelper.hpp"

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
    , myFilterWidth(2)
    , myGaussianAlpha(1)
    , myGaussianExp(0)
    , myRank(0)
    , myHashTypeName("crypto")
    , myHashType(CRYPTO)
    , myIdTypeName("object")
    , myIdType(OBJECT)
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
    args.stripOptions("w:r:i:h:");

    if (args.found('w')) { myFilterWidth = args.fargp('w'); }
    if (args.found('r')) { myRank        = args.fargp('r'); }

    // hash type
    if (args.found('h')) { 
        myHashTypeName = args.argp('h'); 
        if (std::string(myHashTypeName).compare("mantra") == 0)
            myHashType = Automatte_HashType::MANTRA;
         else 
            myHashType = Automatte_HashType::CRYPTO;
    }

    // id type
    if (args.found('i')) { 
        myIdTypeName  = args.argp('i');
        if (std::string(myIdTypeName).compare("asset") == 0)
            myIdType = Automatte_IdType::ASSET;
        else if (std::string(myIdTypeName).compare("object") == 0)
            myIdType = Automatte_IdType::OBJECT;
        else if (std::string(myIdTypeName).compare("material") == 0)
            myIdType = Automatte_IdType::MATERIAL;
        else if (std::string(myIdTypeName).compare("group") == 0)
            myIdType = Automatte_IdType::GROUP;
    }

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
    
    if (myHashType == MANTRA) {
        addSpecialChannel(imager, VRAY_SPECIAL_OPID);
        addSpecialChannel(imager, VRAY_SPECIAL_MATERIALID);
    }

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

    UT_ASSERT(vectorsize == 4);

    const float *const colordata = getSampleData(source, channel);

    #ifndef VEXSAMPLES
        const float *const Object_ids  = (myHashType == MANTRA) ? \
            getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_OPID)) : NULL;
        const float *const Material_ids = (myHashType == MANTRA) ? \
            getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_MATERIALID)) : NULL;
         // Resolution convention R: Asset*, G: Object, B: Material, A: group*.
         // * - not supported yet.
        const int hash_index = (myIdType == OBJECT) ? 1 : 2;
    #endif

    #ifdef VEXSAMPLES

        AutomatteVexCache * vex_cache = get_AutomatteVexCache();
        AutomatteImage    * vex_image = get_AutomatteImage(); // ?
        VEX_Samples       * samples   = nullptr;

        AutomatteVexCache::const_accessor  channel_reader;
        if(vex_cache->find(channel_reader, AUTOMATTE_CHANNEL_HASH)) {
            const VEX_Samples & ref = channel_reader->second;
            samples = const_cast<VEX_Samples*>(&ref);
        } else {
            // TODO: Do something about it;
            UT_ASSERT(false);
        }
        
        SampleBucket * bucket  = nullptr;
        const int thread_id    = SYSgetSTID();
        
        VEX_Samples::accessor bucket_queue;
        if (samples->find(bucket_queue, thread_id)) {
            // TODO: move to reference, we don't need pointer anymore?.
            SampleBucket & bucket_ref = bucket_queue->second.front();
            bucket = static_cast<SampleBucket*>(&bucket_ref); 
        } else {
            // FIXMED ?
            UT_ASSERT(false);
        }

        if (bucket->size() != 0 && bucket->isRegistered() == 0) {
            bucket->registerBucket();
        } 
        

    #endif
    
    // Run over destination pixels
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
            float   gaussianNorm = 0;

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
                        const float gaussianWeight = gaussianFilter(x*1.66667, y*1.66667, myGaussianExp, myGaussianAlpha);
                    
                        #ifdef VEXSAMPLES

                        const float sx = colordata[vectorsize*sourceidx+0]; // G&B are reserved for id and coverage by bellow setup
                        const float sy = colordata[vectorsize*sourceidx+3]; // se we end up with using R&A for NDC coords.
                        const float pxf = sx * bucket->m_resolution[0] * bucket->m_pixelsamples[0];
                        const float pyf = sy * bucket->m_resolution[1] * bucket->m_pixelsamples[1];
                        const size_t subpxi = std::floor(pxf);
                        const size_t subpyi = std::floor(pyf);
                        const size_t index  = subpyi * bucket->m_resolution[0] * bucket->m_pixelsamples[0] + subpxi;

                        if (index < vex_image->size()) {
                            const Sample & vexsample = vex_image->at(index);
                            const float _id = vexsample[3];
                            // FIXME: cov. should be a sum of all samples behind the current one. (Pz>current sample)
                            const float coverage = vexsample[4] * gaussianWeight; 
                            #ifdef HALTON_FALSE_COLORS
                            // borrowed from: https://github.com/MercenariesEngineering/openexrid
                            const float primes[3] = {2,3,5};
                                sample[0] += gaussianWeight * halton(primes[0], _id);
                                sample[1] += gaussianWeight * halton(primes[1], _id);
                                sample[2] += gaussianWeight * halton(primes[2], _id); 
                            #else
                                uint seed  = static_cast<uint>(_id);
                                sample[1] += gaussianWeight * SYSfastRandom(seed);
                                     seed += 2345;
                                sample[2] += gaussianWeight * SYSfastRandom(seed); 
                            #endif

                            gaussianNorm += gaussianWeight;

                            // Sumarize coverage per hash.
                            if (hash_map.find(_id) == hash_map.end()) {
                                hash_map.insert(std::pair<float, float>(_id, coverage));
                            } else {
                                hash_map[_id] += coverage;
                            }
                        } else {
                            // This should not happen at all:
                            sample[0]  = gaussianWeight * 1.f;
                        }

                        #else

                        gaussianNorm += gaussianWeight;
                        // no transparency support (because of precomposed shader samples).
                        const float coverage = 1.f * gaussianWeight; //fixme
                        // This is ugly, fixme
                        const float object_id   = (myHashType == MANTRA) ? \
                            Object_ids[sourceidx] : colordata[vectorsize*sourceidx+hash_index];    // G -> object_id
                        const float material_id = (myHashType == MANTRA) ? \
                            Material_ids[sourceidx] : colordata[vectorsize*sourceidx+hash_index];  // B -> material_id

                        const float _id = (myIdType == OBJECT) ? object_id : material_id; 

                        // 
                        #ifdef HALTON_FALSE_COLORS
                        // borrowed from: https://github.com/MercenariesEngineering/openexrid
                        const float primes[3] = {2,3,5};
                            sample[0] += gaussianWeight * halton(primes[0], _id);
                            sample[1] += gaussianWeight * halton(primes[1], _id);
                            sample[2] += gaussianWeight * halton(primes[2], _id); 
                        #else
                            uint seed  = static_cast<uint>(_id);
                            sample[1] += gaussianWeight * SYSfastRandom(seed);
                                 seed += 2345;
                            sample[2] += gaussianWeight * SYSfastRandom(seed); 
                        #endif
                        
                        // 
                        if (hash_map.find(_id) == hash_map.end()) {
                            hash_map.insert(std::pair<float, float>(_id, coverage));
                        }
                        else {
                            hash_map[_id]   += coverage;
                        } 

                        #endif // end of VEXSAMPLES
                    }
                }
            }
            
            HashMap coverage_map;
            // sort by coverage 
            HashMap::const_iterator it(hash_map.begin());
            for(; it != hash_map.end(); ++it) {
                const float coverage  = it->second;// / gaussianNorm;
                const float object_id = it->first;// ? coverage != 0.0f: 0.f;
                coverage_map.insert(std::pair<float, float>(coverage, object_id));
            }

            HashMap::const_reverse_iterator rit(coverage_map.rbegin());
            if (myRank == 0) {
                // sample[2] = rit->second;
                for (int i = 0; i< vectorsize; ++i, ++destination) {
                    *destination  = sample[i] / gaussianNorm; 
                }
                    
            } else {
        
                const size_t id_offset = (static_cast<size_t>(myRank) - 1) * 2; 
                std::advance(rit, id_offset);

                for (int i = id_offset; i < 2; ++rit,  ++i) {
                    if (rit == coverage_map.rend()) {
                        destination += i*2;
                        break;
                    }
                    destination[0] = rit->second; // object_id
                    destination[1] = rit->first / gaussianNorm; // coverage
                    destination += 2;
                }   
                       
            }
         // end of canonical way 
        }
    }
    // end of destx/desty loop;


    #ifdef VEXSAMPLES 
        // DEBUG_PRINT("Filter thread: %i,\n", thread_id);
        bucket->clear();
        SampleBucket newbucket;
        newbucket.copyInfo(bucket);
        bucket_queue->second.pop();
        bucket_queue->second.push(newbucket);
    #endif

}



