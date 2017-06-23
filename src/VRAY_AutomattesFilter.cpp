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
#include <UT/UT_Thread.h>
#include <UT/UT_PointGrid.h>


#include <iostream>
#include <map>
#include <unordered_map>
#include <memory>

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

    // VEX_Samples * samples = VEX_Samples_get();
    // std::cout << "Filter: " << samples << std::endl;
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

    // It's not technically necessery, but some convention needs to be taken.
    UT_ASSERT(vectorsize == 4);

    // temporarly to test vex->filter passing.
    #ifdef VEXSAMPLES

    std::unique_ptr<UT_PointGrid<UT_Vector3Point>> pixelgrid(nullptr);
    UT_Vector3Array  positions;
    UT_ValArray<int> indices;
    VEX_Samples  * samples = VEX_Samples_get();
    SampleBucket * bucket  = nullptr;

    // std::cout << "Filter: " << samples << "\n";
    const int thread_id = UT_Thread::getMyThreadId();
    VEX_Samples::const_iterator it = samples->find(thread_id);
    bool use_vex_samples = false;

    if (it != samples->end()) {
        const BucketQueue  & queue  = samples->at(thread_id);
        BucketQueue::const_iterator jt = queue.begin();
        for (; jt != queue.end(); ++jt) {
            if (jt->size() != 0) {
                bucket = &(*jt);
                break; 
            }
        }
        // const SampleBucket & bucket = *jit;
        const size_t size = bucket->size();
        UT_Vector3F bucket_min = {0,0,0};
        UT_Vector3F bucket_max = {0,0,0};
        positions.bumpSize(size);
        indices.bumpSize(size);

        for (int i=0; i<size; ++i) {
            const Sample & vexsample = bucket->at(i);
            UT_Vector3 pos = {vexsample[0], vexsample[1], vexsample[2]};
            bucket_min = SYSmin(pos, bucket_min);
            bucket_max = SYSmax(pos, bucket_max);
            positions.append(pos);
            indices.append(i);
        }
        // Enlarage bbox a bit.
        UT_Vector3 bucket_size(bucket_max - bucket_min);
        UT_Vector3 _sigma(bucket_size*.01);
        bucket_size += _sigma;
        bucket_min  +- _sigma;

        // Point Grid structures/objects:
        UT_Vector3Point accessor(positions, indices);
        pixelgrid = std::unique_ptr<UT_PointGrid<UT_Vector3Point>>\
            (new UT_PointGrid<UT_Vector3Point> (accessor));

        // 
        if (pixelgrid->canBuild(destwidth, destheight, 1)) {
            // Build it:
            pixelgrid->build(bucket_min, bucket_size, destwidth+2, destheight+2, 1);
            use_vex_samples = true;
        }
    }
    
    #endif

    const float *const colordata = getSampleData(source, channel);
    const float *const Object_ids  = (myHashType == MANTRA) ? \
        getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_OPID)) : NULL;
    const float *const Material_ids = (myHashType == MANTRA) ? \
        getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_MATERIALID)) : NULL;

     // Resolution convention R: Asset*, G: Object, B: Material, A: group*.
     // * - not supported yet.
    const int hash_index = (myIdType == OBJECT) ? 1 : 2;

    #ifdef VEXSAMPLES 
    UT_Vector3PointQueue *queue;
    queue = pixelgrid->createQueue();
    UT_PointGridIterator<UT_Vector3Point> iter;
    // BucketQueue  & queue  = samples->at(thread_id);
    // SampleBucket & bucket = samples->at(thread_id);
    #endif

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

            #ifdef VEXSAMPLES 
            // temporarly to test vex->filter passing.
            if (use_vex_samples) {
                iter = pixelgrid->getKeysAt(destx, desty, 0, *queue);
                UT_Vector3 color = {0.0f, 0.0f, 0.0f};
                float alpha = 0.f;
                if (iter.entries() != 0) {
                    for (;!iter.atEnd(); iter.advance()) {
                        const size_t idx = iter.getValue();
                        UT_ASSERT(idx < bucket->size());
                        const Sample & vexsample = bucket->at(idx); //something wrong is here.
                        color += UT_Vector3(vexsample[0], vexsample[1], vexsample[2]);
                        alpha += vexsample[4];
                        const int mythreadid = static_cast<int>(vexsample[5]);
                        if (thread_id != mythreadid)
                            DEBUG_PRINT("%i : %i", thread_id, mythreadid);
                    }
                    color /= iter.entries();
                    // alpha /= iter.entries();
                }

                sample[0] = color.x();
                sample[1] = color.y();
                sample[2] = iter.entries()*1.f;//color.z();
                sample[3] = alpha;
            }

            #endif



            HashMap hash_map;
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
                        const float x = (float(sourcex)-0.5f*float(sourcelastx + sourcefirstx))\
                            / float(mySamplesPerPixelX);
                        const float y = (float(sourcey)-0.5f*float(sourcelasty + sourcefirsty))\
                            / float(mySamplesPerPixelY);

                        // TODO: remove magic number
                        const float gaussianWeight = gaussianFilter(x*1.66667, y*1.66667, myGaussianExp, \
                            myGaussianAlpha);
                        gaussianNorm += gaussianWeight;
                        

                        // For now this is our coverage sample, which means 
                        // no transparency support (because of precomposed shader samples).
                        // thus we can assume 1 instead of sampling alphs / opacity.
                        const float coverage = 1.f * gaussianWeight; //fixme

                        // This is ugly, fixme
                        const float object_id   = (myHashType == MANTRA) ? \
                            Object_ids[sourceidx] : colordata[vectorsize*sourceidx+hash_index];    // G -> object_id
                        const float material_id = (myHashType == MANTRA) ? \
                            Material_ids[sourceidx] : colordata[vectorsize*sourceidx+hash_index];  // B -> material_id

                        const float _id = (myIdType == OBJECT) ? object_id : material_id; 

                        // temporarly to test vex->filter passing.
                        if (!use_vex_samples) {
                            uint seed  = static_cast<uint>(_id);
                            sample[1] += gaussianWeight * SYSfastRandom(seed);
                                 seed += 2345;
                            sample[2] += gaussianWeight * SYSfastRandom(seed); 
                        }

                        if (hash_map.find(object_id) == hash_map.end()) {
                            hash_map.insert(std::pair<float, float>(_id, coverage));
                        }
                        else {
                            hash_map[object_id]   += coverage;
                        } 
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
                //sample[2] = rit->second;
                for (int i = 0; i< vectorsize; ++i, ++destination) {
                    *destination  = sample[i];// / gaussianNorm; 
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
        }
    }

    #ifdef VEXSAMPLES 
    // end of destx/desty loop;
    pixelgrid->destroyQueue(queue);
    const int myBucketCounter = VEX_Samples_increamentBucketCounter(thread_id);
    DEBUG_PRINT("Filter thread: %i, bucket count:%i (size: %lu)\n", thread_id, myBucketCounter, bucket->size());
    BucketQueue  & bqueue = samples->at(thread_id);
    BucketQueue::iterator kt = bqueue.begin();
    SampleBucket new_bucket;
    bqueue.insert(kt, new_bucket);
    // bucket.clear();
    #endif
}



