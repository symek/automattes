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
    #if 0
    // debug: check if samples are consistant
    GU_Detail gdp, gdp2;
    long npoints = 0;
    long npoints2 = 0;
    VEX_Samples  * samples = VEX_Samples_get();
    // UT_PointGrid<UT_Vector3Point> pixelgrid(nullptr);
    

    DEBUG_PRINT(" VEX_Samples size: %lu\n", samples->size());
    const int destwidth  = 64;
    const int destheight = 64;

    // const int thread_id = UT_Thread::getMyThreadId();
    const VEX_Samples::const_iterator vexit = samples->begin();
    for (; vexit!= samples->end(); ++vexit) {
        const BucketQueue & queue = vexit->second;
        const BucketQueue::const_iterator queit = queue.begin();
        for (; queit != queue.end(); ++queit) {
            const SampleBucket & bucket = (*queit);

            if (bucket.size() == 0)
                continue;

            const int size = bucket.size();
            UT_Vector3F bucket_min = {FLT_MAX,FLT_MAX,FLT_MAX};
            UT_Vector3F bucket_max = {FLT_MIN,FLT_MIN,FLT_MIN};
            UT_Vector3Array  positions(size);
            UT_ValArray<int> indices(size);

            const SampleBucket::const_iterator buit = bucket.begin();
            for (int i=0; i< size; ++i) {
                npoints++;
                // const size_t size = bucket.size();
                // positions.bumpSize(size);
                // indices.bumpSize(size);
                const Sample & sample = bucket.at(i);
                const UT_Vector3 pos = {sample[0], sample[1], 0.f};
                GA_Offset ptoff = gdp.appendPoint();
                gdp.setPos3(ptoff, pos);
                bucket_min = SYSmin(pos, bucket_min);
                bucket_max = SYSmax(pos, bucket_max);
                positions.append(pos);
                indices.append(i);
            }

            // Enlarage bbox a bit.
            // UT_Vector3 bucket_size(bucket_max - bucket_min);
            // UT_Vector3 _sigma(bucket_size*.1);
            // bucket_size += _sigma;
            // bucket_min  -= _sigma;
            UT_Vector3Point accessor(positions, indices);
            UT_PointGrid<UT_Vector3Point> pixelgrid(accessor);

            // 

            if (pixelgrid.canBuild(destwidth, destheight, 1)) {
                // Build it:
                pixelgrid.build(bucket_min, bucket_size, destwidth, destheight, 1);
                // DEBUG_PRINT("Build pixelgrid: ", 0);
                
            } else {
                DEBUG_PRINT("Can't build. \n", 0);
            }

            UT_Vector3PointQueue *queue;
            queue = pixelgrid.createQueue();
            UT_PointGridIterator<UT_Vector3Point> iter;
            // npoints2++;
            for (int desty=0; desty < destheight; ++desty) {
                for (int destx=0; destx < destwidth; ++destx) {
                    iter = pixelgrid.getKeysAt(destx, desty, 0, *queue);    
                    for (;!iter.atEnd(); iter.advance()) {
                        npoints2 ++;
                        const size_t idx = iter.getValue();
                        // UT_ASSERT(idx < bucket.size());
                        const Sample & sample = bucket.at(idx); 
                        const UT_Vector3 pos = {sample[0], sample[1], sample[2]};
                        GA_Offset ptoff = gdp2.appendPoint();
                        gdp2.setPos3(ptoff, pos);
                    }
                    
                }
            }   

            pixelgrid.destroyQueue(queue);
        }
    }
    DEBUG_PRINT("Saving %lu and %lu points in total.\n", npoints, npoints2);
    gdp.save("/tmp/testauto.bgeo", 0);
    gdp2.save("/tmp/testauto2.bgeo", 0);

    #endif


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

void VRAY_AutomatteFilter::updateSourceBoundingBox(
    const int & destwidth, 
    const int & destheight,
    const int & sourcewidth,
    const int & sourceheight,
    const int & destxoffsetinsource,
    const int & destyoffsetinsource,
    const int & vectorsize,
    const float * colordata, 
    UT_BoundingBox * sourceBbox) const
{
    UT_Vector3 source_min = {FLT_MAX, FLT_MAX, FLT_MAX};
    UT_Vector3 source_max = {FLT_MIN, FLT_MIN, FLT_MIN};
    for (int desty = 0; desty < destheight; ++desty) 
    {
        for (int destx = 0; destx < destwidth; ++destx) 
        {
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

            for (int sourcey = sourcefirstry; sourcey <= sourcelastry; ++sourcey)
            {
                for (int sourcex = sourcefirstrx; sourcex <= sourcelastrx; ++sourcex) 
                {
                    const int sourceidx = sourcex + sourcewidth*sourcey;
                    if(sourcex >= sourcefirstox && sourcex <= sourcelastox &&\
                      sourcey >= sourcefirstoy && sourcey <= sourcelastoy) {
                        const float sx = colordata[vectorsize*sourceidx+0]; // G&B are reserved for id and coverage by bellow setup
                        const float sy = colordata[vectorsize*sourceidx+3]; // se we end up with using R&A for NDC coords.
                        const UT_Vector3 position = {sx, sy, 0.f};
                        source_min = SYSmin(source_min, position);
                        source_max = SYSmax(source_max, position);

                    }
                }
            }
        }
    }

    source_min.z()  = -.01;
    source_max.z()  =  .01;
    sourceBbox->initBounds(source_min, source_max);
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
    const float *const Object_ids  = (myHashType == MANTRA) ? \
        getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_OPID)) : NULL;
    const float *const Material_ids = (myHashType == MANTRA) ? \
        getSampleData(source, getSpecialChannelIdx(imager, VRAY_SPECIAL_MATERIALID)) : NULL;

     // Resolution convention R: Asset*, G: Object, B: Material, A: group*.
     // * - not supported yet.
    const int hash_index = (myIdType == OBJECT) ? 1 : 2;

    #ifdef VEXSAMPLES
    const int sourcetodestwidth  = sourcewidth;// / mySamplesPerPixelX;
    const int sourcetodestheight = sourceheight;// / mySamplesPerPixelY;

    // std::unique_ptr<UT_PointGrid<UT_Vector3Point>> pixelgrid(nullptr);
    UT_Vector3Array  positions;
    UT_ValArray<int> indices;
    VEX_Samples  * samples = VEX_Samples_get();
    SampleBucket * bucket  = nullptr;

    int offset = 0;
    int bucket_threadid = 0;
    int fullbuckets = 0;
    int foundDeepSamples = 0;
    int bucketgridsize = 0;
    int bucketsFoundInStore = 0;

    // find first non empty bucket. this is wrong way
    // buckets are empty for background, needs to find a way ot find them.
    // const VEX_Samples::const_iterator it = samples->find(thread_id);
    const int thread_id = SYSgetSTID();
    const int myBucketCounter = VEX_Samples_increamentBucketCounter(thread_id);

    UT_BoundingBox sourcebbox;
    updateSourceBoundingBox(destwidth, destheight, sourcewidth, sourceheight, 
        destxoffsetinsource, destyoffsetinsource, vectorsize, colordata, &sourcebbox);

    // const int q = VEX_getBucket(thread_id, bucket, offset);
    // std::cout << "q: " << bucket << "\n";
    UT_ASSERT(samples->find(thread_id) != samples->end());
    {
        const BucketQueue  & queue  = samples->at(thread_id);
        BucketQueue::const_iterator jt = queue.begin();
        for (; jt != queue.end(); ++jt, ++offset) {
            if (jt->size() != 0) {
                bucket = &(*jt);
                break; 
            }
        }
    }
    if (offset == 0) {
        bucketgridsize = bucket->updateBoundingBox(0.0f, 0.0f, 0.01f);
    } else {
        // find lets try to find source:
        // const float xmin = sourcebbox.minvec().x();
        // const float ymin = sourcebbox.minvec().y();
        // const float xmax = sourcebbox.maxvec().x();
        // const float ymax = sourcebbox.maxvec().y();
        const UT_Vector3 min = sourcebbox.minvec();
        const UT_Vector3 max = sourcebbox.maxvec();
        
        // const SampleBucket * oldbucket  = bucket;
        bucketsFoundInStore = bucket->findBucket(min, max, bucket);

        // if (bucketsFoundInStore != 0)
        //     DEBUG_PRINT("Suspected buckets: %i\n",bucketsFoundInStore);

        // if (oldbucket != bucket)
        //     DEBUG_PRINT("%i,\n", bucket);
        
    }
    

    const size_t bucket_size = bucket->size() + bucket->getNeighbourSize();
    positions.bumpSize(bucket_size);
    indices.bumpSize(bucket_size);

    for (int i=0; i<bucket_size; ++i) {
        const Sample & vexsample = bucket->at(i);
        const UT_Vector3 pos = {vexsample[0], vexsample[1], 0.f}; // we ommit Pz, to flatten grid.
        positions.append(pos);
        indices.append(i);
    }

    // Point Grid structures/objects (we could use persistand grid per bucket once we move to bucket class)
    UT_Vector3Point accessor(positions, indices);
    UT_PointGrid<UT_Vector3Point> pixelgrid(accessor);
    // pixelgrid = std::unique_ptr<UT_PointGrid<UT_Vector3Point>> (new UT_PointGrid<UT_Vector3Point> (accessor));

    // 
    if (pixelgrid.canBuild(sourcetodestwidth, sourcetodestheight, 1)) {
        pixelgrid.build(bucket->getBBox()->minvec(), bucket->getBBox()->size(),\
         sourcetodestwidth, sourcetodestheight, 1);
    }

    UT_Vector3PointQueue *queue;
    queue = pixelgrid.createQueue();
    UT_PointGridIterator<UT_Vector3Point> iter;

    #endif
    

    #if 0
    UT_StackBuffer<float> sample(vectorsize);
    for (int i = 0; i < vectorsize; ++i)
        sample[i] = 0.f;

    const float insidemin = sourcebbox.isInside(bucket->getBBox()->minvec());
    const float insidemax = sourcebbox.isInside(bucket->getBBox()->maxvec());
    const float distmin = sourcebbox.minvec().distance2(bucket->getBBox()->minvec());
    const float distmax = sourcebbox.maxvec().distance2(bucket->getBBox()->maxvec()); 

    // sourcebbox.isInside(*(bucket->getBBox()));
    sample[0] = distmin;// 1.f * insidemin * insidemax; //bucket->getBBox()->minvec().x();
    sample[1] = distmax;//bucket->getBBox()->minvec().y();
    sample[2] = sourcebbox.minvec().x();
    sample[3] = sourcebbox.minvec().y();

    for (int desty = 0; desty < destheight; ++desty) {
        for (int destx = 0; destx < destwidth; ++destx) {
             for (int i = 0; i< vectorsize; ++i, ++destination)
                    *destination  = sample[i];
        }
    }


    #else
    // Run over destination pixels
    for (int desty = 0; desty < destheight; ++desty) 
    {
        // const int pixelgridoffsety = (destyoffsetinsource/mySamplesPerPixelY) + desty;

        for (int destx = 0; destx < destwidth; ++destx)
        {
            // const int pixelgridoffsetx = (destxoffsetinsource/mySamplesPerPixelX) + destx;
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
                        const float gaussianWeight = gaussianFilter(x*1.66667, y*1.66667, myGaussianExp, myGaussianAlpha);
                    
                        #ifdef VEXSAMPLES

                        const float sx = colordata[vectorsize*sourceidx+0]; // G&B are reserved for id and coverage by bellow setup
                        const float sy = colordata[vectorsize*sourceidx+3]; // se we end up with using R&A for NDC coords.
                        const UT_Vector3 position = {sx, sy, 0.f};
                        const float radius = 0.000001f;


                        // int idx, idy, idz;
                        // const bool found_voxel = pixelgrid.posToIndex(position, idx, idy, idz, true);
                        // iter = pixelgrid.getKeysAt(idx, idy, idz, *queue);
                        iter = pixelgrid.findCloseKeys(position, *queue, radius);
                        const int entries = SYSmax((float)iter.entries(), 1.f);
                        foundDeepSamples += (iter.entries() - 1);
                        const float repEntries = 1.f/(float)entries; 
                        gaussianNorm += (gaussianWeight*entries);

                        // TMP pseudo color to check offset:
                        if (iter.entries() == 0) {
                            sample[0] += float(offset) * gaussianWeight;
                        } else {
                            for (;!iter.atEnd(); iter.advance()) {
                                const size_t idx = iter.getValue();
                                UT_ASSERT(idx < bucket_size);
                                const Sample & vexsample = bucket->at(idx);
                                const float _id =  vexsample[3];
                                // FIXME: cov. should be a sum of all samples behind the current one. (Pz>current sample)
                                const float coverage = vexsample[4] * gaussianWeight; 
                                uint seed  = static_cast<uint>(_id);
                                sample[1] += gaussianWeight * SYSfastRandom(seed);
                                     seed += 2345;
                                sample[2] += gaussianWeight * SYSfastRandom(seed); 


                                if (hash_map.find(_id) == hash_map.end()) {
                                    hash_map.insert(std::pair<float, float>(_id, coverage));
                                }
                                else {
                                    hash_map[_id] += coverage;
                                }
                            }
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
                        uint seed  = static_cast<uint>(_id);
                        sample[1] += gaussianWeight * SYSfastRandom(seed);
                             seed += 2345;
                        sample[2] += gaussianWeight * SYSfastRandom(seed); 
                        
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

    #endif

    #ifdef VEXSAMPLES 
    // end of destx/desty loop;
    pixelgrid.destroyQueue(queue);
    DEBUG_PRINT("Filter thread: %i, bucket count:%i (size: %lu) (offset: %i), (dim: %i, %i), (deep: %i), (bucketgrid: %i), (neighbours: %i)\n", \
        thread_id, myBucketCounter, bucket_size, offset, destwidth, destheight, foundDeepSamples, bucketgridsize, bucket->getNeighbourSize());
    // BucketQueue  & bqueue = samples->at(thread_id);
    // BucketQueue::iterator kt = bqueue.begin();
    // SampleBucket new_bucket;
    // bqueue.insert(kt, new_bucket);
    // bucket.clear();
    bucket->clearNeighbours();
    VEX_Samples_insertBucket(thread_id);
    #endif


}



