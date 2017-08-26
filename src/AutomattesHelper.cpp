#include <iostream>
#include <map>
#include <mutex>
#include <cstring>
#include <functional>
#include <memory>
#include <queue>


#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <UT/UT_PointGrid.h>
#include <tbb/concurrent_vector.h>
#include "AutomattesHelper.hpp"

#if defined(NO_MUTEX_IN_BUCKETVECTOR) || defined(TBB_VEX_STORE)
#include <tbb/concurrent_hash_map.h>
#endif

namespace HA_HDK {

// used only for creating per thread storage
static std::mutex automattes_mutex;
static std::mutex automattes_mutex2;
// our main storage
static VEX_Samples vexsamples;
// 
// static BucketGrid bucketGridX;
// static BucketGrid bucketGridY;
//
static BucketCounter vexBucketCounter;
static BucketCounter vrayBucketCounter;
// static VEX_SampleClass vexsamplesC;
static BucketSize bucketSize = {0,0};
static bool bucketSizeSet = 0;
static ut_thread_id_t mainThreadId = 0;
static const size_t BucketQueueCapacity = 1024;
static BucketVector bucketVector;

static AutomatteVexCache atm_vex_cache;


int VEX_Samples_create(const int& thread_id)
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    const ut_thread_id_t currentMainThreadId = UT_Thread::getMainThreadId();

    if (currentMainThreadId != mainThreadId) {
        vexsamples.clear();
        bucketVector.clear();
        mainThreadId = currentMainThreadId;
    }

    #ifdef NO_MUTEX_IN_BUCKETVECTOR
    VEX_Samples::const_accessor ra;
    if (vexsamples.find(ra, thread_id))
        return thread_id;
    #else
    VEX_Samples::const_iterator it = vexsamples.find(thread_id);
    if(it != vexsamples.end())
        return thread_id;
    #endif

    BucketQueue queue;
    queue.reserve(BucketQueueCapacity);//?
    bucketVector.reserve(BucketQueueCapacity*256);
    SampleBucket bucket;    
    vexsamples.insert(std::pair<int, BucketQueue>(thread_id, queue));
    #ifdef NO_MUTEX_IN_BUCKETVECTOR
    VEX_Samples::accessor wa;
    vexsamples.find(wa, thread_id);
    wa->second.push_back(bucket);
    #else
    vexsamples[thread_id].push_back(bucket);
    #endif
    //debug
    vexBucketCounter.insert(std::pair<int, int>(thread_id, 0));
    vrayBucketCounter.insert(std::pair<int, int>(thread_id, 0));

    return thread_id;
} 

int VEX_Samples_insert(const int& thread_id, const Sample& sample)
{
    #ifdef NO_MUTEX_IN_BUCKETVECTOR
    VEX_Samples::const_accessor ra;
    const bool result = vexsamples.find(ra, thread_id);
    UT_ASSERT(result);
    BucketQueue::const_iterator jt = ra->second.begin();
    #else
    VEX_Samples::const_iterator it = vexsamples.find(thread_id);
    UT_ASSERT(it != vexsamples.end());
    BucketQueue::iterator jt = vexsamples[thread_id].begin();
    #endif

    jt->push_back(sample);
    const size_t size = jt->size();
    if (size == 1) {
        vexBucketCounter[thread_id] += 1;
        std::cout << "VEX    thread: " << thread_id << ", bucket count:" \
            << vexBucketCounter[thread_id] << "\n";  
    }
    return static_cast<int>(size); //thread_id;
}

void VEX_Samples_insertBucket(const int & thread_id)
{
    #ifdef NO_MUTEX_IN_BUCKETVECTOR
    VEX_Samples::accessor wa;
    const bool result = vexsamples.find(wa, thread_id);
    UT_ASSERT(result);
    BucketQueue & bqueue = wa->second;
    #else
    std::lock_guard<std::mutex> guard(automattes_mutex);
    BucketQueue & bqueue = vexsamples.at(thread_id);
    #endif
    BucketQueue::iterator kt = bqueue.begin();
    SampleBucket new_bucket;
    bqueue.insert(kt, new_bucket);
}

VEX_Samples * VEX_Samples_get() 
{
    return &vexsamples;
}

int VEX_Samples_increamentBucketCounter(const int& thread_id)
{
    vrayBucketCounter[thread_id] += 1;
    return vrayBucketCounter[thread_id];
}

BucketSize * VEX_getBucketSize() {
    return &bucketSize;
}

void VEX_setBucketSize(int x, int y) {
    std::lock_guard<std::mutex> guard(automattes_mutex);
    if (bucketSize[0] != 0 || bucketSize[1] != 0) 
        return;
    
    bucketSize[0] = x;
    bucketSize[1] = y;
}

int VEX_bucketSizeSet() { return bucketSizeSet; }

const Sample & SampleBucket::at(const int & index) const 
{
    const int size = mySamples.size();
    // const int idx = SYSmin(index, size-1);
    // return mySamples.at(idx);
    if (index < size) {
        return mySamples.at(index);
    } else {
        const int neighbour_size = myNeighbours.size();
        const int idx = SYSclamp(index-size, 0, neighbour_size-1);
        return myNeighbours.at(idx);
    }
}

const size_t SampleBucket::getNeighbourSize() const noexcept
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    size_t size = 0;
    return myNeighbours.size();
    // std::vector<SampleBucket*>::const_iterator it = myNeighbours.begin();
    // for (; it!= myNeighbours.end(); ++it) {
    //     size += (*it)->size();
    // }
    // return size;
}

void SampleBucket::clearNeighbours() noexcept
{ 
    std::lock_guard<std::mutex> guard(automattes_mutex);
    myNeighbours.clear(); 
}

void SampleBucket::clear() noexcept
{ 
    SampleBucketV tmp;
    mySamples.swap(tmp); 
}

void SampleBucket::updateBoundingBox(const float & expx, const float & expy, const float & expz) 
{
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    const size_t size = mySamples.size();
    // this should never happen;
    if (size == 0)
        return -1;

    UT_Vector3 bucket_min = { FLT_MAX,  FLT_MAX,  FLT_MAX};
    UT_Vector3 bucket_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for(int i=0; i < size; ++i) {
        const Sample sample = mySamples.at(i);
        const UT_Vector3 position = {sample[0], sample[1], 0.f};
        bucket_min = SYSmin(bucket_min, position);
        bucket_max = SYSmax(bucket_max, position);
    }
    myBbox.initBounds(bucket_min, bucket_max);
    myBbox.expandBounds(expx, expy, expz);
}

size_t SampleBucket::registerBucket() 
{
    // update BucketGrid;
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    // temporarly heavily inefficient just to prove the point.
    bucketVector.push_back(*this);
    myRegisteredFlag  = 1;
    return bucketVector.size();

    // const coord_t xmin = myBbox.minvec().x();
    // const coord_t ymin = myBbox.minvec().y();
    // BucketGrid::const_iterator xit = bucketGridX.find(xmin);
    // BucketGrid::const_iterator yit = bucketGridY.find(ymin);

    // if (xit == bucketGridX.end()) {
    //     // BucketLine line;
    //     // line.push_back(std::pair<coord_t, SampleBucket*>(xmax, this));
    //     bucketGridX.insert(std::pair<coord_t, SampleBucket*>(xmin, this));
    // } else {
    //     BucketLine line = it->second;
    //     BucketLine::const_iterator jt = line.find(xmax);
    //     if (jt == line.end()) {
    //         line.insert(std::pair<coord_t, SampleBucket*>(xmax, this));
    //     }
    // }
}

int SampleBucket::fillBucket(const UT_Vector3 & min, const UT_Vector3 & max, SampleBucket * bucket) 
{
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    int counter = 0;
    const UT_Vector3 dir(max - min);
    const UT_Vector3 max2(min.x(), max.y(), max.z());
    const UT_Vector3 min2(max.x(), min.y(), min.z());
    const BucketVector::const_iterator it = bucketVector.begin();
    for(; it!=bucketVector.end(); ++it) {
        // const SampleBucket * store = static_cast<SampleBucket*>(*it);
        const SampleBucket & store = *it;
        const UT_BoundingBox * bbox = store.getBBox();
        // if (bbox->isInside(min) || bbox->isInside(max) \
            // || bbox->isInside(min2) || bbox->isInside(max2)) {
        if (bbox->isLineInside(min, dir)) {
            counter++;
            size_t size = store.size();
            for(size_t i=0; i < size; ++i) {
                const Sample & vexsample = store.at(i);
                mySamples.push_back(vexsample);
            }
        }
    }

    return counter;
}


int create_vex_storage(const std::string & channel_name, const int & thread_id, \
    const std::vector<float> & res, const std::vector<float> & samples)
{
    AutomatteVexCache::accessor channel_writer;
    VEX_SamplesQ::accessor       bucket_queue_writer;
    std::hash<std::string>      hasher;

    // It's a bug, right? We need int32_t for VEX, size_t won't fit
    // we could take murmur hasher...
    const size_t  channel_hash_64 = hasher(channel_name);
    const int32_t channel_hash_32 = static_cast<int32_t>(channel_hash_64);
    std::cout << channel_hash_32 << " ";

    if (!atm_vex_cache.find(channel_writer, channel_hash_32)) 
    {
        VEX_SamplesQ vex_samples;
        BucketQueueQ queue;
        SampleBucket bucket;    
        // queue.reserve(BucketQueueCapacity);
        queue.push(bucket);
        vex_samples.insert(bucket_queue_writer, std::pair<int, BucketQueueQ>(thread_id, queue));
        atm_vex_cache.insert(channel_writer, std::pair<int32_t, VEX_SamplesQ>(channel_hash_32, vex_samples));
        return channel_hash_32;

    } else {
        VEX_SamplesQ & channel = channel_writer->second;
        channel_writer.release();
        if (!channel.find(bucket_queue_writer, thread_id)) {
            BucketQueueQ queue;
            SampleBucket bucket;    
            // queue.reserve(BucketQueueCapacity);
            queue.push(bucket);
            channel.insert(bucket_queue_writer, std::pair<int, BucketQueueQ>(thread_id, queue));  
        }
        return channel_hash_32; 
    }

    return -1;
}


int insert_vex_sample(const int32_t & handle, const int & thread_id, const Sample & sample)
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    AutomatteVexCache::accessor channel_reader;
    const bool channel_found = atm_vex_cache.find(channel_reader, handle);
    UT_ASSERT(channel_found);
    VEX_SamplesQ::accessor bucket_reader;
    const bool thread_found = channel_reader->second.find(bucket_reader, thread_id);
    UT_ASSERT(thread_found);
    // BucketQueue::const_iterator it = bucket_reader->second.begin();
    SampleBucket & it = bucket_reader->second.front();
    it.push_back(sample);
    const size_t size = it.size();
    return size;
}


AutomatteVexCache * get_AutomatteVexCache()
{ 
    return &atm_vex_cache; 
}


#if 0
 // const BucketGrid::const_iterator it = bucketGrid.begin();
 //    for (; it!= bucketGrid.end(); ++it) {
 //        if (it->first <= ymax) {
 //            const BucketLine & line = it->second;
 //            BucketLine::const_iterator jt = line.begin();
 //            for (; jt != line.end(); ++jt) {
 //                if(jt->first <= xmax) {
 //                    // std::cout << "selecting new bucket\n";
 //                    jt++;
 //                    if (jt == line.end())
 //                        jt--;
 //                    bucket = static_cast<SampleBucket*>(jt->second);
 //                    break;
 //                }
 //            } 
            
 //        }
 //    }
}


void SampleBucket::findBucket(const float & xmin, const float & ymin, 
    const float & xmax, const float & ymax, SampleBucket * bucket) const 
{
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    // const BucketGrid::const_iterator it = bucketGrid.begin();
 //    for (; it!= bucketGrid.end(); ++it) {
 //        if (it->first <= ymax) {
 //            const BucketLine & line = it->second;
 //            BucketLine::const_iterator jt = line.begin();
 //            for (; jt != line.end(); ++jt) {
 //                if(jt->first <= xmax) {
 //                    // std::cout << "selecting new bucket\n";
 //                    jt++;
 //                    if (jt == line.end())
 //                        jt--;
 //                    bucket = static_cast<SampleBucket*>(jt->second);
 //                    break;
 //                }
 //            } 
            
 //        }
 //    }
}

int VEX_getBucket(const int thread_id, SampleBucket * bucket, int & offset)
{
    // #ifdef NO_MUTEX_IN_BUCKETVECTOR
    // VEX_Samples::accessor wa;
    // const bool result = vexsamples.find(wa, thread_id);
    // UT_ASSERT(result);

    // #else
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    // UT_ASSERT(vexsamples.find(thread_id) != vexsamples.end());
    // const BucketQueue  & queue  = vexsamples.at(thread_id);
    // #endif
    // BucketQueue::const_iterator jt = queue.begin();
    // for (; jt != queue.end(); ++jt, ++offset) {
    //     if (jt->size() != 0) {
    //         bucket = &(*jt);
    //         break; 
    //     }
    // }

    // return offset;
}

int 
VEX_SampleClass::open_channel(const std::string & channel)
{
    VEX_Channels::const_iterator it = myChannels.find(channel);
    if (it == myChannels.end()) {
        const int index = myChannels.size();
        VEX_Samples samples;
        myChannels.insert(std::pair<std::string, VEX_Samples>(channel, samples));
        return index;
    } else {
        return -1; // TODO: should return an index of a channel;
    }
}

int 
VEX_SampleClass::insert_queue(const std::string & channel, const int & thread_id)
{
    // std::lock_guard<std::mutex> guard(automattes_mutex2);
    VEX_Channels::const_iterator chit = myChannels.find(channel);
    UT_ASSERT(chit != myChannels.end());

    if (chit == myChannels.end()) {
        const int result = open_channel(channel);
        UT_ASSERT(result != -1);
        if (result == -1)
            return -1;
    }

    VEX_Samples & vexsamples  = myChannels.at(channel);
    VEX_Samples::iterator sit = vexsamples.find(thread_id);

    if (sit == vexsamples.end()) {
        BucketQueue queue;
        vexsamples.insert(std::pair<int, BucketQueue>(thread_id, queue));
        return thread_id;
    } else {
        return 0;
    }

}



int
VEX_SampleClass::insert_bucket(const std::string & channel, const int & thread_id)
{
    VEX_Channels::const_iterator chit = myChannels.find(channel);
    UT_ASSERT(chit != myChannels.end());
    VEX_Samples vexsamples = chit->second;
    VEX_Samples::const_iterator sit = vexsamples.find(thread_id);
    UT_ASSERT(sit != vexsamples.end());
    BucketQueue  & bqueue = vexsamples.at(thread_id);
    BucketQueue::iterator kt = bqueue.begin();
    SampleBucket new_bucket;
    bqueue.insert(kt, new_bucket);
    return bqueue.size();

 }

// int
// VEX_SampleClass::insert_sample(const int &thread_id, const std::string& channel, const Sample& sample)
// {
//  VEX_Channels::const_iterator chit = myChannels.find(channel);
//  VEX_Samples samples = chit->second;
//  VEX_Samples::const_iterator sit = samples.find(thread_id);

// }
#endif

} // end of HA_HDK