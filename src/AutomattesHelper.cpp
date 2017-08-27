#include <iostream>
#include <map>
#include <mutex>
#include <cstring>
#include <functional>
#include <memory>
#include <queue>
// #include <tuple>


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

// make it atomic;
static ut_thread_id_t mainThreadId = 0;

static const size_t BucketQueueCapacity = 1024;
static const size_t max_opacity_samples = 1;

/// AutomatteVexCache stores buckets per thread per channel
static AutomatteVexCache atm_vex_cache;
// This is output grid composed from buckets in pixel filter.
static AutomatteImage    atm_image;


const Sample & SampleBucket::at(const int & index) const 
{
    const int size = mySamples.size();
    UT_ASSERT (index < size) ;
    return mySamples.at(index);
}

const size_t SampleBucket::getNeighbourSize() const noexcept
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    size_t size = 0;
    return myNeighbours.size();
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


void SampleBucket::copyInfo(const SampleBucket * bucket) noexcept
{
    m_resolution[0] = bucket->m_resolution[0];
    m_resolution[1] = bucket->m_resolution[1];
    m_pixelsamples[0] = bucket->m_pixelsamples[0];
    m_pixelsamples[1] = bucket->m_pixelsamples[1];
}


void SampleBucket::copyInfo(const std::vector<int> & res, const std::vector<int> & samples) noexcept
{
    m_resolution[0] = res.at(0);
    m_resolution[1] = res.at(1);
    m_pixelsamples[0] = samples.at(0);
    m_pixelsamples[1] = samples.at(1);
}

size_t SampleBucket::registerBucket() 
{
    const size_t size = mySamples.size();
    for(int i=0; i < size; ++i) {
        const Sample vexsample = mySamples.at(i);
        const float pxf = vexsample[0] * m_resolution[0] * m_pixelsamples[0];
        const float pyf = vexsample[1] * m_resolution[1] * m_pixelsamples[1];
        const size_t  subpxi = std::floor(pxf);
        const size_t  subpyi = std::floor(pyf);
        const size_t index = subpyi * m_resolution[0] * m_pixelsamples[0] + subpxi;
        if(index < atm_image.size()) {
            atm_image.at(index) = vexsample;
        } else {
            // place empty sample here?
        }
    }

    myRegisteredFlag = 1;
    return atm_image.size();
}


int create_vex_storage(const std::string & channel_name, const int & thread_id, \
    const std::vector<int> & res, const std::vector<int> & samples)
{

    std::lock_guard<std::mutex> guard(automattes_mutex);
    const ut_thread_id_t currentMainThreadId = UT_Thread::getMainThreadId();

    if (atm_image.capacity() == 0)
        atm_image.resize(res[0]*res[1]*samples[0]*samples[1]*max_opacity_samples);

    if (currentMainThreadId != mainThreadId) {
        atm_vex_cache.clear();
        AutomatteImage tmp;
        tmp.resize(res[0]*res[1]*samples[0]*samples[1]*max_opacity_samples);
        atm_image.swap(tmp);
        mainThreadId = currentMainThreadId;
    }

    AutomatteVexCache::accessor channel_writer;
    VEX_SamplesQ::accessor       bucket_queue_writer;
    std::hash<std::string>      hasher;

    // It's a bug, right? We need int32_t for VEX, size_t won't fit
    // we could take murmur hasher...
    const size_t  channel_hash_64 = hasher(channel_name);
    const int32_t channel_hash_32 = static_cast<int32_t>(channel_hash_64);

    if (!atm_vex_cache.find(channel_writer, channel_hash_32)) 
    {
        VEX_SamplesQ vex_samples;
        BucketQueueQ queue;
        SampleBucket bucket; 
        bucket.copyInfo(res, samples);
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
            bucket.copyInfo(res, samples);   
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
    VEX_SamplesQ::accessor      bucket_reader;
    atm_vex_cache.find(channel_reader, handle);
    channel_reader->second.find(bucket_reader, thread_id);
    SampleBucket & it = bucket_reader->second.front();
    it.push_back(sample);
    return it.size();
}


AutomatteVexCache * get_AutomatteVexCache()
{ 
    return &atm_vex_cache; 
}

AutomatteImage * get_AutomatteImage()
{
    return &atm_image;
}

} // end of HA_HDK