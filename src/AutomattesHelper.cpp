#include <iostream>
#include <map>
#include <mutex>
#include <cstring>
#include <functional>
#include <memory>
#include <queue>
#include <atomic>

#define TBB_CACHE_ALIGNED_ALLOC
#define TBB_PREVIEW_MEMORY_POOL 1


#ifdef TBB_PREVIEW_MEMORY_POOL
#include "tbb/memory_pool.h"
#include "tbb/scalable_allocator.h"
#endif

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>

#ifdef TBB_CACHE_ALIGNED_ALLOC
#include "tbb/cache_aligned_allocator.h"
#include "tbb/scalable_allocator.h"
#endif

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <UT/UT_PointGrid.h>
#include <IMG/IMG_DeepShadow.h>
#include "AutomattesHelper.hpp"



namespace HA_HDK {

// used only for creating per thread storage
static std::mutex automattes_mutex;

// make it atomic;
static std::atomic<ut_thread_id_t> mainThreadId;

static const size_t BucketQueueCapacity = 1024;
static const size_t max_opacity_samples = 1;

/// AutomatteVexCache stores buckets per thread per channel
static AutomatteVexCache atm_vex_cache;
// This is output grid composed from buckets in pixel filter.
static AutomatteImage atm_image;
static ImageInfo      atm_image_info;

// Another experiment
#ifdef USE_DEEP_MAP
static IMG_DeepShadow dsm;
#endif

#ifdef TBB_PREVIEW_MEMORY_POOL
static tbb::memory_pool<tbb::scalable_allocator<SampleVector>> sample_memory_pool;
#endif

SampleBucket::SampleBucket()
{
    #ifdef TBB_PREVIEW_MEMORY_POOL
    m_samples =  new SampleVector(sample_vector_allocator_t(sample_memory_pool));
    #endif
}

const Sample & SampleBucket::at(const int & index) const 
{
    UT_ASSERT (index < m_samples->size()) ;
    return m_samples->at(index);
}


void SampleBucket::clear() noexcept
{ 

    m_samples->clear(); 
    myRegisteredFlag = 0;
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
    #ifdef USE_DEEP_MAP
    IMG_DeepPixelWriter writer(dsm);
    #endif

    const size_t size = m_samples->size();
    for(int i=0; i < size; ++i) {
        const Sample vexsample = m_samples->at(i);
        const float pxf = vexsample[0] * atm_image_info.gridresx;
        const float pyf = vexsample[1] * atm_image_info.gridresy;
        const int subpxi = std::floor(pxf) + atm_image_info.image_margin;
        const int subpyi = std::floor(pyf) + atm_image_info.image_margin;
        const int index  = subpyi *atm_image_info.gridresx + subpxi;

        if(index < atm_image.size() && index > 0) {
            Sample & pixel = atm_image.at(index);
            pixel = Sample(vexsample);
            // for (uint i = 0; i < vexsample.size(); ++i)
                // pixel.push_back(vexsample[i]);
        } else {
            // place empty sample here?
            DEBUG_PRINT("index not found: %i, sub_pix: (%i, %i), ndc: (%f, %f)\n", index, subpxi, subpyi, \
                vexsample[0], vexsample[1]);
        }

        #ifdef USE_DEEP_MAP
        writer.open(subpxi, subpyi);
        const float * sample = reinterpret_cast<const float*>(&vexsample);
        writer.write(0.0f, sample, 6, 0x00, 0);
        #endif

    }
    clear();
    #ifdef USE_DEEP_MAP
    writer.close();
    #endif
    myRegisteredFlag = 1;
    return atm_image.size();
}

bool ImageInfo::update_size(const std::vector<int> & res, 
                            const std::vector<int> & samples) noexcept
{
    // if ( res[0] != m_resolution[0] || res[1] != m_resolution[1] ||
        // samples[0] != m_samples[0] || samples[1] != m_samples[1] ) {
        gridresx = res[0] * samples[0] + 2 * image_margin * samples[0];
        gridresy = res[1] * samples[1] + 2 * image_margin * samples[1];
        image_size = gridresx * gridresy * max_samples;
        m_resolution = res;
        m_samples    = samples;
        // return true;
    // } 

    return false;
}

void close_vex_storage()
{
    #ifdef USE_DEEP_MAP
    dsm.close();
    #endif
}


int create_vex_storage(const std::string & channel_name, const int & thread_id, \
    const std::vector<int> & res, const std::vector<int> & samples)
{

    // std::lock_guard<std::mutex> guard(automattes_mutex);
    const ut_thread_id_t currentMainThreadId = UT_Thread::getMainThreadId();
    if (atm_image_info.image_size == 0) {
        atm_image_info.update_size(res, samples);
        mainThreadId = currentMainThreadId;
    }

    if (atm_image.capacity() == 0)
    {
        atm_image.resize(atm_image_info.image_size);
         mainThreadId = currentMainThreadId;

    }

    if (currentMainThreadId != mainThreadId &&\
        atm_image.size() != atm_image_info.image_size) {
       
        atm_vex_cache.clear();
        AutomatteImage tmp;
        atm_image.resize(atm_image_info.image_size);
        atm_image.swap(tmp);

        mainThreadId = currentMainThreadId;
    }


    // if (atm_image_info.image_size == 0) {
    //     atm_image_info.update_size(res, samples);
    //     atm_image.resize(atm_image_info.image_size);  
    //     mainThreadId = currentMainThreadId;
    // }

    // if (currentMainThreadId != mainThreadId && \
    //     atm_image_info.image_size != atm_image.size()) {
    //     atm_vex_cache.clear();
    //     atm_image.resize(atm_image_info.image_size);
    //     mainThreadId = currentMainThreadId;
    // }
     
    AutomatteVexCache::accessor channel_writer;
    VEX_Samples::accessor       bucket_queue_writer;
    std::hash<std::string>      hasher;

    // It's a bug, right? We need int32_t for VEX, size_t won't fit
    // we could take murmur hasher...
    const size_t  channel_hash_64 = hasher(channel_name);
    const int32_t channel_hash_32 = static_cast<int32_t>(channel_hash_64);

    if (!atm_vex_cache.find(channel_writer, channel_hash_32)) 
    {
        VEX_Samples vex_samples;
        BucketQueue queue;
        SampleBucket bucket; 
        bucket.copyInfo(res, samples);
        queue.push(bucket);
        vex_samples.insert(bucket_queue_writer, std::pair<int, BucketQueue>(thread_id, queue));
        atm_vex_cache.insert(channel_writer, std::pair<int32_t, VEX_Samples>(channel_hash_32, vex_samples));
        return channel_hash_32;

    } else {
        VEX_Samples & channel = channel_writer->second;
        channel_writer.release();
        if (!channel.find(bucket_queue_writer, thread_id)) {
            BucketQueue queue;
            SampleBucket bucket; 
            bucket.copyInfo(res, samples);   
            queue.push(bucket);
            channel.insert(bucket_queue_writer, std::pair<int, BucketQueue>(thread_id, queue));  
        }
        return channel_hash_32; 
    }

    return -1;
}


int insert_vex_sample(const int32_t & handle, const int & thread_id, const UT_StackBuffer<float> & sample)
{
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    AutomatteVexCache::const_accessor channel;
    atm_vex_cache.find(channel, handle);
    VEX_Samples::accessor queue;
    channel->second.find(queue, thread_id);
    SampleBucket & bucket = queue->second.front();
    const Sample vex{sample[0], sample[1], sample[2], 
        sample[3], sample[4], sample[5]};
    bucket.push_back(vex);
    return bucket.size();
    // return 1;
}


AutomatteVexCache * get_AutomatteVexCache()
{ 
    return &atm_vex_cache; 
}

AutomatteImage * get_AutomatteImage()
{
    return &atm_image;
}

ImageInfo * get_ImageInfo()
{
    return &atm_image_info; 
}

} // end of HA_HDK