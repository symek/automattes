#include <iostream>
#include <map>
#include <mutex>
#include <cstring>
#include <functional>
#include <memory>
#include <queue>
#include <atomic>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <UT/UT_PointGrid.h>
#include <IMG/IMG_DeepShadow.h>
#include "AutomattesHelper.hpp"



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
static AutomatteImage atm_image;
static ImageInfo      atm_image_info;

// Another experiment
#ifdef USE_DEEP_MAP
static IMG_DeepShadow dsm;
#endif

const Sample & SampleBucket::at(const int & index) const 
{
    UT_ASSERT (index < m_samples.size()) ;
    return m_samples.at(index);
}


void SampleBucket::clear() noexcept
{ 
    std::vector<Sample> tmp;
    m_samples.swap(tmp); 
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

    const size_t size = m_samples.size();
    for(int i=0; i < size; ++i) {
        const Sample vexsample = m_samples.at(i);
        const float pxf = vexsample[0] * atm_image_info.gridresx;
        const float pyf = vexsample[1] * atm_image_info.gridresy;
        const int subpxi = std::floor(pxf) + atm_image_info.image_margin;
        const int subpyi = std::floor(pyf) + atm_image_info.image_margin;
        const int index  = subpyi *atm_image_info.gridresx + subpxi;
        if(index < atm_image.size() && index > 0) {
            atm_image.at(index) = vexsample;
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
    #ifdef USE_DEEP_MAP
    writer.close();
    #endif
    myRegisteredFlag = 1;
    return atm_image.size();
}

inline void ImageInfo::set_image_size(const std::vector<int> & res, 
                                      const std::vector<int> & samples) noexcept
{
    gridresx = res[0] * samples[0] + 2 * image_margin * samples[0];
    gridresy = res[1] * samples[1] + 2 * image_margin * samples[1];
    image_size = gridresx * gridresy * max_samples;
    m_resolution = res;
    m_samples    = samples;
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

    std::lock_guard<std::mutex> guard(automattes_mutex);
    const ut_thread_id_t currentMainThreadId = UT_Thread::getMainThreadId();

    if (atm_image_info.image_size == 0) {
        atm_image_info.set_image_size(res, samples);
         #ifdef USE_DEEP_MAP
        dsm.setOption("compression", "5");
        dsm.setOption("zbias", "0.05");
        dsm.setOption("depth_planes", "zfront,zback");
        dsm.create("/tmp/automatte.rat", res[0], res[1], samples[0], samples[1]);
         #endif
    }

    if (atm_image.capacity() == 0)
    {
        atm_image.resize(atm_image_info.image_size);
        for (int i=0; i<atm_image_info.image_size; ++i) {
            const Sample sample{0,0,0,0,0,0};
            atm_image[i] = sample;
        }
    }

    if (currentMainThreadId != mainThreadId) {
        #ifdef USE_DEEP_MAP
        dsm.setOption("compression", "5");
        dsm.setOption("zbias", "0.05");
        dsm.setOption("depth_planes", "zfront,zback");
        dsm.create("/tmp/automatte.rat", res[0], res[1], samples[0], samples[1]);
        #endif
        atm_vex_cache.clear();
        AutomatteImage tmp;
        tmp.resize(atm_image_info.image_size);
        atm_image.swap(tmp);
        for (int i=0; i<atm_image_info.image_size; ++i) {
            const Sample sample{0,0,0,0,0,0};
            atm_image[i] = sample;
        }
        mainThreadId = currentMainThreadId;
    }

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


int insert_vex_sample(const int32_t & handle, const int & thread_id, const Sample & sample)
{
    // std::lock_guard<std::mutex> guard(automattes_mutex);
    AutomatteVexCache::const_accessor channel;
    atm_vex_cache.find(channel, handle);
    VEX_Samples::accessor queue;
    channel->second.find(queue, thread_id);
    SampleBucket & bucket = queue->second.front();
    bucket.push_back(sample);
    return bucket.size();
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