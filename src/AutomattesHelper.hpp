#pragma once

#define DEBUG
// #define USE_DEEP_MAP


#ifdef DEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif

namespace HA_HDK {

// our fixed sample: x,y,z,id,Af 
#ifdef TBB_CACHE_ALIGNED_ALLOC
typedef std::vector<float, tbb::scalable_allocator<float>> Sample;
#else 
typedef std::vector<float> Sample;   
#endif

#ifdef TBB_PREVIEW_MEMORY_POOL
typedef tbb::memory_pool_allocator<Sample> sample_vector_allocator_t;
typedef std::vector<Sample, sample_vector_allocator_t> SampleVector;
#else
// typedef std::allocator<Sample>  sample_vector_allocator_t;
typedef std::vector<Sample> SampleVector;
#endif

class ImageInfo
{
public:
    ImageInfo() : image_size(0), gridresx(0), gridresy(0), initalized(0) {}
    bool update_size(const std::vector<int> & res, 
                     const std::vector<int> & samples) noexcept;


    std::atomic<size_t> image_size;//?
    std::atomic<int>    gridresx;
    std::atomic<int>    gridresy;
    std::atomic<int>    initalized;

    const size_t        image_margin = 3;
    const size_t        max_samples  = 1;
private:
    std::vector<int>    m_resolution{0,0};
    std::vector<int>    m_samples{0,0};

};

class SampleBucket
{
public:
    SampleBucket();
    const Sample & at(const int & index) const;
    size_t size() const noexcept { return m_samples->size(); }
    int isRegistered() const noexcept { return myRegisteredFlag; } 
    void clear() noexcept;
    void push_back(const Sample & sample) { m_samples->emplace_back(sample); }
    void copyInfo(const SampleBucket *) noexcept;
    void copyInfo(const std::vector<int> &, const std::vector<int> &) noexcept;
    size_t registerBucket();
private:
    // std::vector<Sample> m_samples;
    SampleVector * m_samples;
    int myRegisteredFlag = 0;
public:
    uint m_resolution[2]   = {0,0};
    uint m_pixelsamples[2] = {0,0};
};

typedef std::queue<SampleBucket>                       BucketQueue;
typedef tbb::concurrent_hash_map<int, BucketQueue>     VEX_Samples;
typedef tbb::concurrent_hash_map<int32_t, VEX_Samples> AutomatteVexCache; 
typedef tbb::concurrent_vector<Sample>                 AutomatteImage;
typedef uint32_t shader_id_t;

// VEX access:

int allocate_vex_storage(const shader_id_t &, const int &, \
    const std::vector<int> &, const std::vector<int> &);

int create_vex_storage(const std::string &, const int &, 
                       const std::vector<int> &, 
                       const std::vector<int> &);

void create_vex_storage2(shader_id_t & shader_id);

int insert_vex_sample (const shader_id_t &, const uint32_t &, const UT_StackBuffer<float> &);
void close_vex_storage();

inline void compute_atm_image_size(const std::vector<int> & res, 
                                   const std::vector<int> & samples,
                                   const size_t & margin, size_t & gridresx,
                                   size_t & gridresy, size_t & size) noexcept;
// Pixel filter access:
AutomatteVexCache * get_AutomatteVexCache();
AutomatteImage    * get_AutomatteImage();
ImageInfo         * get_ImageInfo();
// ImageInfo

} // end of HA_HDK Space

