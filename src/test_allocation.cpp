#include <vector>
#include <queue>
#include <iostream>
#include <thread>
#include <atomic>
#include <memory>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "../../microbench/microbench.h"
#include "bench_utils.h"



#define VECTOR_ALLOCATOR TBB_SCALABLE  // STD, TBB_SCALABLE, TBB_CACHE_ALLIGNED 
#define TBB_SCALABLE_ALLOCATOR
#define TBB_PREVIEW_MEMORY_POOL 1


#include "tbb/memory_pool.h"
#include "tbb/scalable_allocator.h"
#include "tbb/cache_aligned_allocator.h"


constexpr size_t SAMPLE_DEEP = 6;//should be at least 6
constexpr size_t SAMPLE_SIZE = 6;
constexpr size_t nsamples    = 10*10*SAMPLE_DEEP;
constexpr size_t BUCKET_WIDTH= 64;
constexpr size_t BUCKET_SIZE = BUCKET_WIDTH*BUCKET_WIDTH*nsamples;
constexpr size_t BUCKET_MIN  = BUCKET_SIZE; // the ratio of SIZE we start randomly to grow buckets
constexpr size_t _OVERFLOW   = 1.f; // How much to grow buckets from starting point (see above)
constexpr size_t NTHREADS    = 8;
constexpr size_t RESX        = 1920;
constexpr size_t RESY        = 1080;
constexpr size_t NBUCKETS    = (RESX/BUCKET_WIDTH) * (RESY/BUCKET_WIDTH) / NTHREADS; // avarage number of buckets per thread
constexpr size_t HISTORY     = 4; // how many buckets we store per thread before prunning them out.




/////////////////////////1/////////////////////////////

struct Sample { float data[SAMPLE_SIZE]{0,1,2,3,4,5}; };
typedef std::vector<Sample> StdBucket;
void allocate_std_sample_vector()
{
    std::queue<StdBucket> bucketQueue;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        StdBucket bucket;
        bucket.reserve(BUCKET_MIN);
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            Sample sample;
            bucket.emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 
        bucketQueue.push(std::move(bucket));
    }
}

//////////////////////////2///////////////////////////////

typedef std::vector<float, tbb::scalable_allocator<float>> TBBAllocSample;
void allocate_tbb_sample_vector()
{
    //std::vector<TBBAllocSample> bucket;
    std::queue<std::unique_ptr<std::vector<TBBAllocSample>>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<std::vector<TBBAllocSample>> 
            bucket(new std::vector<TBBAllocSample>());
        bucket->reserve(BUCKET_MIN);
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            TBBAllocSample sample{0,1,2,3,4,5};
            bucket->emplace_back(sample);
        }
        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////3//////////////////////////////////////
typedef  std::vector<Sample, tbb::scalable_allocator<Sample>> 
BucketVectorTBBAllocSampleVector;
void allocate_tbb_bucket_vector()
{
    // std::vector<>
    std::queue<std::unique_ptr<BucketVectorTBBAllocSampleVector>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketVectorTBBAllocSampleVector> 
            bucket(new BucketVectorTBBAllocSampleVector());
        bucket->reserve(BUCKET_MIN);
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            Sample sample;
            bucket->emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}

////////////////////////////4/////////////////////////////////////

typedef tbb::memory_pool_allocator<float> sample_vector_allocator_t;
typedef std::vector<float, sample_vector_allocator_t> MemPoolSample;
typedef std::vector<MemPoolSample> BucketVectorTBBPoolSampleVector;
static  tbb::memory_pool<tbb::scalable_allocator<MemPoolSample>> 
sample_vector_memory_pool;

void allocate_tbb_sample_vector_pool()
{
    std::queue<std::unique_ptr<BucketVectorTBBPoolSampleVector>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN+ ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketVectorTBBPoolSampleVector>  bucket(new BucketVectorTBBPoolSampleVector());
        bucket->reserve(BUCKET_MIN);
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            MemPoolSample sample(6, 0, sample_vector_allocator_t(sample_vector_memory_pool));
            bucket->emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}


////////////////////////////5/////////////////////////////////////

struct Bucket {
public:
    Bucket() { 
        m_samples = new float[m_size]();
    }
    ~Bucket() { deallocate(); }
    void emplace_back(const float * sample) {
        if (pointer < m_size) {
            std::memcpy(&m_samples[pointer], sample, sizeof(Sample));
            pointer+=SAMPLE_SIZE;

        } else {
            const size_t old_size = m_size;
            m_size = old_size + BUCKET_MIN*SAMPLE_SIZE;
            float * tmp = new float[m_size]();
            std::memcpy(tmp, m_samples, sizeof(float)*old_size);
            delete[] m_samples;
            m_samples = tmp;
            emplace_back(sample); } } 
    void clear()      { pointer = 0; }
    void deallocate() { 
        clear(); 
        m_size = 0; 
        if (m_samples != nullptr)
            delete[] m_samples;
        m_samples = nullptr; 
    }
private:
    float  * m_samples = nullptr;
    size_t pointer = 0;
    size_t m_size  =  BUCKET_MIN*SAMPLE_SIZE;
};

void allocate_struct_bucket_pool_CppNewDelete()
{
    
    std::queue<std::unique_ptr<Bucket>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<Bucket> bucket(new Bucket());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket->emplace_back((float*)sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 
        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////6//////////////////////////////////
class TBB_Bucket;
static tbb::memory_pool<tbb::scalable_allocator<float>> 
tbb_bucket_mem_pool;

struct TBB_Bucket {
public:
    TBB_Bucket() { 
        m_size    = BUCKET_MIN*SAMPLE_SIZE;
        m_samples = (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*m_size);
    }
    ~TBB_Bucket() { deallocate(); }
    inline void emplace_back(const float * sample) noexcept {
        if (pointer < m_size) {
            std::memcpy(&m_samples[pointer], sample, sizeof(float)*SAMPLE_SIZE);
            pointer+=SAMPLE_SIZE;
        } else {
            const float old_size = m_size;
            m_size = old_size + BUCKET_MIN*SAMPLE_SIZE;
            float * tmp = (float*) tbb_bucket_mem_pool.malloc(2*sizeof(float)*m_size);
            std::memcpy(tmp, m_samples, sizeof(float)*old_size);
            tbb_bucket_mem_pool.free(m_samples);
            m_samples = tmp;
            emplace_back(sample); } }

    void clear()      { pointer = 0; }
    void deallocate() { 
        clear(); 
        if (m_samples != nullptr)
            tbb_bucket_mem_pool.free(m_samples);
        m_size = 0;
        m_samples = nullptr;
    }

private:
    float  * m_samples;
    size_t pointer = 0;
    size_t m_size  = 0;
};

void allocate_struct_bucket_pool_TbbMallocFree()
{
    std::queue<std::unique_ptr<TBB_Bucket>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<TBB_Bucket> bucket(new TBB_Bucket());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket->emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}


///////////////////////////////7//////////////////////////////////

struct TBB_Bucket_Sliced {
public:
    TBB_Bucket_Sliced() { 
        float * slice = (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
        m_samples.push_back(slice);
        m_capacity = m_slice_size;
    }
    ~TBB_Bucket_Sliced() { deallocate(); }

    inline void emplace_back(const float * sample) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            const uint _item = m_current_item % m_slice_size;
            float * buffer = m_samples.at(slice_pointer);
            std::memcpy(&buffer[_item], sample, sizeof(float)*SAMPLE_SIZE);
        } else {
            float * _slice = (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
            std::memcpy(_slice, sample, sizeof(float)*SAMPLE_SIZE);
            m_samples.push_back(_slice);
            m_capacity *= 2;
        }

        m_current_item += SAMPLE_SIZE;
    }

    void clear()      { m_current_item = 0; }
    void deallocate() { 
        clear(); m_capacity = 0;
        std::vector<float*>::iterator it = m_samples.begin();
        for (; it!=m_samples.end(); ++it)  {
            tbb_bucket_mem_pool.free(*it); 
        }
        m_samples.clear();
    }
private:
    std::vector<float*> m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = BUCKET_MIN*SAMPLE_SIZE;
    size_t m_capacity     = 0;
};


void allocate_struct_bucket_pool_TbbMallocFree_Slices()
{

    std::queue<std::unique_ptr<TBB_Bucket_Sliced>> bucketQueue;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<TBB_Bucket_Sliced> bucket(new TBB_Bucket_Sliced());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket->emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////8//////////////////////////////////
typedef std::vector<float, sample_vector_allocator_t> MemPoolSample;
typedef tbb::memory_pool_allocator<float> sample_vector_allocator_t;

struct TBBA_Bucket_Sliced_Std {
public:
    #if   VECTOR_ALLOCATOR == STD
    typedef std::vector<Sample>
    #elif VECTOR_ALLOCATOR == TBB_SCALABLE
    typedef std::vector<Sample, tbb::scalable_allocator<Sample>> 
    #elif VECTOR_ALLOCATOR == TBB_CACHE_ALLIGNED
    typedef std::vector<Sample, tbb::cache_aligned_allocator<Sample>>
    #elif VECTOR_ALLOCATOR == TBB_MEMORY_POOL
    typedef tbb::memory_pool_allocator<float> sample_vector_allocator_t;
    typedef std::vector<float, sample_vector_allocator_t> MemPoolSample;
    #endif
    Std_Slice;
    
    TBBA_Bucket_Sliced_Std() { 
        Std_Slice slice;
        slice.reserve(m_slice_size);
        m_samples.emplace_back(std::move(slice));
        m_capacity = m_slice_size;
    }

    inline void emplace_back(const Sample & sample) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            m_samples[slice_pointer].emplace_back(sample);
        } else {
            Std_Slice slice;
            slice.reserve(m_slice_size);
            slice.emplace_back(sample);
            m_samples.emplace_back(std::move(slice));
            m_capacity *= 2;
        }
        m_current_item++;
    }

    void clear()      { m_current_item = 0; m_samples.clear(); }
    void deallocate() { }
private:
    std::vector<Std_Slice> m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = BUCKET_MIN;
    size_t m_capacity     = 0;
};


void allocate_struct_bucket_pool_Std_TbbAlloc_Slices()
{

    std::queue<std::unique_ptr<TBBA_Bucket_Sliced_Std>> bucketQueue;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<TBBA_Bucket_Sliced_Std> bucket(new TBBA_Bucket_Sliced_Std());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const Sample sample;
            bucket->emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}

/////////////////////////////////////////////////////////////////
//                      END OF BENCHMARK FUNCS                 //
/////////////////////////////////////////////////////////////////


template<void (*BenchFunc)()>
void benchmark_runner()
{
    std::vector<std::thread> threads;
    int memory = -1;
    for(size_t t = 0; t < NTHREADS; ++t) {
            threads.push_back(std::thread(BenchFunc));
            memory = std::max(memory, getPhysicalMem());
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Max Physical Memory: %iMB, ", memory);
}

/////////////////////////////////////////////////////////////////

#define SINGLE_TEST
#define TESTS 8
#ifdef SINGLE_TEST
#define SIGN ==
#else
#define SIGN >=
#endif

using namespace moodycamel;

int main()
{
    const int repetition = 1;

     #if TESTS SIGN 1
    // std::vector<float[6]>
    printf("std::vector<Sample> // Sample -> Struct {float[6]}  : %.1fms (buckets in queue: %i ) \n", 
        microbench(&benchmark_runner<allocate_std_sample_vector>, 1, repetition), HISTORY);
    #endif

    #if  TESTS SIGN 2
    // std::vector<std::vector<float, tbb::scalable_alloc>>
    printf("std::vector<std::vector<float, tbb::scalable_alloc>>: %.1fms (buckets in queue: %i )\n", 
        microbench(&benchmark_runner<allocate_tbb_sample_vector>, 1, repetition), HISTORY);
    #endif

    #if TESTS SIGN 3
    // std::vector<Sample, std::vector<tbb::scalable_allocator>>
    printf("std::vector<Sample,std::vector<tbb::scalable_alloc>>: %.1fms (buckets in queue: %i ) \n", 
        microbench(&benchmark_runner<allocate_tbb_bucket_vector>, 1, repetition), HISTORY);
    #endif

    #if TESTS SIGN 4
    // std::vector<Sample, bucket_allocation >
    printf("std::vector<std::vector<float,tbb:memory_pool_alloc>: %.1fms (buckets in queue: %i ) \n", 
        microbench(&benchmark_runner<allocate_tbb_sample_vector_pool>, 1, repetition), HISTORY);
    #endif

    #if TESTS SIGN 5
     // Bucket <Sample, new/delete >
    printf("Struct Bucket[float[6]] -> new[all samples]/delete> : %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_CppNewDelete>, 1, repetition), HISTORY);
    #endif 

    #if TESTS SIGN 6
     // Bucket <Sample, tbb:mempool >
    printf("Struct TBB_Bucket[float[6]], tbb:malloc/tbb::free > : %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_TbbMallocFree>, 1, repetition), HISTORY);
    #endif

    #if TESTS SIGN 7
    // Bucket <Sample, tbb:mempool >
    printf("Struct TBB_Bucket_Sliced[float[6]], tbb:malloc...>  : %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_TbbMallocFree_Slices>, 1, repetition), HISTORY);
    #endif

     #if TESTS SIGN 8
    // Bucket <Sample, tbb:mempool >
    printf("Struct TBB_Bucket_VecTbbA[vec<Sample>,tbb::scal..]> : %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_Std_TbbAlloc_Slices>, 1, repetition), HISTORY);
    #endif

    return 0;

}