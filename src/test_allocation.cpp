#include <vector>
#include <queue>
#include <iostream>
#include <thread>
#include <atomic>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "../../microbench/microbench.h"
#include "bench_utils.h"



// #define TBB_CACHE_ALIGNED_ALLOC
#define TBB_SCALABLE_ALLOCATOR
#define TBB_PREVIEW_MEMORY_POOL 1


#include "tbb/memory_pool.h"
#include "tbb/scalable_allocator.h"


constexpr size_t SAMPLE_DEEP = 6;//should be at least 6
constexpr size_t SAMPLE_SIZE = 6;
constexpr size_t nsamples    = 10*10*SAMPLE_DEEP;
constexpr size_t BUCKET_SIZE = 64*64*nsamples;
constexpr size_t NTHREADS    = 32;
constexpr size_t NBUCKETS    = (1920/64) * (1080/64) / NTHREADS; // avarage number of buckets per thread
constexpr size_t _OVERFLOW   = .7f; // .5 means num samples == bucket size. .6 means 10% more samples, so buckets will have to reallocate.
constexpr size_t HISTORY     = 4;




/////////////////////////1/////////////////////////////

struct Sample { float data[SAMPLE_SIZE]{0,1,2,3,4,5}; };
typedef std::vector<Sample> StdBucket;
void allocate_std_sample_vector()
{
    std::queue<StdBucket> bucketQueue;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        StdBucket bucket;
        bucket.reserve(BUCKET_SIZE);
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            Sample sample;
            bucket.emplace_back(sample);
        }

        if (bucketQueue.size() >= HISTORY) {
            bucketQueue.pop();
        } 
        bucketQueue.push(bucket);
    }
}

//////////////////////////2///////////////////////////////

typedef std::vector<float, tbb::scalable_allocator<float>> TBBAllocSample;
void allocate_tbb_sample_vector()
{
    std::vector<TBBAllocSample> bucket;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            TBBAllocSample sample{0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }
}

///////////////////////////3//////////////////////////////////////

void allocate_tbb_bucket_vector()
{
    // std::vector<>
    std::vector<Sample, tbb::scalable_allocator<Sample>> bucket;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            Sample sample; //{0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }
}

////////////////////////////4/////////////////////////////////////

typedef tbb::memory_pool_allocator<float> sample_vector_allocator_t;
typedef std::vector<float, sample_vector_allocator_t> MemPoolSample;
static  tbb::memory_pool<tbb::scalable_allocator<MemPoolSample>> 
sample_vector_memory_pool;

void allocate_tbb_sample_vector_pool()
{
    // MemPoolSample * sample = nullptr; 
    std::vector<MemPoolSample> bucket;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            MemPoolSample sample(6, 0, sample_vector_allocator_t(sample_vector_memory_pool));
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }
}


////////////////////////////5/////////////////////////////////////

struct Bucket {
public:
    Bucket() { m_samples = new float[BUCKET_SIZE*SAMPLE_SIZE](); }
    void emplace_back(const float * sample) {
        if (pointer < m_size) {
            std::memcpy(&m_samples[pointer], sample, sizeof(Sample));
            pointer+=SAMPLE_SIZE;
        } else {
            float * tmp = new float[2*BUCKET_SIZE*SAMPLE_SIZE]();
            std::memcpy(tmp, m_samples, sizeof(float)*m_size);
            delete[] m_samples;
            m_samples = tmp;
            m_size = BUCKET_SIZE*SAMPLE_SIZE*2;
            emplace_back(sample); } } 
    void clear()      { pointer = 0; }
    void deallocate() { clear(); m_size = 0; delete[] m_samples; }
private:
    float  * m_samples;
    size_t pointer = 0;
    size_t m_size  = 0;
};

void allocate_struct_bucket_pool_CppNewDelete()
{
    Bucket bucket;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket.emplace_back((float*)sample);
        }
        bucket.clear();
        // bucket_memory_pool.free(bucket.m_samples);
    }
    bucket.deallocate();
}

///////////////////////////////6//////////////////////////////////
class TBB_Bucket;
static tbb::memory_pool<tbb::scalable_allocator<float>> 
tbb_bucket_mem_pool;

struct TBB_Bucket {
public:
    TBB_Bucket() { 
        m_samples = (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*BUCKET_SIZE*SAMPLE_SIZE);
        m_size    = BUCKET_SIZE*SAMPLE_SIZE; }
    inline void emplace_back(const float * sample) noexcept {
        if (pointer < m_size) {
            std::memcpy(&m_samples[pointer], sample, sizeof(float)*SAMPLE_SIZE);
            pointer+=SAMPLE_SIZE;
        } else {
            float * tmp = (float*) tbb_bucket_mem_pool.malloc(2*sizeof(float)*BUCKET_SIZE*SAMPLE_SIZE);
            std::memcpy(tmp, m_samples, sizeof(float)*m_size);
            tbb_bucket_mem_pool.free(m_samples);
            m_samples = tmp;
            m_size = BUCKET_SIZE*SAMPLE_SIZE*2;
            emplace_back(sample); } }

    void clear()      { pointer = 0; }
    void deallocate() { clear(); m_size = 0; tbb_bucket_mem_pool.free(m_samples); }
private:
    float  * m_samples;
    size_t pointer = 0;
    size_t m_size  = 0;
};

void allocate_struct_bucket_pool_TbbMallocFree()
{
    TBB_Bucket bucket;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }

    bucket.deallocate();
}


///////////////////////////////7//////////////////////////////////

struct TBB_Bucket_Sliced {
public:
    TBB_Bucket_Sliced() { 
        float * slice = (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
        m_samples.push_back(slice);
        m_capacity = m_slice_size;
    }

    inline void emplace_back(const float * sample) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            float * buffer = m_samples.at(slice_pointer);
            std::memcpy(buffer, sample, sizeof(float)*SAMPLE_SIZE);
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
        for (; it!=m_samples.end(); ++it)  
            tbb_bucket_mem_pool.free(*it); 
    }
private:
    std::vector<float*> m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = BUCKET_SIZE*SAMPLE_SIZE;
    size_t m_capacity     = 0;
};

void allocate_struct_bucket_pool_TbbMallocFree_Slices()
{
    TBB_Bucket_Sliced bucket;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_SIZE/2 + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }

    bucket.deallocate();
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
#define TESTS 1
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
    printf("std::vector<Sample> // Sample -> Struct {float[6]}  : %.1fms\n", 
        microbench(&benchmark_runner<allocate_std_sample_vector>, 1, repetition));
    #endif

    #if TESTS SIGN 2
    // std::vector<std::vector<float, tbb::scalable_alloc>>
    printf("std::vector<std::vector<float, tbb::scalable_alloc>>: %.1fms\n", 
        microbench(&benchmark_runner<allocate_tbb_sample_vector>, 1, repetition));
    #endif

    #if TESTS SIGN 3
    // std::vector<Sample, std::vector<tbb::scalable_allocator>>
    printf("std::vector<Sample,std::vector<tbb::scalable_alloc>>: %.1fms\n", 
        microbench(&benchmark_runner<allocate_tbb_bucket_vector>, 1, repetition));
    #endif

    #if TESTS SIGN 4
    // std::vector<Sample, bucket_allocation >
    printf("std::vector<std::vector<float,tbb:memory_pool_alloc>: %.1fms\n", 
        microbench(&benchmark_runner<allocate_tbb_sample_vector_pool>, 1, repetition));
    #endif

    #if TESTS SIGN 5
     // Bucket <Sample, new/delete >
    printf("Struct Bucket[float[6]] -> new[all samples]/delete> : %.1fms\n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_CppNewDelete>, 1, repetition));
    #endif 

    #if TESTS SIGN 6
     // Bucket <Sample, tbb:mempool >
    printf("Struct TBB_Bucket[float[6]], tbb:malloc/tbb::free > : %.1fms\n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_TbbMallocFree>, 1, repetition));
    #endif

    #if TESTS SIGN 7
    // Bucket <Sample, tbb:mempool >
    printf("Struct TBB_Bucket_Sliced[float[6]], tbb:malloc...>  : %.1fms\n",
        microbench(&benchmark_runner<allocate_struct_bucket_pool_TbbMallocFree_Slices>, 1, repetition));
    #endif


    return 0;

}