#include <vector>
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


constexpr size_t sample_deep = 6;//should be at least 6
constexpr size_t sample_size = 6;
constexpr size_t nsamples    = 10*10*sample_deep;
constexpr size_t bucket_size = 64*64*nsamples;
constexpr size_t nthreads    = 32;
constexpr size_t nbuckets    = 1920/64;


struct Sample { float data[sample_size]{0,1,2,3,4,5}; };
typedef std::vector<float, tbb::scalable_allocator<float>> TBBAllocSample;

typedef tbb::memory_pool_allocator<float> sample_vector_allocator_t;
typedef std::vector<float, sample_vector_allocator_t> MemPoolSample;
typedef tbb::memory_pool_allocator<Sample> sample_bucket_allocator_t;



static tbb::memory_pool<tbb::scalable_allocator<MemPoolSample>> sample_memory_pool;
static tbb::memory_pool<tbb::scalable_allocator<std::vector<Sample>>> bucket_memory_pool;
class TBB_Bucket;
static tbb::memory_pool<tbb::scalable_allocator<float>> tbb_bucket_mem_pool;


struct Bucket {
public:
    Bucket() { m_samples = new float[bucket_size](); }
    Bucket(const size_t size) { 
        m_samples = new float[size]();
        m_size = size;
    } 
    void emplace_back(const float * sample) {
        std::memcpy(m_samples, sample, sizeof(Sample));
        pointer+=sample_size;
    }
    void clear() {
        pointer = 0;
        // delete m_samples;
        // m_samples = nullptr;
    }
private:
    float  * m_samples;
    size_t pointer = 0;
    size_t m_size  = 0;
};


struct TBB_Bucket {
public:
    TBB_Bucket() { 
        m_samples =  (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*bucket_size);
        m_size    = bucket_size;
    }

    TBB_Bucket(const size_t size) { 
        m_samples =  (float*) tbb_bucket_mem_pool.malloc(sizeof(float)*size);
        m_size    = size;
    } 

    inline void emplace_back(const float * sample) noexcept {
        // if (pointer < m_size) {
            std::memcpy(&m_samples[pointer], sample, sizeof(float)*sample_size);
            pointer = pointer+5;
        // } 
            // else {
            // float * tmp = tbb_bucket_mem_pool.malloc(sizeof(float)*bucket_size);
        // }
    }

    void clear()      { pointer = 0; }
    void deallocate() { clear(); m_size = 0; tbb_bucket_mem_pool.free(m_samples); }
private:
    float  * m_samples;
    size_t pointer = 0;
    size_t m_size  = 0;
};



// static tbb::memory_pool<tbb::scalable_allocator<std::vector<Sample>>> struct_bucket_memory_pool;

//////////////////////////////////////////////////////
void allocate_std_vector()
{
    std::vector<Sample> bucket;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            Sample sample;
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }
}
void run_std_vector_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_std_vector));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}
/////////////////////////////////////////////////////////

void allocate_tbb_vector()
{
    std::vector<TBBAllocSample> bucket;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            TBBAllocSample sample{0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }
}
void run_tbb_vector_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_tbb_vector));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}

/////////////////////////////////////////////////////////////////

void allocate_tbb_bucket()
{
    std::vector<Sample, tbb::scalable_allocator<Sample>> bucket;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            Sample sample;//{0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
    }
}
void run_tbb_bucket_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_tbb_bucket));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}

/////////////////////////////////////////////////////////////////

void allocate_tbb_pool()
{
    MemPoolSample * sample = nullptr; 
    std::vector<MemPoolSample*, tbb::scalable_allocator<MemPoolSample>> bucket;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            sample = new MemPoolSample(sample_vector_allocator_t(sample_memory_pool));
            bucket.emplace_back(sample);
        }
        // bucket.clear();
         std::vector<MemPoolSample*, tbb::scalable_allocator<MemPoolSample>>::iterator it;
         it = bucket.begin();
         // for (; it != bucket.end(); ++it)
            // sample_memory_pool.free(*it);
        bucket.clear();
        // sample_memory_pool.recycle();
    }
}
void run_tbb_pool_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_tbb_pool));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}

/////////////////////////////////////////////////////////////////
// This is busted.


void allocate_tbb_bucket_pool()
{
    std::vector<Sample, sample_bucket_allocator_t> * bucket = nullptr;
    Sample sample;
    bucket = new std::vector<Sample, sample_bucket_allocator_t>
    (bucket_size, sample_bucket_allocator_t(bucket_memory_pool));

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            Sample sample;
            bucket->emplace_back(sample);
        }
        bucket->clear();
        bucket_memory_pool.free(bucket);
    }
}
void run_tbb_bucket_pool_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_tbb_bucket_pool));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}

/////////////////////////////////////////////////////////////////



void allocate_Cppnew_struct_bucket_pool()
{
    Bucket bucket;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket.emplace_back((float*)sample);
        }
        bucket.clear();
        // bucket_memory_pool.free(bucket.m_samples);
    }
}
void run_Cppnew_struct_bucket_pool_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_Cppnew_struct_bucket_pool));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}

/////////////////////////////////////////////////////////////////



void allocate_tbbpool_struct_bucket_pool()
{
    TBB_Bucket bucket;

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_bucket_size = bucket_size/2 + ceil(bucket_size * r *.5);
  
    for(size_t i =0; i < nbuckets; ++i) {
        for(size_t p=0; p < rnd_bucket_size; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket.emplace_back(sample);
        }
        bucket.clear();
        // bucket_memory_pool.free(bucket.m_samples);
    }

    bucket.deallocate();
}
void run_tbbpool_struct_bucket_pool_bench()
{
    std::vector<std::thread> threads;
    for(size_t t = 0; t < nthreads; ++t) {
            threads.push_back(std::thread(allocate_tbbpool_struct_bucket_pool));
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    int memory = -1;
    for(; it != threads.end(); ++it) {
        memory = std::max(memory, getPhysicalMem());
        it->join();
    }
    printf("Physical mem: %iMB\n", memory);
}

/////////////////////////////////////////////////////////////////




using namespace moodycamel;

int main()
{
    const int repetition = 1;
     #if 0
    // std::vector<float[6]>
    printf("std::vector<Struct[float[6]] : %.1fms\n", 
        microbench(&run_std_vector_bench, 1, repetition));
    #endif

    #if 0
    // std::vector<float, std::vector<tbb::scalable_alloc>>
    printf("std::vector<std::vector<tbb:alloc>> : %.1fms\n", 
        microbench(&run_tbb_vector_bench, 1, repetition));
    #endif

    #if 0
    // std::vector<Sample, std::vector<tbb::scalable_alloc>>
    printf("std::vector<Sample, tbb::alloc> : %.1fms\n", 
        microbench(&run_tbb_bucket_bench, 1, repetition));
    #endif

    #if 0
    // std::vector<Sample, bucket_allocation >
    printf("std::vector<std::vector<float>, tbb:memory_pool_allocator> : %.1fms\n", 
        microbench(&run_tbb_bucket_pool_bench, 1, repetition));
    #endif

    #if 0
    // std::vector<MemPoolsample, sample_vector_allocator_t >
    printf("std::vector<MemPoolsample, sample_vector_allocator_t >>> : %.1fms\n", 
        microbench(&run_tbb_pool_bench, 1, repetition));
    #endif

    #if 1
     // Bucket <Sample, new/delete >
    printf("Bucket <Bucket[float[6]], new / delete > : %.1fms\n",
        microbench(&run_Cppnew_struct_bucket_pool_bench, 1, repetition));
    #endif 

    #if 1
     // Bucket <Sample, tbb:mempool >
    printf("Bucket <Bucket[float[6]], tbb:mempoo > : %.1fms\n",
        microbench(&run_tbbpool_struct_bucket_pool_bench, 1, repetition));
    #endif

    return 0;

}