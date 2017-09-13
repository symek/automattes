#include <vector>
#include <queue>
#include <iostream>
#include <thread>
#include <atomic>
#include <memory>
#include <algorithm>
#include <chrono>
#include <time.h>
#include <math.h>
#include <cstdio>
#include <cstring>
#include <x86intrin.h>
#include "../../microbench/microbench.h"
#include "bench_utils.h"



#define VECTOR_ALLOCATOR TBB_CACHE_ALLIGNED  // STD, TBB_SCALABLE, TBB_CACHE_ALLIGNED 
#define TBB_SCALABLE_ALLOCATOR
#define TBB_PREVIEW_MEMORY_POOL 1


#include "tbb/memory_pool.h"
#include "tbb/scalable_allocator.h"
#include "tbb/cache_aligned_allocator.h"


using namespace BENCHMARK;

/////////////////////////1/////////////////////////////

struct Sample { float data[SAMPLE_SIZE]{0,1,2,3}; };
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
        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        } 
        bucketQueue.push(std::move(bucket));
    }
}

//////////////////////////2///////////////////////////////

typedef std::vector<float, tbb::scalable_allocator<float>> TBBAllocSample;
void allocate_tbb_sample_vector()
{
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
        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        } 
        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////3//////////////////////////////////////
typedef  std::vector<Sample, tbb::scalable_allocator<Sample>> 
BucketStdVecOfSamples_TBBScalableAlloc;
void allocate_tbb_bucket_vector()
{
    std::queue<std::unique_ptr<BucketStdVecOfSamples_TBBScalableAlloc>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketStdVecOfSamples_TBBScalableAlloc> 
            bucket(new BucketStdVecOfSamples_TBBScalableAlloc());
        bucket->reserve(BUCKET_MIN);
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            Sample sample;
            bucket->emplace_back(sample);
        }
        if (bucketQueue.size() >= QUEUE_SIZE) {
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
        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        } 
        bucketQueue.push(std::move(bucket));
    }
}


////////////////////////////5/////////////////////////////////////

struct BucketRawArray {
public:
    BucketRawArray() { m_samples = new float[m_size](); }
    ~BucketRawArray() { deallocate(); }
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

void allocate_BuckeRawArray()
{
    std::queue<std::unique_ptr<BucketRawArray>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
  
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketRawArray> bucket(new BucketRawArray());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket->emplace_back((float*)sample);
        }

        if (bucketQueue.size() >= QUEUE_SIZE) {
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

        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}


///////////////////////////////7//////////////////////////////////

struct BucketSliced_TBBMalloc {
public:
    BucketSliced_TBBMalloc() { 
        float * slice = (float*)tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
        m_samples.emplace_back(slice);
        m_capacity = m_slice_size;
    }
    ~BucketSliced_TBBMalloc() { deallocate(); }

    inline void emplace_back(const float * sample) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            const uint _item = m_current_item % m_slice_size;
            float * slice = m_samples.at(slice_pointer);
            std::memcpy(&slice[_item], sample, sizeof(float)*SAMPLE_SIZE);
        } else {
            float * slice = (float*)tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
            std::memcpy(slice, sample, sizeof(float)*SAMPLE_SIZE);
            m_samples.emplace_back(slice);
            m_capacity *= 2;
        }
        m_current_item += SAMPLE_SIZE;
    }

    void clear()      { /*m_samples.clear()*/ ; m_current_item = 0; }
    void deallocate() { 
        clear(); 
        m_capacity = 0;
        std::vector<float*>::iterator it = m_samples.begin();
        for (; it!=m_samples.end(); ++it)  {
            tbb_bucket_mem_pool.free(*it); 
        }
        m_samples.clear();
    }
private:
    std::vector<float*>     m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = BUCKET_MIN*SAMPLE_SIZE;
    size_t m_capacity     = 0;
};


void allocate_bucketSliced_TBBMalloc()
{
    std::queue<std::unique_ptr<BucketSliced_TBBMalloc>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketSliced_TBBMalloc> bucket(new BucketSliced_TBBMalloc());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const float sample[6] = {0,1,2,3,4,5};
            bucket->emplace_back(sample);}

        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////8//////////////////////////////////

struct Bucket_Sliced_Std {
public:
    typedef std::vector<Sample> BucketSlice;
    Bucket_Sliced_Std() { 
            BucketSlice slice;
            slice.reserve(m_slice_size);
            m_samples.emplace_back(std::move(slice));
            m_capacity = m_slice_size;
    }

    inline void emplace_back(const Sample & sample) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            m_samples[slice_pointer].push_back(sample);
        } else {
            BucketSlice slice;
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
    std::vector<BucketSlice> m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = BUCKET_MIN;
    size_t m_capacity     = 0;
};

void allocate_bucket_Slices_Std()
{
    std::queue<std::unique_ptr<Bucket_Sliced_Std>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<Bucket_Sliced_Std> bucket(new Bucket_Sliced_Std());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const Sample sample;
            bucket->emplace_back(sample); }
        if (bucketQueue.size() >= QUEUE_SIZE)
            bucketQueue.pop();
        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////9//////////////////////////////////

struct Bucket_Sliced_TBBScalable : public Bucket_Sliced_Std {
public:
    typedef std::vector<Sample, tbb::scalable_allocator<Sample>>  BucketSlice;
};

void allocate_bucket_Std_Slices_TBBScalabe()
{
    std::queue<std::unique_ptr<Bucket_Sliced_TBBScalable>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<Bucket_Sliced_TBBScalable> 
        bucket(new Bucket_Sliced_TBBScalable());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const Sample sample;
            bucket->emplace_back(sample); }
        if (bucketQueue.size() >= QUEUE_SIZE)
            bucketQueue.pop();
        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////10//////////////////////////////////

struct Bucket_Sliced_TBBCacheAligned : public Bucket_Sliced_Std {
public:
    typedef std::vector<Sample, tbb::cache_aligned_allocator<Sample>> BucketSlice;
};

void allocate_bucket_Std_Slices_TBBCacheAligned()
{
    std::queue<std::unique_ptr<Bucket_Sliced_TBBCacheAligned>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<Bucket_Sliced_TBBCacheAligned> 
        bucket(new Bucket_Sliced_TBBCacheAligned());
        for(size_t p=0; p < rnd_BUCKET_SIZE; ++p) {
            const Sample sample;
            bucket->emplace_back(sample); }
        if (bucketQueue.size() >= QUEUE_SIZE)
            bucketQueue.pop();
        bucketQueue.push(std::move(bucket));
    }
}


///////////////////////////////11//////////////////////////////////


struct BucketSliced_TBBMalloc_SIMD {
public:
    BucketSliced_TBBMalloc_SIMD() { 
        float * slice = (float*)tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
        m_samples.emplace_back(slice);
        m_capacity = m_slice_size;
    }
    ~BucketSliced_TBBMalloc_SIMD() { deallocate(); }

    inline void emplace_back(const float * sample) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            const uint _item = m_current_item % m_slice_size;
            float * slice = m_samples.at(slice_pointer);
            std::memcpy(&slice[_item], sample, sizeof(float)*SAMPLE_SIZE);
        } else {
            float * slice = (float*)tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
            std::memcpy(slice, sample, sizeof(float)*SAMPLE_SIZE);
            m_samples.emplace_back(slice);
            m_capacity *= 2;
        }
        m_current_item += SAMPLE_SIZE;
    }

    inline void emplace_back_samples(const float * samples) noexcept {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            const uint _item = m_current_item % m_slice_size;
            const size_t i = SAMPLE_SIZE;
            float * slice = m_samples[slice_pointer];
            SIMD::simd_copy_stream((const float *)samples, /*4*samples*/ SAMPLE_SIZE*4, 4, &slice[_item]);
        } else {
            const size_t i = SAMPLE_SIZE;
            float * slice = (float*)tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size);
            SIMD::simd_copy_stream((const float *)samples, /*4*samples*/ SAMPLE_SIZE*4, 4, slice);
            m_samples.emplace_back(slice);
            m_capacity *= 2;
        }
        m_current_item += SAMPLE_SIZE*4;
    }

    void clear()      {m_current_item = 0; }
    void deallocate() { 
        clear(); 
        m_capacity = 0;
        std::vector<float*>::iterator it = m_samples.begin();
        for (; it!=m_samples.end(); ++it)  {
            tbb_bucket_mem_pool.free(*it); 
        }
        m_samples.clear();
    }
private:
    const size_t m_simd_size = 4;
    std::vector<float*>     m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = BUCKET_MIN*SAMPLE_SIZE;
    size_t m_capacity     = 0;
};


void allocate_bucketSliced_TBBMalloc_SIMD()
{
    std::queue<std::unique_ptr<BucketSliced_TBBMalloc_SIMD>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketSliced_TBBMalloc_SIMD> bucket(new BucketSliced_TBBMalloc_SIMD());
        for(size_t p=0; p < rnd_BUCKET_SIZE; p+=4) {
            const float sample[SAMPLE_SIZE*4] = \
            {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
            bucket->emplace_back_samples(sample);
        }

        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        } 

        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////12//////////////////////////////////


struct BucketSliced_TBBMalloc_SIMDAligned {
public:
    BucketSliced_TBBMalloc_SIMDAligned() { 
        void * ptr = tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size*4*SAMPLE_SIZE);
        if(ptr) {
            SIMD::align(alignof(double), sizeof(float), ptr, m_slice_size*4*SAMPLE_SIZE);
        }

        if(ptr) {
            float * slice = (float*)ptr;
            m_samples.emplace_back(slice);
            m_capacity = m_slice_size;
        }
        
    }
    ~BucketSliced_TBBMalloc_SIMDAligned() { deallocate(); }

    inline void emplace_back_samples(const float * samples) {
        if (m_current_item < m_capacity) {
            const uint slice_pointer = (const uint) \
            std::ceil(m_current_item / m_slice_size);
            const uint _item = (m_current_item % m_slice_size)*4*SAMPLE_SIZE;
            float * slice = m_samples[slice_pointer];
            SIMD::simd_copy_stream((const float *)samples, /*4*samples*/ 16, 4, &slice[_item]);
        } else {
            void * ptr = tbb_bucket_mem_pool.malloc(sizeof(float)*m_slice_size*4*SAMPLE_SIZE);
            if(ptr) {
                SIMD::align(alignof(double), sizeof(float), ptr, m_slice_size*4*SAMPLE_SIZE);
            }
            if (ptr) {
                float * slice = (float*)ptr;
                SIMD::simd_copy_stream((const float *)samples, /*4*samples*/ 16, 4, slice);
                m_samples.emplace_back(slice);
                m_capacity *= 2;
            }
        }
        m_current_item += 1;
    }

    void clear()      {m_current_item = 0; }
    void deallocate() { 
        clear(); 
        m_capacity = 0;
        std::vector<float*>::iterator it = m_samples.begin();
        for (; it!=m_samples.end(); ++it)  {
            tbb_bucket_mem_pool.free(*it); 
        }
        m_samples.clear();
    }
private:
    const size_t m_simd_size = 4;
    std::vector<float*>     m_samples;
    size_t m_current_item = 0;
    size_t m_slice_size   = 512; // number of sample packs (=4*4) in cache_level_1?
    size_t m_capacity     = 0;
};


void allocate_bucketSliced_TBBMalloc_SIMDAligned()
{
    static const float sample_array[16] __attribute__((aligned(64))) = \
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
   
    std::queue<std::unique_ptr<BucketSliced_TBBMalloc_SIMDAligned>> bucketQueue;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    uint rnd_BUCKET_SIZE = BUCKET_MIN + ceil(BUCKET_SIZE * r *_OVERFLOW);
    for(size_t i =0; i < NBUCKETS; ++i) {
        std::unique_ptr<BucketSliced_TBBMalloc_SIMDAligned> 
            bucket(new BucketSliced_TBBMalloc_SIMDAligned());
        for(size_t p=0; p < rnd_BUCKET_SIZE; p+=4) {
            bucket->emplace_back_samples(sample_array);
        }
        if (bucketQueue.size() >= QUEUE_SIZE) {
            bucketQueue.pop();
        }
        bucketQueue.push(std::move(bucket));
    }
}

///////////////////////////////13,14//////////////////////////////////

template < size_t simd, void (*CopyFunction) 
    (const float * source, 
     const size_t size, 
     const size_t stride, 
     float * dest) >

void copy_RawArray()
{
    const size_t stride = simd * SAMPLE_SIZE;
    const size_t tail   = BUCKET_SIZE % stride;
    size_t   array_size = 0;
    if (!tail) {
        array_size = BUCKET_SIZE;
    } else {
        array_size = BUCKET_SIZE + (stride - tail); 
    }
    float * bucket = new (_mm_malloc(sizeof(float)*array_size, alignof(float))) \
        float[array_size];

    for(size_t b = 0; b < NBUCKETS; ++b) {
        for(size_t p = 0; p < array_size; p+=stride) {
            float sample[stride]; // stride is knows at compile time...
            for (size_t s=0; s<stride; ++s) {
                sample[s] = s*1.f;//*xorshf96();// rand() is terrible slow at multithread...;
            } 
            CopyFunction((const float*) sample, stride, simd, (float*) &bucket[p]);
        }
    }
    delete[] bucket;
}



/////////////////////////////////////////////////////////////////
//                      END OF BENCHMARK FUNCS                 //
/////////////////////////////////////////////////////////////////


template<void (*BenchFunc)()>
void benchmark_runner()
{
    static int memory  =-1;
    static int times   = 0;

    std::vector<std::thread> threads;
    for(size_t t = 0; t < NTHREADS; ++t) {
        threads.push_back(std::thread(BenchFunc));
        const int inst_mem = getPhysicalMem();
        memory = std::max(memory, inst_mem); 
        }

    std::vector<std::thread>::iterator it;
    it = threads.begin();
    for(; it != threads.end(); ++it) {
        const int inst_mem = getPhysicalMem();
        memory = std::max(memory, inst_mem);
        it->join();
    }
    times++;
    if (times == REPETITION)
        printf("Max Memory: %iMB, ", memory);
}

/////////////////////////////////////////////////////////////////

#define SINGLE_TEST
#define TESTS 12    
#ifdef SINGLE_TEST
#define SIGN ==
#else
#define SIGN >=
#endif

using namespace moodycamel;

int main()
{

    //  #if TESTS SIGN 1
    // // std::vector<float[6]>
    // printf("std::vector<Sample> // Sample -> Struct {float[6]}  : %.1fms (buckets in queue: %i ) \n", 
    //     microbench(&benchmark_runner<allocate_std_sample_vector>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    // #if  TESTS SIGN 2
    // // std::vector<std::vector<float, tbb::scalable_alloc>>
    // printf("std::vector<std::vector<float, tbb::scalable_alloc>>: %.1fms (buckets in queue: %i )\n", 
    //     microbench(&benchmark_runner<allocate_tbb_sample_vector>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    // #if TESTS SIGN 3
    // // std::vector<Sample, std::vector<tbb::scalable_allocator>>
    // printf("std::vector<Sample,std::vector<tbb::scalable_alloc>>: %.1fms (buckets in queue: %i ) \n", 
    //     microbench(&benchmark_runner<allocate_tbb_bucket_vector>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    // #if TESTS SIGN 4
    // // std::vector<Sample, bucket_allocation >
    // printf("std::vector<std::vector<float,tbb:memory_pool_alloc>: %.1fms (buckets in queue: %i ) \n", 
    //     microbench(&benchmark_runner<allocate_tbb_sample_vector_pool>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    // #if TESTS SIGN 5
    //  // Bucket <Sample, new/delete >
    // printf("Struct BucketRawArray[float[6]] -> new/delete>      : %.1fms (buckets in queue: %i ) \n",
    //     microbench(&benchmark_runner<allocate_BuckeRawArray>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif 

    // #if TESTS SIGN 6
    // Bucket <Sample, tbb:mempool >
    // printf("Struct TBB_Bucket[float[6]], tbb:malloc/tbb::free > : %.1fms (buckets in queue: %i ) \n",
    //    microbench(&benchmark_runner<allocate_struct_bucket_pool_TbbMallocFree>, 1, REPETITION), QUEUE_SIZE);
    //std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    // #if TESTS SIGN 7
    // Bucket <Sample, tbb:mempool >
    // printf("Struct BucketSliced_TBBMalloc[float[6]],tbb:mallo>  : %.1fms (buckets in queue: %i ) \n",
    //     microbench(&benchmark_runner<allocate_bucketSliced_TBBMalloc>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    #if TESTS SIGN 8
    // Bucket <Sample, tbb:mempool >
    printf("Struct BucketSliced_StdVec[vec<Sample>,std::alloc.]>: %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_bucket_Slices_Std>, 1, REPETITION), QUEUE_SIZE);
    std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    #endif

    //  #if TESTS SIGN 9
    // // Bucket <Sample, tbb:mempool >
    // printf("Struct BucketSliced_StdVec[vec<Sample>,tbb::scal..]>: %.1fms (buckets in queue: %i ) \n",
    //     microbench(&benchmark_runner<allocate_bucket_Std_Slices_TBBScalabe>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    //  #if TESTS SIGN 10
    // // Bucket <Sample, tbb:mempool >
    // printf("Struct BucketSliced_StdVec[vec<Sample>,tbb::cache.]>: %.1fms (buckets in queue: %i ) \n",
    //     microbench(&benchmark_runner<allocate_bucket_Std_Slices_TBBCacheAligned>, 1, REPETITION), QUEUE_SIZE);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

     #if TESTS SIGN 11
    // Bucket <Sample, tbb:mempool >
    printf("Struct BucketSliced_TBBMalloc_SIMD[float*,tbb:mallo>: %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_bucketSliced_TBBMalloc_SIMD>, 1, REPETITION), QUEUE_SIZE);
    std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    #endif


   #if TESTS SIGN 12
    // Bucket <Sample, tbb:mempool >
    printf("Struct BucketSliced_TBBMalloc_SIMD[float*,tbb:align>: %.1fms (buckets in queue: %i ) \n",
        microbench(&benchmark_runner<allocate_bucketSliced_TBBMalloc_SIMDAligned>, 1, REPETITION), QUEUE_SIZE);
    std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    #endif


    // Copy tests:
    // #if TESTS SIGN 13
    // // std::memcpy
    // printf(" sample float[] -> bucket float[] std::memcpy>      : %.1fms(stride:(x*Sample): %i ) \n",
    //     microbench(&benchmark_runner<copy_RawArray<1, SIMD::std_copy>>, 1, REPETITION), 1);

    // printf(" sample float[] -> bucket float[] std::memcpy>      : %.1fms(stride:(x*Sample): %i ) \n",
    //     microbench(&benchmark_runner<copy_RawArray<4, SIMD::std_copy>>, 1, REPETITION), 4);
    // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif

    // #if TESTS SIGN 14
    // printf(" sample float[] -> bucket float[] _mm*_store_ps>    : %.1fms(stride:(x*Sample): %i ) \n",
    //     microbench(&benchmark_runner<copy_RawArray<4, SIMD::simd_copy>>, 1, REPETITION), 4);
    // // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // printf(" sample float[] -> bucket float[] _mm*_stream_ps>   : %.1fms(stride:(x*Sample): %i ) \n",
    //     microbench(&benchmark_runner<copy_RawArray<4, SIMD::simd_copy_stream>>, 1, REPETITION), 4);
    // // std::this_thread::sleep_for(std::chrono::milliseconds(COOL_ALLOC));
    // #endif


    return 0;

}
