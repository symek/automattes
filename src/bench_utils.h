#pragma once
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#define DEBUG
#ifdef DEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif


constexpr size_t SAMPLE_DEEP = 3;//should be at least 6
constexpr size_t SAMPLE_SIZE = 4;
constexpr size_t nsamples    = 10*10*SAMPLE_DEEP;
constexpr size_t BUCKET_WIDTH= 64;
constexpr size_t BUCKET_SIZE = BUCKET_WIDTH*BUCKET_WIDTH*nsamples;
constexpr size_t BUCKET_MIN  = BUCKET_SIZE; // the ratio of SIZE we start randomly to grow buckets
constexpr size_t _OVERFLOW   = 2.f; // How much to grow buckets from starting point (see above)
constexpr size_t NTHREADS    = 32;
constexpr size_t RESX        = 2048;
constexpr size_t RESY        = 1556;
constexpr size_t NBUCKETS    = (RESX/BUCKET_WIDTH) * (RESY/BUCKET_WIDTH) / NTHREADS; // avarage number of buckets per thread
constexpr size_t QUEUE_SIZE  = 6; // how many buckets we store per thread before prunning them out.
constexpr size_t COOL_ALLOC  = 5 * 1000; // we need this to let allocator do deallocation.
constexpr size_t REPETITION  = 5;

namespace BENCHMARK
{

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVirtualMem(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result / 1024;
}

int getPhysicalMem(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result / 1024;
}


static unsigned long x_seed=123456789, y_seed=362436069, z_seed=521288629;

inline unsigned long xorshf96(void) 
{   //period 2^96-1
    unsigned long t;
    x_seed ^= x_seed << 16;
    x_seed ^= x_seed >> 5;
    x_seed ^= x_seed << 1;

   t = x_seed;
   x_seed = y_seed;
   y_seed = z_seed;
   z_seed = t ^ x_seed ^ y_seed;

  return z_seed;
}

namespace SIMD
{

    inline void *align( std::size_t alignment, std::size_t size,
                        void *&ptr, std::size_t &space ) {
        std::uintptr_t pn = reinterpret_cast< std::uintptr_t >( ptr );
        std::uintptr_t aligned = ( pn + alignment - 1 ) & - alignment;
        std::size_t padding = aligned - pn;
        if ( space < size + padding ) return nullptr;
        space -= padding;
        return ptr = reinterpret_cast< void * >( aligned );
    }

    inline void std_copy(const float * source, const size_t size,  const size_t stride, float * dest)
    {
        const size_t step = SAMPLE_SIZE*stride;
        std::memcpy(dest, source, sizeof(float)*step);
    }

    inline void simd_copy(const float * source, const size_t size, const size_t stride, float * dest) 
    {   
        // I assume stride == 4 for now:
        // 6 floats per sample * 4 samples = 24 floats = (3 * __m256) ||  (6*__m128) stores.
        #ifdef __AVX__ 
        __m256 * varray = (__m256 *) source;
        _mm256_store_ps(dest, *varray);
                 varray = (__m256 *) &source[8];
        _mm256_store_ps(&dest[8], *varray);
                 // varray = (__m256 *) &source[16];
        // _mm256_store_ps(&dest[16], *varray);
        #else
            #if defined (__SSE2__) || defined(__SSE3__) || defined(__SSE4__) 
                __m128 * varray = (__m128 *) source;
                _mm_store_ps(dest, *varray);
                         varray = (__m128 *) &source[4];
                _mm_store_ps(&dest[4], *varray);
                         varray = (__m128 *) &source[8];
                _mm_store_ps(&dest[8], *varray);
                         varray = (__m128 *) &source[12];
                _mm_store_ps(&dest[12], *varray);
                         // varray = (__m128 *) &source[16];
                // _mm_store_ps(&dest[16], *varray);
                         // varray = (__m128 *) &source[20];
                // _mm_store_ps(&dest[20], *varray);
            #else
                const size_t step = SAMPLE_SIZE*stride;
                std::memcpy(dest, source, sizeof(float)*step);
            #endif
        #endif
    }

    inline void simd_copy_stream(const float * source, const size_t size, const size_t stride, float * dest) 
    {
         #ifdef __AVX__ 
        __m256 * varray = (__m256 *) source;
        _mm256_stream_ps(dest, *varray);
                 varray = (__m256 *) &source[8];
        _mm256_stream_ps(&dest[8], *varray);
                 // varray = (__m256 *) &source[16];
        // _mm256_stream_ps(&dest[16], *varray);
        #else
            #if defined (__SSE2__) || defined(__SSE3__) || defined(__SSE4__) 
                __m128 * varray = (__m128 *) source;
                _mm_stream_ps(dest, *varray);
                         varray = (__m128 *) &source[4];
                _mm_stream_ps(&dest[4], *varray);
                         varray = (__m128 *) &source[8];
                _mm_stream_ps(&dest[8], *varray);
                         varray = (__m128 *) &source[12];
                _mm_stream_ps(&dest[12], *varray);
                         // varray = (__m128 *) &source[16];
                // _mm_stream_ps(&dest[16], *varray);
                         // varray = (__m128 *) &source[20];
                // _mm_stream_ps(&dest[20], *varray);
            #else
                const size_t step = SAMPLE_SIZE*stride;
                std::memcpy(dest, source, sizeof(float)*step);
            #endif
        #endif
    }
}

} // end of BENCHMARK namespace
