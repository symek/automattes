#include <iostream>
#include <type_traits>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>
#include <x86intrin.h>
#include <vector>
#include <memory>

#include "tbb/aligned_space.h"
#include "tbb/cache_aligned_allocator.h"
#include "tbb/scalable_allocator.h"

#define TBB_PREVIEW_MEMORY_POOL 1
#include "tbb/memory_pool.h"

namespace BENCHMARK  
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
} // end of BENCHMARK namespace

template<class T, std::size_t N>
class static_vector
{
    // properly aligned uninitialized storage for N T's
    typename std::aligned_storage<sizeof(T), alignof(T)>::type data[N];
    std::size_t m_size = 0;
 
public:
    // Create an object in aligned storage
    template<typename ...Args> void emplace_back(Args&&... args) 
    {
        if( m_size >= N ) // possible error handling
            throw std::bad_alloc{};
        new(data+m_size) T(std::forward<Args>(args)...);
        ++m_size;
    }
 
    // Access an object in aligned storage
    const T& operator[](std::size_t pos) const 
    {
        return *reinterpret_cast<const T*>(data+pos);
    }
 
    // Delete objects from aligned storage
    ~static_vector() 
    {
        for(std::size_t pos = 0; pos < m_size; ++pos) {
            reinterpret_cast<T*>(data+pos)->~T();
        }
    }
};
// little embarassingn but whatever...
template <typename T>
T intPtr(void * ptr)
{
    T result;
    char b[20];
    sprintf(b, "%p", ptr);
    std::string s(b);
    std::stringstream ss;
    ss << std::hex << s;
    ss >> result;

    return result; //!ss.fail();
}
 
static tbb::memory_pool<tbb::scalable_allocator<float>> my_pool;

int main()
{
    constexpr size_t size = 16;
    static_vector<float, size> p1;


    float source[16] __attribute__((aligned(64)));
    std::cout << "Source alignment: " << intPtr<long>((void*)source) % 64 << '\n';
    
    // float dest[size] __attribute__((aligned(64)));
    // std::vector<float, tbb::cache_aligned_allocator<float>> destv(size);
    // float *dest = (float*)&destv;
    // tbb::cache_aligned_allocator<float> allocator;
    // float *dest = allocator.allocate(8192);

    size_t maxsize = 8192;

    float * dest = nullptr;
    void  * ptr  = my_pool.malloc(8192);
    BENCHMARK::align(alignof(int), sizeof(float), ptr, maxsize); 
    if (ptr) {
        std::cout << "Dest   maxsize: " << maxsize <<  '\n';
        dest = (float*)ptr;
    }   

    const long dest_ptr = intPtr<long>((void*)dest);
    std::cout << "Dest   alignment: " << dest_ptr % 64 << '\n';


    for (int i = 0; i < size; ++i)
      source[i] = i;

   
    for(int i=0; i<8192; i+=8) {
        // const int idx = i*8;
        _mm256_stream_ps(&dest[i],   *(__m256*)source);
        i+= 8;
        _mm256_stream_ps(&dest[i], *(__m256*)&source[8]);
    }

    for (int i = 0; i < 512; ++i)
        std::cout << dest[i] << ", ";
    std::cout << "\n";


    my_pool.free((void*)ptr);

    //  long x, y, z1, z2, w;

    // x  = intPtr<long>((void*)&p1);
    // y  = intPtr<long>(p2);
    // z1 = intPtr<long>((void*)&p3[0]);
    // z2 = intPtr<long>((void*)&p3[1]);
    // w  = intPtr<long>((void*)p4);

    // std::cout\
    // << "static_vector: " << x % 32 << "\n"\
    // << "aligned_alloc: " << y % 64 << "\n"\
    // << "__attrib__() : " << (double)z1 / 32 << "\n"\
    // << "aligned_space: " << w << "\n";
    // // << "diff         : " << z2 - z1 << "\n"\
    // delete[] (float*)p2;

    return 0;
    // tbb::aligned_space<float, size> space;
    // float * p4 = space.begin();

}