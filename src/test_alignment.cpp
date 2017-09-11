#include <iostream>
#include <type_traits>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>
 
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
bool address_as_int(void * ptr, T& result)
{
    char b[20];
    sprintf(b, "%p", ptr);
    std::string s(b);
    std::stringstream ss;
    ss << std::hex << s;
    ss >> result;

    return !ss.fail();
}
 
int main()
{
    static_vector<double, 10> p1;
    void * p2 = aligned_alloc(alignof(double), 10*sizeof(float));
    static float p3[10] __attribute__((aligned(alignof(double))));// = {1.0, 2.0, 3.0, 4.0};
    
    int x, y, z, w;

   address_as_int<int>((void*)&p1[0], x);
   address_as_int<int>(p2, y);
   address_as_int<int>((void*)&p3[1], z);

   std::cout\
   << "static_vector: " << (float)x / 32 << "\n"\
   << "aligned_alloc: " << (float)y / 64 << "\n"\
   << "__attrib__() : " << (float)z / 32 << "\n";
   delete[] (float*)p2;

   return 0;

}