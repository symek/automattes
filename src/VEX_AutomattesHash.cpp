/*

 */

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <VEX/VEX_VexOp.h>

#include <cstring>

#include "MurmurHash3.h"

namespace HA_HDK {

// From Cryptomatte specification[1]
float hash_to_float(uint32_t hash)
{
    uint32_t mantissa = hash & (( 1 << 23) - 1);
    uint32_t exponent = (hash >> 23) & ((1 << 8) - 1);
    exponent = std::max(exponent, (uint32_t) 1);
    exponent = std::min(exponent, (uint32_t) 254);
    exponent = exponent << 23;
    uint32_t sign = (hash >> 31);
    sign = sign << 31;
    uint32_t float_bits = sign | exponent | mantissa;
    float f;
    std::memcpy(&f, &float_bits, 4);
    return f;
}

static void
murmurhash3(int argc,  void *argv[], void *data)
{
    float *result     = static_cast<float*>(argv[0]);
    const char * name = static_cast<const char*>(argv[1]);
    uint32_t m3hash = 0;
    MurmurHash3_x86_32(name, std::strlen(name), 0, &m3hash);
    *result = hash_to_float(m3hash);
}

}

//
// Installation function
//
using namespace HA_HDK;
void
newVEXOp(void *)
{
    new VEX_VexOp("murmurhash3@&FS",  // Signature
        murmurhash3,      // Evaluator
        VEX_ALL_CONTEXT,    // Context mask
        NULL,           // init function
        NULL,           // cleanup function
        VEX_OPTIMIZE_2, // Optimization level
        true);    
    
}
