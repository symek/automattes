/*
 Partaly taken from HDK examples:
 http://www.sidefx.com/docs/hdk15.5/_v_e_x_2_v_e_x__example_8_c-example.html

 Contains snippets from PsyOp
 [1] - Jonah Friedman, Andrew C. Jones, Fully automatic ID mattes with support for motion blur and transparency.
 */
#include <functional>
#include <memory>
#include <map>
#include <tbb/concurrent_vector.h> 

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <VEX/VEX_VexOp.h>
#include <UT/UT_PointGrid.h>
#include <UT/UT_BoundingBox.h>

#include <cstring>

#include "MurmurHash3.h"
#include "AutomattesHelper.hpp"
 
#ifdef CONCURRENT_HASH_MAP
#include <tbb/concurrent_hash_map.h>
#endif

namespace HA_HDK {

static const uint32_t background = 2287214504; // precomputed from  MurmurHash3_x86_32("_ray_fog_object_internal_xyzzy", ...);

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
murmurhash3F(int argc,  void *argv[], void *data)
{
    float *result     = static_cast<float*>(argv[0]);
    const char * name = static_cast<const char*>(argv[1]);
    uint32_t m3hash = 0;
    MurmurHash3_x86_32(name, std::strlen(name), 0, &m3hash);
    *result = (m3hash != background) ? hash_to_float(m3hash) : 0.f;
}

static void
murmurhash3I(int argc,  void *argv[], void *data)
{
    uint32_t *result  = static_cast<uint32_t*>(argv[0]);
    const char * name = static_cast<const char*>(argv[1]);
    uint32_t m3hash = 0;
    MurmurHash3_x86_32(name, std::strlen(name), 0, &m3hash);
    *result = (m3hash != background) ? m3hash : 0;
}

static void vex_store_open(int argc, void *argv[], void *data)
{
    int            *result    = (int*)            argv[0];
    const char     *channel   = (const char*)     argv[1];

    // const int thread_id = UT_Thread::getMyThreadId();
    const int thread_id = SYSgetSTID();
    result[0] = VEX_Samples_create(thread_id);

}

static void vex_store_save(int argc, void *argv[], void *data)
{
          VEXint   *result = (      VEXint* )  argv[0];
    const VEXint   *handle = (const VEXint* )  argv[1];
    const VEXvec3  *P      = (const VEXvec3*)  argv[2];
    const VEXfloat *id     = (const VEXfloat*) argv[3];
    const VEXfloat *Af     = (const VEXfloat*) argv[4];

    // const int thread_id = UT_Thread::getMyThreadId();
    const int thread_id = SYSgetSTID();
    const float thf = static_cast<float>(thread_id);
    Sample sample = {P->x(), P->y(), P->z(), *id, *Af, thf};
    *result = VEX_Samples_insert(*handle, sample);
}

}// end of HA_HDK namespace



//
// Installation function
//
using namespace HA_HDK;
void
newVEXOp(void *)
{
    new VEX_VexOp("murmurhash3@&FS",  // Signature
        murmurhash3F,      // Evaluator
        VEX_ALL_CONTEXT,    // Context mask
        NULL,           // init function
        NULL,           // cleanup function
        VEX_OPTIMIZE_2, // Optimization level
        true);    

    new VEX_VexOp("murmurhash3@&IS",  // Signature
        murmurhash3I,      // Evaluator
        VEX_ALL_CONTEXT,    // Context mask
        NULL,           // init function
        NULL,           // cleanup function
        VEX_OPTIMIZE_2, // Optimization level
        true);

    new VEX_VexOp("vexstoreopen@&IS",  // Signature
        vex_store_open,      // Evaluator
        VEX_ALL_CONTEXT,    // Context mask
        NULL,           // init function
        NULL,           // cleanup function
        VEX_OPTIMIZE_2 // Optimization level
        );
     new VEX_VexOp("vexstoresave@&IIVFF",  // Signature
        vex_store_save,      // Evaluator
        VEX_ALL_CONTEXT,    // Context mask
        NULL,           // init function
        NULL,           // cleanup function
        VEX_OPTIMIZE_2, // Optimization level
        true);           
    
}
