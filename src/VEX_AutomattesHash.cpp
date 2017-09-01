/*
 Partaly taken from HDK examples:
 http://www.sidefx.com/docs/hdk15.5/_v_e_x_2_v_e_x__example_8_c-example.html

 Contains snippets from PsyOp
 [1] - Jonah Friedman, Andrew C. Jones, Fully automatic ID mattes with support for motion blur and transparency.
 */
#include <functional>
#include <memory>
#include <map>
#include <queue>
#include <mutex>
#include <atomic>
#include <tbb/concurrent_vector.h> 

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <VEX/VEX_VexOp.h>
#include <UT/UT_PointGrid.h>
#include <UT/UT_BoundingBox.h>

#include <cstring>

#include "MurmurHash3.h"
#include "AutomattesHelper.hpp"
#include <tbb/concurrent_hash_map.h>


namespace HA_HDK {
class ShaderStore {

    typedef uint32_t shader_id_t;
    typedef tbb::concurrent_hash_map
    <uint32_t, shader_id_t> ThreadAccountant;
    std::mutex register_lock;

public:
    ShaderStore(): current_cache(-1) {}
    shader_id_t * register_shader(uint32_t thread_id) {
        std::lock_guard<std::mutex> guard(register_lock);
        ThreadAccountant::accessor threads_writer;
        if (current_cache == -1) {
            m_threads.insert(threads_writer, std::pair<uint32_t, shader_id_t>(thread_id, 0));
            create_vex_storage2(threads_writer->second);
            current_cache = threads_writer->second;
            std::cout << "Creating new shader : " << current_cache << "\n";
            return &(threads_writer->second);
        } else {
            if (!m_threads.find(threads_writer, thread_id)) {
                m_threads.insert(threads_writer, 
                    std::pair<uint32_t, shader_id_t>(thread_id, 0));
            } else {
                threads_writer->second++;   
                if (threads_writer->second > current_cache) {
                    std::cout << "This is second shader? : " << current_cache << "\n";
                    create_vex_storage2(threads_writer->second);
                    current_cache = threads_writer->second;
                }
            }
                
            std::cout << "Returning current one: " << current_cache << "\n";
            return &(threads_writer->second);
        }
    }
private:
    std::atomic<shader_id_t> 
    current_cache;
    ThreadAccountant m_threads;
};

static const uint32_t background = 2287214504; // precomputed from  MurmurHash3_x86_32("_ray_fog_object_internal_xyzzy", ...);

static ShaderStore shader_store;

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


static void * automatte_open_init()
{
    const uint32_t thread_id     = SYSgetSTID();
    const uint32_t numProcessors = UT_Thread::getNumProcessors ();
    std::cout << "automatte_open_init: " << numProcessors << ": " << thread_id << "is main: " \
    << UT_Thread::isMainThread() << "\n";
    void * shader_id_ptr = static_cast<void*>(shader_store.register_shader(thread_id));
    return shader_id_ptr;
}


static void automatte_open(int argc, void *argv[], void *data)
{
    int            *result  = (int*)            argv[0];
    const char     *channel = (const char*)     argv[1];
    const VEXvec3  *res     = (const VEXvec3*)  argv[2];
    const VEXvec3  *samples = (const VEXvec3*)  argv[3];

    const uint32_t thread_id  = SYSgetSTID();
    const uint32_t shader_id  = *static_cast<uint32_t*>(data);
    const std::vector<int> image_res    {(int)res->x(), (int)res->y()};
    const std::vector<int> pix_samples  {(int)samples->x(), (int)samples->y()};
    const std::string      channel_name (channel);
    int allocated = allocate_vex_storage(shader_id, thread_id, image_res, pix_samples);
     // std::cout << "automatte_open shader_id: " << shader_id << " allocated: " << allocated << "\n";

    if (allocated) {
        result[0] = shader_id+1; // from index to bool 
    } else {
        result[0] = 0;
    }
}

static void automatte_write(int argc, void *argv[], void *data)
{
          VEXint   *result = (      VEXint* )  argv[0];
    const VEXint   *handle = (const VEXint* )  argv[1];
    const VEXvec3  *P      = (const VEXvec3*)  argv[2];
    const VEXfloat *id     = (const VEXfloat*) argv[3];
    const VEXfloat *Af     = (const VEXfloat*) argv[4];

    const uint32_t thread_id = SYSgetSTID();
    const float thf = static_cast<float>(thread_id);

    UT_StackBuffer<float> sample(6);
    sample[0] = P->x(); sample[1] = P->y(); sample[2] = P->z();
    sample[3] = *id;    sample[4] = *Af;    sample[5] = thf;
    // const Sample sample{P->x(), P->y(), P->z(), *id, *Af, thf};
    // bool to index...
    const shader_id_t shader = (*handle) - 1;
    int size = insert_vex_sample(shader, thread_id, sample);
    *result  = static_cast<uint32_t>(size);
}

static void automatte_close(void *data)
{
    close_vex_storage();
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
    new VEX_VexOp("automatte_open@&ISVV",  // Signature
        automatte_open,      // Evaluator
        VEX_SHADING_CONTEXT,    // Context mask
        automatte_open_init,// init function
        NULL,           // cleanup function
        VEX_OPTIMIZE_2, // Optimization level
        true); 
    new VEX_VexOp("automatte_write@&IIVFF",  // Signature
        automatte_write,      // Evaluator
        VEX_SHADING_CONTEXT,    // Context mask
        NULL,           // init function
        automatte_close,// cleanup function
        VEX_OPTIMIZE_2, // Optimization level
        true);  
}
