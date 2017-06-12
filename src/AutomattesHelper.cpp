#include <iostream>
#include <map>
#include <mutex>
#include <cstring>
#include <functional>
#include <memory>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <VEX/VEX_VexOp.h>

#include "AutomattesHelper.hpp"


namespace HA_HDK {

static std::mutex automattes_mutex;
static std::map<int, OpacitySamples> AutoSampleStore;


const int create_shader(const char* filename, const VEXfloat uscale, 
                        const VEXfloat vscale, const VEXfloat intensity, 
                        const VEXint realuv)
{
    std::lock_guard<std::mutex> guard(automattes_mutex);


    std::map<int, AutoSampleStore>::const_iterator it;

    it = AutoSampleStore.find(hash);
    if (it == AutoSampleStore.end()) {
        OpacitySamples samples;
    	
        // shader.intersection = &tlIntersectionData(); //make_unique<tlIntersectionData>();
        // IrawanStore.insert(std::pair<int, IrawanInstance>(hash, shader));

    //     return hash;
    // } else {
    //     return it->first;
    }
}

} // end of HA_HDK