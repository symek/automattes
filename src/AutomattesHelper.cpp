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

// used only for creating per thread storage
static std::mutex automattes_mutex;
// our main storage
static VEX_Samples vexsamples;

int VEX_Samples_create(const int& thread_id)
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    VEX_Samples::const_iterator it = vexsamples.find(thread_id);
    if(it == vexsamples.end()) {
    	SampleBucket bucket;	
    	vexsamples.insert(std::pair<int, SampleBucket>(thread_id, bucket));
    }
    // std::cout << "Helper: " << &vexsamples << std::endl;
    return thread_id;
} 

int VEX_Samples_insert(const int& thread_id, const OpacitySample& sample)
{
	VEX_Samples::const_iterator it = vexsamples.find(thread_id);
	UT_ASSERT(it != vexsamples.end());
	vexsamples[thread_id].push_back(sample);
	return vexsamples[thread_id].size(); //thread_id;
}

VEX_Samples * VEX_Samples_get() {
	return &vexsamples;
}



} // end of HA_HDK