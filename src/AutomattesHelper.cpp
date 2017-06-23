#include <iostream>
#include <map>
#include <mutex>
#include <cstring>
#include <functional>
#include <memory>


#include <UT/UT_DSOVersion.h>
#include <UT/UT_Thread.h>
#include <UT/UT_PointGrid.h>

#include "AutomattesHelper.hpp"


namespace HA_HDK {

// used only for creating per thread storage
static std::mutex automattes_mutex;
// our main storage
static VEX_Samples vexsamples;
static BucketCounter vexBucketCounter;
static BucketCounter vrayBucketCounter;
static VEX_SampleClass vexsamplesC;

int VEX_Samples_create(const int& thread_id)
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    VEX_Samples::const_iterator it = vexsamples.find(thread_id);
    if(it == vexsamples.end()) {
        BucketQueue queue;
    	SampleBucket bucket;	
    	vexsamples.insert(std::pair<int, BucketQueue>(thread_id, queue));
        vexsamples[thread_id].push_back(bucket);
        //debug
        vexBucketCounter.insert(std::pair<int, int>(thread_id, 0));
        vrayBucketCounter.insert(std::pair<int, int>(thread_id, 0));
    }
    return thread_id;
} 

int VEX_Samples_insert(const int& thread_id, const Sample& sample)
{
	VEX_Samples::const_iterator it = vexsamples.find(thread_id);
	UT_ASSERT(it != vexsamples.end());
    BucketQueue::iterator jt = vexsamples[thread_id].begin();
	jt->push_back(sample);
    const size_t size = jt->size();
    if (size == 1) {
        vexBucketCounter[thread_id] += 1;
        std::cout << "VEX    thread: " << thread_id << ", bucket count:" \
            << vexBucketCounter[thread_id] << "\n";  
    }
	return static_cast<int>(size); //thread_id;
}

VEX_Samples * VEX_Samples_get() 
{
	return &vexsamples;
}

int VEX_Samples_increamentBucketCounter(const int& thread_id)
{
    vrayBucketCounter[thread_id] += 1;
    return vrayBucketCounter[thread_id];
}

#if 0
int 
VEX_SampleClass::create_channel(const std::string channel)
{
	VEX_Channels::const_iterator it = myChannels.find(channel);
	if (it == myChannels.end()) {
		const int index = myChannels.size();
		VEX_Samples samples;
		myChannels.insert(std::pair<std::string, VEX_Samples>(channel, samples));
		return index;
	} else {
		return -1;
	}
}

int
VEX_SampleClass::create_bucket(const int &thread_id, const std::string &channel)
{
	VEX_Channels::const_iterator chit = myChannels.find(channel);
	UT_ASSERT(chit != myChannels.end());
	VEX_Samples samples = chit->second;
	VEX_Samples::const_iterator sit = samples.find(thread_id);
	if (sit == samples.end()) {
		SampleBucket bucket;
		samples.insert(std::pair<int, SampleBucket>(thread_id, bucket));
	} 

	return thread_id;
}

int
VEX_SampleClass::insert_sample(const int &thread_id, const std::string& channel, const Sample& sample)
{
	VEX_Channels::const_iterator chit = myChannels.find(channel);
	VEX_Samples samples = chit->second;
	VEX_Samples::const_iterator sit = samples.find(thread_id);

}
#endif

} // end of HA_HDK