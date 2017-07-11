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
static std::mutex automattes_mutex2;
// our main storage
static VEX_Samples vexsamples;
static BucketGrid bucketGrid;
static BucketCounter vexBucketCounter;
static BucketCounter vrayBucketCounter;
// static VEX_SampleClass vexsamplesC;
static BucketSize bucketSize = {0,0};
static bool bucketSizeSet = 0;
static ut_thread_id_t mainThreadId = 0;


int VEX_Samples_create(const int& thread_id)
{
    std::lock_guard<std::mutex> guard(automattes_mutex);
    const ut_thread_id_t currentMainThreadId = UT_Thread::getMainThreadId();

    if (currentMainThreadId != mainThreadId) {
        vexsamples.clear();
        mainThreadId = currentMainThreadId;
    }

    VEX_Samples::const_iterator it;
    it = vexsamples.find(thread_id);

    if(it != vexsamples.end())
        return thread_id;

    BucketQueue queue;
    SampleBucket bucket;    
    vexsamples.insert(std::pair<int, BucketQueue>(thread_id, queue));
    vexsamples[thread_id].push_back(bucket);
    //debug
    vexBucketCounter.insert(std::pair<int, int>(thread_id, 0));
    vrayBucketCounter.insert(std::pair<int, int>(thread_id, 0));

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

void VEX_Samples_insertBucket(const int & thread_id)
{
    std::lock_guard<std::mutex> guard(automattes_mutex2);
    BucketQueue  & bqueue = vexsamples.at(thread_id);
    BucketQueue::iterator kt = bqueue.begin();
    SampleBucket new_bucket;
    bqueue.insert(kt, new_bucket);
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

BucketSize * VEX_getBucketSize() {
	return &bucketSize;
}

void VEX_setBucketSize(int x, int y) {
	std::lock_guard<std::mutex> guard(automattes_mutex2);
	if (bucketSize[0] != 0 || bucketSize[1] != 0) 
		return;
	
	bucketSize[0] = x;
	bucketSize[1] = y;
}

int VEX_bucketSizeSet() { return bucketSizeSet; }


int SampleBucket::updateBoundingBox(const float & expx, const float & expy, const float & expz) 
{
    UT_Vector3 bucket_min = {FLT_MAX, FLT_MAX, FLT_MAX};
    UT_Vector3 bucket_max = {FLT_MIN, FLT_MIN, FLT_MIN};
    const size_t size = mySamples.size();
    for(int i=0; i < size; ++i) {
        const Sample sample = mySamples.at(i);
        const UT_Vector3 position = {sample[0], sample[1], 0.f};
        bucket_min = SYSmin(bucket_min, position);
        bucket_max = SYSmax(bucket_max, position);
    }
    myBbox.initBounds(bucket_min, bucket_max);
    myBbox.expandBounds(expx, expy, expz);
    registerBucket();
    myFinishedFlag = 1;
    return bucketGrid.size();
}

void SampleBucket::registerBucket() const 
{
    // update BucketGrid;
    std::lock_guard<std::mutex> guard(automattes_mutex2);
    const Xmax xmax = myBbox.maxvec().x();
    const Ymax ymax = myBbox.maxvec().y();
    BucketGrid::const_iterator it = bucketGrid.find(ymax);
    if (it == bucketGrid.end()) {
        BucketLine line;
        line.insert(std::pair<Xmax, SampleBucket*>(xmax, this));
        bucketGrid.insert(std::pair<Ymax, BucketLine>(ymax, line));
    } else {
        BucketLine line = it->second;
        BucketLine::const_iterator jt = line.find(xmax);
        if (jt == line.end()) {
            line.insert(std::pair<Xmax, SampleBucket*>(xmax, this));
        }
    }
}

// void SampleBucket::findBuckets(const UT_BoudingBox & bbox) const 
// {
// 	std::lock_guard<std::mutex> guard(automattes_mutex2);
	

// }

int VEX_getBucket(const int & thread_id, SampleBucket * bucket, int & offset)
{
    UT_ASSERT(vexsamples.find(thread_id) != vexsamples.end());
    const BucketQueue & queue = vexsamples.at(thread_id);
    BucketQueue::const_iterator jt = queue.begin();
    for (; jt != queue.end(); ++jt, ++offset) {
        if (jt->size() != 0) 
            bucket = &(*jt); 
            return 1;
    }

    return 0;
}

#if 0
int 
VEX_SampleClass::open_channel(const std::string & channel)
{
	VEX_Channels::const_iterator it = myChannels.find(channel);
	if (it == myChannels.end()) {
		const int index = myChannels.size();
		VEX_Samples samples;
		myChannels.insert(std::pair<std::string, VEX_Samples>(channel, samples));
		return index;
	} else {
		return -1; // TODO: should return an index of a channel;
	}
}

int 
VEX_SampleClass::insert_queue(const std::string & channel, const int & thread_id)
{
    // std::lock_guard<std::mutex> guard(automattes_mutex2);
    VEX_Channels::const_iterator chit = myChannels.find(channel);
    UT_ASSERT(chit != myChannels.end());

    if (chit == myChannels.end()) {
        const int result = open_channel(channel);
        UT_ASSERT(result != -1);
        if (result == -1)
            return -1;
    }

    VEX_Samples & vexsamples  = myChannels.at(channel);
    VEX_Samples::iterator sit = vexsamples.find(thread_id);

    if (sit == vexsamples.end()) {
        BucketQueue queue;
        vexsamples.insert(std::pair<int, BucketQueue>(thread_id, queue));
        return thread_id;
    } else {
        return 0;
    }

}



int
VEX_SampleClass::insert_bucket(const std::string & channel, const int & thread_id)
{
	VEX_Channels::const_iterator chit = myChannels.find(channel);
	UT_ASSERT(chit != myChannels.end());
	VEX_Samples vexsamples = chit->second;
	VEX_Samples::const_iterator sit = vexsamples.find(thread_id);
    UT_ASSERT(sit != vexsamples.end());
    BucketQueue  & bqueue = vexsamples.at(thread_id);
    BucketQueue::iterator kt = bqueue.begin();
    SampleBucket new_bucket;
    bqueue.insert(kt, new_bucket);
    return bqueue.size();

 }

// int
// VEX_SampleClass::insert_sample(const int &thread_id, const std::string& channel, const Sample& sample)
// {
// 	VEX_Channels::const_iterator chit = myChannels.find(channel);
// 	VEX_Samples samples = chit->second;
// 	VEX_Samples::const_iterator sit = samples.find(thread_id);

// }
#endif

} // end of HA_HDK