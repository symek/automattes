#pragma once

#ifndef __AutomattesHelper__
#define __AutomattesHelper__

namespace HA_HDK {

// our fixed sample: x,y,z,id,Af 
typedef std::vector<float> Sample; 
// vector of samples per thread (reused by many buckets)
typedef std::vector<Sample> SampleBucketV;



class SampleBucket
{
public:
	size_t size() const { return mySamples.size(); }
	const Sample & at(const int & index) const { return mySamples.at(index);}
	void push_back(const Sample & sample) { mySamples.push_back(sample); }
	void updateBoundingBox(const float & expx, const float & expy, const float & expz) 
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

    	// bucket_min.z() = expz * -1.0f;
        // bucket_max.z() = expz;
    	myBbox.initBounds(bucket_min, bucket_max);
    	myBbox.expandBounds(expx, expy, expz);
	}
	const UT_BoundingBox * getBBox() const { return &myBbox; }
private:
	SampleBucketV mySamples;
	UT_BoundingBox myBbox;
};

typedef std::vector<SampleBucket>BucketQueue;
typedef std::map<int, BucketQueue> VEX_Samples;
typedef std::map<std::string, VEX_Samples> VEX_Channels;

// storage per channel.
typedef std::map<int, int> BucketCounter;
// main storage container.
typedef std::array<int, 2> BucketSize;


// pointgrid stuff useful for buliding filter side accesor.
typedef  UT_PointGridVector3ArrayAccessor<int, int> UT_Vector3Point;
typedef  UT_PointGrid<UT_Vector3Point>::queuetype UT_Vector3PointQueue;

// function exposed on vex side (temporarily instead of proper class)
int VEX_Samples_create(const int&);
int VEX_Samples_insert(const int&, const Sample&);
void VEX_Samples_insertBucket(const int&);
VEX_Samples * VEX_Samples_get();
int VEX_Samples_increamentBucketCounter(const int&);
BucketSize * VEX_getBucketSize();
void VEX_setBucketSize(int x, int y);
int VEX_bucketSizeSet();

} // end of HA_HDK Space

#endif

// // this keeps main container.
// struct VEX_SampleClass
// {
//     int open_channel(const std::string&);
//     int insert_queue(const std::string&, const int&);
//     int insert_bucket(const std::string&, const int&);
//     int insert_sample(const std::string&, const int&, const int&, const Sample&);
//     int get_bucket_offset(const std::string&, const int&) const;

// private:

//     VEX_Channels myChannels;
//     int resx;
//     int resy;

// };