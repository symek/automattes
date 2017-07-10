#pragma once

#ifndef __AutomattesHelper__
#define __AutomattesHelper__

namespace HA_HDK {
// our fixed sample: x,y,z,id,Af 
typedef std::vector<float> Sample; 
// vector of samples per thread (reused by many buckets)
typedef std::vector<Sample> SampleBucket;
// storage per channel.
typedef std::vector<SampleBucket>BucketQueue;
typedef std::map<int, BucketQueue> VEX_Samples;
typedef std::map<int, int> BucketCounter;
// main storage container.
typedef std::map<std::string, VEX_Samples> VEX_Channels;
typedef std::array<int, 2> BucketSize;

// tmp
static constexpr const int& size_at_least = 16*16*3*3;



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

// this keeps main container.
struct VEX_SampleClass
{
    int create_channel(const std::string);
    int create_bucket(const int&, const std::string&);
    int insert_sample(const int&, const std::string&, const Sample&);

private:

    VEX_Channels myChannels;
    int resx;
    int resy;

};

} // end of HA_HDK Space

#endif