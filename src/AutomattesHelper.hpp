#pragma once

#ifndef __AutomattesHelper__
#define __AutomattesHelper__

#define CONCURRENT_HASH_MAP

namespace HA_HDK {

// our fixed sample: x,y,z,id,Af 
typedef std::vector<float> Sample; 
// vector of samples per thread (reused by many buckets)
typedef std::vector<Sample> SampleBucketV;



class SampleBucket
{
public:
    const size_t size() const noexcept { return mySamples.size(); }
    const size_t getNeighbourSize() const noexcept ;
    const Sample & at(const int & index) const;
    const UT_BoundingBox * getBBox() const noexcept { return &myBbox; }
    const SampleBucketV & getMySamples() const noexcept { return mySamples; }
    const int isRegistered() const noexcept { return myRegisteredFlag; } 
    void clearNeighbours() noexcept;
    void clear() noexcept;
    void push_back(const Sample & sample) { mySamples.push_back(sample); }
    void updateBoundingBox(const float &, const float &, const float &);
    size_t registerBucket();
    int  fillBucket(const UT_Vector3 &, const UT_Vector3 &, SampleBucket *);
    void findBucket(const float &, const float &, 
        const float &, const float &, SampleBucket *) const;
private:
    SampleBucketV mySamples;
    UT_BoundingBox myBbox;
    int myRegisteredFlag = 0;
    SampleBucketV myNeighbours;
    size_t myNeighbourSize = 0;
};

typedef std::vector<SampleBucket>BucketQueue;
#ifdef CONCURRENT_HASH_MAP
typedef tbb::concurrent_hash_map<int, BucketQueue> VEX_Samples;
#else
typedef std::map<int, BucketQueue> VEX_Samples;
#endif
typedef std::map<std::string, VEX_Samples> VEX_Channels;

// storage per channel.
typedef std::map<int, int> BucketCounter;
// main storage container.
typedef std::array<int, 2> BucketSize;

//
typedef float coord_t;
// typedef std::vector<SampleBucket*>    BucketVector;
typedef tbb::concurrent_vector<SampleBucket>    BucketVector;
// typedef std::map<coord_t, SampleBucket*> BucketGrid;


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
int VEX_getBucket(const int, SampleBucket *, int &);

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