#pragma once

#ifndef __AutomattesHelper__
#define __AutomattesHelper__

#define NO_MUTEX_IN_BUCKETVECTOR
#define TBB_VEX_STORE
#define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif

namespace HA_HDK {

// our fixed sample: x,y,z,id,Af 
typedef std::vector<float> Sample;
// temporal
// struct Sample 
// {
//     Sample(const float x, const float y, const float z, \
//         const float id, const float Af, const float thread) noexcept
//     {
//         storage[0] = x; storage[1] = y; storage[2] = z;
//         storage[3] = id; storage[4] = Af; storage[5] = thread;
//     }
//     float operator[](size_t i) const { return storage[i]; }
// private: 
//     float storage[6]; 
// };

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
    size_t copyToAutomatteImage();
    void copyInfo(const SampleBucket *);
    void copyInfo(const std::vector<int> &, const std::vector<int> &);
    int  fillBucket(const UT_Vector3 &, const UT_Vector3 &, SampleBucket *);
    void findBucket(const float &, const float &, 
        const float &, const float &, SampleBucket *) const;

private:
    SampleBucketV mySamples;
    UT_BoundingBox myBbox;
    int myRegisteredFlag = 0;
    SampleBucketV myNeighbours;
    size_t myNeighbourSize = 0;
public:
    uint m_resolution[2] = {0,0};
    uint m_pixelsamples[2] = {0,0};
};

typedef std::vector<SampleBucket>BucketQueue;
typedef std::queue<SampleBucket> BucketQueueQ;
#ifdef NO_MUTEX_IN_BUCKETVECTOR
typedef tbb::concurrent_hash_map<int, BucketQueue> VEX_Samples;
typedef tbb::concurrent_hash_map<int, BucketQueueQ> VEX_SamplesQ;
#else
typedef std::map<int, BucketQueue> VEX_Samples;
#endif
typedef std::map<std::string, VEX_Samples> VEX_Channels;
#if defined(TBB_VEX_STORE) && defined(NO_MUTEX_IN_BUCKETVECTOR)
typedef tbb::concurrent_hash_map<int32_t, VEX_SamplesQ> AutomatteVexCache; 
#endif
// storage per channel.
typedef std::map<int, int> BucketCounter;
// main storage container.
typedef std::array<int, 2> BucketSize;

//
typedef float coord_t;
// typedef std::vector<SampleBucket*>    BucketVector;
typedef tbb::concurrent_vector<SampleBucket>    BucketVector;
typedef tbb::concurrent_vector<Sample>          AutomatteImage;
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
// new stuff
int create_vex_storage(const std::string &, const int &, const std::vector<int> &,\
        const std::vector<int> &);
int insert_vex_sample(const int32_t &, const int &, const Sample&);
AutomatteVexCache * get_AutomatteVexCache();
AutomatteImage    * get_AutomatteImage();

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