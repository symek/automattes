#pragma once

#ifndef __AutomattesHelper__
#define __AutomattesHelper__

namespace HA_HDK {
// our fixed sample: x,y,z,id,Af 
typedef std::vector<float> Sample; 
// vector of samples per thread (reused by many buckets)
typedef std::vector<Sample> SampleBucket;
// storage per channel.
typedef std::map<int, SampleBucket> VEX_Samples;
// main storage container.
typedef std::map<std::string, VEX_Samples> VEX_Channels;
// tmp
static constexpr const int& size_at_least = 16*16*3*3;

// pointgrid stuff useful for buliding filter side accesor.
typedef  UT_PointGridVector3ArrayAccessor<int, int> UT_Vector3Point;
typedef  UT_PointGrid<UT_Vector3Point>::queuetype UT_Vector3PointQueue;

// function exposed on vex side (temporarily instead of proper class)
int VEX_Samples_create(const int&);
int VEX_Samples_insert(const int&, const Sample&);
VEX_Samples * VEX_Samples_get();

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