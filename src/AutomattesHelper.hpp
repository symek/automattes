#pragma once

#ifndef __AutomattesHelper__
#define __AutomattesHelper__

namespace HA_HDK {
// our fixed sample: x,y,z,id,Af 
typedef std::vector<float> OpacitySample; 
// vector of samples per thread (reused by many buckets)
typedef std::vector<OpacitySample> SampleBucket;
// main storage
typedef std::map<int, SampleBucket> VEX_Samples;
// tmp
static constexpr const int& size_at_least = 16*16*3*3;
// function exposed on vex side
int VEX_Samples_create(const int&);
int VEX_Samples_insert(const int&, const OpacitySample&);
VEX_Samples * VEX_Samples_get();

} // end of HA_HDK Space

#endif