#pragma once
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#define DEBUG
#ifdef DEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif

namespace BENCHMARK
{

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVirtualMem(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result / 1024;
}

int getPhysicalMem(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result / 1024;
}


static unsigned long x_seed=123456789, y_seed=362436069, z_seed=521288629;

inline unsigned long xorshf96(void) 
{   //period 2^96-1
    unsigned long t;
    x_seed ^= x_seed << 16;
    x_seed ^= x_seed >> 5;
    x_seed ^= x_seed << 1;

   t = x_seed;
   x_seed = y_seed;
   y_seed = z_seed;
   z_seed = t ^ x_seed ^ y_seed;

  return z_seed;
}

} // end of BENCHMARK namespace