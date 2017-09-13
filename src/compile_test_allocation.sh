g++ -O3 -mavx -std=c++11 -I$HT/include -I../../microbench/ -L/usr/lib64 -ltbb -ltbbmalloc -o test_allocation test_allocation.cpp ../../microbench/systemtime.cpp 
# $HFS/dsolib/libjemalloc.so.1
# hcustom -i . -s -ltbbmalloc test_allocation.cpp ../../microbench/systemtime.cpp
