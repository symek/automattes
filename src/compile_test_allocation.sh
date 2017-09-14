g++ -O3 -mavx -g -std=c++11 -I$HT/include -I../../microbench/ -L/$HFS/dsolib -ltbb -ltbbmalloc -o test_allocation test_allocation.cpp ../../microbench/systemtime.cpp 
# $HFS/dsolib/libjemalloc.so.1
# hcustom -i . -s -ltbbmalloc test_allocation.cpp ../../microbench/systemtime.cpp
