# Minimal Automattes Makefile
INSTDIR = $(HIH)
SOURCES =  ./src/MurmurHash3.cpp ./src/VRAY_AutomattesFilter.cpp 
OPTIMIZER = -O3
DSONAME = VRAY_AutomattesFilter.so
# Include HDK Makefile.
include $(HT)/makefiles/Makefile.gnu

