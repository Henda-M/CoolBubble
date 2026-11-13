# CoolBubble makefile
#
# Needs a C++17 compiler with OpenMP, FFTW3 and the HDF5 C++ library.
#
#   make                          build ./coolbubble
#   make HDF5_DIR=/path/to/hdf5   use an HDF5 installation outside the system paths
#   make CXX=h5c++                use the HDF5 compiler wrapper (common on clusters)
#   make clean                    remove the build files
#
# Settings for your own machine can go in a file `make.local` next to this one
# (not tracked by git), for example:
#   HDF5_DIR = /home/me/hdf5
#   FFTW_DIR = /opt/fftw

-include make.local

CXX      ?= g++
CXXFLAGS ?= -O2          # optimisation and warning flags, e.g. make CXXFLAGS="-O3 -march=native"

# flags the code needs in any case
REQFLAGS = -std=c++17 -fopenmp -MMD -MP -DPYTHON_DIR='"$(CURDIR)/python"'
LIBS     = -lhdf5_cpp -lhdf5 -lfftw3 -lz -lm

ifdef HDF5_DIR
REQFLAGS += -I$(HDF5_DIR)/include
LIBS     += -L$(HDF5_DIR)/lib -Wl,-rpath,$(HDF5_DIR)/lib
endif
ifdef FFTW_DIR
REQFLAGS += -I$(FFTW_DIR)/include
LIBS     += -L$(FFTW_DIR)/lib -Wl,-rpath,$(FFTW_DIR)/lib
endif

TARGET = coolbubble
SRCS   = src/main.cpp src/bubble.cpp src/gridfuncs.cpp
OBJS   = $(SRCS:src/%.cpp=build/%.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(REQFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) $(REQFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	rm -rf build $(TARGET)

.PHONY: clean

-include $(OBJS:.o=.d)
