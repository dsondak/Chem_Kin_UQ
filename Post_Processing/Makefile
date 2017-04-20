ANTIOCH_DIR := /workspace/dsondak/antioch_build

INC_PATHS := \
	-I. \
	-I$(HDF5_DIR)/include \
	-I$(EIGEN_INCLUDE) \
        -I$(ANTIOCH_DIR)/include

LIBS := \
	-L$(ANTIOCH_DIR)/lib -lantioch \
        -L$(HDF5_DIR)/lib -lhdf5

CXX = g++
CXXFLAGS += -Wall -O3 -c -std=c++11

default: calc_rhs

.SUFFIXES: .o .C

all:	calc_rhs

clean:
	rm -f *~
	rm -f *.o
	rm -f calc_rhs

calc_rhs: calc_rhs.o write_data.o truth_data.o
	$(CXX) calc_rhs.o \
	       write_data.o \
	       truth_data.o \
	       -o calc_rhs $(LIBS)

%.o: %.c
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
