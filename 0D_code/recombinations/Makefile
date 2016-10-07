#ANTIOCH_DIR := /h1/dsondak/bin/antioch_new
ANTIOCH_DIR := /workspace/dsondak/antioch_build

INC_PATHS := \
	-I. \
	-I$(GSL_DIR)/include \
	-I$(GRVY_DIR)/include \
	-I$(HDF5_DIR)/include \
	-I$(EIGEN_INCLUDE) \
	-I$(MPI_DIR)/ \
        -I$(ANTIOCH_DIR)/include

LIBS := \
	-L$(ANTIOCH_DIR)/lib -lantioch \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas \
	-L$(GRVY_DIR)/lib -lgrvy \
	-L$(MPI_DIR)/lib \
	-L$(HDF5_DIR)/lib -lhdf5 

CXX = g++
#CXX = mpic++ #g++
CXXFLAGS += -O3 -g -Wall -c -std=c++0x

default: forward_model

.SUFFIXES: .o .C

all:	forward_model

clean:
	rm -f *~
	rm -f *.o
	rm -f forward_model

forward_model: forward_model.o model.o reaction_info.o problem_size.o write_data.o
	$(CXX) forward_model.o \
	       model.o \
	       reaction_info.o \
               problem_size.o \
	       write_data.o \
	       -o forward_model $(LIBS)

%.o: %.c
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
