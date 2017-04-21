ANTIOCH_DIR := /workspace/dsondak/antioch_build
CVODE_DIR := /workspace/dsondak/sundials/install
NLOPT_DIR := /workspace/dsondak/nlopt/install

INC_PATHS := \
	-I. \
        -I$(CVODE_DIR)/include \
	-I$(HDF5_DIR)/include \
	-I$(EIGEN_INCLUDE) \
        -I$(ANTIOCH_DIR)/include \
        -I$(NLOPT_DIR)/include

LIBS := \
	-L$(ANTIOCH_DIR)/lib -lantioch \
        -L$(HDF5_DIR)/lib -lhdf5 \
	-L$(CVODE_DIR)/lib -lsundials_cvode -lsundials_nvecserial -lm /usr/lib64/librt.so \
        -L$(NLOPT_DIR)/lib -lnlopt

LINKFLAGS = -Wl,-rpath,/workspace/dsondak/sundials/install/lib

CXX = g++
#CXXFLAGS += -ggdb -Wall -c -std=c++0x
CXXFLAGS += -Wall -O3 -c -std=c++11

default: kinetics

.SUFFIXES: .o .C

all:	kinetics

clean:
	rm -f *~
	rm -f *.o
	rm -f kinetics

kinetics: kinetics.o reaction_info.o problem_size.o write_data.o
	$(CXX) kinetics.o \
	       reaction_info.o \
               problem_size.o \
	       write_data.o \
	       -o kinetics $(LIBS) ${LINKFLAGS}

%.o: %.c
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
