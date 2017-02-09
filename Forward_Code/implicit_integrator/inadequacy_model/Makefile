ANTIOCH_DIR := /workspace/dsondak/antioch_build
#CVODE_DIR := /workspace/dsondak/sundials/install
CVODE_DIR := /workspace/dsondak/sundials_gdb/install_gdb

INC_PATHS := \
	-I. \
        -I$(CVODE_DIR)/include \
	-I$(HDF5_DIR)/include \
	-I$(EIGEN_INCLUDE) \
        -I$(ANTIOCH_DIR)/include

LIBS := \
	-L$(ANTIOCH_DIR)/lib -lantioch \
        -L$(HDF5_DIR)/lib -lhdf5 \
	-L$(CVODE_DIR)/lib -lsundials_cvode -lsundials_nvecserial -lm /usr/lib64/librt.so

LINKFLAGS = -Wl,-rpath,/workspace/dsondak/sundials_gdb/install_gdb/lib

CXX = g++
#CXXFLAGS += -g -O0 -Wall -c -std=c++0x
CXXFLAGS += -O3 -Wall -c -std=c++11

default: cvRoberts_dns_eigs

.SUFFIXES: .o .C

all:	cvRoberts_dns_eigs

clean:
	rm -f *~
	rm -f *.o
	rm -f cvRoberts_dns_eigs

cvRoberts_dns_eigs: cvRoberts_dns_eigs.o reaction_info.o problem_size.o write_data.o
	$(CXX) cvRoberts_dns_eigs.o \
	       reaction_info.o \
               problem_size.o \
	       write_data.o \
	       -o cvRoberts_dns_eigs $(LIBS) ${LINKFLAGS}

%.o: %.c
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
