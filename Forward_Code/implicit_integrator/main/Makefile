ANTIOCH_DIR := /n/home01/dsondak/software/antioch
CVODE_DIR := /n/home01/dsondak/sundials/install
NLOPT_DIR := /n/home01/dsondak/software/nlopt

INC_PATHS := \
	-I. \
        -I$(CVODE_DIR)/include \
        -I$(ANTIOCH_DIR)/include \
        -I$(NLOPT_DIR)/include \
	-I$(EIGEN_INCLUDE) \
	-I$(HDF5_DIR)/include

LIBS := -L$(ANTIOCH_DIR)/lib -lantioch \
	-L$(CVODE_DIR)/lib -lsundials_cvode -lsundials_nvecserial -lm /usr/lib64/librt.so \
        -L$(NLOPT_DIR)/lib -lnlopt \
        -L$(HDF5_DIR)/lib -lhdf5

#LINKFLAGS = -Wl,-rpath,/n/home01/dsondak/software/antioch/lib
LINKFLAGS = -Wl,-rpath,/n/home01/dsondak/software/antioch/lib,-rpath,/n/home01/dsondak/sundials/install/lib

CXX = g++
#CXXFLAGS += -Wall -O3 -c -std=c++0x
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
