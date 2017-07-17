QUESO_DIR := /h1/dsondak/bin/queso
ANTIOCH_DIR := /workspace/dsondak/antioch_build
CVODE_DIR := /workspace/dsondak/sundials/install
NLOPT_DIR := /workspace/dsondak/nlopt/install

INC_PATHS := \
	-I. \
	-I$(QUESO_DIR)/include \
	-I$(BOOST_DIR)/include \
	-I$(GSL_DIR)/include \
	-I$(GRVY_DIR)/include \
	-I$(HDF5_DIR)/include \
	-I$(GLPK_DIR)/include \
	-I$(MPI_DIR)/ \
	-I$(EIGEN_INCLUDE) \
        -I$(ANTIOCH_DIR)/include \
        -I$(CVODE_DIR)/include \
        -I$(NLOPT_DIR)/include

LIBS := \
	-L$(QUESO_DIR)/lib -lqueso \
	-L$(ANTIOCH_DIR)/lib -lantioch \
	-L$(BOOST_DIR)/lib -lboost_program_options \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas\
	-L$(GRVY_DIR)/lib -lgrvy \
	-L$(MPI_DIR)/lib \
	-L$(HDF5_DIR)/lib -lhdf5 \
	-L$(CVODE_DIR)/lib -lsundials_cvode -lsundials_nvecserial -lm /usr/lib64/librt.so \
        -L$(NLOPT_DIR)/lib -lnlopt

LINKFLAGS = -Wl,-rpath,/workspace/dsondak/sundials/install/lib

CXX = mpic++
CXXFLAGS += -O3 -Wall -c -std=c++11

default: ip-catchall

.SUFFIXES: .o .c

all:	ip-catchall

clean:
	rm -f *~
	rm -f *.o
	rm -f ip-catchall

ip-catchall: ip_catchall.o compute.o problem_size.o LinearInterpolation.o likelihood.o reaction_info.o  kinetics_solve.o inadequacy_model.o write_data.o truth_data.o chemistryVectorRV.o chemistryJointPdf.o
	$(CXX) ip_catchall.o \
	       compute.o \
               problem_size.o \
               LinearInterpolation.o \
	       likelihood.o \
               reaction_info.o \
               kinetics_solve.o \
               inadequacy_model.o \
               write_data.o \
               truth_data.o \
               chemistryVectorRV.o \
	       chemistryJointPdf.o \
	       -o ip-catchall $(LIBS) ${LINKFLAGS}

%.o: %.c
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
