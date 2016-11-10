QUESO_DIR := /h1/dsondak/bin/queso
ANTIOCH_DIR := /workspace/dsondak/antioch_build

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

LIBS := \
	-L$(QUESO_DIR)/lib -lqueso \
	-L$(ANTIOCH_DIR)/lib -lantioch \
	-L$(BOOST_DIR)/lib -lboost_program_options \
	-L$(GSL_DIR)/lib -lgsl -lgslcblas\
	-L$(GRVY_DIR)/lib -lgrvy \
	-L$(MPI_DIR)/lib \
	-L$(HDF5_DIR)/lib -lhdf5 

CXX = mpic++
CXXFLAGS += -O3 -g -Wall -c -std=c++0x

default: ip-catchall

.SUFFIXES: .o .c

all:	ip-catchall

clean:
	rm -f *~
	rm -f *.o
	rm -f ip-catchall

ip-catchall: ip_catchall.o compute.o problem_size.o likelihood.o reaction_info.o  model.o inadequacy_model.o write_data.o truth_data.o chemistryVectorRV.o chemistryJointPdf.o
	$(CXX) ip_catchall.o \
	       compute.o \
               problem_size.o \
	       likelihood.o \
               reaction_info.o \
               model.o \
               inadequacy_model.o \
               write_data.o \
               truth_data.o \
               chemistryVectorRV.o \
	       chemistryJointPdf.o \
	       -o ip-catchall $(LIBS)

%.o: %.c
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
