NEIGHBOR = tree.o stringOpt.o neighbor.o distmatrix.o

BINS = neighbor 

FC = gfortran
CC = gcc

CXX = g++
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CXX = g++-4.2
endif

FFLAGS  = -O3
CFLAG = -O3  
LIBS = -lz -lnetcdf -lnetcdf_c++ 

all: $(BINS)

neighbor: $(NEIGHBOR)
	$(CXX) $(NEIGHBOR) -o $@ $(LIBS) $(CFLAGS)	

#############################################################

clean:
	rm -f *.o

cleanbin:
	rm -f $(BINS)

cleanall: clean cleanbin

##############################################################
.SUFFIXES: .f90 .f .o .mod .F .c .cpp

.c.o:
	$(CC)  $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

.f.o:
	$(FC)  $(FFLAGS) -c $<

.f90.o:
	$(FC)  $(FFLAGS) -c $<
