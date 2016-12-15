NEIGHBOR = tree.o stringOpt.o neighbor.o distmatrix.o neighborJoint.o
NJPLUS   = tree.o stringOpt.o neighbor.o distmatrix.o neighborJointPlus.o

BINS = neighbor njplus

CC = gcc -std=c11
CXX = g++ -std=c++11

LIBS = -lz -lnetcdf -lnetcdf_c++ -L/usr/local/lib 
ifeq ($(DEF), DEBUG)
        CFLAGS = -g -Wall -Wextra -I/usr/local/include
else
	ifeq ($(DEF), PROF)
	        CFLAGS = -O3 -pg -I/usr/local/include
	else
		CFLAGS = -O3 -I/usr/local/include
	endif
endif

all: $(BINS)

neighbor: $(NEIGHBOR)
	$(CXX) $(NEIGHBOR) -o $@ $(LIBS) $(CFLAGS)	

njplus: $(NJPLUS)
	$(CXX) $(NJPLUS) -o $@ $(LIBS) $(CFLAGS)	

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
