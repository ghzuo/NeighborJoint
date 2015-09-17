NEIGHBOR = tree.o stringOpt.o neighbor.o distmatrix.o neighborJoint.o

BINS = neighbor 

CC = gcc -std=c11
CXX = g++ -std=c++11

CFLAGS = -O3 -I/usr/local/include
LIBS = -lz -lnetcdf -lnetcdf_c++ -L/usr/local/lib 

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
