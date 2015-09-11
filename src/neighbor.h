#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <limits>
//#include <omp.h>

#include "stringOpt.h"
#include "global.h"
#include "tree.h"
#include "distmatrix.h"
using namespace std;

// read arguments
struct Args{
    string program, distfile, listfile, outfile;
    bool netcdf;
    
    Args(int, char**);
    void usage();
};

typedef vector<Node*> StarTree;

// select the genomes to build the tree
void selectLeafs(const Mdist&, const string&, vector<Node*>&);

// build the tree by neighbor joint algorithm
Node* neighborJoint(Mdist&, const vector<Node*>&);

// the distance from the star point
void lenStar(const StarTree&, const Mdist&);

// reset the distance of the nearest neighbor
void njnearest(const Mdist&, StarTree&, StarTree::iterator&, StarTree::iterator&);

// joint the two neighbors

// reset the distance of the nearest neighbor
void joint(Mdist&, StarTree&, StarTree::iterator&, StarTree::iterator&);
void recjoint(Mdist&, StarTree&, StarTree::iterator&, StarTree::iterator&);

#endif
