#ifndef DISTMATRIX_H
#define DISTMATRIX_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <netcdfcpp.h>

#include "tree.h"
#include "global.h"

using namespace std;

class Mdist{
    size_t ng;
    vector<double> dist;

    void _setdist(size_t, size_t, double);
    double _getdist(size_t, size_t) const;

public:
    Mdist();

    double getdist(size_t, size_t) const;
    void setdist(size_t, size_t, double);

    void push_back(const vector<double>&);
    void erase(size_t);

    void resize(size_t);
    void extend(size_t);

    size_t size() const;
    size_t msize() const;
    size_t capacity() const;
};

// read and write the tranditional infile for distance matrix and name list
void readmtx(const string&, Mdist&, vector<string>&);
void readmtxnc(const string&, Mdist&, vector<string>&);
void readmtxtxt(const string&, Mdist&, vector<string>&);
void writemtx(const string&, const Mdist&, const vector<string>&);
void writemtxnc(const string&, const Mdist&, const vector<string>&);
void writemtxtxt(const string&, const Mdist&, const vector<string>&);

// readjust the distance matrix and their name by the index list
void adjustmtx(const vector<size_t>&, Mdist&, vector<string>&);
#endif

#ifndef NCERR
#define NCERR
static const int NC_ERR = 10;
#endif

