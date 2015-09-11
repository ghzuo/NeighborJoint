#include "distmatrix.h"

Mdist::Mdist():ng(0){};

double Mdist::getdist(size_t i, size_t j) const{
    if(i<j)
	return _getdist(i,j);
    return _getdist(j,i);
};

double Mdist::_getdist(size_t i, size_t j) const{
    if(j>=ng){
	cerr << "Error: the index out of matrix" << endl;
	exit(4);
    }

    return dist[i+(j+1)*j/2];
};

void Mdist::setdist(size_t i, size_t j, double d){
    if(i<j)
	return _setdist(i,j, d);
    return _setdist(j,i, d);
};

void Mdist::_setdist(size_t i, size_t j, double d){
    if(j>=ng){
	cerr << "Error: the index out of matrix" << endl;
	exit(4);
    }

    dist[i+j*(j+1)/2] = d;
};

void Mdist::push_back(const vector<double>& dd){
    if(dd.size() != ng+1){
	cerr << "Error: the number of the append vector is unmatched!!" << endl;
	exit(4);
    }

    ++ng;
    foreach(double d, dd)
	dist.push_back(d);    
};

void Mdist::erase(size_t n){
    if(n>=ng){
	cerr << "Error: the index out of matrix in delete" << endl;
	exit(4);
    }
    
    vector<double>::iterator iter=dist.begin();
    for(size_t i=ng-1; i>n; --i)
	dist.erase(iter+(i*(i+1)/2+n));

    dist.erase(iter+(n*(n+1)/2), iter+((n+1)*(n+2)/2)-1);
    
};

void Mdist::extend(size_t n){
    resize(n+ng);
};

void Mdist::resize(size_t n){
    ng = n;
    dist.resize(n*(n+1)/2);
};

size_t Mdist::msize() const{
    return ng*(ng+1)/2;
};

size_t Mdist::size() const{
    return ng;
};

size_t Mdist::capacity() const{
    return dist.capacity()*sizeof(double) + sizeof(ng);
};

// read the tranditional infile for the distance matrix
void readmtx(const string& file, Mdist& dm, vector<string>& name){
    if(getsuffix(file) == "nc")
	readmtxnc(file, dm, name);
    else
	readmtxtxt(file, dm, name);
};

// read the tranditional infile for the distance matrix
void readmtxtxt(const string& file, Mdist& dm, vector<string>& name){
    ifstream dd(file.c_str());
    if(!dd.is_open()){
	cerr << "Error opening file: " << file << endl;
	exit(1);
    }

    double ng; dd >> ng;
    for(size_t i=0; i<ng; ++i){
	string str; dd >> str; name.push_back(str);
	vector<double> vd;
	for(size_t j=0; j<ng; ++j){
	    double d; dd >> d;
	    if(i>=j) vd.push_back(d);
	}
	dm.push_back(vd);
    }
};

// read the netcdf file for the distance matrix
void readmtxnc(const string& file, Mdist& dm, vector<string>& name){
    NcFile mtxFile(file.c_str(), NcFile::ReadOnly);
    if (!mtxFile.is_valid()){
	cout << "Couldn't open file!\n";
	cerr <<  NC_ERR << endl;
    }

    //get the dim information

    size_t lenWord=mtxFile.get_dim("word")->size();
    size_t ngenome=mtxFile.get_dim("genome")->size();
    size_t nmtx=mtxFile.get_dim("matrix")->size();
    
    // read the space name
    NcVar *spname=mtxFile.get_var("spname");
    for(int i=0; i<ngenome; ++i){
	long counts[2] = {1, lenWord};
	char aname[lenWord];
	spname->set_cur(i,0);
	spname->get(aname, counts);
	name.push_back(aname);
    }

    // read the distance matrix
    float *dist = new float[nmtx];
    NcVar* distance=mtxFile.get_var("distance");
    distance->get(dist, nmtx);
    dm.resize(ngenome);
    int index(0);
    for(int i=0; i<ngenome; ++i)
	for(int j=0; j<=i; ++j)
	    dm.setdist(i,j,dist[index++]);    
    mtxFile.close();
    delete [] dist;
};


// select the distance matrix file type
void writemtx(const string& file, const Mdist& dm, const vector<string>& name){
    if(getsuffix(file) == "nc")
	writemtxnc(file, dm, name);
    else
	writemtxtxt(file, dm, name);
};


// write the tranditional infile for the distance matrix
void writemtxtxt(const string& file, const Mdist& dm, const vector<string>& name){
    ofstream dd(file.c_str());
    if(!dd.is_open()){
	cerr << "Error opening file: " << file << endl;
	exit(1);
    }
	
    size_t ng = name.size();

    dd << ng << endl;
    for(size_t i=0; i<ng; ++i){
	dd << name[i] << " ";
	for(size_t j=0; j<ng; ++j)
	    dd << left << fixed << setprecision(10) 
	       << setw(13) << dm.getdist(i, j);
	dd << endl;
    }
}

// write the tranditional infile for the distance matrix
void writemtxnc(const string& file, const Mdist& dm, const vector<string>& name){
    NcFile mtxFile(file.c_str(), NcFile::Replace);
    if (!mtxFile.is_valid()){
	cout << "Couldn't open file!\n";
	cerr <<  NC_ERR << endl;
    }

    //output the space name
    size_t wlength(0);
    foreach(const string& str, name)
	if(str.size()>wlength)
	    wlength = str.size();
    ++wlength;

    NcDim* lenWord=mtxFile.add_dim("word", wlength);
    NcDim* ngenome=mtxFile.add_dim("genome", name.size());
    NcVar *spname=mtxFile.add_var("spname", ncChar, ngenome, lenWord);
    for(long i=0; i<name.size(); ++i){
    	long counts[2] = {1,name[i].size()+1};
    	spname->set_cur(i,0);
    	spname->put(name[i].c_str(),counts);
    }

    // output the distance matrix
    NcDim* mtxdist=mtxFile.add_dim("matrix", dm.msize());
    NcVar *distance=mtxFile.add_var("distance", ncFloat, mtxdist);
    float *vdm = new float[dm.msize()];
    size_t index(0);
    for(int i=0; i<name.size(); ++i)
    	for(int j=0; j<=i; ++j)
    	    vdm[index++] = dm.getdist(i,j);
    distance->put(vdm, dm.msize());
    delete [] vdm;

    mtxFile.close();
}

// readjust the distance matrix and their name by the index list
void adjustmtx(const vector<size_t>& ndx, Mdist& dm, vector<string>& name){
    Mdist ndm;
    vector<string> nm;

    ndm.resize(ndx.size());
    nm.resize(ndx.size());
    for(int i=0; i<ndx.size(); ++i){
	nm[i] = name[ndx[i]];
	for(int j=0; j<=i; ++j)
	    ndm.setdist(i,j,dm.getdist(ndx[i],ndx[j]));
    }

    dm = ndm;
    name = nm;
}
