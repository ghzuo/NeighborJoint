#include "neighbor.h"

int main(int argc, char* argv[]){
    // set a timer
    boost::timer mytimer;
    
    // get the input arguments
    Args myargs(argc,argv);

    // init the distance matrix and name list
    Mdist dm;
    vector<string> names;
    if(myargs.netcdf)
	readmtxnc(myargs.distfile, dm, names);
    else
	readmtx(myargs.distfile, dm, names);
    
    // made the star tree by listfile
    // if no, all items in distance matrix are used
    vector<Node*> nodes;
    selectLeafs(dm, myargs.listfile, nodes);

    // do the NJ algorithm and return the NJ tree
    Node* aTree = neighborJoint(dm, nodes);

    //named the leafs of the tree
    foreach(Node* nd, nodes)
	(*nd).name = names[(*nd).id];

    //set the outgroup
    Node* outgroup = *nodes.begin();
    aTree = (*aTree).resetroot(outgroup);

    //delete the outgrouup
    (*aTree).children.erase(outgroup);

    // output the Tree
    ofstream nwk(myargs.outfile.c_str());
    (*aTree).outnwk(nwk);
    nwk.close();

    cerr << "*** Time Elapsed: " << mytimer.elapsed() << "s" << endl;
}

Args::Args(int argc, char** argv):distfile("infile"),listfile(""),
				  outfile("NJtree.nwk"),netcdf(false){

    program = argv[0];
    string wkdir("");

    char ch;
    while ((ch = getopt(argc, argv, "i:d:o:D:Ch")) != -1){
      switch (ch){
      case 'i': 
	  listfile = optarg; break;
      case 'd':
	  distfile = optarg; break;
      case 'o':
	  outfile = optarg; break;
      case 'D':
	  wkdir = optarg; break;
      case 'C':
	  netcdf = true; break;
      case 'h':
	  usage();
      case '?':
	  usage();
      }
    }

    if(wkdir != ""){
	addsuffix(wkdir, '/');
	distfile = wkdir + distfile;
	outfile = wkdir + outfile;
    }
}

void Args::usage(){
    cerr << "\nProgram Usage (VERSION: " << GHZ_VERSION << "): \n\n" 
	 << program  <<"\n" 
	 <<" [ -d infile ]        input distance matrix, defaut: infile\n"
	 <<" [ -o NJtree.nwk ]    output newick tree, defaut: NJtree.nwk\n"
	 <<" [ -i <listfile> ]    the selection index list of the distance matrix\n"
	 <<"                      if no defined, whole distance matrix are used\n"
	 <<" [ -C ]               use the netcdf input format, default false\n"
	 <<" [ -h ]               disply this information\n"
	 << endl;
    exit(1);
}

// select the leafs to of the output tree
void selectLeafs(const Mdist& dm, const string& listfile, vector<Node*>& nodes){
    if(listfile.empty()){
	for(size_t i=0; i<dm.size(); ++i){
	    Node* np = new Node;
	    (*np).id = i;
	    nodes.push_back(np);
	}
    }else{
	ifstream inf(listfile.c_str());
	if(!inf.is_open()){
	    cerr << "Error opening file: " << listfile << endl;
	    exit(3);
	}
	
	while(inf.good()){
	    Node* np = new Node;
	    inf >> (*np).id;
	    nodes.push_back(np);
	}
    }
}


/// For neighborJoint
Node* neighborJoint(Mdist& dm, const StarTree& leafs){
    // the function to neighbor joint the start tree
    // the distance matrix and leafs list are copy from the main program
    // so it will be changed in the main program

    StarTree nodes;
    foreach(Node* nd, leafs)
	nodes.push_back(nd);

    // get the new length of the star tree
    lenStar(nodes,dm);

    while(nodes.size() > 3){
	// get two nearest items and set length and distance matrix
	StarTree::iterator itx, ity;
	njnearest(dm, nodes, itx, ity);

	// joint the two nearest neighbors
	recjoint(dm, nodes, itx, ity);
	// joint(dm, nodes, itx, ity);

    }

    // get the root of the tree which include three branches
    (*nodes[0]).length= 0.5*((*nodes[0]).length - dm.getdist((*nodes[1]).id,(*nodes[2]).id));
    (*nodes[1]).length= 0.5*((*nodes[1]).length - dm.getdist((*nodes[2]).id,(*nodes[0]).id));
    (*nodes[2]).length= 0.5*((*nodes[2]).length - dm.getdist((*nodes[0]).id,(*nodes[1]).id));

    Node* aTree = new Node;
    foreach(Node* nd, nodes)
	(*aTree).addChild(*nd);

    return aTree;
}

void lenStar(const StarTree& vn, const Mdist& dm){
    foreach(Node* np, vn){
	(*np).length = 0.0;
	foreach(Node* nd, vn)
	    (*np).length +=dm.getdist((*np).id, (*nd).id);
    }
};

void njnearest(const Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){

    // get the nearest neighbor
    double minddxy = numeric_limits<double>::max();
    double m2star = nodes.size()-2;

    vector<double> length;
    foreach(Node* nd, nodes)
    	length.push_back((*nd).length/m2star);

    vector<double>::iterator iterA = length.begin();
    StarTree::iterator ita = nodes.begin();
    for( ; ita != nodes.end(); ++ita){
	StarTree::iterator itb = ita+1;
	vector<double>::iterator iterB=iterA;
    	for( ; itb != nodes.end(); ++itb){
    	    double ddxy = dm.getdist((**ita).id, (**itb).id) - *iterA - *(++iterB);
    	    if(ddxy < minddxy){
    		itx = ita;
    		ity = itb;
    		minddxy = ddxy;
    	    }
    	}
	++iterA;
    }

    // for(StarTree::iterator ita=nodes.begin(); ita != nodes.end(); ++ita){
    // 	double dx = (*(*ita)).length;
    // 	for(StarTree::iterator itb=ita+1; itb != nodes.end(); ++itb){
    // 	    double dy = (*(*itb)).length;
    // 	    double ddxy = m2star*dm.getdist((*(*ita)).id, (*(*itb)).id) - dx - dy;
    // 	    if(ddxy < minddxy){
    // 		itx = ita; 
    // 		ity = itb;
    // 		minddxy = ddxy;
    // 	    }
    // 	}
    // }
};



void joint(Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){

    // the new parent node
    Node* nz = new Node;
    
    // move the two nearest nodes to the new node
    Node* ny = *ity;
    nodes.erase(ity);
    (*nz).addChild(*ny);
	
    Node* nx = *itx;
    nodes.erase(itx);
    (*nz).addChild(*nx);

    // reset the length branch
    double dxdy = ((*nx).length - (*ny).length)/nodes.size();
    double dxy  = dm.getdist((*nx).id, (*ny).id);
    (*nx).length = 0.5*(dxy + dxdy);
    (*ny).length = 0.5*(dxy - dxdy);

    // reset the distrance matrix
    // use the nx index to the index of nz in distrance matrix
    (*nz).id = (*nx).id;
    foreach(Node* nu, nodes){
	double dux = dm.getdist((*nx).id, (*nu).id);
	double duy = dm.getdist((*ny).id, (*nu).id);
	double duz = 0.5*(dux + duy - dxy);
	dm.setdist((*nz).id, (*nu).id, duz);
    }

    // add the new node to the star tree
    nodes.push_back(nz);
    
    // renew the length of branch
    lenStar(nodes, dm);
};

// 
void recjoint(Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){
    // the new parent node
    Node* nz = new Node;

    // add the two nodes to the new node
    Node* ny = *ity;
    nodes.erase(ity);
    (*nz).addChild(*ny);
	
    Node* nx = *itx;
    nodes.erase(itx);
    (*nz).addChild(*nx);

    // reset the length branch
    double m2star = nodes.size();
    double dx = (*nx).length;
    double dy = (*ny).length;
    double dxdy = (dx - dy)/m2star;
    double dxy  = dm.getdist((*nx).id, (*ny).id);
    (*nx).length = 0.5*(dxy + dxdy);
    (*ny).length = 0.5*(dxy - dxdy);
    (*nz).length = 0.5*(dx + dy - dxy*(m2star+2));
    
    // set the distrance between node z and other nodes
    // use the nx index to the index of nz in distrance matrix
    (*nz).id = (*nx).id;
    //(*nz).length = 0.0;
    foreach(Node* nu, nodes){
	double dux = dm.getdist((*nx).id, (*nu).id);
	double duy = dm.getdist((*ny).id, (*nu).id);
	double duz = 0.5*(dux + duy - dxy);
	dm.setdist((*nz).id, (*nu).id, duz);
	(*nu).length += (duz - dux - duy);
	// (*nz).length += duz;
    }

    // add the new node to the star tree
    nodes.push_back(nz);
}

