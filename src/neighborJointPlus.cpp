#include "neighborJoint.h"

/// For neighborJoint
Node* neighborJoint(Mdist& dm){
    // the function to neighbor joint the start tree
    // the distance matrix and leafs list are copy from the main program
    // so it will be changed in the main program

    StarTree nodes(dm.size());
    for(size_t i=0; i<dm.size(); ++i){
        nodes[i] = new Node(i, dm.getname(i));
    }
    Node* outgrp = nodes[0];

    // get the new length of the star tree
    lenStar(nodes, dm);

    while(nodes.size() > 3){
	// get two nearest items and set length and distance matrix
	StarTree::iterator itx, ity;
        njnearestPlus(dm, nodes, itx, ity);

	// joint the two nearest neighbors
	recjoint(dm, nodes, itx, ity);
    }

    // get the root of the tree which include three branches
    (*nodes[0]).length= 0.5*((*nodes[0]).length - dm.getdist((*nodes[1]).id,(*nodes[2]).id));
    (*nodes[1]).length= 0.5*((*nodes[1]).length - dm.getdist((*nodes[2]).id,(*nodes[0]).id));
    (*nodes[2]).length= 0.5*((*nodes[2]).length - dm.getdist((*nodes[0]).id,(*nodes[1]).id));

    Node* aTree = new Node;
    for(auto &nd : nodes)
	(*aTree).addChild(nd);

    aTree = (*aTree).resetroot(outgrp);
    
    return aTree;
}

void lenStar(const StarTree& vn, const Mdist& dm){
    for(auto &np : vn){
	(*np).length = 0.0;
	for(auto &nd : vn)
	    (*np).length += dm.getdist((*np).id, (*nd).id);
    }
};

void njnearestPlus(const Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){

    // get the nearest neighbor
    double minddxy(numeric_limits<double>::max());
    size_t nNode(nodes.size());
    double m2star(nNode);
    m2star -= 2.0;
    
    vector<double> length(nNode);
    vector<size_t> id(nNode);
    for(int i=0; i<nNode; ++i){
	length[i] = (*(nodes[i])).length/m2star;
        id[i] = (*(nodes[i])).id;
    }

    int nx, ny;
    for(int i=1; i<nNode; ++i){
        size_t iId = id[i];
        double iLen = length[i];
        for(int j=0; j<i; ++j){
    	    double ddxy = dm.getdist(iId, id[j]) - iLen - length[j];
    	    if(ddxy < minddxy){
                nx = j; ny = i;
    		minddxy = ddxy;
    	    }
    	}
    }

    // get the two iterator
    itx = nodes.begin() + nx;
    ity = nodes.begin() + ny;
};

void recjoint(Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){

    // the new parent node
    Node* nz = new Node;

    // add the two nodes to the new node
    Node* ny = *ity;
    nodes.erase(ity);
    (*nz).addChild(ny);
	
    Node* nx = *itx;
    nodes.erase(itx);
    (*nz).addChild(nx);

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
    for(auto &nu : nodes){
	double dux = dm.getdist((*nx).id, (*nu).id);
	double duy = dm.getdist((*ny).id, (*nu).id);
	double duz = 0.5*(dux + duy - dxy);
	dm.setdist((*nz).id, (*nu).id, duz);
	(*nu).length = (*nu).length - dux - duy + duz;
    }

    // add the new node to the star tree
    nodes.emplace_back(nz);
}

