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
	njnearest(dm, nodes, itx, ity);

	// joint the two nearest neighbors
	joint(dm, nodes, itx, ity);        
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
    double m2star  = vn.size()-2;
    for(auto &np : vn){
	(*np).length = 0.0;
	for(auto &nd : vn)
	    (*np).length += dm.getdist((*np).id, (*nd).id);
        (*np).length /= m2star;
    }
};

void njnearest(const Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){

    // get the nearest neighbor
    double minddxy = numeric_limits<double>::max();
    
    StarTree::iterator ita = nodes.begin(), itb;
    for( ; ita != nodes.end(); ++ita){
    	for(itb = ita + 1; itb != nodes.end(); ++itb){
    	    double ddxy = dm.getdist((**ita).id, (**itb).id) - (**itb).length - (**ita).length;
    	    if(ddxy < minddxy){
    		itx = ita;
    		ity = itb;
    		minddxy = ddxy;
    	    }
    	}
    }
};


void joint(Mdist& dm, StarTree& nodes, StarTree::iterator& itx, StarTree::iterator& ity){

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
    double dxdy = ((*nx).length - (*ny).length);
    double dxy  = dm.getdist((*nx).id, (*ny).id);
    (*nx).length = 0.5*(dxy + dxdy);
    (*ny).length = 0.5*(dxy - dxdy);

    // reset the distrance matrix
    // use the nx index to the index of nz in distrance matrix
    (*nz).id = (*nx).id;
    for(auto &nu : nodes){
	double dux = dm.getdist((*nx).id, (*nu).id);
	double duy = dm.getdist((*ny).id, (*nu).id);
	dm.setdist((*nz).id, (*nu).id, 0.5*(dux+duy-dxy));
    }
	 
    // add the new node to the star tree
    nodes.emplace_back(nz);
    
    // renew the length of branch
    lenStar(nodes, dm);
};


