#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <tr1/unordered_map>
#include <math.h>

#include "stringOpt.h"
#include "global.h"
using namespace std;

struct Node{
    string name;
    size_t id;
    double length;
    Node* parent;
    set<Node*> children;

    size_t taxSize, nleaf, nxleaf;
    bool todefine;

    Node();

    void clear();
    void addChild(Node&);
    void getDescendants(vector<Node*>&);
    void getLeafs(vector<Node*>&);
    bool isLeaf();
    
    size_t checkTodefine();
    void getDefineLeafs(vector<Node*>&);

    void _getPrediction(string&);
    void outPrediction(ostream&);

    Node* resetroot(const string&);
    Node* resetroot(Node*);

    void _outnwk(ostream&);
    void outnwk(ostream&);
    void _innwk(istream&);
    void innwk(istream&);

    void outjson(ostream&);
    void outjsonAbbr(ostream&);
    void renewId(const tr1::unordered_map<string,size_t>&);
};

#endif
