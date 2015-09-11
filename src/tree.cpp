#include "tree.h"

/*******************************************************************/
/********* Member Functions For Node Class *************************/
/*******************************************************************/
Node::Node():parent(NULL),name(""),length(NAN),nxleaf(0),
	     id(0),taxSize(1),nleaf(1),todefine(false){};

bool Node::isLeaf(){
    return children.empty();
};

void Node::addChild(Node& nd){
    children.insert(&nd);
    nd.parent = this;
};

void Node::clear(){
    vector<Node*> nodes;
    getDescendants(nodes);
    vector<Node*>::iterator iter=nodes.begin();
    vector<Node*>::iterator iterEnd=nodes.end();
    for(;iter!=iterEnd; ++iter)
	delete *iter;
    children.clear();
}

void Node::getDescendants(vector<Node*>& nds){
    if(!children.empty()){
	set<Node*>::const_iterator iter = children.begin();
	set<Node*>::const_iterator iterEnd = children.end();
	for(;iter!=iterEnd; ++iter){
	    nds.push_back(*iter);
	    (*(*iter)).getDescendants(nds);
	}
    }
};

void Node::getLeafs(vector<Node*>& nodes){
    if(children.empty()){
	nodes.push_back(this);
    }else{
	set<Node*>::const_iterator iter = children.begin();
	set<Node*>::const_iterator iterEnd = children.end();
	for(;iter!=iterEnd; ++iter)
	    (*(*iter)).getLeafs(nodes);
    }
};

Node* Node::resetroot(const string& str){
    vector<Node*> leafs;
    getLeafs(leafs);
    Node* og = NULL;
    foreach(Node* nd, leafs)
	if((*nd).name == str){
	    og = nd; break;
	}

    if(og == NULL)
	return og;
    return resetroot(og);
}

Node* Node::resetroot(Node* np){
    vector<Node*> nlist;
    while((*np).parent != NULL){
	nlist.push_back((*np).parent);
	np = (*np).parent;
    }

    for(vector<Node*>::reverse_iterator iter = nlist.rbegin(); 
	iter != nlist.rend()-1; ++iter){
	vector<Node*>::reverse_iterator next = iter+1;
	(**iter).length = (**next).length;
	(**next).length = NAN;
	(**iter).children.erase(*next);
	(**next).parent = NULL;
	(**next).addChild(**iter);
    }

    return *(nlist.begin());
};

void Node::_outnwk(ostream& os){
    if(!children.empty()){
	os << "(";
	set<Node*>::const_iterator iter = children.begin();
	set<Node*>::const_iterator iterEnd = children.end();
	(*(*iter))._outnwk(os);
	for(++iter; iter!=iterEnd; ++iter){
	    os << ",";
	    (*(*iter))._outnwk(os);
	}
	os << ")";
    }

    if(name.find(' ')==std::string::npos)
	os << name ;
    else
	os << '"' << name << '"';

    if(!isnan(length))
	os << ":" << fixed << setprecision(5) << length;
};

void Node::outnwk(ostream& os){
    _outnwk(os);
    os << ";" << endl;
};

void Node::_innwk(istream& is){
    Node* np = new Node;
    addChild(*np);
    string* sp = &((*np).name);
    string nstr;

    while(is.good()){
	char c = is.get();
	if(c == ','){
	    if(!nstr.empty()){
		(*np).length = str2double(nstr);
		nstr.clear();		
	    }
	    
	    np = new Node;
	    sp = &((*np).name);
	    addChild(*np);
	}else if(c == ')'){
	    if(!nstr.empty()){
		(*np).length = str2double(nstr);
		nstr.clear();
	    }
	    break;
	}else if(c == '('){
	    (*np)._innwk(is);
	}else if(c == ':'){
	    sp = &nstr;
	}else if(c!='"' && c!='\n' && c!='\t' && c!='\'') {
	    (*sp).push_back(c);
	}
    }
};

void Node::innwk(istream& is){

    char c = is.get();
    while(c != ';'){
	if(c =='('){
	    _innwk(is);
	}else if(c!='"' && c!='\n' && c!='\t' && c!='\'') {
	    name.push_back(c);
	}
	
	if(is.good()){
	    c = is.get();
	}else{
	    cerr << "some wrong in the nwk file" << endl;
	    exit(1);
	}
    }
};

size_t Node::checkTodefine(){
    vector<Node*> leafs;
    getLeafs(leafs);
    for(vector<Node*>::iterator it=leafs.begin(); it != leafs.end(); ++it)
	if((*(*it)).name.find("Unclassified<")!=string::npos){
	    (*(*it)).todefine = true;
	    ++nxleaf;
	}
    return nxleaf;
};

void Node::getDefineLeafs(vector<Node*>& nodes){
    if(children.empty()){
	if(!todefine)
	    nodes.push_back(this);
    }else{
	set<Node*>::const_iterator iter = children.begin();
	set<Node*>::const_iterator iterEnd = children.end();
	for(;iter!=iterEnd; ++iter)
	    (*(*iter)).getDefineLeafs(nodes);
    }
};

void Node::outjson(ostream& os){

    os << "{";
    if(!name.empty()){
	os << '"' << "name" << "\":\"" << name << "{" << nleaf;
	if(nleaf < taxSize)
	    os << "/" << taxSize;
	os <<"}\"";
    }

    if(todefine)
	os << ',' << '"' << "upload" << "\":\"" << "true" <<'"';

    if(nleaf == taxSize)
	os << ',' << '"' << "ntype" << "\":\"" << "Coincide" <<'"';

    if(!isnan(length))
	os << ',' << '"' << "length" << "\":\"" << fixed << setprecision(5) << length <<'"';

    if(!children.empty()){
	os << ",\"children\":[";
	
	set<Node*>::const_iterator iter = children.begin();
	set<Node*>::const_iterator iterEnd = children.end();
	(*(*iter)).outjson(os);
	for(++iter; iter!=iterEnd; ++iter){
	    os << ",";
	    (*(*iter)).outjson(os);
	}
	os << "]";
    }
    os << "}";
};

void Node::outjsonAbbr(ostream& os){

    os << "{";
    if(!name.empty()){
	os << '"' << "n" << "\":\"" << name << "{" << nleaf;
	if(nleaf < taxSize)
	    os << "/" << taxSize;
	os <<"}\"";
    }

    if(todefine)
	os << ',' << '"' << "u" << "\":\"" << "1" <<'"';

    if(nleaf == taxSize)
	os << ',' << '"' << "t" << "\":\"" << "C" <<'"';

    // if(!isnan(length))
    // 	os << ',' << '"' << "l" << "\":\"" << fixed << setprecision(5) << length <<'"';

    if(!children.empty()){
	os << ",\"c\":[";
	
	set<Node*>::const_iterator iter = children.begin();
	set<Node*>::const_iterator iterEnd = children.end();
	(*(*iter)).outjsonAbbr(os);
	for(++iter; iter!=iterEnd; ++iter){
	    os << ",";
	    (*(*iter)).outjsonAbbr(os);
	}
	os << "]";
    }
    os << "}";
};

void Node::renewId(const tr1::unordered_map<string,size_t>& mgi){
    if(children.empty()){
	id = mgi.find(name)->second;
    }else{
	set<Node*>::iterator iter = children.begin();
	set<Node*>::iterator iterEnd = children.end();
	for(;iter!=iterEnd; ++iter)
	    (*(*iter)).renewId(mgi);
    }
};

void Node::outPrediction(ostream& os){
    vector<Node*> leafs;
    getLeafs(leafs);
    for(vector<Node*>::iterator it=leafs.begin(); it != leafs.end(); ++it)
	if((*(*it)).todefine){
	    string nm = (*(*it)).name;
	    string p; (*(*it))._getPrediction(p);
	    os << nm.substr(nm.find_last_of('>')+1) << "\t" << p << endl;
	}
};

void Node::_getPrediction(string& p){
    if(parent == NULL){
	p = name;
    }else{
	if((*parent).nleaf == (*parent).taxSize){
	    p = (*parent).name;
	}else{
	    (*parent)._getPrediction(p);
	}
    }
};

