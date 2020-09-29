
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <map>
using namespace std;

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "morphfuncs.h"


map<char *, vector<morph_sample *> *> mpcache;

char *make_key(int ct, int cn)
{
	char *buff = new char[24];
	sprintf(buff, "%d_%d", ct, cn);
	return buff;
}


// Recursive helper function to calculate total surface area
//  distal to a given cable.

void calc_dist_surf_area(int celltype, int cellnum, cable *start, double *totsurf)
{
	node *n1 = nd(start->node1a, start->node1b, start->node1c);
	node *n2 = nd(start->node2a, start->node2b, start->node2c);

	double len = sqrt((n2->xloc - n1->xloc)*(n2->xloc - n1->xloc) + 
			(n2->yloc - n1->yloc)*(n2->yloc - n1->yloc) + 
			(n2->zloc - n1->zloc)*(n2->zloc - n1->zloc));

	*totsurf += PI*start->dia*len;
	//find cables whose proximal nodes are this dendrite's distal node
	vector<cable *> *dends = find_cables_by_endpoint(celltype, cellnum, -1, start->node1c);
	for (int k = 0; k < dends->size(); k++) {
		calc_dist_surf_area(celltype, cellnum, (*dends)[k], totsurf);
	}
	delete dends;
}


// A recursive helper function for get_morph_props, which iterates
// distally from a given start cable and computes parameters such as
// distance from soma, and distal surface area.

void iterate_tree(int celltype, int cellnum, cable *start,
		double dist, double edist, double cdm, int order, int system,
		vector<morph_sample *> *props)
{
	node *n1, *n2;
	double mismatch, len, lam, sa, totsa;
	n1 = nd(start->node1a, start->node1b, start->node1c);
	n2 = nd(start->node2a, start->node2b, start->node2c);

	len = sqrt((n2->xloc - n1->xloc)*(n2->xloc - n1->xloc) + 
			(n2->yloc - n1->yloc)*(n2->yloc - n1->yloc) + 
			(n2->zloc - n1->zloc)*(n2->zloc - n1->zloc));
	dist += len;
	sa = PI*start->dia*len;

	//space_const = sqrt( (d/4) * (Rm/Ri) ) Biophysics of Computation, Ch. 2
	lam = sqrt((2500*drm*start->dia)/dri);

	edist += len / lam;

	totsa = 0;
	calc_dist_surf_area(celltype, cellnum, start, &totsa);
	totsa -= sa;

	(*props)[start->node1c] = new morph_sample;
	(*props)[start->node1c]->node2soma_dist = dist;
	(*props)[start->node1c]->cdm = cdm;
	(*props)[start->node1c]->dist_surf_area = totsa;

	if (strcmp(start->elabl, regname[DEND])  == 0 ||
	    strcmp(start->elabl, regname[DENDP]) == 0) {
		(*props)[start->node1c]->system = system;
	} else {
		(*props)[start->node1c]->system = 0;  
	}

	//(*props)[start->node1c]->system = system;
	(*props)[start->node1c]->order = order;
	(*props)[start->node1c]->dia = start->dia;
	(*props)[start->node1c]->surf_area = sa;
	(*props)[start->node1c]->space_const = lam;
	(*props)[start->node1c]->node2soma_edist = edist;

	//find cables whose proximal nodes are this dendrite's distal node
	vector<cable *> *dends = find_cables_by_endpoint(celltype, cellnum, -1, start->node1c);
	if (dends->size() > 1) order++;
	for (int k = 0; k < dends->size(); k++) {
		mismatch = abs(start->dia - (*dends)[k]->dia);
		iterate_tree(celltype, cellnum, (*dends)[k], dist, edist, cdm+mismatch, order, system, props);
	}
	delete dends;
}


vector <morph_sample*> *get_morph_props(int celltype, int cellnum)
{
	int ct, cn, n, numnodes = 0;
	node *npnt;

	//if this function has been called for the same cell before,
	//returned the cached result, otherwise compute the distances
	//and store in cache
	char *key = make_key(celltype, cellnum);
	vector<morph_sample *> *samps = mpcache[key];
	delete key;

	if (samps == NULL) {
		//count the number of nodes
		for (npnt=nodepnt; npnt=foreach(npnt, celltype, cellnum, -1, &ct, &cn, &n); npnt=npnt->next) numnodes++;

		if (numnodes == 0) return NULL;

		sphere *somaElem = (sphere *) foreach(elempnt, SPHERE, celltype, cellnum, soma);
		double somarad = somaElem->dia / 2;

		//initialize vector to hold distances
		samps = new vector<morph_sample *>(numnodes);
		(*samps)[0] = new morph_sample;
		(*samps)[0]->node2soma_dist = 0.0;     
		(*samps)[0]->node2soma_edist = 0.0;     
		(*samps)[0]->cdm = 0;
		(*samps)[0]->order = 0;
		(*samps)[0]->system = 0;
		(*samps)[0]->dia = somaElem->dia;
		(*samps)[0]->surf_area = 4*PI*somarad*somarad;
		(*samps)[0]->space_const = sqrt((2500*drm*somaElem->dia)/dri); 

		//find primary cables (hillock and primary dendrites)
		vector<cable *> *pcbls = find_cables_by_endpoint(celltype, cellnum, -1, soma);

		int syscnt = 0;
		int sys;

		//recursively iterate distally from soma to compute properties
		for (int k = 0; k < pcbls->size(); k++) {

			if ( strcmp((*pcbls)[k]->elabl, regname[DEND])  == 0 ||
			     strcmp((*pcbls)[k]->elabl, regname[DENDP]) == 0) {
				sys = ++syscnt;
			} else {
				sys = 0;
			}      

			iterate_tree(celltype, cellnum, (*pcbls)[k],
					somarad, 0, 0, 1, sys, samps);
		}
		delete pcbls;

		//store properties in cache
		char *key = make_key(celltype, cellnum);
		mpcache[key] = samps;
		delete key;
	}

	return samps;
}


int is_branch_point(int celltype, int cellnum, int n)
{
	int ret = 0;
	//find cables whose proximal node is n
	vector<cable *> *cbls = find_cables_by_endpoint(celltype, cellnum, -1, n);
	if (cbls->size() > 1) ret = 1;
	delete cbls;
	return ret;
}


vector <vector<double>*> *get_node2node_dists(int celltype, int cellnum)
{
	int ct, cn, n, numnodes = 0;
	int k, j, p, q;
	node *npnt;

	//count nodes
	for (npnt=nodepnt; npnt=foreach(npnt, celltype, cellnum, -1, &ct, &cn, &n); npnt=npnt->next) numnodes++;

	//initialize vector of node2node distances to all -1
	vector< vector<double> *> *node2node_dists = new vector< vector<double> *>(numnodes);
	for (k = 0; k < numnodes; k++) {
		(*node2node_dists)[k] = new vector<double>(numnodes);
		for (j = 0; j < numnodes; j++) (*(*node2node_dists)[k])[j] = -1;
	}

	double dist = 0;
	vector<int> *bps1 = new vector<int>;
	vector<int> *bps2 = new vector<int>;
	//compute uncalculated distances for each node
	for (k = 0; k < numnodes; k++) {

		dist = 0;
		bps1->clear();

		double len = 0;
		vector<cable *> *cbls;    
		node *pnode, *dnode;
		int nextnode = k;
		if ( (*(*node2node_dists)[k])[k] == -1)    (*(*node2node_dists)[k])[k] = 0;

		//travel from node k to the soma, recording distances and
		//branch points along the way
		while (nextnode != soma) {

			if (is_branch_point(celltype, cellnum, nextnode)) bps1->push_back(nextnode);

			//find cable whose distal node is "nextnode"
			cbls = find_cables_by_endpoint(celltype, cellnum, nextnode, -1);
			if (cbls->size() > 0) {	

				dnode = nd((*cbls)[0]->node1a, (*cbls)[0]->node1b, (*cbls)[0]->node1c);
				pnode = nd((*cbls)[0]->node2a, (*cbls)[0]->node2b, (*cbls)[0]->node2c);

				len = sqrt((pnode->xloc - dnode->xloc)*(pnode->xloc - dnode->xloc) + 
						(pnode->yloc - dnode->yloc)*(pnode->yloc - dnode->yloc) + 
						(pnode->zloc - dnode->zloc)*(pnode->zloc - dnode->zloc));

				dist += len;

				nextnode = (*cbls)[0]->node2c;
				if ( (*(*node2node_dists)[k])[nextnode] == -1) {
					//fprintf(stderr, "# distance between %d and %d: %g\n", k, nextnode, dist);
					(*(*node2node_dists)[k])[nextnode] = dist;
					(*(*node2node_dists)[nextnode])[k] = dist;
				}

			} else {
				fprintf(stderr, "# iteration 1: could not find cable with distal node %d!\n", k);
			}
			bps1->push_back(soma);
			delete cbls;      
		}

		//compute node2node distances between k and all other nodes by
		//comparing branch points
		for (j = 0; j < numnodes; j++) {

			if ( (*(*node2node_dists)[k])[j] == -1) {

				//now travel from node j to the soma, recording distances and
				//branch points along the way, and reusing some variables
				bps2->clear();
				nextnode = j;
				dist = 0;
				len = 0;
				while (nextnode != soma) {

					if (is_branch_point(celltype, cellnum, nextnode)) bps2->push_back(nextnode);

					//find cable whose distal node is "nextnode"
					cbls = find_cables_by_endpoint(celltype, cellnum, nextnode, -1);
					if (cbls->size() > 0) {	

						dnode = nd((*cbls)[0]->node1a, (*cbls)[0]->node1b, (*cbls)[0]->node1c);
						pnode = nd((*cbls)[0]->node2a, (*cbls)[0]->node2b, (*cbls)[0]->node2c);

						len = sqrt((pnode->xloc - dnode->xloc)*(pnode->xloc - dnode->xloc) + 
								(pnode->yloc - dnode->yloc)*(pnode->yloc - dnode->yloc) + 
								(pnode->zloc - dnode->zloc)*(pnode->zloc - dnode->zloc));

						dist += len;

						nextnode = (*cbls)[0]->node2c;
						if ( (*(*node2node_dists)[j])[nextnode] == -1) {
							//fprintf(stderr, "# distance between %d and %d: %g\n", j, nextnode, dist);
							(*(*node2node_dists)[j])[nextnode] = dist;
							(*(*node2node_dists)[nextnode])[j] = dist;
						}

					} else {
						fprintf(stderr, "# iteration 2: could not find cable with distal node %d!\n", k);
					}
					bps2->push_back(soma);	  
					delete cbls;
				}

				//find common branch point and compute node2node distance
				for (p = 0; p < bps1->size(); p++) {
					for (q = 0; q < bps2->size(); q++) {
						if ((*bps2)[q] == (*bps1)[p]) {
							if ( (*(*node2node_dists)[k])[j] == -1) {
								dist = (*(*node2node_dists)[k])[(*bps2)[q]] + (*(*node2node_dists)[j])[(*bps2)[q]];
								(*(*node2node_dists)[k])[j] = dist;
								(*(*node2node_dists)[j])[k] = dist;
								//fprintf(stderr, "# distance between %d and %d: %g\n", k, j, dist);
							}
						}	    
					}
				}	
			}  
		}    

	}

	return node2node_dists;
}

vector<cable *> *find_cables_by_endpoint(int celltype, int cellnum, int distNode, int proxNode)
{
	vector<cable *> *cbls = new vector<cable *>;    
	elem *epnt;

	for (epnt=elempnt; epnt=foreach(epnt, CABLE, celltype, cellnum); epnt=epnt->next) {
		if ((distNode != -1) && (proxNode == -1)
				&& (epnt->node1c == distNode))  cbls->push_back((cable *) epnt);
		if ((distNode == -1) && (proxNode != -1)
				&& (epnt->node2c == proxNode))  cbls->push_back((cable *) epnt);    
		if ((distNode != -1) && (proxNode != -1)
				&& (epnt->node1c == distNode) && (epnt->node2c == proxNode)) cbls->push_back((cable *) epnt);
	}

	return cbls;
}
