#include <vector>
using namespace std;

struct morph_sample
{
  double cdm;              //cumulative diameter mismatch
  double node2soma_dist;   //distance to soma
  double node2soma_edist;  //electrotonic distance to soma
  int    order;            //somafugal order number
  int    system;           //dendritic system
  double dia;              
  double surf_area;        
  double space_const;      //computed space constant
  double dist_surf_area;   //distal dendritic surface area
};

/*
  Get a list of morphological properties sampled from the cell
  specified. Element j of the vector is a set of properties
  measured for node j, contained in a morph_sample struct.
*/
vector <morph_sample*> *get_morph_props(int celltype, int cellnum);


/*
  Locate cable(s) with endpoint nodes "distNode" and "proxNode. "distNode" is the
  distal node, and "proxNode" is the proximal node, either of which can be -1 (unspecified).
  If both nodes are specified, the vector returned should only contain one cable, if
  no nodes are found the returned vector will be empty. Returns the element # of the cable.
*/
vector <cable*> *find_cables_by_endpoint(int celltype, int cellnum, int distNode, int proxNode);


/*
  Compute the node-to-node distance for each node in the
  dendritic tree of the specified cell. This is a very
  expensive function and takes a long time to run. Returns a
  two-dimensional vector where element ij is the distance in
  microns between node i and node j.
*/
vector <vector<double>*> *get_node2node_dists(int celltype, int cellnum);
