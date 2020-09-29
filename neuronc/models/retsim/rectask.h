#include "ncfuncs.h"
#include "scheduler.h"

#include <map>
#include <vector>
#include <string>
using namespace std;

extern "C" {
#include <stdio.h>
}

//type definitions for function pointers
typedef double (*node_func_ptr)(node *);
typedef double (*elem_func_ptr)(elem *);

/*
  Holds information about a node or element
  with an associated "recording function", i.e.
  a function of a node* or elem* which computes
  returns some scalar property.
*/
struct recpair {
  int isnode;

  node *n;
  node_func_ptr nfunc;

  elem *e;
  elem_func_ptr efunc;
};

//string to function library for nodes (global)
map<string, node_func_ptr> nfunclib;

//string to function library for elements (global)
map<string, elem_func_ptr> efunclib;

//initialize function pointer maps (global), called
//by constructor for rectask if not the maps are
//not initialized
void init_func_pnt_maps();

//flag for whether or not the function pointer maps
//have been initialized (global)
int func_pnt_maps_flag = 0;


/*
  A task designed to be called by the scheduler
  which writes properties of specified nodes and/or
  elements to a file, in columnar format (one row
  per time point).
*/
class rectask : public task
{
 protected:

  //delimiter
  char delim;

  //output file name  
  const char *fname;

  //file being written to
  FILE *fp;

  //nodes/elements to record from
  vector<recpair> rpairs;  

 public:

  rectask(const char *filename, char delimiter);

  //record property "prop" from node n
  void add(node *n, const char *prop);

  //record property computed by "nfunc" from node n
  void add(node *n, node_func_ptr nfunc);

  //record property "prop" from element e
  void add(elem *e, const char *prop);

  //record property computed by "efunc" from element e
  void add(elem *e, elem_func_ptr efunc);

  //write recorded properties to a file
  void execute();
};
