#include "rectask.h"

#include <cstring>
using namespace std;


void init_func_pnt_maps()
{
  nfunclib[string("v")] = &v;  //voltage at a node
  nfunclib[string("i")] = &i;  //current at a node

  func_pnt_maps_flag = 1;
}


rectask::rectask(const char *filename, char delimiter)
{
  fname = filename;
  delim = delimiter;

  fp = fopen(fname, "w");

  if (!func_pnt_maps_flag) init_func_pnt_maps();
}

void rectask::add(node *n, const char *prop)
{
  recpair rp;
  rp.isnode = 1;
  rp.n = n;

  string sprop(prop);

  map<string, node_func_ptr>::iterator iter = nfunclib.find(sprop);
  if(iter == nfunclib.end()) {
    fprintf(stderr, "[rectask::add()] Function for property '%s' not found!\n", prop);
    return;
  } else {
    rp.nfunc = nfunclib[sprop];      
    rpairs.push_back(rp);
  }
}

void rectask::add(elem *e, const char *prop)
{
  recpair rp;
  rp.isnode = 0;
  rp.e = e;

  string sprop(prop);

  map<string, elem_func_ptr>::iterator iter = efunclib.find(sprop);
  if(iter == efunclib.end()) {
    fprintf(stderr, "[rectask::add()] Function for property '%s' not found!\n", prop);
    return;
  } else {
    rp.efunc = efunclib[sprop];      
    rpairs.push_back(rp);
  }
}

void rectask::add(node *n, node_func_ptr nfunc)
{
  recpair rp;
  rp.isnode = 1;
  rp.n = n;
  rp.nfunc = nfunc;    

  rpairs.push_back(rp);
}

void rectask::add(elem *e, elem_func_ptr efunc)
{
  recpair rp;
  rp.isnode = 0;
  rp.e = e;
  rp.efunc = efunc;    

  rpairs.push_back(rp);
}

void rectask::execute()
{
  int k;
  double val;
  for (k = 0; k < rpairs.size(); k++) {
    
    if (rpairs[k].isnode)  val = rpairs[k].nfunc(rpairs[k].n);
    else                   val = rpairs[k].efunc(rpairs[k].e);    

    if (k != (rpairs.size()-1)) fprintf(fp, "%g%c", val, delim);  
    else                        fprintf(fp, "%g", val);
  }
  fprintf(fp, "\n");
}
