/* Experiment cell_vclamp */
/*  for nc script retsim.cc */

#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "ncio.h"
#include "stimfuncs.h"

#define PRINTER(name) printer((char *)#name, (name))
void printer(char *name, double value) { fprintf(stderr, "%s = %G\t", name, value); }
void printer(char *name, int value)    { fprintf(stderr, "%s = %i\t", name, value); }

// Input variables for stimulus.
const char *stimfile;
const char *stimtype;

// Internal variables for stimulus.
std::vector<double> stim_time;
std::vector<double> stim;

// BC parameters.
int is_OFF_bp;
int n_c_bp;
int bp_type_idx;
int gc_type_idx;

// Local parameters.
declare_params_here

// Parameters that can be passed in python.
double set_ploti;
double set_timinc;
double set_stiminc;
double set_srtimestep;
double set_dscaeg; 

int set_pigm;
double set_cone_maxcond;
double cone_soma_z;
const char *rectype;
		
int print_comps;
const char *compfile;

double light_min = 35e4;
double light_rng = 70e4;

/*------------------------------------------------------*/
void defparams(void) 
{ 
  // Stimulus parameters.
  setptr("is_OFF_bp",	     &is_OFF_bp);
  setptr("n_c_bp",         &n_c_bp);          n_c_bp = 1;
  
  setptr("stimfile",	     &stimfile);
  setptr("stimtype",	     &stimtype);      stimtype = "clamp";
  
  setptr("set_ploti",      &set_ploti);
  setptr("set_timinc",     &set_timinc);
  setptr("set_stiminc",    &set_stiminc);
  setptr("set_srtimestep", &set_srtimestep);
  setptr("set_dscaeg",     &set_dscaeg);
  
  setptr("rectype",        &rectype);
  setptr("cone_soma_z",    &cone_soma_z);
  setptr("print_comps",    &print_comps);
  setptr("compfile",       &compfile);  
  
  // Local parameters.
  init_params_here
}

/*------------------------------------------------------*/
void setparams(void)
{ 
  SOMA = R_3;	
  
  // Make cone type.
  make_ct(xcone);
  setn(xcone,MORPH,5); // see makcel.cc
  setn(xcone,NCOLOR,RCOLOR);
  if (!notinit(cone_soma_z)) setn(xcone,SOMAZ,cone_soma_z);
  
  // Find BC type.
  if (is_OFF_bp) {
    fprintf(stderr, "# Simulate OFF BC\n");
    bp_type_idx = hbp1;
    gc_type_idx = gcaoff;
  } else {
    fprintf(stderr, "# Simulate ON BC\n");
    bp_type_idx = dbp1;
    gc_type_idx = gca;
  }
  
  // Make BC type.
  make_ct(bp_type_idx);
  set_ncel(bp_type_idx,1); // Set number of cells to one.
  setn(bp_type_idx,MORPH,0);
  setn(bp_type_idx,NCOLOR,RCOLOR);
  setn(bp_type_idx,AXARBT,BRANCHED);
  
  // Don't make connections yet.
  make_cone_dbp1 = 0;
  make_cone_hbp1 = 0;
 
  // Create cells.
  checkcellout(xcone);
  fprintf(stderr, "# Make cone\n");
  #define CONEARR 30 
  conexarr  = (double *)emalloc(CONEARR*sizeof(double));
  coneyarr  = (double *)emalloc(CONEARR*sizeof(double));
  conetharr = (double *)emalloc(CONEARR*sizeof(double));
  conenarr  = (int *) emalloc(CONEARR*sizeof(int));
  conexarr[0] = 0;    coneyarr[0] = 0;    conetharr[0] = 0;    conenarr[0] = 1;
  n_cones = 1;
  
  if (bp_type_idx == dbp1) {
    checkcellout(dbp1);
    fprintf(stderr, "# Make dbp1\n");
    #define DBP1ARR 30
    dbp1xarr  = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1yarr  = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1tharr = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1narr  = (int *) emalloc(DBP1ARR*sizeof(int));
    dbp1xarr[0] = 0;    dbp1yarr[0] = 0;    dbp1tharr[0] = 0;    dbp1narr[0] = 1;
    n_dbp1 = 1;
  } else {
    checkcellout(hbp1);
    fprintf(stderr, "# Make hbp1\n");
    #define HBP1ARR 30
    hbp1xarr  = (double *)emalloc(HBP1ARR*sizeof(double));
    hbp1yarr  = (double *)emalloc(HBP1ARR*sizeof(double));
    hbp1tharr = (double *)emalloc(HBP1ARR*sizeof(double));
    hbp1narr  = (int *) emalloc(HBP1ARR*sizeof(int));
    hbp1xarr[0] =  0;    hbp1yarr[0] = 0;    hbp1tharr[0] = 0;    hbp1narr[0] = 1;
    n_hbp1 = 1;
  }
 
  // Print compartment file.
  if (notinit(print_comps)) print_comps = 0;
  if (print_comps) {
    fprintf(stderr, "# print_comps to %s\n", compfile);
    if (print_comps) comp_file(compfile);
  }
 
  // Set general parameters.
  vnoise   = 0;  		// =1 -> vesicle noise
  pnoise   = 0;  		// =1 -> photon noise
  dnoise   = 0;  		// =1 -> continuous dark noise
  Chnoise  = 0;  		// =1 -> channel noise
  
  // Calcium parameters.
  dscavg   = 1E+06; // Calcium sensitivity @ vesicle release  
  dicafrac = 1;
  if (notinit(set_dscaeg)) set_dscaeg = 1;
  dscaeg = set_dscaeg;
  
  // Cell temperature.
  tempcel  = 37;
  
  // Ion concentrations from Franke 2017.
  dnao = 151.5e-3;
  dko  = 2.5e-3;  
  dclo = 133.5e-3; 
  dcao = 2e-3;

  use_ghki = 1;

  // Reversal potentials, ions.
  vna = +65e-3;
  vk  = -89e-3;
  vcl = -70e-3;
  
  // Cone parameters.
  cone_timec   = 0.2;		
  cone_loopg   = 0.0;
  cone_pigm    = 14;
  cone_maxcond = 2e-10;
}

/*------------------------------------------------------*/
void addcells() 
{
  fprintf(stderr, "# Add GCs.\n");
  make_gc_comps(bp_type_idx, 1, R8, gc_type_idx);

  // Connect cones.
  int contact_nodes_dbp1 [] = {};
  int contact_nodes_hbp1 [] = {};
  
  fprintf(stderr, "# Add cone connections ");
  for(int i=0;i<sizeof(contact_nodes_dbp1)/sizeof(int);i++) {
    if (i==0) fprintf(stderr, "ON BC: ");
    fprintf(stderr, "%d (%d), ", i, n_c_bp);
    for (int j=0; j<n_c_bp; j++) {
      connect_synapse(xcone,1,1,dbp1, 1,contact_nodes_dbp1[i]);
    }
  }
  for(int i=0;i<sizeof(contact_nodes_hbp1)/sizeof(int);i++) {
    if (i==0) fprintf(stderr, "OFF BC: ");
    fprintf(stderr, "%d (%d), ", i, n_c_bp);
    for (int j=0; j<n_c_bp; j++) {
      connect_synapse(xcone,1,1,hbp1, 1,contact_nodes_hbp1[i]);
    }
  }
  fprintf(stderr, "\n");
}

/*------------------------------------------------------*/
void runexpt(void)
{ 
  elem *epnt, *s1, *s2;
  synapse *cur_syn;
  double Vmin, Vmax, fmin, fmax, Imin, Imax, maxca;
  int colr, cell_n, node_n;
  
  // Initialize if necessary. 
  timinc     = (notinit(set_timinc))     ? 2e-5 : set_timinc;
  stiminc    = (notinit(set_stiminc))    ? 1e-4 : set_stiminc;
  srtimestep = (notinit(set_srtimestep)) ? 1e-4 : set_srtimestep;
  ploti      = (notinit(set_ploti))      ? 1e-3 : set_ploti;  
  if (notinit(setxmin)) setxmin = 0;
  
  // ////////////////////////////////////////////////
  // Initialize plotting.
  // ////////////////////////////////////////////////
  
  int plot_light = 0;
  int plot_cone = 0;
  int plot_cone_cascade = 0;
  int plot_bc_vm_soma = 0;
  int plot_bc_i_soma = 0;
  int plot_bc_vm = 0;
  int plot_bc_ca = 0;
  int plot_bc_rate = 0;
  
  plot_light = 1;
  if (strcmp(rectype, "all") == 0) {
    fprintf(stderr, "# Rec-Type: all\n");
    plot_cone = 1;
    plot_cone_cascade = 1;
    plot_bc_vm_soma = 1;
    plot_bc_i_soma = 1;
    plot_bc_vm = 1;
    plot_bc_ca = 1;
    plot_bc_rate = 1;
  } else if (strcmp(rectype, "test") == 0) {
    fprintf(stderr, "# Rec-Type: test\n");
    plot_cone = 1;
    plot_bc_vm_soma = 1;
    plot_bc_i_soma = 1;
    plot_bc_vm = 1;
    plot_bc_rate = 1;
    plot_bc_ca = 1;
  } else if (strcmp(rectype, "optimize") == 0) {
    fprintf(stderr, "# Rec-Type: optimize\n");
    plot_bc_vm_soma = 1;
    plot_bc_rate = 1;
  } else if ((strcmp(rectype, "cones") == 0) || strcmp(rectype, "Cones") == 0) {
    fprintf(stderr, "# Rec-Type: cones\n");
    plot_cone = 1;
  } else {
    fprintf(stderr, "# No Rec-Type given!\n");
  }

  // ////////////////////////////////////////////////
  // Plot.
  // ////////////////////////////////////////////////
  
  if (plot_light) {
    // Plot absorbed light at central cones.
    plot_l_nod(xcone, cell_n=1, node_n=0, Vmin = 0, Vmax = 10000, colr=white, "Light", 22, 0.5);
  }
  
  if (plot_cone) {
    // Plot cone.
    plot_vm_nod(xcone, cell_n=1, node_n=0, Vmin = -1, Vmax = 1, colr=blue, "Cone Vm Soma", 21, 1.0);
    plot_vm_nod(xcone, cell_n=1, node_n=1, Vmin = -1, Vmax = 1, colr=blue, "Cone Vm Axon Terminal", 21, 1.0);
    plot_vm_nod(xcone, cell_n=1, node_n=2, Vmin = -1, Vmax = 1, colr=blue, "Cone Vm Axon", 21, 1.0);
    
    // Plot cone glutamate release.
    plot_synrate(findsyn(xcone, cell_n=1, node_n=1, bp_type_idx), fmin=0, fmax=400,  cyan,  20, "Cone",0.3);
  }

  if (plot_bc_vm_soma) {
    // Plot Vm of BC soma.
    plot_vm_nod(bp_type_idx, cell_n=1, soma, Vmin = -1, Vmax = 1, colr=blue, "BC Vm Soma", 21, 1.0);
  }
  if (plot_bc_i_soma) {
    plot_i_nod(bp_type_idx, cell_n=1, soma, -5e-9, 5e-9, colr=blue, "Isoma", 21, 1.0);
  }
  if (plot_bc_vm) {
    // Plot Vm of several nodes of BC.
    int rec_nodes []  = {10, 107, 152, 200, 300, 365, 583, 601, 645, 693, 804, 1060, 1079, 1039};
    
    for (int node_n=0 ; node_n<sizeof(rec_nodes)/sizeof(int); node_n++) {
      plot_vm_nod(bp_type_idx, cell_n=1, rec_nodes[node_n], Vmin = -1, Vmax = 1, colr=blue, "BC Vm", 21, 1.0);
    }
  }
  
  if (plot_bc_ca) {
    for (node_n = 0; node_n < 5000; node_n++) {
      cur_syn = findsyn(bp_type_idx, cell_n=1, node_n, gc_type_idx);
      if (cur_syn != NULL) {
        plot_ca_syn(cur_syn,  1.0e-6, blue, "BC Ca", 6,0.5);
      }
    }
  }
  
  if (plot_bc_rate) {
    // Plot BC glutamate release.
    for (node_n = 0; node_n < 5000; node_n++) {
      cur_syn = findsyn(bp_type_idx, cell_n=1, node_n, gc_type_idx);
      if (cur_syn != NULL) {
        plot_synrate(cur_syn, fmin=0,fmax=400, colr=blue, 8, "BC", 0.5);
      }
    }
  }
  
  if (plot_cone_cascade == 1) {
    // Plot details of cone cascade.
    for (epnt=elempnt; epnt=foreach(epnt, CONE); epnt=epnt->next) {
      s1 = epnt;
    }
    plot(G, 0, s1->elnum);
    
    for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE); epnt=epnt->next) {
      s2 = epnt;
    }
    for (int ii = 1; ii <= 10; ii++) {
       // Plot isomerized rhodopsin
      plot(G, ii, s2->elnum);
    }
  }
  
  // ////////////////////////////////////////////////
  // Get stimulus.
  // ////////////////////////////////////////////////
  
  std::string line;
  std::ifstream stim_file(stimfile);
  
  if (stim_file.is_open()) {
    while (getline(stim_file,line)) {
      int idx_comma = line.find(',');
      stim_time.push_back(atof(line.substr(0, idx_comma).c_str()));
      stim.push_back(atof(line.substr(idx_comma+1, line.size()-(idx_comma+1)).c_str()));
    }
    stim_file.close();
  }
  
  // ////////////////////////////////////////////////
  // Run.
  // ////////////////////////////////////////////////
  
  fprintf(stderr, "# Run stimulus ... \n");
  
  double V0    =  30e-3;
  double V0_dt1 = 0.01;
  double V0_dt2 = 0.1;
  double V0_dt3 = 0.3;
  
  PRINTER(bp_vst);
  PRINTER(V0);
  PRINTER(V0_dt1);
  PRINTER(V0_dt2);
  PRINTER(V0_dt3);
  fprintf(stderr, "\n");
  
  fprintf(stderr, "# Run step 1 ... ");
  step(V0_dt1);
  fprintf(stderr, "Done!\n");
  for (int i = 0; i < 1400; i++) {
    if (nde(bp_type_idx,1,i) != NULL) {
      if (strcmp(stimtype, "ramp") == 0) {
        ramp_v (ndn(bp_type_idx,1,i), bp_vst, V0, simtime, V0_dt2, stiminc);
      } else if (strcmp(stimtype, "tworamps") == 0) {
        ramp_v (ndn(bp_type_idx,1,i), bp_vst, V0, simtime, V0_dt2, stiminc);
        step(V0_dt2*1.1);
        ramp_v (ndn(bp_type_idx,1,i), bp_vst, V0, simtime, V0_dt2, stiminc);
      } else {
        vclamp (ndn(bp_type_idx,1,i), V0, simtime, V0_dt2);
      }
    }
  }
  step(V0_dt2 + V0_dt3);
  
  fprintf(stderr, "# Successfully terminated!\n");
}
