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
bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }

// Define what to record.
const char *rectype;

// Input variables for stimulus.
const char *stimtype;
const char *stimfile;

std::vector<double> light_stim;
std::vector<double> light_stim_time;
double light_min = 35e4;
double light_rng = 70e4;

int ex_comp_node_offset = 5000;
double stimdur;

// Compartment variables.
std::vector<int> comps;
std::vector<int> comps_ct;
std::vector<int> comps_cn;
std::vector<int> comps_nd;

// BC parameters.
int make_dbp1;
int make_hbp1;
double cone_soma_z;
int merge_cones;
int make_cones;
int n_c_bp;
double set_dscaeg; 
double set_tempcel; 
double set_c_vst;

// Equilibrium state.
int save_eq;   // Save equilibrium state?
int load_eq;   // Load equilibrium state?
int run_predur_only;
const char *eqfile;

// Parameters that can be passed in python.
double predur;
double rec_predur;
double pre_load_predur;
double post_load_predur;
double set_ploti;
double set_timinc;
double predur_timinc = 1e-4;
double set_stiminc;
double set_srtimestep;
double post_load_predur_timinc;

int print_comps;
const char *compfile;

// Local parameters.

// Connect cones.
int contact_nodes_dbp1 [] = {};
int contact_nodes_hbp1 [] = {};

// Cone positions
double cone_x_dbp1 [] = {};
double cone_y_dbp1 [] = {};

double cone_x_hbp1 [] = {};
double cone_y_hbp1 [] = {};

// BC positions
double dbp1_x = 0;
double dbp1_y = 0;
double hbp1_x = 0;
double hbp1_y = 0;

/*------------------------------------------------------*/
void defparams(void) 
{ 
  setptr("make_dbp1",	     &make_dbp1);       make_dbp1 = 1;
  setptr("make_hbp1",	     &make_hbp1);       make_hbp1 = 0;
  setptr("n_c_bp",         &n_c_bp);          n_c_bp = 2;
  setptr("cone_soma_z",    &cone_soma_z);     cone_soma_z = 20;
  setptr("make_cones",     &make_cones);      make_cones = 1;
  setptr("merge_cones",    &merge_cones);     merge_cones = 1;
  
  // Stimulus and recording parameters.
  setptr("stimfile", &stimfile);
  
  setptr("rectype",  &rectype);
  setptr("stimtype", &stimtype);
  setptr("stimdur",  &stimdur);
  
  // Time parameters.
  setptr("predur",           &predur);
  setptr("rec_predur",       &rec_predur);
  setptr("pre_load_predur",  &pre_load_predur);
  setptr("post_load_predur", &post_load_predur);
  setptr("post_load_predur_timinc",  &post_load_predur_timinc);
  setptr("set_ploti",        &set_ploti);
  setptr("set_timinc",       &set_timinc);
  setptr("predur_timinc",    &predur_timinc);
  setptr("set_stiminc",      &set_stiminc);
  setptr("set_srtimestep",   &set_srtimestep);
  setptr("set_dscaeg",       &set_dscaeg);
  setptr("set_tempcel",      &set_tempcel);
  setptr("set_c_vst",        &set_c_vst);
  
  setptr("save_eq",          &save_eq);         save_eq = 0;
  setptr("load_eq",          &load_eq);         load_eq = 0;
  
  setptr("eqfile", &eqfile);
  
  setptr("run_predur_only",  &run_predur_only); run_predur_only = 0;
  
  // Compartment parameters.
  setptr("print_comps",    &print_comps); print_comps = 1;
  setptr("compfile",       &compfile); compfile = "temp_compfile.n";
  
  // Local parameters.
}

/*------------------------------------------------------*/
void setparams(void)
{ 
  SOMA = R_3;	
  
  if (make_cones) {
    make_ct(xcone);
    setn(xcone,MORPH,5); // see makcel.cc
    setn(xcone,NCOLOR,RCOLOR);
    if (!notinit(cone_soma_z)) setn(xcone,SOMAZ,cone_soma_z);
  }
  
  // Make BC types.
  if (make_dbp1) {
    make_ct(dbp1);
    set_ncel(dbp1,1); // Set number of cells to one.
    setn(dbp1,MORPH,0);
    setn(dbp1,NCOLOR,RCOLOR);
    setn(dbp1,AXARBT,BRANCHED);
  }
  
  if (make_hbp1) {
    make_ct(hbp1);
    set_ncel(hbp1,1); // Set number of cells to one.
    setn(hbp1,MORPH,0);
    setn(hbp1,NCOLOR,RCOLOR);
    setn(hbp1,AXARBT,BRANCHED);
  }
  
  if (notinit(eqfile)) {
    if (make_dbp1 and !make_hbp1) {
      eqfile = "dbp1.eq";
    } else if (make_hbp1 and !make_dbp1) {
      eqfile = "hbp1.eq";
    } else {
      eqfile = "temp.eq";
    }
  }
  
  // Don't make connections yet.
  make_cone_dbp1 = 0;
  make_cone_hbp1 = 0;
  
  // Create cones?
  if (make_cones) {
    checkcellout(xcone);
    fprintf(stderr, "# Make cones\n");
    #define CONEARR 30 
    conexarr  = (double *)emalloc(CONEARR*sizeof(double));
    coneyarr  = (double *)emalloc(CONEARR*sizeof(double));
    conetharr = (double *)emalloc(CONEARR*sizeof(double));
    conenarr  = (int *) emalloc(CONEARR*sizeof(int));
    
    if (merge_cones) {
      conexarr[0] = 0;    coneyarr[0] = 0;    conetharr[0] = 0;    conenarr[0] = 1;
      n_cones = 1;
    } else {
      int n_cones_on  = sizeof(contact_nodes_dbp1)/sizeof(int);
      int n_cones_off = sizeof(contact_nodes_hbp1)/sizeof(int);
      
      for(int i=0;i<n_cones_on;i++) {
        conexarr[i] = cone_x_dbp1[i];
        coneyarr[i] = cone_y_dbp1[i];
        conetharr[i] = 0;
        conenarr[i] = i+1;
      }
      for(int i=0;i<n_cones_off;i++) {
        conexarr[n_cones_on+i] = cone_x_hbp1[i];
        coneyarr[n_cones_on+i] = cone_y_hbp1[i];
        conetharr[n_cones_on+i] = 0;
        conenarr[n_cones_on+i] = n_cones_on+i+1;
      }
      n_cones = n_cones_on + n_cones_off;
    }
  }
  
  if (make_dbp1) {
    checkcellout(dbp1);
    fprintf(stderr, "# Make dbp1\n");
    #define DBP1ARR 30
    dbp1xarr  = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1yarr  = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1tharr = (double *)emalloc(DBP1ARR*sizeof(double));
    dbp1narr  = (int *) emalloc(DBP1ARR*sizeof(int));
    dbp1xarr[0] = dbp1_x;    dbp1yarr[0] = dbp1_y;    dbp1tharr[0] = 0;    dbp1narr[0] = 1;
    n_dbp1 = 1;
  }
  if (make_hbp1) {
    checkcellout(hbp1);
    fprintf(stderr, "# Make hbp1\n");
    #define HBP1ARR 30
    hbp1xarr  = (double *)emalloc(HBP1ARR*sizeof(double));
    hbp1yarr  = (double *)emalloc(HBP1ARR*sizeof(double));
    hbp1tharr = (double *)emalloc(HBP1ARR*sizeof(double));
    hbp1narr  = (int *) emalloc(HBP1ARR*sizeof(int));
    hbp1xarr[0] =  hbp1_x;    hbp1yarr[0] = hbp1_y;    hbp1tharr[0] = 0;    hbp1narr[0] = 1;
    n_hbp1 = 1;
  }
 
  // Print compartment file.
  if (print_comps) {
    fprintf(stderr, "# print_comps to %s\n", compfile);
    comp_file(compfile);
  }
 
  // Set general parameters.
  vnoise   = 0;  		// =1 -> vesicle noise
  pnoise   = 0;  		// =1 -> photon noise
  dnoise   = 0;  		// =1 -> continuous dark noise
  Chnoise  = 0;  		// =1 -> channel noise
  
  // Calcium parameters.
  dscavg   = 1E+06; // Calcium sensitivity @ vesicle release  
  dicafrac = 1;
  
  // Cell temperature.
  if (notinit(set_tempcel)) set_tempcel = 37;
  tempcel  = set_tempcel;
  
  // Ion concentrations from Franke 2017.
  dnao = 151.5e-3;
  dko  = 2.5e-3;  
  dclo = 133.5e-3; 
  dcao = 2e-3;
  if (notinit(set_dscaeg)) set_dscaeg = 2;
  dscaeg = set_dscaeg;

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
void add_cone_connections(void) {
  node *npnt;

  fprintf(stderr, "# Add cone connections\n");
  for(int i=0;i<sizeof(contact_nodes_dbp1)/sizeof(int);i++) {
    if (i==0) fprintf(stderr, "# DBP1: ");
    if (!merge_cones) {
      fprintf(stderr, "%d (%d), ", i, n_c_bp);
      for (int j=0; j<n_c_bp; j++) {
        connect_synapse(xcone,i+1,1,dbp1, 1,contact_nodes_dbp1[i]);
      }
    } else {
      npnt=ndn(dbp1,1,contact_nodes_dbp1[i]);
      fprintf(stderr, "%d (%d), ", i, n_c_bp);
      for (int j=0; j<n_c_bp; j++) {
        connect_synapse(xcone,1,1,dbp1,1,contact_nodes_dbp1[i]);
      }
    }
  }
  fprintf(stderr,"\n");
  for(int i=0;i<sizeof(contact_nodes_hbp1)/sizeof(int);i++) {
    if (i==0) fprintf(stderr, "# HBP1: ");
    if (!merge_cones) {
      fprintf(stderr, "%d (%d), ", i, n_c_bp);
      for (int j=0; j<n_c_bp; j++) {
        connect_synapse(xcone,sizeof(contact_nodes_dbp1)/sizeof(int)+i+1,1,hbp1,1,contact_nodes_hbp1[i]);
      }
    } else {
      npnt=ndn(hbp1,1,contact_nodes_hbp1[i]);
      fprintf(stderr, "%d (%d), ", i, n_c_bp);
      for (int j=0; j<n_c_bp; j++) {
        connect_synapse(xcone,1,1,hbp1,1,contact_nodes_hbp1[i]);
      }
    }
  }
  fprintf(stderr,"\n");
}

/*------------------------------------------------------*/
void addcells() 
{
  fprintf(stderr, "# Add GCs.\n");
  if (make_dbp1) {
    make_gc_comps(dbp1, 1, R8, gca);
  }
  if (make_hbp1) {
    make_gc_comps(hbp1, 1, R8, gcaoff);
  }
  
  if (make_cones) add_cone_connections();
}

/*------------------------------------------------------*/
void read_comp_file(void) {
  
  fprintf(stderr, "# Reading compartment file: %s\n", compfile);
  
  std::ifstream cmp_file(compfile);
  std::string line; 
  
  int idx_x;
  int idx1;
  int idx2;
  int line_comp;
  int line_ct;
  int line_cn;
  int line_nd;
  
  // Read line by line. Extract information.
  if (cmp_file.is_open()) {
    while (std::getline(cmp_file, line)) {
      
      std::string::iterator new_end = std::unique(line.begin(), line.end(), BothAreSpaces);
      line.erase(new_end, line.end()); 

      idx1 = line.find("comp", 0)+5;
      idx2 = line.find(" ", idx1);
      line_comp = atoi(line.substr(idx1,idx2-idx1).c_str());
      
      comps.push_back(line_comp);
      
      idx_x = line.find("(", 0);

      if (idx_x!=std::string::npos) {
        
        
        idx1 = line.find("(", 0) + 1;
        idx2 = line.find(",", idx1);
        line_ct = atoi(line.substr(idx1,idx2-idx1).c_str());
        
        idx1 = idx2+1;
        idx2 = line.find(",", idx1);
        line_cn = atoi(line.substr(idx1,idx2-idx1).c_str());
        
        idx1 = idx2+1;
        idx2 = line.find(")", idx1);
        line_nd = atoi(line.substr(idx1,idx2-idx1).c_str());
        
        // Store in vectors.
        comps_ct.push_back(line_ct);
        comps_cn.push_back(line_cn);
        comps_nd.push_back(line_nd);
        
      } else {
        
        // Use previous values.
        fprintf(stderr, "# Comp %d has no node. Use previous cell %d and node -1!\n", line_comp, line_cn);
        
        comps_ct.push_back(line_ct);
        comps_cn.push_back(line_cn);
        comps_nd.push_back(-1);
      }
    }
    cmp_file.close();
  } else {
    fprintf(stderr, "# Could not open comp file. Print it first!\n");
  }

  // Count number of compartments.
  int n_cone_comps   = std::count(comps_ct.begin(), comps_ct.end(), xcone);
  int n_dbp1_comps   = std::count(comps_ct.begin(), comps_ct.end(), dbp1);
  int n_hbp1_comps   = std::count(comps_ct.begin(), comps_ct.end(), hbp1);
  int n_gca_comps    = std::count(comps_ct.begin(), comps_ct.end(), gca);
  int n_gcaoff_comps = std::count(comps_ct.begin(), comps_ct.end(), gcaoff);
  int n_aii_comps    = std::count(comps_ct.begin(), comps_ct.end(), aii);
  
  // Print number of compartments.
  if (n_cone_comps   > 0) fprintf(stderr,"# Number of Cone comps: %d\n", n_cone_comps); 
  if (n_dbp1_comps   > 0) fprintf(stderr,"# Number of DBP1 comps: %d\n", n_dbp1_comps); 
  if (n_hbp1_comps   > 0) fprintf(stderr,"# Number of HBP1 comps: %d\n", n_hbp1_comps); 
  if (n_gca_comps    > 0) fprintf(stderr,"# Number of GCA comps: %d\n", n_gca_comps);   
  if (n_gcaoff_comps > 0) fprintf(stderr,"# Number of GCAOFF comps: %d\n", n_gcaoff_comps);
  if (n_aii_comps    > 0) fprintf(stderr,"# Number of AII comps: %d\n", n_aii_comps);
}

/*------------------------------------------------------*/
void create_external_compartments(void) {
  fprintf(stderr, "# Creating %d compartments.\n", (int)comps.size());
  // Create spheres for every external compartment. Add offset not to interfere with original nodes.
  sphere *s;
  for (int i=0; i<(int)comps.size(); i++) {
    // Parameters: Node and diameter, diameter is arbitrary.
	  s = make_sphere(nd(comps_ct[i],comps_cn[i],ex_comp_node_offset+i),1);
  }
} 
  

/*------------------------------------------------------*/
void connect_external_compartments(void) {

  comp *comp_ex;
  comp *comp_in;

  fprintf(stderr, "# Connect external compartments ...\n");

  int comp_nd = comps_nd[0];
  int i = 0;

  do {
     if (i==0) {
       comp_in = ndn(comps_ct[i], comps_cn[i], comps_nd[i])->comptr;
     } else {
       comp_in = comp_in->next;
     }
    
    // Get external compartment.
    comp_ex = ndn(comps_ct[i], comps_cn[i], ex_comp_node_offset+i)->comptr;
      
    // Connect.
    if (i < 7)  fprintf(stderr, "nd (%d, %d, %d): comp %d <-> external nd %d\n", comps_ct[i], comps_cn[i], comps_nd[i], comp_in->num, ex_comp_node_offset+i);
    if (i == 7) fprintf(stderr, "... \n");
    addcomp_extern(comp_ex, comp_in);
    
    // Next i.
    i++;
    
    if (i == (int)comps.size()) break;
    
  } while(comp_in->next);
  
  fprintf(stderr, "# Done!\n");
}
  
/*------------------------------------------------------*/
void runexpt(void)
{ 
  elem *epnt, *s1, *s2;
  synapse *cur_syn;
  double Vmin, Vmax, fmin, fmax, Imin, Imax, maxca;
  int colr, cell_n, node_n;
  
  double x,y;
  double start;
  double wavel = 527; // wavelength nm. peak of rabbit mcone see wave.cc
  double mask = 0; // no masking
  int spotdia = 300;
  
  // Initialize if necessary. 
  timinc     = (notinit(set_timinc))     ? 1e-5 : set_timinc;
  stiminc    = (notinit(set_stiminc))    ? 1e-4 : set_stiminc;
  srtimestep = (notinit(set_srtimestep)) ? 1e-4 : set_srtimestep;
  ploti      = (notinit(set_ploti))      ? 1e-3 : set_ploti;
  
  if (notinit(rec_predur)) rec_predur = 0.0;
  if (notinit(setxmin)) setxmin = rec_predur;

  // ////////////////////////////////////////////////
  // Print.
  // ////////////////////////////////////////////////

  PRINTER(timinc);
  PRINTER(predur_timinc);
  PRINTER(stiminc);
  PRINTER(srtimestep);
  PRINTER(ploti);
  fprintf(stderr, "\n");

  fprintf(stderr, "# Stimfile: %s\n", stimfile);
  fprintf(stderr, "# dscaeg: %G\n", dscaeg);

  // ////////////////////////////////////////////////
  // Initialize Stimulus.
  // ////////////////////////////////////////////////
  
  // Light Stimulus.
  if (notinit(stimtype)) fprintf(stderr, "# ERROR: Initialize stimtype\n");
  
  if (strcmp(stimtype, "Light") == 0) {
    fprintf(stderr, "# Stimulus Type: Light\n");
    std::string line;
    std::ifstream read_stimfile(stimfile);
    
    if (read_stimfile.is_open()) {
      while (getline(read_stimfile,line)) {
        int idx_comma = line.find(',');
        light_stim_time.push_back(atof(line.substr(0, idx_comma).c_str()));
        light_stim.push_back(atof(line.substr(idx_comma+1, line.size()-(idx_comma+1)).c_str()));
      }
      read_stimfile.close();
    }
    
    fprintf(stderr, "# Stimulus Time: [%G, %G, %G]\n", -predur, light_stim_time[0], light_stim_time[light_stim_time.size()-1]);
    fprintf(stderr, "# Stimulus:      [%G, %G]\n\n", *min_element(light_stim.begin(), light_stim.end()), *max_element(light_stim.begin(), light_stim.end()));
    
    for (int i = 0; i < (light_stim.size()-1); i++) {
      stim_spot(spotdia, x=0, y=0, light_rng*light_stim[i], light_stim_time[i], light_stim_time[i+1]-light_stim_time[i], wavel=wavel, mask=mask);
    }
    stimdur = light_stim_time[light_stim_time.size()-1];
    
  // Vext Stimulus.
  } else if (strcmp(stimtype, "Vext") == 0) {
    fprintf(stderr, "# Stimulus Type: Vext\n");
    stim_file(stimfile);
    read_comp_file();
    create_external_compartments();
  } else {
    fprintf(stderr, "# Warning: No stimulus type given.\n");
  }

  // ////////////////////////////////////////////////
  // Initialize plotting.
  // ////////////////////////////////////////////////
  
  int plot_light = 0;
  int plot_cone = 0;
  int plot_cone_V = 0;
  int plot_cone_cascade = 0;
  int plot_bc_vm_soma = 0;
  int plot_bc_v_soma = 0;
  int plot_bc_i_soma = 0;
  int plot_heatmap_Ca = 0;
  int plot_heatmap_Vm = 0;
  int plot_heatmap_V = 0;
  int plot_bc_ca = 0;
  int plot_bc_rate = 0;
  int plot_bc_syn = 0;
  
  int gc_idx;
  int bp_idx;
  
  if (make_dbp1) {
    bp_idx = dbp1;
    gc_idx = gca;
  } else if (make_hbp1) {
    bp_idx = hbp1;
    gc_idx = gcaoff;
  }
  
  if (strcmp(stimtype, "Light") == 0) plot_light = 1;
     
  if (notinit(rectype)) fprintf(stderr, "# ERROR: Initialize rectype\n");
     
  if (strcmp(rectype, "all") == 0) {
    fprintf(stderr, "# Rec-Type: all\n");
    plot_cone = 1;
    //plot_cone_cascade = 1;
    plot_bc_vm_soma = 1;
    plot_bc_i_soma = 1;
    plot_bc_ca = 1;
    plot_bc_rate = 1;
  } else if (strcmp(rectype, "heatmap") == 0) {
    plot_heatmap_Vm = 1;
    plot_heatmap_V = 1;
    plot_heatmap_Ca = 1;
  } else if (strcmp(rectype, "heatmap_vm") == 0) {
    plot_heatmap_Vm = 1;
  } else if (strcmp(rectype, "heatmap_vm_ca") == 0) {
    plot_heatmap_Vm = 1;
    plot_heatmap_Ca = 1;
  } else if (strcmp(rectype, "heatmap_v") == 0) {
    plot_heatmap_V = 1;
  } else if (strcmp(rectype, "heatmap_ca") == 0) {
    plot_heatmap_Ca = 1;
  } else if (strcmp(rectype, "synapses") == 0) {
    plot_bc_syn = 1;
  } else if (strcmp(rectype, "test") == 0) {
    fprintf(stderr, "# Rec-Type: test\n");
    plot_cone = 1;
    plot_cone_V = 1;
    plot_bc_vm_soma = 1;
    plot_bc_v_soma = 1;
    plot_bc_i_soma = 1;
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
    for (int i = 1; i <= n_cones; i++) {
      plot_vm_nod(xcone, cell_n=i, node_n=0, Vmin = -1, Vmax = 1, colr=blue, ("Cone Vm Soma" + std::to_string(i)).c_str(), 21, 1.0);
      plot_vm_nod(xcone, cell_n=i, node_n=1, Vmin = -1, Vmax = 1, colr=blue, ("Cone Vm Axon Terminal"  + std::to_string(i)).c_str(), 21, 1.0);
      plot_vm_nod(xcone, cell_n=i, node_n=2, Vmin = -1, Vmax = 1, colr=blue, ("Cone Vm Axon"  + std::to_string(i)).c_str(), 21, 1.0);
      
      // Plot cone glutamate release.
      plot_synrate(findsyn(xcone, cell_n=1, node_n=1, bp_idx), fmin=0, fmax=400,  cyan,  20, "Cone",0.3);
    }
  }

  if (plot_cone_V) {
    for (int i = 1; i <= n_cones; i++) {
      plot_v_nod(xcone, cell_n=i, node_n=0, Vmin = -1, Vmax = 1, colr=blue, ("Cone V Soma" + std::to_string(i)).c_str(), 21, 1.0);
      plot_v_nod(xcone, cell_n=i, node_n=1, Vmin = -1, Vmax = 1, colr=blue, ("Cone V Axon Terminal"  + std::to_string(i)).c_str(), 21, 1.0);
      plot_v_nod(xcone, cell_n=i, node_n=2, Vmin = -1, Vmax = 1, colr=blue, ("Cone V Axon"  + std::to_string(i)).c_str(), 21, 1.0);
    }
  }
  
  if (plot_bc_vm_soma) {
    // Plot Vm of BC soma.
    plot_vm_nod(bp_idx, cell_n=1, soma, Vmin = -1, Vmax = 1, colr=blue, "BC Vm Soma", 21, 1.0);
  }
  if (plot_bc_v_soma) {
    // Plot Vm of BC soma.
    plot_v_nod(bp_idx, cell_n=1, soma, Vmin = -1, Vmax = 1, colr=blue, "BC V Soma", 21, 1.0);
  }
  if (plot_bc_i_soma) {
    plot_i_nod(bp_idx, cell_n=1, soma, -5e-9, 5e-9, colr=blue, "Isoma", 21, 1.0);
  }
  
  if (plot_heatmap_Vm | plot_heatmap_V | plot_heatmap_Ca) {
    
    // Plot BC glutamate release.
    for (node_n = 0; node_n < 5000; node_n++) {
      cur_syn = findsyn(bp_idx, cell_n=1, node_n, gc_idx);
      if (cur_syn != NULL) {
        plot_synrate(cur_syn, fmin=0, fmax=400, colr=blue, 8, ("BC ("+ std::to_string(node_n) + ")").c_str(), 0.5);
      }
    }
    
    if ((int)comps.size() == 0) read_comp_file();

    comp *comp_in;
    int i;
    
    int node_nr;
    
    // Get V first.
    i = 0;
    do {
      if (i==0) {
        comp_in = ndn(comps_ct[i], comps_cn[i], comps_nd[i])->comptr;
      } else {
        comp_in = comp_in->next;
      }

      if (comps_ct[i] == bp_idx) { // | comps_ct[i] == xcone
        node_nr = comps_nd[i];
        if (node_nr == -1) {
          fprintf(stderr, "# Warning. Compartment %d has no node. Use node 2.\n", comps[i]);
          node_nr = 2;
        }
        
        if (plot_heatmap_Vm) {
          plot_vm_nod(comps_ct[i], comps_cn[i], node_nr, Vmin = -1, Vmax = 1, colr=blue, ("Vm ("+ std::to_string(comps[i]) + ")").c_str(), 21, 1.0);
        }
        if (plot_heatmap_V) {
          plot_v_nod(comps_ct[i], comps_cn[i], node_nr,  Vmin = -1, Vmax = 1, colr=blue, ("V ("+ std::to_string(comps[i]) + ")").c_str(), 21, 1.0);
        }
      }
      i++;      
    } while(i <= (int)comps.size());
    
    // Ca if there is still capacity.
    i = 0;
    do {
      if (i==0) {
        comp_in = ndn(comps_ct[i], comps_cn[i], comps_nd[i])->comptr;
      } else {
        comp_in = comp_in->next;
      }


      int exclude_comps[] = {-1};

      bool comp_found_idx = std::find(std::begin(exclude_comps), std::end(exclude_comps), comps[i]) != std::end(exclude_comps);

      if (!comp_found_idx) { // comp not in 
        if (comps_ct[i] == bp_idx) { // | comps_ct[i] == xcone
          node_nr = comps_nd[i];
          if (node_nr == -1) {
            fprintf(stderr, "# Warning. Compartment %d has no node. Use node 2.\n", comps[i]);
            node_nr = 2;
          }
          
          if (plot_heatmap_Ca) {
            plot_ca_nod(comps_ct[i], comps_cn[i], node_nr,  1.0e-6, blue, ("Ca ("+ std::to_string(comps[i]) + ")").c_str(), 6,0.5);
          }
        }
      }
      i++;      
    } while(i <= (int)comps.size());  
  }
  
  if (plot_bc_ca) {
    for (node_n = 0; node_n < 3000; node_n++) {
      cur_syn = findsyn(bp_idx, cell_n=1, node_n, gc_idx);
      if (cur_syn != NULL) {
        plot_ca_syn(cur_syn,  1.0e-6, blue, "BC Ca", 6,0.5);
      }
    }
  }
  
  if (plot_bc_syn) {
    for (node_n = 0; node_n < 3000; node_n++) {
      cur_syn = findsyn(bp_idx, cell_n=1, node_n, gc_idx);
      if (cur_syn != NULL) {
        plot_synrate(cur_syn, fmin=0, fmax=400, colr=blue, 8, ("BC ("+ std::to_string(node_n) + ")").c_str(), 0.5);
        plot_vm_nod(bp_idx, cell_n=1, node_n, Vmin = -1, Vmax = 1, colr=blue, ("Vm ("+ std::to_string(node_n) + ")").c_str(), 21, 1.0);
        plot_v_nod(bp_idx, cell_n=1, node_n, Vmin = -1, Vmax = 1, colr=blue, ("V ("+ std::to_string(node_n) + ")").c_str(), 21, 1.0);
        plot_ca_nod(bp_idx, cell_n=1, node_n,  1.0e-6, blue, ("Ca ("+ std::to_string(node_n) + ")").c_str(), 6,0.5);
      }
    }
  }

  if (plot_bc_rate) {
    // Plot BC glutamate release.
    for (node_n = 0; node_n < 3000; node_n++) {
      cur_syn = findsyn(bp_idx, cell_n=1, node_n, gc_idx);
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
  // Run.
  // ////////////////////////////////////////////////

  if (notinit(predur)) fprintf(stderr, "# ERROR: Initialize predur\n");
  
  // Time step before loading if loading.
  
  double sim_timinc = timinc;
  
  if (load_eq) {
    
    if (notinit(pre_load_predur)) pre_load_predur = sim_timinc;
    if (notinit(post_load_predur)) post_load_predur = sim_timinc;
    if (notinit(post_load_predur_timinc)) post_load_predur_timinc = sim_timinc;
    
    fprintf(stderr, "# Run pre load predur: %G\n", pre_load_predur);
    simtime = 0 - pre_load_predur - post_load_predur;
    step (pre_load_predur);
    fprintf(stderr, "# Load equilibrium state from %s.\n", eqfile);
    restoremodel(eqfile);
    if (strcmp(stimtype, "Light") == 0) {
      fprintf(stderr, "# Activate background light, t=%G.\n", simtime);
      stim_backgr(light_min, wavel=wavel, mask=mask, 0, start=simtime);
    }
    fprintf(stderr, "# Run post load predur: %G\n", post_load_predur);
    
    timinc = post_load_predur_timinc;
    step (post_load_predur);
    timinc = sim_timinc;
  } else {
    if (predur < 0.0) predur = 0.0;
    fprintf(stderr, "# predur=%G, t=%G\n", predur, simtime);

    simtime = 0 - predur;
    
    if (strcmp(stimtype, "Light") == 0) {
      fprintf(stderr, "# Activate background light, t=%G.\n", simtime);
      stim_backgr(light_min, wavel=wavel, mask=mask, 0, start=simtime);
    }
    
    // Wait a little.
    if (predur > 0.0) {      
    
      double small_step = min(0.2, predur*0.1);
      
      timinc = sim_timinc;
      synaptau = 1.0;
      fprintf(stderr, "# Run predur with (synaptau=%G, t=%G, dt=%G).\n", synaptau, simtime, timinc);
      step (small_step);
      
      timinc = sim_timinc;
      synaptau = 0.01;
      fprintf(stderr, "# Run predur with (synaptau=%G, t=%G, dt=%G).\n", synaptau, simtime, timinc);
      step (small_step);
      
      timinc = predur_timinc;
      synaptau = 0.01;
      fprintf(stderr, "# Run predur with (synaptau=%G, t=%G, dt=%G).\n", synaptau, simtime, timinc);
      step (predur*0.5);

      timinc = predur_timinc;
      synaptau = 1.0;
      fprintf(stderr, "# Run predur with (synaptau=%G, t=%G, dt=%G).\n", synaptau, simtime, timinc);
      step (-simtime-small_step);

      timinc = sim_timinc;
      synaptau = 1.0;
      fprintf(stderr, "# Run predur with (synaptau=%G, t=%G, dt=%G).\n", synaptau, simtime, timinc);
      step (small_step);
      
    }
    
    // Save equilibrium state?
    if (save_eq) {
      fprintf(stderr, "# Save equilibrium state to %s.\n", eqfile);
      savemodel(eqfile);
    }
  }

  // Step to zero, e.g. for rounding errors.
  timinc = sim_timinc;
  simtime = 0.0;
  synaptau = 1.0;

  // Create external compartments for Vext stimulation.
  if (strcmp(stimtype, "Vext") == 0) connect_external_compartments();

  if (rec_predur > 0.0) {
    fprintf(stderr, "# Run rec_predur=%G, t=%G, dt=%G\n", rec_predur, simtime, timinc);
    step(min(rec_predur, stimdur));
    stimdur -= rec_predur;
  }
  
  if (!run_predur_only) {
    // Run stimulation.
    if (notinit(stimdur)) fprintf(stderr, "# ERROR: Initialize stimdur\n");
    if (stimdur > 0.0) {
      fprintf(stderr, "# Run stimdur=%G, t=%G, dt=%G\n", stimdur, simtime, timinc);
      step (stimdur);
    }
  }
    
  fprintf(stderr, "# Successfully terminated!\n");

}
