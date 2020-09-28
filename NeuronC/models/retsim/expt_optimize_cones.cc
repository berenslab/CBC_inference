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

// Internal variables for stimulus.
std::vector<double> stim_time;
std::vector<double> stim;

// BC parameters.
const char *set_bp_type;
int dbp1;
int gc_type_idx;

// Local parameters.

// Parameters that can be passed in python.
double set_ploti;
double set_timinc;
double set_stiminc;
double set_srtimestep;
double predur_timinc = 1e-4;

int set_pigm;
double set_cone_maxcond;
double predur;
double cone_soma_z;
const char *rectype;
double set_dscaeg; 
		
int print_comps;
const char *compfile;

double light_min = 35e4;
double light_rng = 70e4;

/*------------------------------------------------------*/
void defparams(void) 
{ 
  // Stimulus parameters.
  setptr("stimfile",	     &stimfile);
  setptr("set_bp_type",	   &set_bp_type);
  
  setptr("set_ploti",      &set_ploti);
  setptr("set_timinc",     &set_timinc);
  setptr("predur_timinc",  &predur_timinc);
  setptr("set_stiminc",    &set_stiminc);
  setptr("set_srtimestep", &set_srtimestep);
  setptr("set_dscaeg",     &set_dscaeg);
  
  setptr("rectype",        &rectype);
  setptr("predur",         &predur);
  setptr("cone_soma_z",    &cone_soma_z);
  setptr("print_comps",    &print_comps);
  setptr("compfile",       &compfile);  
  
  // Local parameters.
}

/*------------------------------------------------------*/
void setparams(void)
{ 
  SOMA = R_3;

  // Make cone.
  make_ct(xcone);
  setn(xcone,MORPH,5); // see makcel.cc
  setn(xcone,BIOPHYS,1);
  setn(xcone,NCOLOR,RCOLOR);
  if (!notinit(cone_soma_z)) setn(xcone,SOMAZ,cone_soma_z);
  
  // Make BC type.
  make_ct(dbp1);
  setn(dbp1,MORPH,0);
  setn(dbp1,NCOLOR,RCOLOR);
  setn(dbp1,AXARBT,BRANCHED);
  
  // Print compartment file.
  if (notinit(print_comps)) print_comps = 0;
  if (print_comps) {
    fprintf(stderr, "# print_comps to %s\n", compfile);
    if (print_comps) comp_file(compfile);
  }
  
  checkcellout(xcone);  
  fprintf(stderr, "# Make cones\n");
  
  #define CONEARR 30 
  conexarr  = (double *)emalloc(CONEARR*sizeof(double));
  coneyarr  = (double *)emalloc(CONEARR*sizeof(double));
  conetharr = (double *)emalloc(CONEARR*sizeof(double));
  conenarr  = (int *) emalloc(CONEARR*sizeof(int));
  
  n_cones = 1;
  conexarr[0] = 0; coneyarr[0] = 0; conetharr[0] = 0; conenarr[0] = 1;
  
  checkcellout(dbp1);
  fprintf(stderr, "# Make dbp1\n");
  
  #define DBP1ARR 30
  dbp1xarr  = (double *)emalloc(DBP1ARR*sizeof(double));
  dbp1yarr  = (double *)emalloc(DBP1ARR*sizeof(double));
  dbp1tharr = (double *)emalloc(DBP1ARR*sizeof(double));
  dbp1narr  = (int *) emalloc(DBP1ARR*sizeof(int));
  dbp1xarr[0] = 0; dbp1yarr[0] = 0; dbp1tharr[0] = 0; dbp1narr[0] = 1;
  n_dbp1 = 1;

  // Set general parameters.
  vnoise   = 0;  		// =1 -> vesicle noise
  pnoise   = 0;  		// =1 -> photon noise
  dnoise   = 0;  		// =1 -> continuous dark noise
  Chnoise  = 0;  		// =1 -> channel noise
  
  // Calcium parameters.
  dscavg   = 1E+06; // Calcium sensitivity @ vesicle release  
  dicafrac = 1;
  
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
  
  if (notinit(set_dscaeg)) set_dscaeg = 1;
  dscaeg = set_dscaeg;
}

/*------------------------------------------------------*/
void setdens(void)
// set density and conductance parameters
{
	int cn;
  ndens[xcone][cn=1] = 0;
}

/*------------------------------------------------------*/
void runexpt(void)
{ 

  fprintf(stderr, "# dscaeg: %G\n", dscaeg);

  // Initialize if necessary. 
  timinc        = (notinit(set_timinc))        ? 1e-5 : set_timinc;
  stiminc       = (notinit(set_stiminc))       ? 1e-4 : set_stiminc;
  srtimestep    = (notinit(set_srtimestep))    ? 1e-4 : set_srtimestep;
  ploti         = (notinit(set_ploti))         ? 1e-3 : set_ploti;  

  elem *epnt, *s1, *s2;
  synapse *cur_syn;
  node *npnt;
  double Vmin, Vmax, fmin, fmax, Imin, Imax, maxca;
  int colr, cell_n, node_n;
  
  if (notinit(setxmin)) setxmin = 0;
  
  // ////////////////////////////////////////////////
  // Initialize plotting.
  // ////////////////////////////////////////////////
  
  fprintf(stderr, "# Init recording.\n");
  
  int plot_light = 0;
  int plot_cone_Vm = 0;
  int plot_cone_Vm_soma = 0;
  int plot_cone_rate = 0;
  int plot_cone_cascade = 0;
  int plot_cone_ca = 0;
  
  plot_light = 1;
  if (strcmp(rectype, "all") == 0) {
    fprintf(stderr, "# Rec-Type: all\n");
    plot_cone_Vm = 1;
    plot_cone_rate = 1;
    plot_cone_ca = 1;
    plot_cone_cascade = 1;
  } else if (strcmp(rectype, "test") == 0) {
    fprintf(stderr, "# Rec-Type: test\n");
    plot_cone_Vm = 1;
    plot_cone_rate = 1;
    plot_cone_ca = 1;
  } else if (strcmp(rectype, "optimize") == 0) {
    fprintf(stderr, "# Rec-Type: optimize\n");
    plot_cone_Vm_soma = 1;
    plot_cone_rate = 1;
  } else if (strcmp(rectype, "cascade") == 0) {
    fprintf(stderr, "# Rec-Type: optimize\n");
    plot_cone_cascade = 1;
  } else {
    fprintf(stderr, "# No Rec-Type given!\n");
  }
  
  if (plot_light) {
    // Plot absorbed light at central cones.
    plot_l_nod(xcone, cell_n=1, node_n=0, Vmin = 0, Vmax = 10000, colr=white, "Light", 22, 0.5);
  }
  
  if (plot_cone_Vm) {
    // Plot cone.
    for (node_n = 0; node_n < 100; node_n++) {
      npnt = nde(xcone, 1, node_n);
      if (npnt != NULL) {
        plot_vm_nod(npnt, Vmin = -1, Vmax = 1, colr=blue, ("Vm " + std::to_string(node_n)).c_str(), 21, 1.0);
        plot_v_nod(npnt, Vmin = -1, Vmax = 1, colr=blue, ("V " + std::to_string(node_n)).c_str(),  21, 1.0);
      }
    }
  }
  
  if (plot_cone_Vm_soma) {
    // Plot cone.
    plot_vm_nod(xcone, 1, 0, Vmin = -1, Vmax = 1, colr=blue, "Vm Soma", 21, 1.0);
  }
  
  if (plot_cone_rate) {   
    // Plot cone glutamate release.
    plot_synrate(findsyn(xcone, cell_n=1, node_n=1, dbp1), fmin=0, fmax=400,  cyan,  20, "Cone",0.3);
  }
  
  if (plot_cone_ca) {
    plot_ca_nod(xcone, cell_n=1, node_n=1, 1.0e-6, blue, "Cone Ca", 6,0.5);
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
    for (int ii = 1; ii <= 11; ii++) {
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
  // Print parameters.
  // ////////////////////////////////////////////////
  
  fprintf(stderr, "\n### Stimulus parameters\n");
  fprintf(stderr, "Time:     [%G, %G, %G]\n", -predur, stim_time[0], stim_time[stim_time.size()-1]);
  fprintf(stderr, "Stimulus: [%G, %G]\n", *min_element(stim.begin(), stim.end()), *max_element(stim.begin(), stim.end()));
  
  // Add printer command here.
  fprintf(stderr, "\n");
  
  PRINTER(cone_timec);
  PRINTER(cone_loopg);
  fprintf(stderr, "\n");
  
  PRINTER(dki);
  PRINTER(dnai);
  PRINTER(dcli);
  PRINTER(dcashell);
  fprintf(stderr, "\n");
  
  PRINTER(vna);
  PRINTER(vk);
  PRINTER(vcl);
  fprintf(stderr, "\n");
  
  PRINTER(timinc);
  PRINTER(predur_timinc);
  PRINTER(stiminc);
  PRINTER(srtimestep);
  PRINTER(ploti);
  fprintf(stderr, "\n");
  
  // ////////////////////////////////////////////////
  // Run.
  // ////////////////////////////////////////////////

  double x,y;
  double start;
  double stimdur;
  double wavel = 527; // wavelength nm. peak of rabbit mcone see wave.cc
  
  double mask = 0; // no masking
  int spotdia = 300;
    
  for (int i = 0; i < (stim.size()-1); i++) {
    stim_spot(spotdia, x=0, y=0, light_rng*stim[i], stim_time[i], stim_time[i+1]-stim_time[i], wavel=wavel, mask=mask);
  }
  stimdur = stim_time[stim_time.size()-1];
  
  // Wait a little.
  if (notinit(predur)) predur = 0.3;
  if (predur < 0) predur = 0;
  synaptau = 1;
  simtime = 0 - predur;
  stim_backgr(light_min, wavel=wavel, mask=mask, 0, start=simtime);
  
  double sim_timinc = timinc;
  
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
  
  timinc = sim_timinc;
  synaptau = 1.0;
  simtime = 0.0;

  // Run stimulation.  
  fprintf(stderr, "# Run stimdur: %G\n", stimdur);  
  step (stimdur);
  fprintf(stderr, "# Successfully terminated!\n");

}
