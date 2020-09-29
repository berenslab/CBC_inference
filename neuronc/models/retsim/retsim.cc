/* Simulation of circuit that creates ganglion cell
   receptive field.  Contains array of cones,
   array of cone horizontal cells,
   cone bipolars and ganglion cell. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ncio.h"
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.cc"
#include "ncinit.h"
#include "setexpt.h"


/**********************************************/

int main(int argc, char **argv)
{
   int i, n, ct;
   int morph;
   double retsim_version;
   double scal, nscal;
   double dzmax=LARGENUM, dzmin=-LARGENUM;
   FILE *filout;

 ncinit(argc,argv);
 filout = stderr;
 retsim_version = 1.824;

 timinc	 = 1e-4;
 ploti	 = 1e-3;
 crit	 = 1e-8;
 tempcel = 35;
 relax	 = .1;
 setplsep = 1;

 rseed=13745811;

			/* default membrane params in densfile */
 drm	= 50000;
 gc_rm  = 50000;
 gc_vs  = -0.060;
 gc_vr  = -0.065;
 vcl	= -0.065;
 dvs	= -0.070;	/* default vstart */
 dri	= 200;

 // dcavoff = 0;
 calcnernst = 0;
 dpcak	= 0;		/* permeability of K channels to Ca */
 dpcana	= 0;		/* permeability of Na channels to Ca */
 dpcaampa = 0;		/* permeability of AMPA channels to Ca */

 /* experiments */

 expt = "help";  /* Show "expt" values on command line or here */

 /* dnoise    = 0.2;	/* dark continuous noise in rod */
 dnoise    = 1;		/* dark continuous noise in rod */
 pnoise    = 1;
 vnoise	   = 1;		/* =1 -> synaptic vesicle noise */
 cnoise    = 0;		/* =1 -> synaptic channel noise */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 mxrot     = 0;		/* X rotation */
 myrot     = 0;		/* Y rotation */

/*---------------- define params ---------------------*/

setexpt();		/* find and initialize experiment library functions */

(*defparams)();		/* Define params for the experiment */
			/*  defined in the experiment file */
setvar();		/* set variables from command line */

/*-------------------------------------------------------*/

/* anatomy files */


if (notinit(nvalfile))        nvalfile       = (char *)"nval.n";     /* default neuron params */
if (notinit(chanparamsfile))  chanparamsfile = (char *)"chanparams";
if (notinit(confdir))         confdir        = (char *)"runconf";    /* default run configuration directory */
 
if (notinit(dbp1_file))  dbp1_file   = (char *)"morph_bp";  	/* bipolar cell anatomy */
if (notinit(dbp2_file))  dbp2_file   = (char *)"morph_bp";  	/* bipolar cell anatomy */
if (notinit(dbp3_file))  dbp3_file   = (char *)"morph_bp";  	/* bipolar cell anatomy */
if (notinit(dbp4_file))  dbp4_file   = (char *)"morph_bp";  	/* bipolar cell anatomy */
if (notinit(hbp1_file))  hbp1_file   = (char *)"morph_bp";  	/* bipolar cell anatomy */
if (notinit(hbp2_file))  hbp2_file   = (char *)"morph_bp";  	/* bipolar cell anatomy */
if (notinit(am_file))    am_file     = (char *)"morph_wfamac1";	/* amac cell anatomy */
if (notinit(am2_file))   am2_file    = (char *)"morph_wfamac1";	/* amac cell anatomy */
if (notinit(am3_file))   am3_file    = (char *)"morph_wfamac1";	/* amac cell anatomy */
if (notinit(am4_file))   am4_file    = (char *)"morph_wfamac1";	/* amac cell anatomy */
if (notinit(amh_file))   amh_file    = (char *)"morph_wfamac2";	/* amac cell anatomy */
if (notinit(amh2_file))  amh2_file   = (char *)"morph_wfamac2";	/* amac cell anatomy */
if (notinit(sbac_file))  sbac_file   = (char *)"morph_sb1";  	/* sb   cell anatomy */
if (notinit(rbp_file))   rbp_file    = (char *)"morph_rbp";  	/* rbp  cell anatomy */
if (notinit(aii_file))   aii_file    = (char *)"morph_aii";  	/* aii  cell anatomy */
if (notinit(a17_file))   a17_file    = (char *)"morph_a17";  	/* a17  cell anatomy */
if (notinit(gca_file))   gca_file    = (char *)"morph_beta8b"; 	/* gc   cell anatomy */
if (notinit(gcb_file))   gcb_file    = (char *)"morph_beta8b"; 	/* gc   cell anatomy */
if (notinit(dsgc_file))  dsgc_file   = (char *)"morph_ds1e";   	/* dsgc cell anatomy */
if (notinit(gcaoff_file))gcaoff_file = (char *)"morph_beta8b";  /* gcoff cell anatomy */
if (notinit(gcboff_file))gcboff_file = (char *)"morph_beta8b";  /* gcoff cell anatomy */

if (notinit(def_densfile)) def_densfile=(char *)"dens_default.n"; /* def biophys */

/*----------------------------------------------------*/

if (strcmp(expt,"help")==0) {
  fprintf (stderr,"Retina simulation, version %s\n",print_version(retsim_version));
  fprintf (stderr,"\n");
  fprintf (stderr,"Usage: retsim [options] [-v] file.n [ | vid ] \n");
  fprintf (stderr," options:  -d 1     (display morphology; see 'nc -h')\n");
  fprintf (stderr,"           -R       (display in ray trace mode)\n");
  fprintf (stderr,"           -v       (display in video mode)\n");
  fprintf (stderr,"           --expt surf_area    \n");
  fprintf (stderr,"                | sb_cc_stim   | sb_vc_sclamp\n");
  fprintf (stderr,"                | sb_cc_sine   | dsgc_cc_stim\n");
  fprintf (stderr,"                | dsgc_dbp_stim| dsgc_dbp_bar\n");
  fprintf (stderr,"                | rbp_flash    | aii_flash\n");
  fprintf (stderr,"                | dbp_flash    | hbp_flash \n");
  fprintf (stderr,"                | gc_dbp_flash \n");
  fprintf (stderr,"                | gcoff_hbp_flash \n");
  fprintf (stderr,"                | gc_Rin       | disp_only\n");
  fprintf (stderr,"           --confdir s       (runconf)\n");
  fprintf (stderr,"           --gca_file s      (beta8)\n");
  fprintf (stderr,"           --gcb_file s      (beta8)\n");
  fprintf (stderr,"           --sbac_file s     (sb1)\n");
  fprintf (stderr,"           --dsgc_densfile s (dens_dsgc.n)\n");
  fprintf (stderr,"           --gca_densfile s  (dens_gca.n)\n");
  fprintf (stderr,"           --sbac_densfile s (dens_sbac.n)\n");
  fprintf (stderr,"           --am_densfile s   (dens_am.n)\n");
  fprintf (stderr,"           --make_dbp1  n    (1)\n");
  fprintf (stderr,"           --make_hbp1  n    (1)\n");
  fprintf (stderr,"           --make_rbp   n    (0)\n");
  fprintf (stderr,"           --make_aii   n    (0)\n");
  fprintf (stderr,"           --make_a17   n    (0)\n");
  fprintf (stderr,"           --make_gca   n    (1)\n");
  fprintf (stderr,"           --make_gcb   n    (0)\n");
  fprintf (stderr,"           --am_morph   n    (0)\n");
  fprintf (stderr,"           --sb_morph   n    (2)\n");
  fprintf (stderr,"           --gca_morph  n    (0)\n");
  fprintf (stderr,"           --gca_file   n    (0)\n");
  fprintf (stderr,"           --gcb_file   n    (0)\n");
  fprintf (stderr,"           --stimdia n (stimulus diameter (default=30um))\n");
  exit(0);
};

/*---------------- set build params ---------------------*/

initneurvals();		/* initialize neuron params from nval.n */

(*setparams)();		/* Set or possibly change the neuron params */
			/*  defined in the experiment file */

save_synapses();        /* file for remembering synapses to save time */
restore_synapses();     /* file for restoring synapses from a previous run */

setmorphfiles();	/* set morphology filenames */

initsynconn();		/* initialize synaptic connection table */
                        /* after possible modifications in setparams() */

/*-----------------------------------------------------*/
 
if (!notinit(arrsiz)) setarrsiz = 1;
else                  setarrsiz = 0;

if (!notinit(xarrsiz)) setxarrsiz = 1;
else                   setxarrsiz = 0;

/*------------------- make ------------------------*/

/* Set "make" params consistent with numbers. */
/* User may turn off a celltype by setting number "n_xx" to zero */

 if (n_cones == 0) make_cones = 0;
 if (n_rods  == 0) make_rods  = 0;
 if (n_ha    == 0) make_ha    = 0;
 if (n_hb    == 0) make_hb    = 0;
 if (n_hbat  == 0) make_hbat  = 0;
 if (n_dbp1  == 0) make_dbp1  = 0;
 if (n_dbp2  == 0) make_dbp2  = 0;
 if (n_dbp3  == 0) make_dbp3  = 0;
 if (n_dbp4  == 0) make_dbp4  = 0;
 if (n_hbp1  == 0) make_hbp1  = 0;
 if (n_hbp2  == 0) make_hbp2  = 0;
 if (n_rbp   == 0) make_rbp   = 0;
 if (n_aii   == 0) make_aii   = 0;
 if (n_a17   == 0) make_a17   = 0;
 if (n_am    == 0) make_am    = 0;
 if (n_am2   == 0) make_am2   = 0;
 if (n_am3   == 0) make_am3   = 0;
 if (n_am4   == 0) make_am4   = 0;
 if (n_amh   == 0) make_amh   = 0;
 if (n_amh2  == 0) make_amh2  = 0;
 if (n_ams   == 0) make_ams   = 0;
 if (n_amhs  == 0) make_amhs  = 0;
 if (n_sbac  == 0) make_sbac  = 0;
 if (n_gca   == 0) make_gca   = 0;
 if (n_gcb   == 0) make_gcb   = 0;
 if (n_dsgc  == 0) make_dsgc  = 0;
 if (n_gcaoff == 0) make_gcaoff= 0;
 if (n_gcboff == 0) make_gcboff= 0;

/* set "make" params for automatic processing from command line */

  setn (xcone,MAKE,make_cones);
  setn (xrod, MAKE,make_rods);
  setn (ha,   MAKE,make_ha);
  setn (hb,   MAKE,make_hb);
  setn (hbat, MAKE,make_hbat);
  setn (dbp1, MAKE,make_dbp1);
  setn (dbp2, MAKE,make_dbp2);
  setn (dbp3, MAKE,make_dbp3);
  setn (dbp4, MAKE,make_dbp4);
  setn (hbp1, MAKE,make_hbp1);
  setn (hbp2, MAKE,make_hbp2);
  setn (rbp,  MAKE,make_rbp);
  setn (aii,  MAKE,make_aii);
  setn (a17,  MAKE,make_a17);
  setn (sbac, MAKE,make_sbac);
  setn (am,   MAKE,make_am);
  setn (am2,  MAKE,make_am2);
  setn (am3,  MAKE,make_am3);
  setn (am4,  MAKE,make_am4);
  setn (amh,  MAKE,make_amh);
  setn (amh2, MAKE,make_amh2);
  setn (ams,  MAKE,make_ams);
  setn (amhs, MAKE,make_amhs);
  setn (gca,  MAKE,make_gca);
  setn (gcb,  MAKE,make_gcb);
  setn (dsgc, MAKE,make_dsgc);
  setn (gcaoff,MAKE,make_gcaoff);
  setn (gcboff,MAKE,make_gcaoff);


  if (!notinit(sb_morph))   setn (sbac,  MORPH,sb_morph);
  if (!notinit(rbp_morph))  setn (rbp,   MORPH,rbp_morph);
  if (!notinit(aii_morph))  setn (aii,   MORPH,aii_morph);
  if (!notinit(a17_morph))  setn (a17,   MORPH,a17_morph);
  if (!notinit(am_morph))   setn (am,    MORPH,am_morph);
  if (!notinit(am2_morph))  setn (am2,   MORPH,am2_morph);
  if (!notinit(am3_morph))  setn (am3,   MORPH,am3_morph);
  if (!notinit(am4_morph))  setn (am4,   MORPH,am4_morph);
  if (!notinit(amh_morph))  setn (amh,   MORPH,amh_morph);
  if (!notinit(amh2_morph)) setn (amh2,  MORPH,amh2_morph);
  if (!notinit(gca_morph))  setn (gca,   MORPH,gca_morph);
  if (!notinit(gcb_morph))  setn (gcb,   MORPH,gcb_morph);
  if (!notinit(dsgc_morph)) setn (dsgc,  MORPH,dsgc_morph);
  if (!notinit(gcaoff_morph)) setn (gcaoff,MORPH,gcaoff_morph);
  if (!notinit(gcboff_morph)) setn (gcboff,MORPH,gcboff_morph);
  

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

if (!notinit(gca_biophys))   setn(GCA,BIOPHYS,gca_biophys);
if (!notinit(gcb_biophys))   setn(GCB,BIOPHYS,gcb_biophys);
if (!notinit(dsgc_biophys))  setn(DSGC,BIOPHYS,dsgc_biophys);
if (!notinit(gcaoff_biophys))setn(GCAOFF,BIOPHYS,gcaoff_biophys);
if (!notinit(gcboff_biophys))setn(GCBOFF,BIOPHYS,gcboff_biophys);
if (!notinit(am_biophys))    setn(AM,BIOPHYS,am_biophys);
if (!notinit(am2_biophys))   setn(AM,BIOPHYS,am2_biophys);
if (!notinit(am3_biophys))   setn(AM,BIOPHYS,am3_biophys);
if (!notinit(am4_biophys))   setn(AM,BIOPHYS,am4_biophys);
if (!notinit(amh_biophys))   setn(AMH,BIOPHYS,amh_biophys);
if (!notinit(amh2_biophys))  setn(AMH2,BIOPHYS,amh2_biophys);
if (!notinit(rbp_biophys))   setn(RBP,BIOPHYS,rbp_biophys);
if (!notinit(aii_biophys))   setn(AII,BIOPHYS,aii_biophys);
if (!notinit(a17_biophys))   setn(A17,BIOPHYS,a17_biophys);
if (!notinit(sb_biophys))    setn(SBAC,BIOPHYS,sb_biophys);


findceldens();		/* read and set default channel densities from density files */

(*setdens)();		/* User defined changes for cell channel densities. */
			/*  defined in the experiment file */

modchandens();		/* possibly modify channel densities */
			/* automatic algorithm */

if (anysetn(BIOPHYS)) {	  /* reduce time step for membr. channels */
  timinc = 5e-6;
  ploti  = 5e-5;
  if (anysetn(CHNOISE)) {	  /* reduce time step for membr. channels */
    timinc = 1e-6;
  };
};

if (!vnoise) {
  for (ct=0; ct<nceltypes; ct++) {
    // if ((ct==dbp1))
    for (n=1; n<=NCONNO; n++) 
      setsv(ct,SVNOISE,n,0);
  };
};

if (!cnoise) {
  for (ct=0; ct<nceltypes; ct++) {
    for (n=1; n<=NCONNO; n++) 
      setsv(ct,SCNOISE,n,0);
  };
};

if (notinit(Chnoise)) Chnoise = 0;

if (Chnoise) {          /* turn on all channel noise */
  if (notinit(timinc)) timinc  = 1e-6;

  if (notinit(na1nois))  na1nois  = 1;
  if (notinit(na2nois))  na2nois  = 1;
  if (notinit(na3nois))  na3nois  = 1;
  if (notinit(na4nois))  na4nois  = 1;
  if (notinit(na5nois))  na5nois  = 1;
  if (notinit(na6nois))  na6nois  = 1;
  if (notinit(na8nois))  na8nois  = 1;
  if (notinit(k1nois))   k1nois   = 1;
  if (notinit(k3nois))   k3nois   = 1;
  if (notinit(k4nois))   k4nois   = 1;
  if (notinit(k5nois))   k5nois   = 1;
  if (notinit(k6nois))   k6nois   = 1;
  if (notinit(k7nois))   k7nois   = 1;
  if (notinit(hcn1nois))   hcn1nois   = 1;
  if (notinit(hcn2nois))   hcn2nois   = 1;
  if (notinit(hcn3nois))   hcn3nois   = 1;
  if (notinit(hcn4nois))   hcn4nois   = 1;
  if (notinit(kca1nois)) kca1nois = 1;
  if (notinit(kca3nois)) kca3nois = 1;
  if (notinit(kca4nois)) kca4nois = 1;
  if (notinit(kca5nois)) kca5nois = 1;
  if (notinit(kca6nois)) kca6nois = 1;
  if (notinit(ca1nois))  ca1nois  = 1;
  if (notinit(ca3nois))  ca3nois  = 1;
  if (notinit(ca5nois))  ca5nois  = 1;
  if (notinit(ca6nois))  ca6nois  = 1;
  if (notinit(ca7nois))  ca7nois  = 1;
  if (notinit(clcanois)) clcanois = 1;

}
else {
  if (notinit(na1nois))  na1nois  = 0;
  if (notinit(na2nois))  na2nois  = 0;
  if (notinit(na3nois))  na3nois  = 0;
  if (notinit(na4nois))  na4nois  = 0;
  if (notinit(na5nois))  na5nois  = 0;
  if (notinit(na6nois))  na6nois  = 0;
  if (notinit(na8nois))  na8nois  = 0;
  if (notinit(k1nois))   k1nois   = 0;
  if (notinit(k3nois))   k3nois   = 0;
  if (notinit(k4nois))   k4nois   = 0;
  if (notinit(k5nois))   k5nois   = 0;
  if (notinit(k6nois))   k6nois   = 0;
  if (notinit(k7nois))   k7nois   = 0;
  if (notinit(hcn1nois))   hcn1nois   = 0;
  if (notinit(hcn2nois))   hcn2nois   = 0;
  if (notinit(hcn3nois))   hcn3nois   = 0;
  if (notinit(hcn4nois))   hcn4nois   = 0;
  if (notinit(kca1nois)) kca1nois = 0;
  if (notinit(kca3nois)) kca3nois = 0;
  if (notinit(kca4nois)) kca4nois = 0;
  if (notinit(kca5nois)) kca5nois = 0;
  if (notinit(kca6nois)) kca6nois = 0;
  if (notinit(ca1nois))  ca1nois  = 0;
  if (notinit(ca3nois))  ca3nois  = 0;
  if (notinit(ca5nois))  ca5nois  = 0;
  if (notinit(ca6nois))  ca6nois  = 0;
  if (notinit(ca7nois))  ca7nois  = 0;
  if (notinit(clcanois)) clcanois = 0;
};

/*---------------------------------------------------*/

/* print version, input file, date, host and morphology file names */

if (ninfo>=1) {
  if (disp_ray) {
     filout = stdout;
     ncfprintf (filout,"/*\n");
  }
  ncfprintf (filout,"# Retina simulation\n");
  ncfprintf (filout,"#\n");
  ncfprintf (filout,"#   %s version:    %s     \n",  progname, print_version(retsim_version));
  ncfprintf (filout,"#   nc version:        %s     \n", print_version(ncversion));
  ncfprintf (filout,"#   date:              %s     \n", xsystem("date"));
  ncfprintf (filout,"#   machine:           %s     \n", xsystem("hostname -s"));
  ncfprintf (filout,"#   experiment:        %s     \n",  expt);
  printfilenames(filout);
  ncfprintf (filout,"#\n");
  if (disp_ray) ncfprintf (filout,"*/\n");
}

/* print biophysical properties */

printchaninfo();


/*---------------------------------------------------------*/

/* Rationale for order of cell creation:

 1) an expt requires:
   - make gcs
   - make enough amacs to fill gc dend field, update size of array
   - make enough bipolars to fill amac dend field, update size of array
   - if enough bipolars exist to fill hz dend field, add hzs, 
	update size of array, else add extra bipolars, 
	then add hzs, update size of array.
    - make enough photoreceptors to fill bp receptive fields

 2) grow dendrites of hz, am, and bp cells, make synapses.

*/

 if (notinit(maxha))   maxha   = (int)getn(ha,  MAXNUM);
 if (notinit(maxhb))   maxhb   = (int)getn(hb,  MAXNUM);
 if (notinit(maxhbat)) maxhbat = (int)getn(hbat,MAXNUM);
 if (notinit(maxrbp))  maxrbp  = (int)getn(rbp, MAXNUM);
 if (notinit(maxdbp1)) maxdbp1 = (int)getn(dbp1,MAXNUM);
 if (notinit(maxdbp2)) maxdbp2 = (int)getn(dbp2,MAXNUM);
 if (notinit(maxhbp1)) maxhbp1 = (int)getn(hbp1,MAXNUM);
 if (notinit(maxa17))  maxa17  = (int)getn(a17, MAXNUM);
 if (notinit(maxaii))  maxaii  = (int)getn(aii, MAXNUM);
 if (notinit(maxsb))   maxsb   = (int)getn(sbac,MAXNUM);
 if (notinit(maxam))   maxam   = (int)getn(am,  MAXNUM);
 if (notinit(maxam2))  maxam2  = (int)getn(am2, MAXNUM);
 if (notinit(maxam3))  maxam3  = (int)getn(am3, MAXNUM);
 if (notinit(maxam4))  maxam4  = (int)getn(am4, MAXNUM);
 if (notinit(maxamh))  maxamh  = (int)getn(amh, MAXNUM);
 if (notinit(maxamh2)) maxamh2 = (int)getn(amh2, MAXNUM);
 if (notinit(maxams))  maxams  = (int)getn(ams, MAXNUM);
 if (notinit(maxamhs)) maxamhs = (int)getn(amhs,MAXNUM);
 if (notinit(maxgca))  maxgca  = (int)getn(gca, MAXNUM);
 if (notinit(maxgcb))  maxgcb  = (int)getn(gcb, MAXNUM);
 if (notinit(maxdsgc)) maxdsgc = (int)getn(dsgc,MAXNUM);
 if (notinit(maxgcaoff))maxgcaoff = (int)getn(gcaoff,MAXNUM);
 if (notinit(maxgcboff))maxgcaoff = (int)getn(gcboff,MAXNUM);

 maxha   = min(maxha,  (int)getn(ha,  MAXNUM));	/* number of HA cells */
 maxhb   = min(maxhb,  (int)getn(hb,  MAXNUM));	/* number of HB cells */
 maxhbat = min(maxhbat,(int)getn(hbat,MAXNUM));	/* number of HB axon terms */
 maxdbp1 = min(maxdbp1,(int)getn(dbp1,MAXNUM)); /* number of dbp1 cone bipolars  */
 maxdbp2 = min(maxdbp2,(int)getn(dbp2,MAXNUM)); /* number of dbp2 cone bipolars  */
 maxhbp1 = min(maxhbp1,(int)getn(hbp1,MAXNUM)); /* number of hbp1 cone bipolars  */
 maxrbp  = min(maxrbp, (int)getn(rbp, MAXNUM));	/* number of rod bipolars */
 maxa17  = min(maxa17, (int)getn(a17, MAXNUM));	/* number of a17 amacrines */
 maxaii  = min(maxaii, (int)getn(aii, MAXNUM));	/* number of aii amacrines */
 maxsb   = min(maxsb,  (int)getn(sbac,MAXNUM));	/* number of sb amacrines */
 maxam   = min(maxam,  (int)getn(am,  MAXNUM));	/* number of am amacrines */
 maxam2  = min(maxam2, (int)getn(am2, MAXNUM));	/* number of am2 amacrines */
 maxam3  = min(maxam3, (int)getn(am3, MAXNUM));	/* number of am3 amacrines */
 maxam4  = min(maxam4, (int)getn(am4, MAXNUM));	/* number of am4 amacrines */
 maxamh  = min(maxamh, (int)getn(amh, MAXNUM));	/* number of amh amacrines */
 maxamh2 = min(maxamh2,(int)getn(amh2,MAXNUM));	/* number of amh2 amacrines */
 maxams  = min(maxams, (int)getn(ams, MAXNUM));	/* number of ams amacrines */
 maxamhs = min(maxamhs,(int)getn(amhs, MAXNUM));/* number of amhs amacrines */
 maxgca  = min(maxgca, (int)getn(gca, MAXNUM));	/* number of GC cells */
 maxgcb  = min(maxgcb, (int)getn(gcb, MAXNUM));	/* number of GC cells */
 maxdsgc = min(maxdsgc,(int)getn(dsgc,MAXNUM));	/* number of DSGC cells */
 maxgcaoff= min(maxgcaoff,(int)getn(gcaoff,MAXNUM));	/* number of off GC cells */
 maxgcboff= min(maxgcboff,(int)getn(gcboff,MAXNUM));	/* number of off GC cells */

/*---------------------------------------------------------*/

if(!notinit(arrsiz) && notinit(xarrsiz)) {  /* if user specifies arrsiz, make square array */
      xarrsiz = arrsiz;
      yarrsiz = arrsiz;
};

if (!notinit(garrsiz)) {	/* If ganglion cell array size has been set, */
   gxarrsiz = garrsiz;		/*  use it to set gc soma positions. */
   gyarrsiz = garrsiz;
}
else if (notinit(gxarrsiz)) {
    gxarrsiz = xarrsiz;		/* allow rectangular array of gcs */
    gyarrsiz = yarrsiz;
}

if (notinit(arrcentx)) arrcentx = 0;
if (notinit(arrcenty)) arrcenty = 0;

/*-------------------------Ganglion cells------------------------------*/
  
if (ninfo>=2) ncfprintf (stderr,"# Constructing circuit:\n#\n");

if (notinit(update_array)) update_array = 0; 	/* set array size from cells present */

/* First, make large cells (gc, hz cells, amacrines) */

if (make_gca) {
  if (gcxarr==NULL) {
     ngca=setupcells(gca,n_gca,gxarrsiz,gyarrsiz);
     if (update_array) enlarge_array(gca,-1);
  } else {
     ngca=setupcells(gca,n_gca,gcxarr,gcyarr,gctharr);
  }
}
if (make_gcb) {
  if (gcxarr==NULL) {
     ngcb=setupcells(gcb,n_gcb,gxarrsiz,gyarrsiz);
     if (update_array) enlarge_array(gcb,-1);
  } else {
     ngcb=setupcells(gcb,n_gcb,gcxarr,gcyarr);
  }
}

if (make_dsgc) {
  if (gcxarr==NULL) {
     ndsgc=setupcells(dsgc,n_dsgc,gxarrsiz,gyarrsiz);
     if (update_array) enlarge_array(dsgc,-1);
  } else {
     ndsgc=setupcells(dsgc,n_dsgc,gcxarr,gcyarr,gctharr);
  }
};

if (make_gcaoff) {
  if (gcxarr==NULL) {
     ngcaoff=setupcells(gcaoff,n_gcaoff,gxarrsiz,gyarrsiz);
     if (update_array) enlarge_array(gcaoff,-1);
  } else {
     ngcaoff=setupcells(gcaoff,n_gcaoff,gcxarr,gcyarr);
  }
};

if (make_gcboff) {
  if (gcxarr==NULL) {
     ngcboff=setupcells(gcboff,n_gcboff,gxarrsiz,gyarrsiz);
     if (update_array) enlarge_array(gcboff,-1);
  } else {
     ngcboff=setupcells(gcboff,n_gcboff,gcxarr,gcyarr);
  }
};


if (ninfo>=2 && (make_gca || make_gcb || make_dsgc || make_gcaoff || make_gcb || make_gcboff)) {
  if (!notinit(gxarrsiz))
       fprintf(stderr, "# ganglion cells done, garrsiz= x %g y %g um\n", gxarrsiz,gyarrsiz);
  if (!notinit(xarrsiz))     fprintf(stderr, "# arrsiz= x %g y %g um\n", xarrsiz,yarrsiz);
  fprintf(stderr, "# ganglion cells done.\n");
}
/*----------------------Starburst amacrine cells------------------------*/

if (make_sbac) {

  if (sb_morph>0) {
    if (notinit(am_seglen))  am_seglen  = 10;          /* length of seg  */
    if (notinit(am_den_seg)) am_den_seg = xradius/am_seglen;  /* segs/br */
  };

  if (!notinit(sbarrsiz)) {
      sbxarrsiz = sbarrsiz;
      sbyarrsiz = sbarrsiz;
      nsbac = setupcells(sbac,n_sbac,sbxarrsiz,sbyarrsiz);
  }
  else if (!notinit(sbxarrsiz)) {
      nsbac = setupcells(sbac,n_sbac,sbxarrsiz,sbyarrsiz);
  } 
  else if (sbxarr!=NULL) {
       nsbac = setupcells(sbac,n_sbac, sbxarr, sbyarr, sbytharr, sbtharr, sbnarr);
  }
  else
      nsbac = setupcells(sbac,n_sbac);
  // if (ninfo>=2) printf ("# starburst amacrine cells made = %d\n",nsbac);
};

/*----------------------A17 amacrine cells------------------------*/

if (make_a17) {		 		/* A17 cells */
  if (a17xarr!=NULL) {
     if (a17narr!=NULL)  
       na17 = setupcells(a17, n_a17, a17xarr, a17yarr, a17tharr, a17narr);
     else
      na17 = setupcells(a17, n_a17, a17xarr, a17yarr, a17tharr);
  }
  else
      na17=setupcells(a17,n_a17);
}

/*----------------------AM amacrine cells------------------------*/

if (make_am) {		 		/* am cells */
  if (!notinit(amarrsiz)) {
      amxarrsiz = amarrsiz;
      amyarrsiz = amarrsiz;
      nam = setupcells(am,n_am,amxarrsiz,amyarrsiz);
  }
  else if (!notinit(amxarrsiz)) {
      nam = setupcells(am,n_am,amxarrsiz,amyarrsiz);
  } 
  else if (amxarr!=NULL) {
      nam = setupcells(am,n_am, amxarr, amyarr, amtharr);
  }
  else
      nam=setupcells(am,n_am);
};

/*----------------------AM2 amacrine cells------------------------*/

if (make_am2) {		 		/* am2 cells */
  if (!notinit(am2arrsiz)) {
      am2xarrsiz = am2arrsiz;
      am2yarrsiz = am2arrsiz;
      nam2 = setupcells(am2,n_am2,am2xarrsiz,am2yarrsiz);
  }
  else if (!notinit(am2xarrsiz)) {
      nam2 = setupcells(am2,n_am2,am2xarrsiz,am2yarrsiz);
  } 
  else if (am2xarr!=NULL) {
      nam2 = setupcells(am2,n_am2, am2xarr, am2yarr, am2tharr);
  }
  else
      nam2=setupcells(am2,n_am2);
};

/*----------------------AM3 amacrine cells------------------------*/

if (make_am3) {		 		/* am3 cells */
  if (!notinit(am3arrsiz)) {
      am3xarrsiz = am3arrsiz;
      am3yarrsiz = am3arrsiz;
      nam3 = setupcells(am3,n_am3,am3xarrsiz,am3yarrsiz);
  }
  else if (!notinit(am3xarrsiz)) {
      nam3 = setupcells(am3,n_am3,am3xarrsiz,am3yarrsiz);
  } 
  else if (am3xarr!=NULL) {
      nam3 = setupcells(am3,n_am3, am3xarr, am3yarr, am3tharr);
  }
  else
      nam3=setupcells(am3,n_am3);
};

/*----------------------AM4 amacrine cells------------------------*/

if (make_am4) {		 		/* am4 cells */
  if (!notinit(am4arrsiz)) {
      am4xarrsiz = am4arrsiz;
      am4yarrsiz = am4arrsiz;
      nam4 = setupcells(am4,n_am4,am4xarrsiz,am4yarrsiz);
  }
  else if (!notinit(am4xarrsiz)) {
      nam4 = setupcells(am4,n_am4,am4xarrsiz,am4yarrsiz);
  } 
  else if (am4xarr!=NULL) {
      nam4 = setupcells(am4,n_am4, am4xarr, am4yarr, am4tharr);
  }
  else
      nam4=setupcells(am4,n_am4);
};

/*----------------------AMH amacrine cells------------------------*/

if (make_amh) {		 		/* amh cells */
  if (!notinit(amharrsiz)) {
      amhxarrsiz = amharrsiz;
      amhyarrsiz = amharrsiz;
      namh = setupcells(amh,n_amh,amhxarrsiz,amhyarrsiz);
  }
  else if (!notinit(amhxarrsiz)) {
      namh = setupcells(amh,n_amh,amhxarrsiz,amhyarrsiz);
  } 
  else if (amhxarr!=NULL) {
      namh = setupcells(amh,n_amh, amhxarr, amhyarr, amhtharr);
  }
  else
      namh=setupcells(amh,n_amh);
};

/*----------------------AMH2 amacrine cells------------------------*/

if (make_amh2) {	 		/* amh2 cells */
  if (!notinit(amh2arrsiz)) {
      amh2xarrsiz = amh2arrsiz;
      amh2yarrsiz = amh2arrsiz;
      namh2 = setupcells(amh2,n_amh2,amh2xarrsiz,amh2yarrsiz);
  }
  else if (!notinit(amh2xarrsiz)) {
      namh2 = setupcells(amh2,n_amh2,amh2xarrsiz,amh2yarrsiz);
  } 
  else if (amh2xarr!=NULL) {
      namh2 = setupcells(amh2,n_amh2, amh2xarr, amh2yarr, amh2tharr);
  }
  else
      namh2=setupcells(amh2,n_amh2);
};

/*------------------------------------------------------------------------*/
/*----------------------small-field amacrine cells------------------------*/

/*  first, get array size */

if (update_array) enlarge_array(-1,-1);

/*----------------------AMS amacrine cells------------------------*/

if (make_ams) {		 		/* ams cells */
  if (!notinit(amsarrsiz)) {
      amsxarrsiz = amsarrsiz;
      amsyarrsiz = amsarrsiz;
      nams = setupcells(ams,n_ams,amsxarrsiz,amsyarrsiz);
  }
  else if (!notinit(amsxarrsiz)) {
      nams = setupcells(ams,n_ams,amsxarrsiz,amsyarrsiz);
  } 
  else if (amsxarr!=NULL) {
      nams = setupcells(ams,n_ams, amsxarr, amsyarr, amstharr);
  }
  else
      nams=setupcells(ams,n_ams);
};

/*----------------------AMHS amacrine cells------------------------*/

if (make_amhs) {		 		/* amhs cells */
  if (!notinit(amhsarrsiz)) {
      amhsxarrsiz = amhsarrsiz;
      amhsyarrsiz = amhsarrsiz;
      namhs = setupcells(amhs,n_amhs,amhsxarrsiz,amhsyarrsiz);
  }
  else if (!notinit(amhsxarrsiz)) {
      namhs = setupcells(amhs,n_amhs,amhsxarrsiz,amhsyarrsiz);
  } 
  else if (amhsxarr!=NULL) {
      namhs = setupcells(amhs,n_amhs, amhsxarr, amhsyarr, amhstharr);
  }
  else
      namhs=setupcells(amhs,n_amhs);
};

/*----------------------AII amacrine cells------------------------*/

if (make_aii) {		 		/* AII cells */
  if (!notinit(aiiarrsiz)) {
      aiixarrsiz = aiiarrsiz;
      aiiyarrsiz = aiiarrsiz;
      naii = setupcells(aii,n_aii,aiixarrsiz,aiiyarrsiz);
  }
  else if (!notinit(aiixarrsiz)) {
      naii = setupcells(aii,n_aii,aiixarrsiz,aiiyarrsiz);
  } 
  else if (aiixarr!=NULL) {
      naii = setupcells(aii,n_aii, aiixarr, aiiyarr);
  }
  else
      naii=setupcells(aii,n_aii);
};

if ((make_aii || make_sbac || make_am || make_am2 || make_am3 || make_am4 || 
			     make_amh || make_amh2 || make_ams || make_amhs || make_a17) &&
			   (ninfo>=2)) ncfprintf (stderr,"# amacrine cells done\n");

/*------------------------Horizontal cells-----------------------------*/

if (make_hbat) {	/* rod Hz cells */
  nhbat=setupcells(hbat,n_hbat);
};

if (make_ha) {		/* HA cells */
  nha=setupcells(ha,n_ha);
};

if (make_hb) {	 	/* HB cells */
  nhb=setupcells(hb,n_hb);
};

if ((make_ha || make_hb || make_hbat) && (ninfo>=2)) 
	ncfprintf (stderr,"# horizontal cells done\n");

/*---------------------------------------------------------------------*/

/* Before making bps and photoreceptors,
   make array larger so that bp,hz dens can grow out.
   Then fill array with photoreceptors and bipolars
   After dendrites have grown, remove photorecs, bps without connections.
*/
	/* if none of the above cells are made, set arrsiz now */


if(setarrsiz==0) { 

 if(notinit(xarrsiz)) xarrsiz = 10;
 if(notinit(yarrsiz)) yarrsiz = 10;
 
 if (make_ha && make_hb) {/* add a little more than half dend tree diam */
   xarrsiz += max(getn(ha,DTREEDIA),getn(hb,DTREEDIA))*1.5;
   yarrsiz += max(getn(ha,DTREEDIA),getn(hb,DTREEDIA))*1.5;}
 else if (make_ha){
   xarrsiz += getn(ha,DTREEDIA)*1.5;
   yarrsiz += getn(ha,DTREEDIA)*1.5;}
 else if (make_hb){
   xarrsiz += getn(hb,DTREEDIA)*1.5;
   yarrsiz += getn(hb,DTREEDIA)*1.5;}
 else if (make_hbat){
   xarrsiz += getn(hbat,DTREEDIA)*1.5;
   yarrsiz += getn(hbat,DTREEDIA)*1.5;
 }
 if (make_dbp1) {
   axarbdia = getn(dbp1,AXARBDIA);
   if (xarrsiz<axarbdia) xarrsiz=axarbdia;	/* set minimum arrsiz */
   xarrsiz += getn(dbp1,DTREEDIA)*1.5;
   if (yarrsiz<axarbdia) yarrsiz=axarbdia;
   yarrsiz += getn(dbp1,DTREEDIA)*1.5;
 }
 if (make_dbp2) {
   axarbdia = getn(dbp2,AXARBDIA);
   if (xarrsiz<axarbdia) xarrsiz=axarbdia;	/* set minimum arrsiz */
   xarrsiz += getn(dbp2,DTREEDIA)*1.5;
   if (yarrsiz<axarbdia) yarrsiz=axarbdia;
   yarrsiz += getn(dbp2,DTREEDIA)*1.5;
 }
 if (make_dbp3 && !make_dbp2) {
   axarbdia = getn(dbp3,AXARBDIA);
   if (xarrsiz<axarbdia) xarrsiz=axarbdia;	/* set minimum arrsiz */
   xarrsiz += getn(dbp3,DTREEDIA)*1.5;
   if (yarrsiz<axarbdia) yarrsiz=axarbdia;
   yarrsiz += getn(dbp3,DTREEDIA)*1.5;
 }
 if (make_dbp4 && !make_dbp2 && !make_dbp3) {
   axarbdia = getn(dbp4,AXARBDIA);
   if (xarrsiz<axarbdia) xarrsiz=axarbdia;	/* set minimum arrsiz */
   xarrsiz += getn(dbp4,DTREEDIA)*1.5;
   if (yarrsiz<axarbdia) yarrsiz=axarbdia;
   yarrsiz += getn(dbp4,DTREEDIA)*1.5;
 }
 if (make_hbp1) {
   axarbdia = getn(hbp1,AXARBDIA);
   if (xarrsiz<axarbdia) xarrsiz=axarbdia;	/* set minimum arrsiz */
   xarrsiz += getn(hbp1,DTREEDIA)*1.5;
   if (yarrsiz<axarbdia) yarrsiz=axarbdia;
   yarrsiz += getn(hbp1,DTREEDIA)*1.5;
 }
 if (make_rbp) {
   xarrsiz += getn(rbp,DTREEDIA)*1.5;
   yarrsiz += getn(rbp,DTREEDIA)*1.5;
 }

   find_maxmin(-1,-1);
   arrcentx=(xmax+xmin)/2;
   arrcenty=(ymax+ymin)/2;
 
 arrsiz  = max(xarrsiz, yarrsiz);
}
else {
   xarrsiz = arrsiz;
   yarrsiz = arrsiz;
};

if(notinit(arrcentx)) arrcentx=0;
if(notinit(arrcenty)) arrcenty=0;


if (ninfo>=2)
 ncfprintf(stderr,"# xarrsiz = %g; yarrsiz=%g; arrcentx=%g; arrcenty=%g\n",
                          xarrsiz, yarrsiz,arrcentx,arrcenty);
/* */

/*------------------------make rods------------------------------------*/

if (notinit(rod_rect_arr)) rod_rect_arr = 0;

if (make_rods) {
     int i,j,rwid, rodnum;
 
  if (rod_rect_arr) {

    if (notinit(rodspac))
        rodspac   = 1 / sqrt (getn(xrod,DENS)* 1e-6);     	/* rod spacing */

    if (ninfo>=2) ncfprintf (stderr,"# rod spacing: %g\n",rodspac);
    rodarrsiz = int(arrsiz/rodspac);			  /* size of rod array */
    if (ninfo>=2) {
      ncfprintf (stderr,"# rod array size: %d\n",rodarrsiz);
    };
    rodnum = 0;
    rwid = int(rodarrsiz/2);
    for (j= -rwid; j<rwid; j++)
      for (i= -rwid; i<rwid; i++) {
         makcell (xrod, ++rodnum, i*rodspac, j*rodspac,0,0,0);
      };
    nrods = rodnum;
    if (ninfo>=2) ncfprintf (stderr,"# number of rods: %d\n", nrods);

  } else  if (rodxarr!=NULL) {   // look for predefined rod positions

     if (rodnarr!=NULL)  
       nrods = setupcells(xrod, n_rods, rodxarr, rodyarr, rodtharr, rodnarr);
     else
      nrods = setupcells(xrod, n_rods, rodxarr, rodyarr, rodtharr);
  }
  else { /* random rod array, uses either n_rods or arrsiz or rod spacing from nval file */

    if (!notinit(rodarrsiz) || !notinit(rodxarrsiz)) {  // if rod array size defined
        if (!notinit(rodarrsiz)) {
	      rodxarrsiz = rodarrsiz;
	      rodyarrsiz = rodarrsiz;
        }
        if (!notinit(rodarrcentx)) {  			// position for the array
             nrods = setupcells(xrod,n_rods,rodxarrsiz,rodyarrsiz,rodarrcentx,rodarrcenty);
        } else {					// else use general position arrcentx,y def=0
             nrods = setupcells(xrod,n_rods,rodxarrsiz,rodyarrsiz,arrcentx,arrcenty);
        }
     } else {						// otherwise use default array size and center
          nrods = setupcells(xrod,n_rods);
     }
  }

  setn (xrod,NMADE,nrods);

  find_maxmin(xrod,-1);	/* get size of array, update array size */
  if (notinit(arrsiz)) {
    xarrsiz=xmax-xmin;
    yarrsiz=ymax-ymin;
    arrcentx=(xmax+xmin)/2;
    arrcenty=(ymax+ymin)/2;
    arrsiz=max(xarrsiz,yarrsiz);
  } else {
    xarrsiz=arrsiz;
    yarrsiz=arrsiz;
    arrcentx=0;
    arrcenty=0;
  }
}  /* make rods */

/*------------------------make cones---------------------------------- */

if (notinit(cone_rect_arr)) cone_rect_arr = 0;

if (make_cones) {
     int i,j, conenum, rwid;

  if (cone_rect_arr) {				/* make rectangular array */

     if (notinit(conespac))
         conespac = 1 / sqrt (getn(xcone,DENS)* 1e-6);	/* cone spacing */

     if (ninfo>=2) ncfprintf (stderr,"# cone spacing: %.3g um\n",conespac);
     conearrsiz = arrsiz / conespac;		/* size of cone array */
     if (ninfo>=2) {
       ncfprintf (stderr,"# cone square array size: %d\n",conearrsiz);
     };
     conenum = 0;
     rwid = int(conearrsiz/2);
     for (j= -rwid; j<rwid; j++)
        for (i= -rwid; i<rwid; i++) {
          makcell (xcone, ++conenum, i*conespac, j*conespac, 0, 0, 0);
        };
     ncones = conenum;
     if (ninfo>=2) ncfprintf (stderr,"# number of cones: %d\n", ncones);

  } else if (conexarr!=NULL) {   // look for predefined cone positions

     if (conenarr!=NULL)  
       ncones = setupcells(xcone, n_cones, conexarr, coneyarr, conetharr, conenarr);
     else
      ncones = setupcells(xcone, n_cones, conexarr, coneyarr, conetharr);

  } else {	/* random cone array */

    if (!notinit(conearrsiz) || !notinit(conexarrsiz)) {  // if rod array size defined
        if (!notinit(conearrsiz)) {
	      conexarrsiz = conearrsiz;
	      coneyarrsiz = conearrsiz;
        }
        if (!notinit(conearrcentx)) {  			// position for the array
             ncones = setupcells(xcone,n_cones,conexarrsiz,coneyarrsiz,conearrcentx,conearrcenty);
        } else {					// else use general position arrcentx,y def=0
             ncones = setupcells(xcone,n_cones,conexarrsiz,coneyarrsiz,arrcentx,arrcenty);
        }
     } else {						// otherwise use default array size and center
          ncones = setupcells(xcone,n_cones);
     }
 }

  setn (xcone,NMADE,ncones);

  find_maxmin(xcone,-1);	/* get size of array, update arraysiz */
  xarrsiz=xmax-xmin;
  yarrsiz=ymax-ymin;
  if (notinit(arrsiz))
    arrsiz=max(xarrsiz,yarrsiz);

};  /* make cones */

if (ninfo>=2 && (make_cones || make_rods)) {
   if (!notinit(arrsiz)) ncfprintf (stderr,"# photoreceptors done, arrsiz=%g\n",arrsiz);
   else ncfprintf (stderr,"# photoreceptors done, xarrsiz=%g yarrsiz=%g\n",xarrsiz,yarrsiz);
}

/*------------------------make bipolar cells---------------------------*/

if (make_rbp) {		 /* make array of rod bipolars */
  if (rbpxarr!=NULL) 
       nrbp = setupcells(rbp, n_rbp, rbpxarr, rbpyarr, rbptharr, rbpnarr);
  else
      nrbp=setupcells(rbp,n_rbp);
}

if (notinit(make_one_dbp1)) make_one_dbp1 = 0;

if (make_one_dbp1) {	 /* make one cone bipolar at specified location */
    double dbp1_xloc, dbp1_yloc;

   if (notinit(gcdistnod)) gcdistnod = 0;
   make_dbp1 = 1;
   dbp1_xloc = ndn(dsgc,1,gcdistnod) -> xloc;
   dbp1_yloc = ndn(dsgc,1,gcdistnod) -> yloc;
   makcell (dbp1, 1, dbp1_xloc, dbp1_yloc, 0, 0, flip );
   ndbp1 = 1;
   setn (dbp1,MAKE,1);
   setn (dbp1,NMADE,ndbp1);
}
else if (make_dbp1) {		 /* make array of on cone bipolars */
  if (!notinit(dbp1arrsiz)) {
      dbp1xarrsiz = dbp1arrsiz;
      dbp1yarrsiz = dbp1arrsiz;
      ndbp1 = setupcells(dbp1,n_dbp1,dbp1xarrsiz,dbp1yarrsiz);
  }
  else if (!notinit(dbp1xarrsiz)) {
      if (notinit(dbp1arrcentx)) {	/* default center (0,0) */ 
           dbp1arrcentx = 0;
           dbp1arrcenty = 0;
           ndbp1 = setupcells(dbp1,n_dbp1,dbp1xarrsiz,dbp1yarrsiz, 
					  dbp1arrcentx, dbp1arrcenty,dbp1_first_cent);
      } else {
           ndbp1 = setupcells(dbp1,n_dbp1,dbp1xarrsiz,dbp1yarrsiz);
      }
  } 
  else if (dbp1xarr!=NULL) {
      ndbp1 = setupcells(dbp1,n_dbp1, dbp1xarr, dbp1yarr, dbp1ytharr, dbp1tharr, dbp1narr);
  }
  else {			/* array of cells with limits undefined, make first cell in center */
      ndbp1=setupcells(dbp1,n_dbp1,dbp1_first_cent=1);
  }
}

if (make_dbp2) {		 /* make array of on cone bipolars */
  if (!notinit(dbp2arrsiz)) {
      dbp2xarrsiz = dbp2arrsiz;
      dbp2yarrsiz = dbp2arrsiz;
      ndbp2 = setupcells(dbp2,n_dbp2,dbp2xarrsiz,dbp2yarrsiz);
  }
  else if (!notinit(dbp2xarrsiz)) {
      if (notinit(dbp2arrcentx)) {	/* default center (0,0) */ 
           dbp2arrcentx = 0;
           dbp2arrcenty = 0;
           ndbp2 = setupcells(dbp2,n_dbp2,dbp2xarrsiz,dbp2yarrsiz, 
					  dbp2arrcentx, dbp2arrcenty, dbp2_first_cent);
      } else {
           ndbp2 = setupcells(dbp2,n_dbp2,dbp2xarrsiz,dbp2yarrsiz);
      }
  } 
  else if (dbp2xarr!=NULL) {
      ndbp2 = setupcells(dbp2,n_dbp2, dbp2xarr, dbp2yarr, dbp2ytharr, dbp2tharr, dbp2narr);
  }
  else {
      ndbp2=setupcells(dbp2,n_dbp2,dbp2_first_cent=1);
  }
}

if (make_dbp3) {		 /* make array of on cone bipolars */
  if (!notinit(dbp3arrsiz)) {
      dbp3xarrsiz = dbp3arrsiz;
      dbp3yarrsiz = dbp3arrsiz;
      ndbp3 = setupcells(dbp3,n_dbp3,dbp3xarrsiz,dbp3yarrsiz);
  }
  else if (!notinit(dbp3xarrsiz)) {
      if (notinit(dbp3arrcentx)) {	/* default center (0,0) */ 
           dbp3arrcentx = 0;
           dbp3arrcenty = 0;
           ndbp3 = setupcells(dbp3,n_dbp3,dbp3xarrsiz,dbp3yarrsiz, 
					  dbp3arrcentx, dbp3arrcenty, dbp3_first_cent);
      } else {
           ndbp3 = setupcells(dbp3,n_dbp3,dbp3xarrsiz,dbp3yarrsiz);
      }
  } 
  else if (dbp3xarr!=NULL) {
      ndbp3 = setupcells(dbp3,n_dbp3, dbp3xarr, dbp3yarr, dbp3ytharr, dbp3tharr, dbp3narr);
  }
  else {
      ndbp3=setupcells(dbp3,n_dbp3,dbp3_first_cent=1);
  }
}

if (make_dbp4) {		 /* make array of on cone bipolars */
  if (!notinit(dbp4arrsiz)) {
      dbp4xarrsiz = dbp4arrsiz;
      dbp4yarrsiz = dbp4arrsiz;
      ndbp4 = setupcells(dbp4,n_dbp4,dbp4xarrsiz,dbp4yarrsiz);
  }
  else if (!notinit(dbp4xarrsiz)) {
      if (notinit(dbp4arrcentx)) {	/* default center (0,0) */ 
           dbp4arrcentx = 0;
           dbp4arrcenty = 0;
           ndbp4 = setupcells(dbp4,n_dbp4,dbp4xarrsiz,dbp4yarrsiz, 
					  dbp4arrcentx, dbp4arrcenty, dbp4_first_cent);
      } else {
           ndbp4 = setupcells(dbp4,n_dbp4,dbp4xarrsiz,dbp4yarrsiz);
      }
  } 
  else if (dbp4xarr!=NULL) {
      ndbp4 = setupcells(dbp4,n_dbp4, dbp4xarr, dbp4yarr, dbp4ytharr, dbp4tharr, dbp4narr);
  }
  else {
      ndbp4=setupcells(dbp4,n_dbp4,dbp4_first_cent=1);
  }
}

if (make_hbp1) {		 /* make array of off cone bipolars */
  if (!notinit(hbp1arrsiz)) {
      hbp1xarrsiz = hbp1arrsiz;
      hbp1yarrsiz = hbp1arrsiz;
      nhbp1 = setupcells(hbp1,n_hbp1,hbp1xarrsiz,hbp1yarrsiz);
  }
  else if (!notinit(hbp1xarrsiz)) {
      if (notinit(hbp1arrcentx)) {	/* default center (0,0) */ 
           hbp1arrcentx = 0;
           hbp1arrcenty = 0;
           nhbp1 = setupcells(hbp1,n_hbp1,hbp1xarrsiz,hbp1yarrsiz, 
					  hbp1arrcentx, hbp1arrcenty,hbp1_first_cent);
      } else {
           nhbp1 = setupcells(hbp1,n_hbp1,hbp1xarrsiz,hbp1yarrsiz);
      }
  } 
  else if (hbp1xarr!=NULL) {
      nhbp1 = setupcells(hbp1,n_hbp1, hbp1xarr, hbp1yarr, hbp1ytharr, hbp1tharr, hbp1narr);
  }
  else {		/* array of cells with limits undefined, make first cell in center */
      nhbp1=setupcells(hbp1,n_hbp1,hbp1_first_cent=1);
  }
}

if (make_hbp2) {		 /* make array of off cone bipolars */
  if (!notinit(hbp2arrsiz)) {
      hbp2xarrsiz = hbp2arrsiz;
      hbp2yarrsiz = hbp2arrsiz;
      nhbp2 = setupcells(hbp2,n_hbp2,hbp2xarrsiz,hbp2yarrsiz);
  }
  else if (!notinit(hbp2xarrsiz)) {
      if (notinit(hbp2arrcentx)) {	/* default center (0,0) */ 
           hbp2arrcentx = 0;
           hbp2arrcenty = 0;
           nhbp2 = setupcells(hbp2,n_hbp2,hbp2xarrsiz,hbp2yarrsiz, 
					  hbp2arrcentx, hbp2arrcenty,hbp2_first_cent);
      } else {
           nhbp2 = setupcells(hbp2,n_hbp2,hbp2xarrsiz,hbp2yarrsiz);
      }
  } 
  else if (hbp2xarr!=NULL) {
      nhbp2 = setupcells(hbp2,n_hbp2, hbp2xarr, hbp2yarr, hbp2ytharr, hbp2tharr, hbp2narr);
  }
  else  {
      nhbp2=setupcells(hbp2,n_hbp2,hbp2_first_cent=1);
  }
}

if (ninfo>=2 && (make_rbp || make_dbp1 || make_dbp2 || make_dbp3 || make_dbp4 || 
			     make_hbp1 || make_hbp2)) 
   ncfprintf (stderr,"# bipolar cells done.\n");

(*addcells)();		/* User defined additions of cells */
			/*  defined in the experiment file */

/*--------------------------------------------------------------------*/

if (ninfo>=2) {
  ncfprintf (stderr,"# Done making neurons.\n");
  ncfprintf (stderr,"#\n");
};

if (ninfo>=2) {
  for (i=0; i<nceltypes; i++) {
    if (int(getn(i,NMADE))) {
        char cbuf[10];
      sprintf (cbuf,"%s",cname[i]);
      ncfprintf (stderr,"# total %-5s = %d\n",cbuf,(int)getn(i,NMADE));
    }
  };
  ncfprintf (stderr,"#\n");
};

// printhash();

find_maxmin(-1,-1);

xarrsiz=xmax-xmin;
yarrsiz=ymax-ymin;
zarrsiz=zmax-zmin;
arrcentx=(xmax+xmin)/2;
arrcenty=(ymax+ymin)/2;
arrcentz=(zmax+zmin)/2;

/*----------------------- make connections ------------------------------*/
int dyadc,make_rcr;

if (make_cone_cone) connect_types (xcone,xcone); /* */

if (make_rod_rbp)   connect_types (xrod,rbp); /* */

if (make_cone_dbp1)  connect_types (xcone,dbp1); /* */

if (make_cone_dbp2)  connect_types (xcone,dbp2); /* */

if (make_cone_dbp3)  connect_types (xcone,dbp3); /* */

if (make_cone_dbp4)  connect_types (xcone,dbp4); /* */

if (make_cone_hbp1)  connect_types (xcone,hbp1); /* */

if (make_cone_hbp2)  connect_types (xcone,hbp2); /* */

if (make_rod_hbat)  connect_types (xrod,hbat); /* */

if (make_cone_ha)   connect_types (xcone,ha, make_ha_cone, make_ha_dbp1); /* */

if (make_cone_hb)   connect_types (xcone,hb); /* */

if (make_ha_dbp1)   connect_types (ha,dbp1); /* */

if (make_hb_dbp1)   connect_types (hb,dbp1); /* */

if (make_ha_dbp2)   connect_types (ha,dbp2); /* */

if (make_hb_dbp2)   connect_types (hb,dbp2); /* */

if (make_ha_dbp3)   connect_types (ha,dbp3); /* */

if (make_hb_dbp3)   connect_types (hb,dbp3); /* */

if (make_ha_dbp4)   connect_types (ha,dbp4); /* */

if (make_hb_dbp4)   connect_types (hb,dbp4); /* */

if (make_ha_ha)     connect_types (ha,ha); /* */

if (make_hb_hb)     connect_types (hb,hb); /* */

if (make_rbp_aii)   connect_types (rbp,aii); /* */

if (make_rbp_a17)   connect_types (rbp,a17,make_a17_rbp,dyadc=0); /* */

if (make_aii_aii)   connect_types (aii,aii); /* */

if (make_dbp1_aii)   connect_types (dbp1,aii); /* */

if (make_dbp2_aii)   connect_types (dbp2,aii); /* */

if (make_dbp3_aii)   connect_types (dbp3,aii); /* */

if (make_dbp4_aii)   connect_types (dbp4,aii); /* */

if (make_aii_hbp1)   connect_types (aii,hbp1); /* */

if (make_dbp1_am)    connect_types (dbp1,am,make_am_dbp1,dyadc=0); /* make reciprocal inhibitory synapse */

if (make_dbp1_am2)   connect_types (dbp1,am2,make_am2_dbp1,dyadc=0); /* */

if (make_dbp1_am3)   connect_types (dbp1,am3,make_am3_dbp1,dyadc=0); /* */

if (make_dbp1_am4)   connect_types (dbp1,am4,make_am4_dbp1,dyadc=0); /* */

if (make_dbp1_ams)   connect_types (dbp1,ams,make_ams_dbp1,dyadc=0); /* */

if (make_dbp2_am)    connect_types (dbp2,am,make_am_dbp2,dyadc=0); /* */

if (make_dbp2_am2)   connect_types (dbp2,am2,make_am2_dbp2,dyadc=0); /* */

if (make_dbp2_am3)   connect_types (dbp2,am3,make_am3_dbp2,dyadc=0); /* */

if (make_dbp2_am4)   connect_types (dbp2,am4,make_am4_dbp2,dyadc=0); /* */

if (make_dbp2_ams)   connect_types (dbp2,ams,make_ams_dbp2,dyadc=0); /* */

if (make_dbp3_am)    connect_types (dbp3,am,make_am_dbp3,dyadc=0); /* */

if (make_dbp3_am2)   connect_types (dbp3,am2,make_am2_dbp3,dyadc=0); /* */

if (make_dbp3_am3)   connect_types (dbp3,am3,make_am3_dbp3,dyadc=0); /* */

if (make_dbp3_am4)   connect_types (dbp3,am4,make_am4_dbp3,dyadc=0); /* */

if (make_dbp3_ams)   connect_types (dbp3,ams,make_ams_dbp3,dyadc=0); /* */

if (make_dbp4_am)    connect_types (dbp4,am,make_am_dbp4,dyadc=0); /* */

if (make_dbp4_am2)   connect_types (dbp4,am2,make_am2_dbp4,dyadc=0); /* */

if (make_dbp4_am3)   connect_types (dbp4,am3,make_am3_dbp4,dyadc=0); /* */

if (make_dbp4_am4)   connect_types (dbp4,am4,make_am4_dbp4,dyadc=0); /* */

// if (make_am2_dbp1)   connect_types (am2,dbp1); 		/* lateral inhibition, no reciprocal feedback */

if (make_hbp1_amh)   connect_types (hbp1,amh,make_amh_hbp1,dyadc=0); /* */

if (make_hbp1_amh2)  connect_types (hbp1,amh2,make_amh2_hbp1,dyadc=0); /* */

if (make_hbp1_amhs)  connect_types (hbp1,amhs,make_amhs_hbp1,dyadc=0); /* */

if (make_dbp1_sbac)  connect_types (dbp1,sbac,make_sbac_dbp1,dyadc=0); /* */

if (make_dbp2_sbac)  connect_types (dbp2,sbac,make_sbac_dbp2,dyadc=0); /* */

if (make_dbp3_sbac)  connect_types (dbp3,sbac,make_sbac_dbp3,dyadc=0); /* */

if (make_dbp4_sbac)  connect_types (dbp4,sbac,make_sbac_dbp4,dyadc=0); /* */

if (make_dbp1_dbp1)  connect_types (dbp1,dbp1); /* */

if (make_dbp2_dbp2)  connect_types (dbp2,dbp2); /* */

if (make_dbp3_dbp3)  connect_types (dbp3,dbp3); /* */

if (make_dbp4_dbp4)  connect_types (dbp4,dbp4); /* */

if (make_hbp1_hbp1)  connect_types (hbp1,hbp1); /* */

if (make_hbp2_amh)   connect_types (hbp2,amh,make_amh_hbp2,dyadc=0); /* */

if (make_hbp2_amh2)  connect_types (hbp2,amh2,make_amh2_hbp2,dyadc=0); /* */

if (make_dbp1_gca)   connect_types (dbp1,gca); /* */

if (make_dbp2_gca)   connect_types (dbp2,gca); /* */

if (make_dbp3_gca)   connect_types (dbp3,gca); /* */

if (make_dbp4_gca)   connect_types (dbp4,gca); /* */

if (make_dbp1_gcb)   connect_types (dbp1,gcb); /* */

if (make_dbp2_gcb)   connect_types (dbp2,gcb); /* */

if (make_dbp3_gcb)   connect_types (dbp3,gcb); /* */

if (make_dbp4_gcb)   connect_types (dbp4,gcb); /* */

if (make_sbac_dsgc)  connect_types (sbac,dsgc); /* */

if (make_sbac_sbac)  connect_types (sbac,sbac,make_rcr=1,dyadc=0);   /* */

if (make_dbp1_dsgc)  connect_types (dbp1,dsgc); /* */

if (make_dbp2_dsgc)  connect_types (dbp2,dsgc); /* */

if (make_dbp3_dsgc)  connect_types (dbp3,dsgc); /* */

if (make_dbp4_dsgc)  connect_types (dbp4,dsgc); /* */

if (make_hbp1_dsgc)  connect_types (hbp1,dsgc); /* */

if (make_ams_dsgc)   connect_types (ams,dsgc); /* */

if (make_amhs_dsgc)  connect_types (amhs,dsgc); /* */

if (make_amh2_amh)   connect_types (amh2,amh); /* */

if (make_dsgc_dsgc)  connect_types (dsgc,dsgc); /* */

if (make_am_am)      connect_types (am,am); /* */

if (make_am2_am2)    connect_types (am2,am2); /* */

if (make_am3_am3)    connect_types (am3,am3); /* */

if (make_am4_am4)    connect_types (am4,am4); /* */

if (make_am_gca)     connect_types (am,gca); /* */

if (make_am2_gca)    connect_types (am2,gca); /* */

if (make_am3_gca)    connect_types (am3,gca); /* */

if (make_am4_gca)    connect_types (am4,gca); /* */

if (make_ams_gca)    connect_types (ams,gca); /* */

if (make_ams_gcb)    connect_types (ams,gcb); /* */

if (make_ams_gcaoff) connect_types (ams,gcaoff); /* */

if (make_ams_gcboff) connect_types (ams,gcboff); /* */

if (make_hbp1_gcaoff) connect_types (hbp1,gcaoff); /* */

if (make_hbp1_gcboff) connect_types (hbp1,gcboff); /* */

if (make_hbp2_gcaoff) connect_types (hbp2,gcaoff); /* */

if (make_hbp2_gcboff) connect_types (hbp2,gcboff); /* */

/*  -  -  -  -  -  -  -  print connections -  -  -  -  -  -  -  -  -  */

if (ninfo >= 2) {
  if (notinit(print_conns)) print_conns = 0;
  if (print_conns) {
    for (ct=0; ct<nceltypes; ct++) {
      print_connections(ct); 
    }
  }

  if (notinit(print_avg_conns)) print_avg_conns = 1;
  if (print_avg_conns) {
    for (ct=0; ct<nceltypes; ct++) {
      print_avg_connections(ct); 
    }
    ncfprintf (stderr,"# \n");
  }
}

/*----------------------- prune connections -----------------------------*/

if (notinit(remove_nconns)) remove_nconns = 1;

if (remove_nconns) {		/* if remove non-connected neurons */

 /* remove cells without any inputs or outputs */

if (ninfo>=2)
     ncfprintf(stderr,"#\n# Removing neurons that don't connect.\n");

/* Start with bipolar cells */
/* Rationale here is to remove the cells that haven't connected to their
   pre- or post-synaptic cells. First, check small cells that connect to 
   large cells because the number of large cells is limited and this will
   limit the number of small cells. Next, check the cells dependent on those
   removed, proceeding backwards towards the photoreceptors. Ignore 
   gap junctions to cells of same type, but check gap junctions to cells of
   different type.  "made_gc_comps" is selective compartments near one bp cell.
*/

 if (getn(dbp1,NMADE)>0) {		/* check on-bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[dbp1]);
   if (make_cones && make_cone_dbp1) checkcellin(dbp1,xcone);
   if (make_gca && make_dbp1_gca && !made_gc_comps) {
          // if (!make_aii) checkcellout(dbp1,gca);
          checkcellout(dbp1);
	  if (make_aii && make_dbp1_aii) {
		  checkcellout(dbp1,aii,dbp1,gca);
		  checkcellout(aii,dbp1,gca);
	  }
   }
   else if (make_aii && make_dbp1_aii) checkcellout(dbp1,aii);
   else if ((make_am  && make_dbp1_am)  || 
            (make_am2 && make_dbp1_am2) ||
            (make_am3 && make_dbp1_am3) || 
            (make_am4 && make_dbp1_am4)) checkcellout(dbp1);

   if (make_gcb && make_dbp1_gcb) {
          checkcellout(dbp1,gcb);
   }
   if (make_dsgc & make_dbp1_dsgc) {
      if (!make_sbac || !make_dbp1_sbac) checkcellout(dbp1,dsgc); 
      if (make_sbac & make_dbp1_sbac) checkcellout(dbp1,sbac,dsgc);
   }
   else if (make_sbac & make_dbp1_sbac) checkcellout(dbp1,sbac);
   if (ninfo >= 2) {
	if (ncell_erased[dbp1]>0) ncfprintf (stderr,"\n# dbp1s ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[dbp1]);
   }
   setn(dbp1,NMADE,getn(dbp1,NMADE)-ncell_erased[dbp1]);
 }

 if (getn(dbp2,NMADE)>0) {		/* check on-bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[dbp2]);
   if (make_cones && make_cone_dbp2) checkcellin(dbp2,xcone);
   if (make_gca && make_dbp2_gca) {
          checkcellout(dbp2);
	  if (make_aii && make_dbp2_aii) {
		  checkcellout(dbp2,aii,dbp2,gca);
		  checkcellout(aii,dbp2,gca);
	  }
   }
   else if (make_aii && make_dbp2_aii) checkcellout(dbp2,aii);
   else if ((make_am  && make_dbp2_am)  || 
            (make_am2 && make_dbp2_am2) ||
            (make_am3 && make_dbp2_am3) ||
            (make_am4 && make_dbp2_am4)) checkcellout(dbp2);

   if (make_gcb && make_dbp2_gcb) {
          checkcellout(dbp2,gcb);
   }
   if (make_dsgc & make_dbp2_dsgc) {
      if (!make_sbac || !make_dbp2_sbac) checkcellout(dbp2,dsgc); 
      if (make_sbac & make_dbp2_sbac) checkcellout(dbp2,sbac,dsgc);
   }
   else if (make_sbac & make_dbp2_sbac) checkcellout(dbp2,sbac);
   if (ninfo >= 2) {
	if (ncell_erased[dbp2]>0) ncfprintf (stderr,"\n# dbp2s ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[dbp2]);
   }
   setn(dbp2,NMADE,getn(dbp2,NMADE)-ncell_erased[dbp2]);
 }
 
 if (getn(dbp3,NMADE)>0) {		/* check on-bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[dbp3]);
   if (make_cones && make_cone_dbp3) checkcellin(dbp3,xcone);
   if (make_gca && make_dbp3_gca) {
          checkcellout(dbp3);
	  if (make_aii && make_dbp3_aii) {
		  checkcellout(dbp3,aii,dbp3,gca);
		  checkcellout(aii,dbp3,gca);
	  }
   }
   else if (make_aii && make_dbp3_aii) checkcellout(dbp3,aii);
   else if ((make_am  && make_dbp3_am)  || 
            (make_am2 && make_dbp3_am2) ||
            (make_am3 && make_dbp3_am3) ||
            (make_am4 && make_dbp3_am4)) checkcellout(dbp3);

   if (make_gcb && make_dbp3_gcb) {
          checkcellout(dbp3,gcb);
   }
   if (make_dsgc & make_dbp3_dsgc) {
      if (!make_sbac || !make_dbp3_sbac) checkcellout(dbp3,dsgc); 
      if (make_sbac & make_dbp3_sbac) checkcellout(dbp3,sbac,dsgc);
   }
   else if (make_sbac & make_dbp3_sbac) checkcellout(dbp3,sbac);
   if (ninfo >= 2) {
	if (ncell_erased[dbp3]>0) ncfprintf (stderr,"\n# dbp3s ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[dbp3]);
   }
   setn(dbp3,NMADE,getn(dbp3,NMADE)-ncell_erased[dbp3]);
 }
 
 if (getn(dbp4,NMADE)>0) {		/* check on-bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[dbp4]);
   if (make_cones && make_cone_dbp4) checkcellin(dbp4,xcone);
   if (make_gca && make_dbp4_gca) {
          checkcellout(dbp4);
	  if (make_aii && make_dbp4_aii) {
		  checkcellout(dbp4,aii,dbp4,gca);
		  checkcellout(aii,dbp4,gca);
	  }
   }
   else if (make_aii && make_dbp4_aii) checkcellout(dbp4,aii);
   else if ((make_am  && make_dbp4_am)  || 
            (make_am2 && make_dbp4_am2) ||
            (make_am3 && make_dbp4_am3) ||
            (make_am4 && make_dbp4_am4)) checkcellout(dbp4);

   if (make_gcb && make_dbp4_gcb) {
          checkcellout(dbp4,gcb);
   }
   if (make_dsgc & make_dbp4_dsgc) {
      if (!make_sbac || !make_dbp4_sbac) checkcellout(dbp4,dsgc); 
      if (make_sbac & make_dbp4_sbac) checkcellout(dbp4,sbac,dsgc);
   }
   else if (make_sbac & make_dbp4_sbac) checkcellout(dbp4,sbac);
   if (ninfo >= 2) {
	if (ncell_erased[dbp4]>0) ncfprintf (stderr,"\n# dbp4s ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[dbp4]);
   }
   setn(dbp4,NMADE,getn(dbp4,NMADE)-ncell_erased[dbp4]);
 }
 
 if (getn(hbp1,NMADE)>0) {		/* check off-bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[hbp1]);
   if (make_cones && make_cone_hbp1) checkcellin(hbp1);
   if (make_gcaoff && make_hbp1_gcaoff) {
	  if (make_aii && make_aii_hbp1) {
		  checkcellout(aii,hbp1,gcaoff);
	  }
          if (!make_aii) checkcellout(hbp1,gcaoff);
   }
   if (make_amh  && make_hbp1_amh)  checkcellout(hbp1);
   if (make_amh2 && make_hbp1_amh2) checkcellout(hbp1);

   if (make_gcboff && make_hbp1_gcboff) checkcellout(hbp1,gcboff);

   if (make_dsgc && make_hbp1_dsgc) {
      if (!make_sbac || !make_hbp1_sbac) checkcellout(hbp1,dsgc);
      if (make_sbac & make_hbp1_sbac) checkcellout(hbp1,sbac,dsgc);
   }
   else if (make_sbac & make_hbp1_sbac) checkcellout(hbp1,sbac);
   if (ninfo >= 2) {
	if (ncell_erased[hbp1]>0) ncfprintf (stderr,"\n# hbp1s ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[hbp1]);
   }
   setn(hbp1,NMADE,getn(hbp1,NMADE)-ncell_erased[hbp1]);
 }

 if (getn(hbp2,NMADE)>0) {		/* check off-bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[hbp2]);
   if (make_cones && make_cone_hbp2) checkcellin(hbp2);
   if (make_gcaoff && make_hbp2_gcaoff) {
	  if (make_aii && make_aii_hbp2) {
		  checkcellout(aii,hbp2,gcaoff);
	  }
          if (!make_aii) checkcellout(hbp2,gcaoff);
   }
   if (make_amh  && make_hbp2_amh) checkcellout(hbp2);
   if (make_amh2 && make_hbp2_amh2) checkcellout(hbp2);

   if (make_gcboff && make_hbp2_gcboff) checkcellout(hbp2,gcboff);

   if (make_dsgc && make_hbp2_dsgc) {
      if (!make_sbac || !make_hbp2_sbac) checkcellout(hbp2,dsgc);
      if (make_sbac & make_hbp2_sbac) checkcellout(hbp2,sbac,dsgc);
   }
   if (ninfo >= 2) {
	if (ncell_erased[hbp2]>0) ncfprintf (stderr,"\n# hbp2s ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[hbp2]);
   }
   setn(hbp1,NMADE,getn(hbp1,NMADE)-ncell_erased[hbp2]);
 }

/* next check amacrine cells */

 if (getn(ams,NMADE)>0) {		/* check small-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %s's ",cname[ams]);
   if ((make_gcaoff && make_ams_gcaoff) || (make_gca && make_ams_gca) || 
		   (make_dsgc && make_ams_dsgc)) checkcellout(ams);
   if (make_cones) checkcellin(ams);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[ams]);
   setn(ams,NMADE,getn(ams,NMADE)-ncell_erased[ams]);
 }

 if (getn(amhs,NMADE)>0) {		/* check small-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %s's ",cname[amhs]);
   if (ninfo>=2) ncfprintf (stderr,"# amhs ");
   if ((make_gcaoff && make_amhs_gcaoff) || (make_gcboff && make_amhs_gcboff) || 
		   (make_dsgc && make_amhs_dsgc)) checkcellout(amhs);
   if (make_cones) checkcellin(amhs);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[amhs]);
   setn(amhs,NMADE,getn(amhs,NMADE)-ncell_erased[amhs]);
 }

 if (getn(amh,NMADE)>0) {		/* check wide-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[amh]);
   //if (make_dbp1 && make_dbp1_amh) checkcellin(amh);
   if ((make_gcaoff && make_amh_gcaoff) || (make_gcboff && make_amh_gcboff) || 
		   (make_dsgc && make_amh_dsgc)) checkcellout(am);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[amh]);
   setn(amh,NMADE,getn(amh,NMADE)-ncell_erased[amh]);
 }
 if (getn(am,NMADE)>0) {		/* check wide-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[am]);
   //if (make_dbp1 && make_dbp1_am) checkcellin(am);
   if ((make_gcaoff && make_am_gcaoff) || (make_gca && make_am_gca) || 
		   (make_dsgc && make_am_dsgc)) checkcellout(am);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[am]);
   setn(am,NMADE,getn(am,NMADE)-ncell_erased[am]);
 }
 if (getn(am2,NMADE)>0) {		/* check wide-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[am2]);
   if (make_dbp1 && make_dbp1_am2) checkcellin(am2);
   if (make_dbp2 && make_dbp2_am2) checkcellin(am2);
   // if ((make_gcaoff && make_am2_gcaoff) || (make_gca && make_am2_gca) || 
   // 		   (make_dsgc && make_am2_dsgc)) checkcellout(am2);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[am2]);
   setn(am2,NMADE,getn(am2,NMADE)-ncell_erased[am2]);
 }
 if (getn(am3,NMADE)>0) {		/* check wide-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[am3]);
   if (make_dbp1 && make_dbp1_am3) checkcellin(am3);
   if (make_dbp2 && make_dbp2_am3) checkcellin(am3);
   // if ((make_gcaoff && make_am3_gcaoff) || (make_gca && make_am3_gca) || 
   // 		   (make_dsgc && make_am3_dsgc)) checkcellout(am3);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[am3]);
   setn(am3,NMADE,getn(am3,NMADE)-ncell_erased[am3]);
 }
 if (getn(am4,NMADE)>0) {		/* check wide-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[am4]);
   if (make_dbp1 && make_dbp1_am4) checkcellin(am4);
   if (make_dbp2 && make_dbp2_am4) checkcellin(am4);
   // if ((make_gcaoff && make_am4_gcaoff) || (make_gca && make_am4_gca) || 
   // 		   (make_dsgc && make_am4_dsgc)) checkcellout(am4);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[am4]);
   setn(am4,NMADE,getn(am4,NMADE)-ncell_erased[am4]);
 }
 if (getn(amh2,NMADE)>0) {		/* check wide-field amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[amh2]);
   //if (make_dbp2 && make_dbp2_amh2) checkcellin(amh2);
   // if ((make_gcaoff && make_amh2_gcaoff) || (make_gca && make_amh2_gca) || 
   // 		   (make_dsgc && make_amh2_dsgc)) checkcellout(amh2);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[amh2]);
   setn(amh2,NMADE,getn(amh2,NMADE)-ncell_erased[amh2]);
 }

 if (getn(sbac,NMADE)>0) {			/* check starburst connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[sbac]);
   if (make_dbp1 && make_dbp1_sbac) checkcellin(sbac);
   if (make_dsgc && make_sbac_dsgc && !make_sbac_sbac) checkcellout(sbac);
   if (make_dsgc && make_sbac_dsgc &&  make_sbac_sbac) checkcellout(sbac,sbac,dsgc);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[sbac]);
   setn(sbac,NMADE,getn(sbac,NMADE)-ncell_erased[sbac]);
 }

 if (getn(aii,NMADE)>0) {		/* check aii amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[aii]);
   if ((make_dbp1 && make_dbp1_aii) || (make_rbp && make_rbp_aii)) checkcellin(aii);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[aii]);
   setn(aii,NMADE,getn(aii,NMADE)-ncell_erased[aii]);
 }
 if (getn(a17,NMADE)>0) {		/* check a17 amacrine connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[a17]);
   if ((make_rbp && make_rbp_a17)) checkcellin(a17);
   if ((make_rbp && make_a17_rbp)) checkcellout(a17);
   if (ninfo >= 2) ncfprintf (stderr," %d erased\n",ncell_erased[a17]);
   setn(a17,NMADE,getn(a17,NMADE)-ncell_erased[a17]);
 }

/* next check rod bipolars */
 
 if (getn(rbp,NMADE)>0) {		/* check rod bipolar connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[rbp]);
   if (make_aii && make_rbp_aii) checkcellout(rbp);
   if (make_a17 && make_rbp_a17) checkcellout(rbp);
   if (make_a17 && make_a17_rbp)  checkcellin(rbp);
   if (make_rods && make_rod_rbp) checkcellin(rbp);
   if (ninfo >= 2) {
	if (ncell_erased[rbp]>0) ncfprintf (stderr,"\n# rbps ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[rbp]);
   }
   setn(rbp,NMADE,getn(rbp,NMADE)-ncell_erased[rbp]);
 }
 
  /* then check has */

 if (getn(ha,NMADE)>0) {		/* check ha connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[ha]);
   if (make_cones && make_dbp1 && make_cone_ha && make_cone_dbp1) checkcellout(ha,xcone,dbp1);
   if (make_cones && make_hbp1 && make_cone_ha && make_cone_hbp1) checkcellout(ha,xcone,hbp1);
   if (ninfo >= 2) ncfprintf (stderr,"# has erased %d\n",ncell_erased[ha]);
   setn(ha,NMADE,getn(ha,NMADE)-ncell_erased[ha]);
 }
 
  /* then check hbs */

 if (getn(hb,NMADE)>0) {		/* check hb connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[hb]);
   if (make_cones && make_dbp1 && make_cone_hb && make_cone_dbp1) checkcellout(hb,xcone,dbp1);
   if (make_cones && make_hbp1 && make_cone_hb && make_cone_hbp1) checkcellout(hb,xcone,hbp1);
   if (ninfo >= 2) ncfprintf (stderr,"# hbs erased %d\n",ncell_erased[hb]);
   setn(hb,NMADE,getn(hb,NMADE)-ncell_erased[hb]);
 }
 
  /* then go back to check cones */

 if (getn(xcone,NMADE)>0) {		/* check cone connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[xcone]);
   if ((make_dbp1 && make_cone_dbp1) || 
       (make_dbp2 && make_cone_dbp2) || 
       (make_dbp3 && make_cone_dbp3) || 
       (make_dbp4 && make_cone_dbp4) || 
       (make_hbp1 && make_cone_hbp1) || 
       (make_ha && make_cone_ha) || 
       (make_hb && make_cone_hb)) checkcellout(xcone);
   if (ninfo >= 2) {
	if (ncell_erased[xcone]>0) ncfprintf (stderr,"\n# cones ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[xcone]);
   }
   setn(xcone,NMADE,getn(xcone,NMADE)-ncell_erased[xcone]);
 }

 if (getn(xrod,NMADE)>0) {		/* check rod connections */
   if (ninfo>=2) ncfprintf (stderr,"# %ss ",cname[xrod]);
   if (make_rbp && make_rod_rbp) checkcellout(xrod);
   if (ninfo >= 2) {
	if (ncell_erased[xrod]>0) ncfprintf (stderr,"\n# rods ");
   	ncfprintf (stderr,"%d erased\n",ncell_erased[xrod]);
   }
   setn(xrod,NMADE,getn(xrod,NMADE)-ncell_erased[xrod]);
 }

if (ninfo>=2) ncfprintf (stderr,"#\n");

  if (ninfo>=2) {
    for (i=0; i<nceltypes; i++) {
      if (int(getn(i,NMADE))) {
          char cbuf[10];
        sprintf (cbuf,"%s",cname[i]);
        ncfprintf (stderr,"# total %-5s = %d\n",cbuf,(int)getn(i,NMADE));
      }
    }
    ncfprintf (stderr,"#\n");
  }

}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* remove cellparts that have grown outside of the array */

if (notinit(disp_margin)) disp_margin = 0.07;		/* size of margin */

if (setarrsiz) {
  arrcentx = 0;
  arrcenty = 0;
  xmax = arrcentx +  arrsiz / 2 * (1 + 1*disp_margin);
  ymax = arrcenty +  arrsiz / 2 * (1 + 1*disp_margin);
  xmin = arrcentx -  arrsiz / 2 * (1 + 1*disp_margin);
  ymin = arrcenty -  arrsiz / 2 * (1 + 1*disp_margin);
}
else {
  xmax = arrcentx +  xarrsiz / 2 * (1 + 1*disp_margin);
  ymax = arrcenty +  yarrsiz / 2 * (1 + 1*disp_margin);
  xmin = arrcentx -  xarrsiz / 2 * (1 + 1*disp_margin);
  ymin = arrcenty -  yarrsiz / 2 * (1 + 1*disp_margin);
}

if (notinit(limit_array)) if (disp > 0) limit_array = 1;  /* always do when displaying */
if (notinit(limit_array)) limit_array = 0; 		  /* remove cell parts outside of arr */
// limit_array = 0;					  /* always do when displaying */

// if (ninfo>=2) {
//    ncfprintf (stderr,"# Limiting neurons to array size.\n");
// }

if (limit_array) {
    elem *epnt,*tepnt;

  elimit (X,xmax,xmin);				// set up windowing limits
  elimit (Y,ymax,ymin);
  if (!notinit(set_zmax) && !notinit(set_zmin)) {
    elimit (Z,set_zmax,set_zmin); 
  }
  for (epnt=elempnt; epnt=foreach(epnt,ELEMENT); epnt=tepnt) {
      tepnt = epnt->next;
      elimit (epnt);
  } 
}

if (ninfo >= 3) {
  ncfprintf (stderr,"# Final array parameters:\n");
  ncfprintf (stderr,"# xcent %g ycent %g\n",arrcentx,arrcenty);
  ncfprintf (stderr,"# xsize %g ysize %g\n",xarrsiz,yarrsiz);
}

if (ninfo>=2) {
  ncfprintf (stderr,"# Done connecting neurons.\n#\n");
}

if (ninfo >= 2 && remove_nconns) {
  if (notinit(print_avg_conns)) print_avg_conns = 1;
  if (print_avg_conns) {
    for (ct=0; ct<nceltypes; ct++) {
      print_avg_connections(ct); 
    };
    ncfprintf (stderr,"# \n");
  }
}

(*addlabels)();		/* User defined node labels */
			/*  defined in the experiment file */
 
/*--------------------------display------------------------------------*/

if (disp & (DISP|DCOMP|DCONN|DNODE)) {

find_maxmin(-1,-1);

xarrsiz=xmax-xmin;
yarrsiz=ymax-ymin;
zarrsiz=zmax-zmin;
arrcentx=(xmax+xmin)/2;
arrcenty=(ymax+ymin)/2;
arrcentz=(zmax+zmin)/2;

/* add extra margin for displaying scale bar */

//zarrsiz = zarrsiz * (-cos(mxrot/180*PI)-sin(mxrot/180*PI)) *
//		    (cos(myrot/180*PI)+sin(myrot/180*PI));
maxsize = max(xarrsiz,yarrsiz);
maxsize = max(maxsize,zarrsiz);
if (notinit(dispsize)) dispsize = maxsize; 
// dispsize = max(dispsize,maxsize) * (1+disp_margin);
dispsize *=  (1+disp_margin);

//fprintf (stderr,"dispsize %g\n",dispsize);
//fprintf (stderr,"xarrsiz %g yarrsiz %g zarrsiz %g\n",xarrsiz,yarrsiz,zarrsiz);
//fprintf (stderr,"arrcentx %g arrcenty %g\n",arrcentx,arrcenty);

set_disp_rot (mxrot,myrot,0,arrcentx,arrcenty,arrcentz,0,0,0,dispsize); 

if (notinit(disp_calib_len)) {
   disp_calib_len = dispsize/5; 
   disp_calib_len = set_int_val(disp_calib_len);
}

if (disp_calib_len > 0) disp_calib (0.97,0.05,disp_calib_len,dispsize,white);

setcmap(1);

   /* if set, use disp_zmax, disp_zmin for display zmax, zmin */

if (!notinit(disp_zmin) || !notinit(disp_zmax)) {
     dzmin = disp_zmin;
     dzmax = disp_zmax;
}

if (make_rods) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(rod_nscale)) nscal = rod_nscale;
   dispcelltype(xrod,1, nrods, scal=0.9, nscal);
}

if (make_cones) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(cone_nscale)) nscal = cone_nscale;
   dispcelltype(xcone,1,ncones,scal=0.8,nscal);
}

if (make_hbat) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(hbat_nscale)) nscal = hbat_nscale;
   dispcelltype(hbat,1,nhbat,scal=1,nscal);
}

if (make_ha) {
   nscal = -3.15;					/* nodenm => 3 => nd, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(ha_nscale))  nscal = ha_nscale;
   dispcelltype(ha,1,nha,scal=1,nscal);
}

if (make_hb) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(hb_nscale))  nscal = hb_nscale;
   dispcelltype(hb,1,nhb,scal=1,nscal);
}

if (make_rbp) {
   nscal = -2.4;					/* nodenm => 2 => cn, 0.4 => medium */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(rbp_nscale)) nscal = rbp_nscale;
   dispcelltype(rbp,1,nrbp,scal=1,nscal);
}

if (make_dbp1) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(dbp1_nscale)) nscal = dbp1_nscale;
   dispcelltype(dbp1,1,ndbp1,scal=1,nscal);	
}

if (make_dbp2) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(dbp2_nscale)) nscal = dbp2_nscale;
   dispcelltype(dbp2,1,ndbp2,scal=1,nscal);
}

if (make_dbp3) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(dbp3_nscale)) nscal = dbp3_nscale;
   dispcelltype(dbp3,1,ndbp3,scal=1,nscal);
}

if (make_dbp4) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(dbp4_nscale)) nscal = dbp4_nscale;
   dispcelltype(dbp4,1,ndbp4,scal=1,nscal);
}

if (make_hbp1) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(hbp1_nscale)) nscal = hbp1_nscale;
   dispcelltype(hbp1,1,nhbp1,scal=1,nscal);
}

if (make_hbp2) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(hbp2_nscale)) nscal = hbp2_nscale;
   dispcelltype(hbp2,1,nhbp2,scal=1,nscal);
}

if (make_aii) {
   nscal = -3.2;					/* nodenm => 3 => nd, 0.2 => small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(aii_nscale))  nscal = aii_nscale;
   dispcelltype(aii,1,naii,scal=1,nscal);
}

if (make_sbac) {
   nscal = -3.03;					/* nodenm => 3 => nd, 0.03 => very small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(sbac_nscale)) nscal = sbac_nscale;
   //dispcelltype(sbac,1,nsbac,scal=1,nscal);
   dispcelltype(sbac,1,100,scal=1,nscal, dzmax, dzmin);
}

if (make_am) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(am_nscale))  nscal = am_nscale;
   dispcelltype(am,1,nam,scal=1,nscal);
}

if (make_am2) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(am2_nscale)) nscal = am2_nscale;
   dispcelltype(am2,1,nam2,scal=1,nscal);
}

if (make_am3) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(am3_nscale)) nscal = am3_nscale;
   dispcelltype(am3,1,nam3,scal=1,nscal);
}

if (make_am4) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(am4_nscale)) nscal = am4_nscale;
   dispcelltype(am4,1,nam4,scal=1,nscal);
}

if (make_amh) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(amh_nscale)) nscal = amh_nscale;
   dispcelltype(amh,1,namh,scal=1,nscal);
}

if (make_amh2) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(amh2_nscale)) nscal = amh2_nscale;
   dispcelltype(amh2,1,namh2,scal=1,nscal);
}

if (make_ams) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(ams_nscale)) nscal = ams_nscale;
   dispcelltype(ams,1,nams,scal=1,nscal);
}

if (make_amhs) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(amhs_nscale)) nscal = amhs_nscale;
   dispcelltype(amhs,1,namhs,scal=1,nscal);
}

//setn(a17,NCOLOR,0);
if (make_a17) {
   nscal = -2.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(a17_nscale)) nscal = a17_nscale;
   dispcelltype(a17,1,na17,scal=1,nscal);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* if set, use disp_gc_zmax, disp_gc_zmin for ganglion cell display zmax, zmin */

if (!notinit(disp_gc_zmin) || !notinit(disp_gc_zmax)) { 
     dzmin = disp_gc_zmin;
     dzmax = disp_gc_zmax;
}

if (make_gca) {
   nscal = -3.09;					/* nodenm => 3 => nd, 0.09 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(gca_nscale)) nscal = gca_nscale;
   dispcelltype(gca,1,ngca,scal=1,nscal, dzmax, dzmin);
}

if (make_gcb) {
   nscal = -3.09;					/* nodenm => 3 => nd, 0.1 => small */
   if (!notinit(node_scale)) nscal = node_scale;
   if (!notinit(gcb_nscale)) nscal = gcb_nscale;
   dispcelltype(gcb,1,ngcb,scal=1,nscal,dzmax, dzmin);
}

if (make_dsgc) {
   nscal = -3.04;					/* nodenm => 3 => nd, 0.04 => very small */
   if (!notinit(node_scale))  nscal = node_scale;
   if (!notinit(dsgc_nscale)) nscal = dsgc_nscale;
   if (!notinit(disp_dsgc_zmin) || !notinit(disp_dsgc_zmax)) {
	dzmin = disp_dsgc_zmin;
	dzmax = disp_dsgc_zmax;
   }
   // dispcelltype(dsgc,1,ndsgc,scal=1,nscal); 
   dispcelltype(dsgc,1,ndsgc,scal=1,nscal,dzmax,dzmin); 
}

if (make_gcaoff) {
   nscal = -3.1;					/* nodenm => 3 => nd, 0.1 => small */
   if (!notinit(node_scale))    nscal = node_scale;
   if (!notinit(gcaoff_nscale)) nscal = gcaoff_nscale;
   dispcelltype(gcaoff,1,ngcaoff,scal=1,nscal,dzmax,dzmin);
}

if (make_gcboff) {
   nscal = -3.1;					/* nodenm => 2 => cn, 0.1 => small */
   if (!notinit(node_scale))    nscal = node_scale;
   if (!notinit(gcboff_nscale)) nscal = gcboff_nscale;
   dispcelltype(gcboff,1,ngcboff,scal=1,nscal,dzmax,dzmin);
}

// if (info_disp>=1) {
//   print_sb_out_syns();
// }

/* display starburst cells */

/*
 if (make_sbac) {
  double colval;
  int i, q;

  for (i=1; i<=nsbac; i++) {
    colval=getn(sbac,NCOLOR)-1+i;
    
    display(ELEMENT, MATCHING, ndt(sbac,i,-1), (int)colval, 1);
    display(COMPS, MATCHING, ndt(sbac,1,-1), (int)colval, 1);
    display(NODE, MATCHING, ndt(sbac,1,-1), magenta, -3.1);

    // show sbac synapses, syns from each sbcell different color
    if (make_dsgc) display_sb_out_syns(i);
   }
}
/* */

/*------------------------------------------------------------*/

synfuncs_cleanup();		// remove cellin, cellout, soma arrays

if (notinit(plot_freq)) {
  if (make_dsgc && getn(dsgc,BIOPHYS))
    plot_freq = 1;			/* plot spike freq of ganglion cells */
  else 
    plot_freq = 0;		/* otherwise no spikes, so don't plot */
}
 infile = expt;
 draw_plotlabel(dispsize);	/* draw label or expt in upper right hand corner */

 if (disp & 15) exit(0);	/* exit if disp */

} /* end display */

/*------------------------expt-------------------------------------*/

  (*runexpt)(); /* Run the experiment defined in the experiment file */
}


