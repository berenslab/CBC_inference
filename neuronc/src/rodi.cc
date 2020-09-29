/* Module rodi in Program nc */

/* sets up mouse rod from Invergo et al., (2014)

/*  Invergo BM, Dell'Orco D, Montanucci L, Koch KW, Bertranpetit J.
    A comprehensive model of the phototransduction cascade in mouse rod cells.
    Mol Biosyst. 2014 Jun;10(6):1481-9. doi: 10.1039/c3mb70584f
    See supplementary electronic info file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "control.h"
#include "ncsub.h"
#include "ncomp.h"
#include "ndef.h"
#include "drand.h"

// #include "wave.h"
extern recpar *rectypes[TOTREC];
extern double specdens[TOTREC][PIGMSIZ];
extern double pigmlen[TOTREC];

/*------------------------------------*/

void rodiset(int pigm)

/* Set up constants for Invergo et al (2014) mouse rod phototransuction */
/* Parameters taken (Dec 2017) from: https://www.ebi.ac.uk/biomodels-main/BIOMD0000000578 */

/* Called by "initrec()" in ncmak.cc. */

{
    int n;
    double rrate, qrate, sfactor;
    recpari *cpt;
    phvarsi *vpt;

  if (pigm>=NUMREC) pigm = NUMREC-1;
  if (pigm==RODI) {
      cpt = (recpari*) rectypes[pigm];	/* set up cone cascade constants */
      cpt->stype = PINVG;	// invergo receptor type(s)
  }
  else return;

  cpt = (recpari*)rectypes[pigm];  /* set up rodi cascade constants */
				/* last pigm type (pigm=NUMREC) is ksens */
  cpt->ctype = ROD;

  cpt->pigm = pigm;
  cpt->recslow = 1.0;		/* cone is slow by this factor */
  cpt->loopgain = 1;		/* gain for ca loop, sets stability */
  // cpt->vrev  =  0.0;		// reversal pot for TRPM1 channel
  cpt->vrev  =  0.72;		// reversal pot for TRPM1 channel

  cpt->qc = dqrec;		/* Q10 factor for conductance shape */
  cpt->qcond = dqcrec;		/* Q10 factor for conductance cond */
				/* not currently used */

  cpt->qt = 37.0;		/* base temperature for channel rate & cond */
  cpt->oldtempcel = tempcel;
  qrate = exp(log(cpt->qc) * (tempcel - cpt->qt) / 10.0);

  sfactor = 1.1;		/* speedup of cascade & intensity */
  rrate = stiminc * qrate /sfactor;    /* stiminc = timinc for syn/stimuli */

  cpt->pathl = pigmlen[pigm];	/* path length through o.s. */

  /* Signal amplification-related model parameter values and sources */

  cpt->Rtot 	= 1.0e8; 	/* total Rhodopsin */
  cpt->omGt 	= 0.6; 		/* expon rate of decay of Gt affinity for R* */
  cpt->kG1_0 	= 1e-3;		/* binding rate of Gt to unphosph R* */
  cpt->kG2	= 2200;		/* dissoc rate of Rn_Gt complex */
  cpt->kG3	= 8500;		/* dissoc rate of Gdp from Rn_Gt complex */
  cpt->kG4gdp	= 400;		/* assoc rate of Gdp to Rn_Gt complex */
  cpt->kG5gtp	= 3500;		/* assoc rate of Gtp to Rn_Gt complex */
  cpt->kG6	= 8500;		/* dissoc rate of Rn_Ggtp complex */
  cpt->kG7	= 200;		/* dissoc rate of Ggtp into Gbg and Gagtp */
  cpt->kOps	= 6.1172e-13;	/* assoc rate of Ops and Gt */
  cpt->kP1	= 0.05497;	/* binding rate of PDE to Gagtp */
  cpt->kP1rev	= 0;		/* dissoc rate of PDE_Gagtp  */
  cpt->kP2	= 940.7;	/* activ rate of first subunit of PDEg of PDE_Gagtp */
  cpt->kP3	= 1.4983e-9;	/* binding rate of Gagtp to PDE_Gagtp */
  cpt->kP4	= 21.088;	/* activ rate of 2nd PDEg to G_agtp_PDE_Gagtp */

  /* Signal recovery-related model parameter values and sources */

  cpt->kRK1_0	= 0.1724;	/* binding rate of RK to unphosphor R* */
  cpt->om	= 2.5;		/* expon rate of decay of RK affinity for R* with phos */
  cpt->kRK2	= 250;		/* dissoc rate of R* from RK prior to phosphorylation */
  cpt->kRK3atp	= 4000;		/* binding rate of ATP to R*_RK */
  cpt->kRK4	= 250;		/* dissoc rate of R* from R*_RK complex */
  cpt->kArr	= 9.9147e-6;	/* binding rate of Arr to singly-phos R* */
  cpt->kA2	= 0.026;	/* dissoc rate of R* from Arr_R* complex prior to R* inact */
  cpt->mArr	= 9.5475e-6;	/* linear rate of increase Arr aff for R* with phosph */
  cpt->kA3	= 1.1651;	/* dissoc rate of R* from Arr_R* complex following R* inact */
  cpt->kA4	= 2.9965e-7;	/* binding rate of Arr to form homo-oliomers */
  cpt->kA5	= 0.424;	/* dissoc rate of Arr from homo-oliomers */
  cpt->kRrecyc	= 7.0e-4;	/* rate const for R regen from Ops */
  cpt->ktherm	= 0.0238;	/* thermal decay of R* */
  cpt->kGrecyc	   = 2.0;	/* binding rate for Gagdp to Gbg */
  cpt->kGshutoff   = 0.05;	/* rate of Gagtp auto-catalytic GTPase activity */
  cpt->kPDEshutoff = 0.1;	/* rate of PDE-induced spon PDE_Gagtp shutoff */
  cpt->kRGS1	= 4.8182e-5;	/* binding rate of RGS9_1 to PDE_Gagtp */
  cpt->kRGS2	= 98;		/* rate of hydrolysis and dissoc of a PDE from PDE_Gagtp */
  cpt->kRec1	= 0.011;	/* rate of Ca-triggered Rec conf change, uM-1 s-1 */
  cpt->kRec2	= 0.05;		/* rate of Rec conf change */
  cpt->kRec3	= 4.1081e-4;	/* binding rate of Ca_Rec to RK */
  cpt->kRec4	= 0.610084;	/* dissoc rate of RK from Rec_Ca */

  /* Ca2+ and cGMP regulation-related model parameter values and sources */

  cpt->Vcyto	= 0.03916;	/* outer segment cytoplasmic volume, pL */
  cpt->Kc1	= 0.171;	/* EC50 for GCAP1-mediated Ca feedback on GC activity, uM */
  cpt->Kc2	= 0.059;	/* EC50 for GCAP2-mediated Ca feedback on GC activity, uM */
  cpt->m1	= 3.0;		/* Hill coeff for GCAP1 */
  cpt->m2	= 1.5;		/* Hill coeff for GCAP2 */
  cpt->Amax	= 60;		/* max rate of cGMP synthesis, uM s-1 */
  cpt->Bdark	= 3.19;		/* dark rate of cGMP hydrolysis */
  cpt->Bsub	= 2.1826e-3;	/* rate const for one catalytic PDE subunit */
  cpt->fCa	= 0.12;		/* frac of circulating current carried by Ca */
  cpt->Jdark	= 14.87;	/* dark circulating current, pA */
  cpt->cGMPdark = 6.5;		/* dark cGMP concentration, uM */
  cpt->ncg	= 3.8;		/* Hill coeff for opening cGMP-gated channels */
  cpt->gCa	= 981.3558;	/* rate of Ca extrusion Na/Ca K+ pump */
  cpt->Ca_dark	= 0.25;		/* dark Ca concentration, uM */
  cpt->Ca_0	= 0.023;	/* minimum Ca concentration, uM */
  cpt->k1	= 9.37059;	/* binding rate of Ca to cytoplasmic buffers, uM s-1 */
  cpt->k2	= 46.412;	/* dissoc rate of Ca from cytoplasmic buffers */
  cpt->eT	= 400;		/* total Ca buffer molecules conc, uM */
  cpt->gain_spr_mul = 1.0;	/* gain for single photon resp, sets peak ampl, def ~1 pA */
  // cpt->chan_eff = 0.75;	/* eff factor for resp, from pipette leakage, Gross (2012a) erratum */
  cpt->chan_eff = 1;		/* eff factor for resp, from pipette leakage, Gross (2012a) erratum */

  vpt = &cpt->vars;
  				   // starting values of state variables for dark
			 	   // Only cGMP, Ca, and Rect are listed in uM
				   //  others are in absolute n
				   // Equil for 10 sec in dark with tcomp18d
				   

/*
 
  //   for 35 uM Tot Rec   from DellOrco et al (2009)
  //   equilib with tcomp18d1
 
  vpt->Rh      = 9.9e7;            // unphosphorylated rhodopsin
  vpt->R_Gt    = 9.095e5;          // 
  vpt->Gt      = 2.4118e7;         // initial Gt
  vpt->RK      = 1.08371e5;        // initial rhodopsin kinase
  vpt->Arr     = 1.260753e6;       // initial Arrestin 
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 8.2998;           // initial Rec, uM
  vpt->Recr_Ca = 0.3612;           // 
  vpt->Cafree  = 0.2499055;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 26.339;	   // 
  vpt->Cabuff  = 19.2119;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.4993224;        // dark cyclic GMP, uM
  cpt->gain    = 14.28;		   // gain for whole rod
  cpt->gain_spr = 0.425;	   // gain for single photon resp
*/

 
  // from https://www.ebi.ac.uk/biomodels-main/BIOMD0000000578 
  //   Gt Tot = 1e7, Rect Tot = 1e7, with original numbers from website
  //   takes 495 iterations to equilib
  //

/*
  vpt->Rh      = 9.81525e7;        // unphosphorylated rhodopsin
  vpt->Gt      = 8152500;          // initial Gt
  vpt->R_Gt    = 1847500;          // 
  vpt->RK      = 589;		   // initial rhodopsin kinase
  vpt->Arr     = 1.260760e6;       // initial Arrestin 
  vpt->Arr_di  = 1.123300e6;       // arrestin dimer
  vpt->Arr_tetra = 8.91810e5;      // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 9289650;          // initial Rec, mol 
  vpt->Recr_Ca = 510930;           // 
  vpt->Cafree  = 0.250000;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 199420;	   // 
  vpt->Cabuff  = 19.2199;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.4944;           // dark cyclic GMP, uM
  cpt->gain    = 58;		   // gain for whole rod
  cpt->gain_spr = 0.38;		   // gain for single photon resp
*/

  // from https://www.ebi.ac.uk/biomodels-main/BIOMD0000000578 
  //      Gt Tot = 1e7, Rect Tot = 1e7, starting with original numbers
  //      but then equilib with tcomp18d1
  //

/*
  vpt->Rh      = 98168923.8;         // unphosphorylated rhodopsin
  vpt->Gt      = 1847825.24;       // initial Gt
  vpt->R_Gt    = 8152674.75;       // 
  vpt->RK      = 579.907781;       // initial rhodopsin kinase
  vpt->Arr     = 1.26075336e6;     // initial Arrestin 
  vpt->Arr_di = 1.123332e6;        // arrestin dimer
  vpt->Arr_tetra = 8.91795e5;      // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 9289702.89;       // initial Rec, mol 
  vpt->Recr_Ca = 510824.907;       // 
  vpt->Cafree  = 0.249943611;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 199472.198;	   // 
  vpt->Cabuff  = 19.2160611;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.49711103;         // dark cyclic GMP, uM
  cpt->gain    = 66;		   // gain for whole rod
  cpt->gain_spr = 0.38;		   // gain for single photon resp
*/

  //   modified from 35 uM Tot Rec   from DellOrco et al (2009)
  //      to fall faster than invergo, Rect Tot = 1e6
  //   equilib with tcomp18d1
  //   falling edge too slow

/*   
  vpt->Rh      = 98163998.7;         // unphosphorylated rhodopsin
  vpt->Gt      = 1847749.66;       // initial Gt
  vpt->R_Gt    = 8152750.33;       // 
  vpt->RK      = 6822.79595;       // initial rhodopsin kinase
  vpt->Arr     = 1.260753e6;       // initial Arrestin 
  vpt->Arr_di = 1.123335e6;        // arrestin dimer
  vpt->Arr_tetra = 8.91798e5;      // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 764736.737;       // initial Rec, mol 
  vpt->Recr_Ca = 42050.8563;       // 
  vpt->Cafree  = 0.249923037;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 193191.796;	   // 
  vpt->Cabuff  = 19.2144091;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.4981737;        // dark cyclic GMP, uM
  cpt->gain    = 62;		   // gain for whole rod
  cpt->gain_spr = 1.25;		   // gain for single photon resp
*/

  //   modified from invergo (2014)
  //      to fall faster than orig invergo (above), Rect Tot = 6e5
  //   equilib with tcomp18d1
  //
  //   a little too late on peak -- decay time same as Fig 1 invergo

/* 
  vpt->Rh      = 98063295.9;         // unphosphorylated rhodopsin
  vpt->Gt      = 1846204.06;       // initial Gt
  vpt->R_Gt    = 8154295.93;       // 
  vpt->RK      = 12941.2538;       // initial rhodopsin kinase
  vpt->Arr     = 1.260754e6;       // initial Arrestin 
  vpt->Arr_di = 1.123335e6;         // arrestin dimer
  vpt->Arr_tetra = 8.91798e5;       // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 391125.095;       // initial Rec, mol 
  vpt->Recr_Ca = 21503.4809;       // 
  vpt->Cafree  = 0.249902266;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 187385.72;	   // 
  vpt->Cabuff  = 19.2127;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.4992595;        // dark cyclic GMP, uM
  cpt->gain    = 62; 		   // gain for whole rod
  cpt->gain_spr = 1.8;		   // gain for single photon resp (~1 pA), 
  				   //   from Gross et al (2012a), 
 				   //   calib with tcomp18bd
*/

  //   modified from invergo (2014)
  //      to fall faster than orig invergo (above), Rect Tot = 4.5e5
  //   equilib with tcomp18d1
  //
  //   peak just about right (200 ms) -- decay time same as Fig 1 invergo

 
  vpt->Rh      = 92710334.8;         // unphosphorylated rhodopsin
  vpt->Gt      = 1856124.7;       // initial Gt
  vpt->R_Gt    = 8144375.29;       // 
  vpt->RK      = 19151.1818;       // initial rhodopsin kinase
  vpt->Arr     = 1.260754e6;       // initial Arrestin 
  vpt->Arr_di  = 1.123334e6;       // arrestin dimer
  vpt->Arr_tetra = 8.91798e5;      // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 255112.214;       // initial Rec, mol 
  vpt->Recr_Ca = 14025.9885;       // 
  vpt->Cafree  = 0.249904705;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 180875.793;	   // 
  vpt->Cabuff  = 19.2129215;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.49921895;       // dark cyclic GMP, uM
  cpt->gain    = 62; 		   // gain for whole rod
  cpt->gain_spr = 2.3;		   // gain for single photon resp (~1 pA), 
  				   //   from Gross et al (2012a), 
 				   //   calib with tcomp18bd


  //   modified from invergo (2014)
  //      to fall faster than invergo, Rect Tot = 3e5
  //   equilib with tcomp18d1
  //   recovery too fast

 
/*
  vpt->Rh      = 99543993;         // unphosphorylated rhodopsin
  vpt->Gt      = 9634993.01;       // initial Gt
  vpt->R_Gt    = 365506.988;       // 
  vpt->RK      = 34714.6836;       // initial rhodopsin kinase
  vpt->Arr     = 1.260753e6;       // initial Arrestin 
  vpt->Arr_di = 1.123335e6;         // arrestin dimer
  vpt->Arr_tetra = 8.91798e5;       // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 128407.953;       // initial Rec, mol 
  vpt->Recr_Ca = 7059.14528;       // 
  vpt->Cafree  = 0.2499055;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 165012.289;	   // 
  vpt->Cabuff  = 19.2117411;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.49947419;       // dark cyclic GMP, uM
  cpt->gain    = 57; 		   // gain for whole rod
  cpt->gain_spr = 3.0;		   // gain for single photon resp
*/

  //   modified from 35 uM Tot Rec   from DellOrco et al (2009)
  //      to fall faster than invergo, Rect = 1e5
  //   equilib with tcomp18d1

/* 
  vpt->Rh      = 99543993;         // unphosphorylated rhodopsin
  vpt->Gt      = 9634993.01;       // initial Gt
  vpt->R_Gt    = 365506.988;       // 
  vpt->RK      = 119315.2068;      // initial rhodopsin kinase
  vpt->Arr     = 1.260753e6;       // initial Arrestin 
  vpt->Arr_di = 1.123335e6;         // arrestin dimer
  vpt->Arr_tetra = 8.91798e5;       // arrestin tetramer
  vpt->PDE     = 2.0e6;            // 
  vpt->RGS     = 1e5;              // 
  vpt->Rect    = 18310.3433;       // initial Rec, mol 
  vpt->Recr_Ca = 1004.5826;        // 
  vpt->Cafree  = 0.2499055;	   // dark free Ca, uM
  vpt->Recr_Ca_RK = 80710.7923;	   // 
  vpt->Cabuff  = 19.21102;	   // dark Ca buffered, uM
  vpt->cGMP    = 6.49998096;       // dark cyclic GMP, uM
  cpt->gain    = 65; 		   // gain for whole rod
  cpt->gain_spr = 5.6;		   // gain for single photon resp
*/

  for (n=0; n<NPHOS; n++) {
      vpt->Rn[n]	= 0;        // phosphorylated rhodopsin
      vpt->Rn_RKpre[n]  = 0;        // rhodopsin kinase
      vpt->Rn_RKpost[n] = 0;        // rhodopsin kinase
      vpt->Rn_Arr[n]    = 0;        // rhodopsin-arrestin complex
      vpt->Rn_Gt[n]     = 0;        // receptor G complex
      vpt->Rn_G[n]      = 0;        // receptor G complex
      vpt->Rn_Ggtp[n]   = 0;        // receptor G complex
  }
  vpt->Ops = 0e-3;                  // ligand-free receptor 
  vpt->Ops_Gt = 0;                  // receptor Gt complex
  vpt->Ops_G = 0;                   // receptor G complex
  vpt->Ops_Ggtp = 0;                // receptor G complex
  vpt->Ggtp = 0;                    // Ggtp 
  vpt->Gagtp = 0;                   // 
  vpt->Gbg = 0;                     // 
  vpt->PDE_Gagtp  = 0;              // PDE-Gagtp complex
  vpt->PDEa_Gagtp = 0;              // activated PDE-Gagtp complex
  vpt->Gagtp_PDEa_Gagtp = 0;        // activated PDE-Gagtp complex
  vpt->Gagtp_aPDEa_Gagtp = 0;       // activated PDE-Gagtp complex
  vpt->RGS_PDEa_Gagtp = 0;          // activated RGS-PDE-Gagtp complex
  vpt->RGS_Gagtp_aPDEa_Gagtp = 0;   // activated RGS-PDE-Gagtp complex
  vpt->Gagdp      = 0;              // Gagdp 
  vpt->cond       = 0.9996;         // dark conductance normalized 0 to 1

//  fprintf (stderr,"# Rh   tot (1e8) %g\n",vpt->Rh + vpt->R_Gt);
//  fprintf (stderr,"# Gt   tot (1e7) %g\n",vpt->Gt + vpt->R_Gt);
//  fprintf (stderr,"# Rect tot (1e7) %g\n",vpt->Rect + vpt->Recr_Ca + vpt->Recr_Ca_RK);
//  fprintf (stderr,"# RK   tot (2e5) %g\n",vpt->RK + vpt->Recr_Ca_RK);
//  fprintf (stderr,"# Arr  tot %g\n",vpt->Arr + 2*vpt->Arr_di + 4*vpt->Arr_tetra);

}

/*------------------------------------*/

void impl_a (int eqn, double *a, double *ao, double rr, double dt)

 /* solve the equation a <---    (forward rate), increment */
 /*                    a --->    (reverse rate), decrement */
 
{

#define MAXITER 100
#define MAXERR 1e-8
	int i;
	double a_prev;
	double err, rep;

//  fprintf (stderr,"impl_ab_c: %d a %g ao %g\n",eqn,*a,*ao);
  a_prev = *ao;
  err = 1;
  rep = MAXITER;
  rep = 1;
  for (i=0; i<rep && err > MAXERR; i++) {
      *a = *ao / (1.0 + rr * dt);
      err = abs((*a - a_prev) / (a_prev < 1e-10 ? 1.0:a_prev));

      a_prev = *a;
    // fprintf (stderr,"impl_a_b: eqn %d i %d err %g\n",eqn,i,err);
  }
  // *ao = *a;	/* save the updated original for later equations */
  // *bo = *b;
//  fprintf (stderr,"impl_a_b: %d i %d final err %g\n",eqn,i,err);
//  fprintf (stderr,"impl_a_b: %d a %g b %g ao %g bo %g\n",eqn,*a,*b,*ao,*bo);
}

/*------------------------------------*/

void impl_a_b (int eqn, int rep, double *a, double *b, double *ao, double *bo, double fr, double rr, double dt)

 /* solve the equation a --> b (forward rate) */
 /*                    b --> a   (reverse rate) */
 
{

#define MAXITER 100
#define MAXERR 1e-8
	int i;
	double b_prev;
	double diff, err;

//  fprintf (stderr,"impl_ab_c: %d a %g b %g ao %g bo %g\n",eqn,*a,*b,*ao,*bo);
  b_prev = *bo;
  err = 1;
  if (rep==0) rep = MAXITER;
  for (i=0; i<rep && err > MAXERR; i++) {
      *b = (*bo + dt * (*a * fr)) / (1.0 + rr * dt);
      diff = *b - *bo;
      *a = *ao - diff;
      err = abs((*b - b_prev) / (b_prev < 1e-10 ? 1.0:b_prev));

      b_prev = *b;
    // fprintf (stderr,"impl_a_b: eqn %d i %d err %g\n",eqn,i,err);
  }
  // *ao = *a;	/* save the updated original for later equations */
  // *bo = *b;
//  fprintf (stderr,"impl_a_b: %d i %d final err %g\n",eqn,i,err);
//  fprintf (stderr,"impl_a_b: %d a %g b %g ao %g bo %g\n",eqn,*a,*b,*ao,*bo);
}

/*------------------------------------*/

void impl_a_b (int eqn, double *a, double *b, double *ao, double *bo, double fr, double rr, double dt)

{
    impl_a_b (eqn, 0, a, b, ao, bo, fr, rr, dt);
}
/*------------------------------------*/

void impl_ab_c (int eqn, int rep, double *a, double *b, double *c, double *ao, double *bo, double *co, double fr, double rr, double dt)

 /* solve the equation a + b --> c (forward rate) */
 /*                    c -> a + b  (reverse rate) */
 
{

/*
    rpi->R_Gt = (rpi_old->R_Gt + dt * (rpi->Gt * rpi->Rh * cpi->kG1_0 * 1.6)) / (1 + cpi->kG2 * 315 * dt);
    rpi->Rh = (rpi_old->Rh + dt * rpi->R_Gt * cpi->kG2 * 315) / (1 + rpi->Gt * cpi->kG1_0*1.6 * dt);
    rpi->Gt = (rpi_old->Gt + dt * rpi->R_Gt * cpi->kG2 * 315) / (1 + rpi->Rh * cpi->kG1_0*1.6 * dt);
*/
#define MAXITER 100
#define MAXERR 1e-8
	int i;
	double c_prev;
	double diff, err_a, err_b, err;

//  fprintf (stderr,"impl_ab_c: %d a %g b %g c %g ao %g bo %g co %g\n",eqn,*a,*b,*c,*ao,*bo,*co);
  c_prev = *co;
  err = 1;
  if (rep==0) rep = MAXITER;
  for (i=0; i<rep && err > MAXERR; i++) {
      //*a = (*ao + dt * (*c * rr)) / (1.0 + *b * fr * dt); 
      //*b = (*bo + dt * (*c * rr)) / (1.0 + *a * fr * dt); 
      *c = (*co + dt * (*a * *b * fr)) / (1.0 + rr * dt);
      diff = *c - *co;
      *a = *ao - diff;
      *b = *bo - diff;
      // *c = *co - diff;
      err = abs((*c - c_prev) / (c_prev < 1e-10 ? 1.0:c_prev));

      c_prev = *c;
    // if (eqn<0) fprintf (stderr,"impl_ab_c: %d i %d err %g a %g diff %g c %g\n",eqn,i,err,*a,diff,*c);
    // fprintf (stderr,"impl_ab_c: %d i %d err %g\n",eqn,i,err);
  }
  // *ao = *a;	/* save the updated original for later equations */
  // *bo = *b;
  // *co = *c;
  //if (eqn<0) fprintf (stderr,"impl_ab_c: %d i %d final err %g\n",eqn,i,err);
 //  fprintf (stderr,"impl_ab_c: %d i %d final err %g\n",eqn,i,err);
//  fprintf (stderr,"impl_ab_c: %d a %g b %g c %g ao %g bo %g co %g\n",eqn,*a,*b,*c,*ao,*bo,*co);
}

/*------------------------------------*/

void impl_ab_c (int eqn, double *a, double *b, double *c, double *ao, double *bo, double *co, double fr, double rr, double dt)
{
      impl_ab_c (eqn, 0, a, b, c, ao, bo, co, fr, rr, dt);
}

/*------------------------------------*/

double runrodi_orig(photreci *rpnti, double iflux)

/* mouse rod phototransduction from Invergo et al (2014) */

/*  Invergo BM, Dell'Orco D, Montanucci L, Koch KW, Bertranpetit J.
    A comprehensive model of the phototransduction cascade in mouse rod cells.
    Mol Biosyst. 2014 Jun;10(6):1481-9. doi: 10.1039/c3mb70584f
    See supplementary electronic info file.
    Parameters taken (Dec 2017) from: https://www.ebi.ac.uk/biomodels-main/BIOMD0000000578
*/

/*    Called by "runrec()" in ncsub.cc. */

{
	int i,k,n;
	double v, dt, tc, tcx;
	double kA1, cond, chan_cond;
	double Ca2_frac, Icgchan;
#define NRODISIZ 200
	double *d;
        double diff[NRODISIZ];
	recpari *cpi;
	phvarsi *rpi;
	phvarsi *rpo;


     cpi = rpnti->chtyp;
     rpi = &rpnti->vars;	
     // if      (iflux > 1e7) tc = 1000;
     // else if (iflux > 1e6) tc = 800;
     // else if (iflux > 1e5) tc = 250;
     // else if (iflux > 1e4) tc = 150;
     // else if (iflux > 1e3) tc = 50;
     // else if (iflux > 1e2) tc = 20;
     // else tc = 10;
     
     //  without kG1, kG2 equation below 
     
     if      (rpi->Gagtp > 1e3) tc = 1000;
     else if (rpi->Gagtp > 500) tc = 500;
     else if (rpi->Gagtp > 300) tc = 400;
     else if (rpi->Gagtp > 100) tc = 300;
     else if (rpi->Gagtp > 50)  tc = 300;
     else if (rpi->Gagtp > 25)  tc = 200;
     else if (rpi->Gagtp > 10)  tc = 180;
     else if (rpi->Gagtp > 3)   tc = 150;
     else if (rpi->Gagtp > 1)   tc = 120;
     else if (rpi->Gagtp > 0.3) tc = 100;
     else if (rpi->Gagtp > 0.1) tc = 50;
     else if (rpi->Gagtp > 0.01) tc = 20;
     else                       tc = 12;

     tcx = cpi->kG2 * 315 * stiminc;			// prevents instability in kG2 equation below  	
     if (tcx > tc) tc = int(tcx*1.2);

     // if (iflux > 1e3) 
     //      fprintf (stderr,"iflux %g Rh %g\n",iflux,rpi->Rh);

     iflux *= cpi->gain_spr * cpi->gain_spr_mul;	// calibration with tcomp18bd to give SPR of ~1 pA, 
							//  for ~1 pA peak, from Gross et al (2012a)
							// But see Gross et al (2015) -> 0.55 pA
     rpo = rpi;

     dt = stiminc/tc;		// time increment 
    
     for (i=0; i<tc; i++) {

        d = diff;

	d[1] = iflux * rpi->Rh / cpi->Rtot / tc;
	rpi->Rh    = rpo->Rh    - d[1];
	rpi->Rn[0] = rpo->Rn[0] + d[1];

	d[2] = iflux * rpi->R_Gt / cpi->Rtot / tc;
	rpi->R_Gt     = rpo->R_Gt     - d[2];
	rpi->Rn_Gt[0] = rpo->Rn_Gt[0] + d[2];

 // fprintf (stderr,"Rn0 %g iflux %g\n",rpi->Rn[0], iflux);
 //
	for (n=0; n<NPHOS; n++) {
		  d[3+n] = rpi->Rn[n] * rpi->RK * cpi->kRK1_0 * exp (-cpi->om*n) * dt;
		  rpi->Rn[n]        = rpo->Rn[n] - d[3+n]; 
		  rpi->RK           = rpo->RK    - d[3+n];
		  rpi->Rn_RKpre[n]  = rpo->Rn_RKpre[n] + d[3+n];
	}
	for (n=0; n<NPHOS; n++) {
		  d[10+n] = rpi->Rn_RKpre[n] * cpi->kRK2 * dt;
		  rpi->Rn_RKpre[n]  = rpo->Rn_RKpre[n]  - d[10+n];
		  rpi->Rn[n]        = rpo->Rn[n]        + d[10+n]; 
		  rpi->RK           = rpo->RK           + d[10+n];
	}
//  fprintf (stderr,"Rn0 %g RK %g R0_RK %g\n",rpi->Rn[0], rpi->RK, rpi->Rn_RKpre[0]);

	for (n=0; n<NPHOS-1; n++) {
		  d[17+n] = rpi->Rn_RKpre[n] * cpi->kRK3atp * dt;
		  rpi->Rn_RKpre[n]    = rpo->Rn_RKpre[n]    - d[17+n];
		  rpi->Rn_RKpost[n+1] = rpo->Rn_RKpost[n+1] + d[17+n];
        }

	for (n=1; n<NPHOS; n++) {
		  d[24+n] = rpi->Rn_RKpost[n] * cpi->kRK4 * dt;
		  rpi->Rn_RKpost[n]   = rpo->Rn_RKpost[n]   - d[24+n];
		  rpi->RK             = rpo->RK             + d[24+n];
		  rpi->Rn[n]          = rpo->Rn[n]          + d[24+n];
        }

	for (n=1; n<NPHOS; n++) {
		  if (n<=4) kA1 = cpi->kArr + (n-1)*cpi->mArr;
		  else      kA1 = cpi->kArr + 3    *cpi->mArr;
		  d[31+n] = rpi->Rn[n] * rpi->Arr * kA1 * dt;
		  rpi->Rn[n]     = rpo->Rn[n]     - d[31+n];
		  rpi->Arr       = rpo->Arr       - d[31+n];
		  rpi->Rn_Arr[n] = rpo->Rn_Arr[n] + d[31+n];
        }
	for (n=1; n<NPHOS; n++) {
		  d[38+n] = rpi->Rn_Arr[n] * cpi->kA2 * dt;
		  rpi->Rn_Arr[n] = rpo->Rn_Arr[n] - d[38+n];
		  rpi->Rn[n]     = rpo->Rn[n]     + d[38+n];
		  rpi->Arr       = rpo->Arr       + d[38+n];
        }

	for (n=1; n<NPHOS; n++) {
		  d[45+n] = rpi->Rn_Arr[n] * cpi->kA3 * dt;
		  rpi->Rn_Arr[n] = rpo->Rn_Arr[n] - d[45+n];		  
		  rpi->Ops       = rpo->Ops       + d[45+n];
		  rpi->Arr       = rpo->Arr       + d[45+n];
        }

	d[52] = rpi->Arr * rpi->Arr * cpi->kA4 * dt;
	rpi->Arr         = rpo->Arr         - 2 * d[52];
	rpi->Arr_di      = rpo->Arr_di      +     d[52];
	d[53] = rpi->Arr_di * cpi->kA5 * dt;
	rpi->Arr_di      = rpo->Arr_di      -     d[53];
	rpi->Arr         = rpo->Arr         + 2 * d[53];

	d[54] = rpi->Arr_di * rpi->Arr_di * cpi->kA4 * dt;
	rpi->Arr_di      = rpo->Arr_di      - 2 * d[54];
	rpi->Arr_tetra   = rpo->Arr_tetra   +     d[54];
	d[55] = rpi->Arr_tetra * cpi->kA5 * dt;
	rpi->Arr_tetra   = rpo->Arr_tetra   -     d[55];
	rpi->Arr_di      = rpo->Arr_di      + 2 * d[55];

	for (n=0; n<NPHOS; n++) {
	    d[56+n] = rpi->Rn[n] * cpi->ktherm * dt;
	    rpi->Rn[n]    = rpo->Rn[n]    - d[56+n];
	    rpi->Ops      = rpo->Ops      + d[56+n];
	}

	d[63] = rpi->Ops * rpi->Gt * cpi->kOps * dt;
	rpi->Ops          = rpo->Ops          - d[63];
	rpi->Gt           = rpo->Gt           - d[63];
	rpi->Ops_Gt       = rpo->Ops_Gt       + d[63];
	d[64] = rpi->Ops_Gt * cpi->kG2 * dt;
	rpi->Ops_Gt       = rpo->Ops_Gt       - d[64];
	rpi->Ops          = rpo->Ops          + d[64];
	rpi->Gt           = rpo->Gt           + d[64];

	d[65] = rpi->Ops_Gt * cpi->kG3 * dt;
	rpi->Ops_Gt       = rpo->Ops_Gt       - d[65];
	rpi->Ops_G        = rpo->Ops_G        + d[65];
	d[66] = rpi->Ops_G * cpi->kG4gdp * dt;
	rpi->Ops_G        = rpo->Ops_G        - d[66];
	rpi->Ops_Gt       = rpo->Ops_Gt       + d[66];

	d[67] = rpi->Ops_G * cpi->kG5gtp * dt;
	rpi->Ops_G        = rpo->Ops_G        - d[67];
	rpi->Ops_Ggtp     = rpo->Ops_Ggtp     + d[67];
	
	d[68] = rpi->Ops_Ggtp * cpi->kG6 * dt;
	rpi->Ops_Ggtp     = rpo->Ops_Ggtp     - d[68];
	rpi->Ops          = rpo->Ops          + d[68];
	rpi->Ggtp         = rpo->Ggtp         + d[68];

	d[69] = rpi->Ops * cpi->kRrecyc * dt;
	rpi->Ops          = rpo->Ops          - d[69];
	rpi->Rh           = rpo->Rh           + d[69];

		
  // fprintf (stderr,"R_Gt %g tc %g\n",rpi->R_Gt, tc);
 	   d[70] = rpi->Gt * rpi->Rh * cpi->kG1_0 * 1.6 * dt;
	   rpi->Rh           = rpo->Rh           - d[70];
	   rpi->Gt           = rpo->Gt           - d[70];
	   rpi->R_Gt         = rpo->R_Gt         + d[70];
  // fprintf (stderr,"R_Gt %g \n",rpi->R_Gt);

	   d[71] = rpi->R_Gt * cpi->kG2 * 315 * dt;
	   rpi->R_Gt         = rpo->R_Gt         - d[71];
	   rpi->Rh           = rpo->Rh           + d[71];
	   rpi->Gt           = rpo->Gt           + d[71];
  // fprintf (stderr,"R_Gt %g \n\n",rpi->R_Gt);
   //	}  /* */

 	for (n=0; n<NPHOS; n++) {
	    d[72+n] = rpi->Rn[n] * rpi->Gt * cpi->kG1_0 * exp (-cpi->omGt*n) * dt;
	    rpi->Rn[n]     = rpo->Rn[n]     - d[72+n]; 
	    rpi->Gt        = rpo->Gt        - d[72+n];
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  + d[72+n];
	}
	for (n=0; n<NPHOS; n++) {
	    d[79+n] = rpi->Rn_Gt[n] * cpi->kG2 * dt;
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  - d[79+n];
	    rpi->Rn[n]     = rpo->Rn[n]     + d[79+n]; 
	    rpi->Gt        = rpo->Gt        + d[79+n];
	}

	for (n=0; n<NPHOS; n++) {
	    d[86+n] = rpi->Rn_Gt[n] * cpi->kG3 * dt;
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  - d[86+n];
	    rpi->Rn_G[n]   = rpo->Rn_G[n]   + d[86+n];
	}
	for (n=0; n<NPHOS; n++) {
	    d[93+n] = rpi->Rn_G[n] * cpi->kG4gdp * dt;
	    rpi->Rn_G[n]   = rpo->Rn_G[n]   - d[93+n];
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  + d[93+n];
	}

	for (n=0; n<NPHOS; n++) {
	    d[100+n] = rpi->Rn_G[n] * cpi->kG5gtp * dt;
	    rpi->Rn_G[n]    = rpo->Rn_G[n]    - d[100+n];
	    rpi->Rn_Ggtp[n] = rpo->Rn_Ggtp[n] + d[100+n];
	}
	for (n=0; n<NPHOS; n++) {
	    d[107+n] = rpi->Rn_Ggtp[n] * cpi->kG6 * dt;
	    rpi->Rn_Ggtp[n] = rpo->Rn_Ggtp[n] - d[107+n];
	    rpi->Rn[n]      = rpo->Rn[n]      + d[107+n];
	    rpi->Ggtp       = rpo->Ggtp       + d[107+n];
	}

	d[114] = rpi->Ggtp * cpi->kG7 * dt;
	rpi->Ggtp           = rpo->Ggtp           - d[114];
	rpi->Gagtp          = rpo->Gagtp          + d[114];
	rpi->Gbg            = rpo->Gbg            + d[114];
	
	d[115] = rpi->PDE * rpi->Gagtp * cpi->kP1 * dt;
	rpi->PDE            = rpo->PDE            - d[115];
	rpi->Gagtp          = rpo->Gagtp          - d[115];
	rpi->PDE_Gagtp      = rpo->PDE_Gagtp      + d[115];
	d[116] = rpi->PDE_Gagtp * cpi->kP1rev * dt;
	rpi->PDE_Gagtp      = rpo->PDE_Gagtp      - d[116];
	rpi->PDE            = rpo->PDE            + d[116];
	rpi->Gagtp          = rpo->Gagtp          + d[116];

	// dark noise, check with tcomp18bd
	//  set to have s.d. of 1/5 pA, giving SNR of 5 for 1 pA SPR
	
	if (rpnti->dnois>0) {
                 if (rpnti->dstate) drand_setstate (rpnti->dstate);
		 rpi->PDE_Gagtp  += (drand() - 0.5) * 2.0e-1 * rpnti->dnois * cpi->gain_spr_mul;
		 // rpi->PDE_Gagtp  += gasdev() * 4.55e-2 * rpnti->dnois * cpi->gain_spr_mul;
 		 restorstate();  	/* restore default random state */
	}

	d[117] = rpi->PDE_Gagtp * cpi->kP2 * dt;
	rpi->PDE_Gagtp      = rpo->PDE_Gagtp      - d[117];
	rpi->PDEa_Gagtp     = rpo->PDEa_Gagtp     + d[117];

	d[118] = rpi->PDEa_Gagtp * rpi->Gagtp * cpi->kP3 * dt;
	rpi->PDEa_Gagtp       = rpo->PDEa_Gagtp       - d[118];
	rpi->Gagtp            = rpo->Gagtp            - d[118];
	rpi->Gagtp_PDEa_Gagtp = rpo->Gagtp_PDEa_Gagtp + d[118];

	d[119] = rpi->Gagtp_PDEa_Gagtp * cpi->kP4 * dt;
	rpi->Gagtp_PDEa_Gagtp  = rpo->Gagtp_PDEa_Gagtp  - d[119];
	rpi->Gagtp_aPDEa_Gagtp = rpo->Gagtp_aPDEa_Gagtp + d[119];

	d[120] = rpi->RGS * rpi->Gagtp_aPDEa_Gagtp * cpi->kRGS1 * dt;	// 1st kRGS1
	rpi->RGS                   = rpo->RGS                   - d[120];
	rpi->Gagtp_aPDEa_Gagtp     = rpo->Gagtp_aPDEa_Gagtp     - d[120];
	rpi->RGS_Gagtp_aPDEa_Gagtp = rpo->RGS_Gagtp_aPDEa_Gagtp + d[120];

	d[121] = rpi->RGS_Gagtp_aPDEa_Gagtp * cpi->kRGS2 * dt;
	rpi->RGS_Gagtp_aPDEa_Gagtp = rpo->RGS_Gagtp_aPDEa_Gagtp - d[121];
	rpi->PDEa_Gagtp            = rpo->PDEa_Gagtp            + d[121];
	rpi->RGS                   = rpo->RGS                   + d[121];
	rpi->Gagdp                 = rpo->Gagdp                 + d[121];

	d[122] = rpi->RGS * rpi->PDEa_Gagtp * cpi->kRGS1 * dt;		// 2nd kRGS1
	rpi->RGS               = rpo->RGS               - d[122];
	rpi->PDEa_Gagtp        = rpo->PDEa_Gagtp        - d[122];
	rpi->RGS_PDEa_Gagtp    = rpo->RGS_PDEa_Gagtp    + d[122];


	d[123] = rpi->RGS_PDEa_Gagtp * cpi->kRGS2 * dt;
	rpi->RGS_PDEa_Gagtp     = rpo->RGS_PDEa_Gagtp     - d[123];
	rpi->PDE                = rpo->PDE                + d[123];
	rpi->RGS                = rpo->RGS                + d[123];
	rpi->Gagdp              = rpo->Gagdp              + d[123];

	d[124] = rpi->PDEa_Gagtp * cpi->kPDEshutoff * dt;
	rpi->PDEa_Gagtp         = rpo->PDEa_Gagtp         - d[124];
	rpi->PDE                = rpo->PDE                + d[124];
	rpi->Gagdp              = rpo->Gagdp              + d[124];

	d[125] = rpi->Gagtp_aPDEa_Gagtp * cpi->kPDEshutoff * dt;
	rpi->Gagtp_aPDEa_Gagtp  = rpo->Gagtp_aPDEa_Gagtp  - d[125];
	rpi->PDEa_Gagtp         = rpo->PDEa_Gagtp         + d[125];
	rpi->Gagdp              = rpo->Gagdp              + d[125];

	d[126] = rpi->Gagtp * cpi->kGshutoff * dt;
	rpi->Gagtp              = rpo->Gagtp              - d[126];
	rpi->Gagdp              = rpo->Gagdp              + d[126];

	d[127] = rpi->Gagdp * rpi->Gbg * cpi->kGrecyc * dt;
	rpi->Gagdp              = rpo->Gagdp              - d[127];
	rpi->Gbg                = rpo->Gbg                - d[127];
	rpi->Gt                 = rpo->Gt                 + d[127];

	d[128] = rpi->Rect * rpi->Cafree * cpi->kRec1 * dt;
	rpi->Rect               = rpo->Rect               - d[128];
	rpi->Cafree             = rpo->Cafree             - d[128];
	rpi->Recr_Ca            = rpo->Recr_Ca            + d[128];
	d[129] = rpi->Recr_Ca * cpi->kRec2 * dt;
	rpi->Recr_Ca            = rpo->Recr_Ca            - d[129];
	rpi->Rect               = rpo->Rect               + d[129];
	rpi->Cafree             = rpo->Cafree             + d[129];

	d[130] = rpi->Recr_Ca * rpi->RK * cpi->kRec3 * dt;
	rpi->Recr_Ca            = rpo->Recr_Ca            - d[130];
	rpi->RK                 = rpo->RK                 - d[130];
	rpi->Recr_Ca_RK         = rpo->Recr_Ca_RK         + d[130];
	d[131] = rpi->Recr_Ca_RK * cpi->kRec4 * dt;
	rpi->Recr_Ca_RK         = rpo->Recr_Ca_RK         - d[131];
	rpi->Recr_Ca            = rpo->Recr_Ca            + d[131];
	rpi->RK                 = rpo->RK                 + d[131];

	d[132] = rpi->Cafree * (cpi->eT - rpi->Cabuff) * cpi->k1 * dt;
	rpi->Cafree            = rpo->Cafree            - d[132];
	rpi->Cabuff            = rpo->Cabuff            + d[132];
	d[133] = rpi->Cabuff * cpi->k2 * dt;
	rpi->Cabuff            = rpo->Cabuff            - d[133];
	rpi->Cafree            = rpo->Cafree            + d[133];

	d[134] = (rpi->Cafree - cpi->Ca_0) * cpi->gCa * dt;
	rpi->Cafree            = rpo->Cafree            - d[134];

	chan_cond = pow (rpi->cGMP / cpi->cGMPdark, cpi->ncg);

	d[135] = 1e6 * cpi->fCa * cpi->Jdark / 
	    ((2 + cpi->fCa) * Fd * cpi->Vcyto) * chan_cond * dt;
	rpi->Cafree            = rpo->Cafree            + d[135];

	d[136] = (cpi->Amax / (1 + pow(rpi->Cafree/cpi->Kc1, cpi->m1)) +
	     cpi->Amax / (1 + pow(rpi->Cafree/cpi->Kc2, cpi->m2))) * dt;
        rpi->cGMP               = rpo->cGMP               + d[136];

	d[137] = (cpi->Bdark + cpi->Bsub * (rpi->PDEa_Gagtp + 2*rpi->Gagtp_aPDEa_Gagtp) + 
			                                   rpi->Gagtp_PDEa_Gagtp) * rpi->cGMP * dt;        
	rpi->cGMP              = rpo->cGMP              - d[137];

 	 // Ca2_frac = (rpi->Cafree - cpi->Ca_0) / (cpi->Ca_dark - cpi->Ca_0);
	 // Icgchan = 2/(2+cpi->fCa)* pow(rpi->cGMP/cpi->cGMPdark,cpi->ncg) * cpi->Jdark + 
	 //  					cpi->fCa/(cpi->fCa+2) * Ca2_frac * cpi->Jdark;
	 // cond = Icgchan * 1e-2 / 100e-3;
	 
	 cond = chan_cond / cpi->chan_eff;
      }
  // fprintf (stderr,"iflux %-4.2g Rh %-7.3g Rn0 %-7.3g R_Gt %-7.3g RK %-7.3g Arr %-7.3g Ggtp %-7.3g PDE %-7.3g RGS %-7.3g Gagtp %-7.3g cGMP %-7.3g cond %-7.3g tc %g\n",
//  	iflux, rpi->Rh, rpi->Rn[0], rpi->R_Gt, rpi->RK, rpi->Arr, rpi->Ggtp, rpi->PDE, rpi->RGS, rpi->Gagtp, rpi->cGMP, cond, tc);
 // fprintf (stderr,"iflux %g Rn0 %g RK %g Arr %g Ggtp %g Cafree %g Cabuff %g cGMP %g Rh0 %g\n",
 // 	iflux, rpi->Rn[0], rpi->RK, rpi->Arr, rpi->Ggtp, rpi->Cafree,rpi->Cabuff,rpi->cGMP,rpi->Rh);

	return cond;
}

/*------------------------------------*/

double runrodi(photreci *rpnti, double iflux)

/* mouse rod phototransduction from Invergo et al (2014) */

/*  Invergo BM, Dell'Orco D, Montanucci L, Koch KW, Bertranpetit J.
    A comprehensive model of the phototransduction cascade in mouse rod cells.
    Mol Biosyst. 2014 Jun;10(6):1481-9. doi: 10.1039/c3mb70584f
    See supplementary electronic info file.
    Parameters taken (Dec 2017) from: https://www.ebi.ac.uk/biomodels-main/BIOMD0000000578

*/

{
	int i,j,k,n;
	double v, dt, tc, tcx;
	double kA1, cond, chan_cond;
	double Ca2_frac, Icgchan;
#define NRODISIZ 200
	double d;
	recpari *cpi;
	phvarsi *rpi;
	static phvarsi rpi_copy;
	phvarsi *rpc = &rpi_copy;	// copy of original values for updating during integ
	phvarsi *rpd;			// dest for updating vars before separate impl calls
	phvarsi *rpo;			// source of updates for non impl eqns


     cpi = rpnti->chtyp;
     rpi = &rpnti->vars;	
     
     //  without kG1, kG2 equation below 
     
     if      (rpi->Gagtp > 1e4) tc = 1000;
     else if (rpi->Gagtp > 1e3) tc = 500;
     else if (rpi->Gagtp > 500) tc = 500;
     else if (rpi->Gagtp > 300) tc = 400;
     else if (rpi->Gagtp > 100) tc = 300;
     else if (rpi->Gagtp > 50)  tc = 300;
     else if (rpi->Gagtp > 25)  tc = 200;
     else if (rpi->Gagtp > 10)  tc = 180;
     else if (rpi->Gagtp > 3)   tc = 150;
     else if (rpi->Gagtp > 1)   tc = 120;
     else if (rpi->Gagtp > 0.3) tc = 100;
     else if (rpi->Gagtp > 0.1) tc = 50;
     else if (rpi->Gagtp > 0.01) tc = 15;
     else                       tc = 12;

//     tcx = cpi->kG2 * 315 * stiminc;			// prevents instability in kG2 equation below  	
//     if (tcx > tc) tc = int(tcx*1.2);
 
   // fprintf (stderr,"tc %g\n",tc);

     // if (iflux > 1e3) 
     //      fprintf (stderr,"iflux %g Rh %g\n",iflux,rpi->Rh);

     //iflux *= 60 * cpi->gain_spr * cpi->gain_spr_mul;	// calibration with tcomp18bd to give SPR of ~1 pA, 
     iflux *= cpi->gain_spr * cpi->gain_spr_mul;	// calibration with tcomp18bd to give SPR of ~1 pA, 
							//  for ~1 pA peak, from Gross et al (2012a)
							// But see Gross et al (2015) -> 0.55 pA

     dt = stiminc/tc;		// time increment 

      rpo = rpi;		// for orig separate impl calls
     // rpo = rpc;		// for full impl in for loop


     for (i=0; i<tc; i++) {
        *rpc = *rpi;			// copy state variables to incr time step

         rpd = rpc;		// for orig separate impl calls 
        // rpd = rpi;		// for full impl

	v = iflux * rpi->Rh / cpi->Rtot / tc;
	rpi->Rh    = rpo->Rh    - v;
	rpi->Rn[0] = rpo->Rn[0] + v;

	v = iflux * rpi->R_Gt / cpi->Rtot / tc;
	rpi->R_Gt     = rpo->R_Gt     - v;
	rpi->Rn_Gt[0] = rpo->Rn_Gt[0] + v;

 // fprintf (stderr,"Rn0 %g iflux %g\n",rpi->Rn[0], iflux);
 //
 //  fprintf (stderr,"Rn0 %g RK %g R0_RK %g\n",rpi->Rn[0], rpi->RK, rpi->Rn_RKpre[0]);

/*
	for (n=0; n<NPHOS; n++) {
		  d = rpi->Rn[n] * rpi->RK * cpi->kRK1_0 * exp (-cpi->om*n) * dt;
		  rpi->Rn[n]        = rpo->Rn[n]        - d; 
		  rpi->RK           = rpo->RK           - d;
		  rpi->Rn_RKpre[n]  = rpo->Rn_RKpre[n]  + d;
	}
	for (n=0; n<NPHOS; n++) {
		  d = rpi->Rn_RKpre[n] * cpi->kRK2 * dt;
		  rpi->Rn_RKpre[n]  = rpo->Rn_RKpre[n]  - d;
		  rpi->Rn[n]        = rpo->Rn[n]        + d; 
		  rpi->RK           = rpo->RK           + d;
	}

/* */
	
	for (n=0; n<NPHOS; n++) {
	     rpd->Rn[n] = rpi->Rn[n];
	     rpd->RK    = rpi->RK;
	     rpd->Rn_RKpre[n]  = rpi->Rn_RKpre[n];
    	     impl_ab_c (1, &rpi->Rn[n], &rpi->RK, &rpi->Rn_RKpre[n], 
		 	   &rpc->Rn[n], &rpc->RK, &rpc->Rn_RKpre[n], cpi->kRK1_0*exp(-cpi->om*n), cpi->kRK2, dt);
	}
/* */

//  fprintf (stderr,"Rn0 %g RK %g R0_RK %g\n",rpi->Rn[0], rpi->RK, rpi->Rn_RKpre[0]);

	for (n=0; n<NPHOS-1; n++) {
		  d = rpi->Rn_RKpre[n] * cpi->kRK3atp * dt;
		  rpi->Rn_RKpre[n]    = rpo->Rn_RKpre[n]    - d;
		  rpi->Rn_RKpost[n+1] = rpo->Rn_RKpost[n+1] + d;
        }

	for (n=1; n<NPHOS; n++) {
		  d = rpi->Rn_RKpost[n] * cpi->kRK4 * dt;
		  rpi->Rn_RKpost[n]   = rpo->Rn_RKpost[n]  - d;
		  rpi->RK             = rpo->RK            + d;
		  rpi->Rn[n]          = rpo->Rn[n]         + d;
        }

/*
	for (n=1; n<NPHOS; n++) {
		  if (n<=4) kA1 = cpi->kArr + (n-1)*cpi->mArr;
		  else      kA1 = cpi->kArr + 3    *cpi->mArr;
		  d = rpi->Rn[n] * rpi->Arr * kA1 * dt;
		  rpi->Rn[n]     = rpo->Rn[n]     - d;
		  rpi->Arr       = rpo->Arr       - d;
		  rpi->Rn_Arr[n] = rpo->Rn_Arr[n] + d;
        }
	for (n=1; n<NPHOS; n++) {
		  d = rpi->Rn_Arr[n] * cpi->kA2 * dt;
		  rpi->Rn_Arr[n] = rpo->Rn_Arr[n] - d;
		  rpi->Rn[n]     = rpo->Rn[n]     + d;
		  rpi->Arr       = rpo->Arr       + d;
        }
/* */

	for (n=1; n<NPHOS; n++) {
	     rpd->Rn_Arr[n] = rpi->Rn_Arr[n];
	     rpd->Rn[n]     = rpi->Rn[n];
	     rpd->Arr       = rpi->Arr;
	     if (n<=4) kA1 = cpi->kArr + (n-1)*cpi->mArr;
	     else      kA1 = cpi->kArr + 3    *cpi->mArr;
    	     impl_ab_c (2, &rpi->Rn[n], &rpi->Arr, &rpi->Rn_Arr[n], 
		 	   &rpc->Rn[n], &rpc->Arr, &rpc->Rn_Arr[n], kA1, cpi->kA2, dt);
        }
/* */

	for (n=1; n<NPHOS; n++) {
		  d = rpi->Rn_Arr[n] * cpi->kA3 * dt;
		  rpi->Rn_Arr[n] = rpo->Rn_Arr[n] - d;		  
		  rpi->Ops       = rpo->Ops       + d;
		  rpi->Arr       = rpo->Arr       + d;
        }
/*
	d = rpi->Arr * rpi->Arr * cpi->kA4 * dt;
	rpi->Arr         = rpo->Arr         - 2 * d;
	rpi->Arr_di      = rpo->Arr_di      + d;

	d = rpi->Arr_di * cpi->kA5 * dt;
	rpi->Arr_di      = rpo->Arr_di      - d;
	rpi->Arr         = rpo->Arr         + 2 * d;
*/

	rpd->Arr         = rpi->Arr;
	rpd->Arr_di      = rpi->Arr_di;
    	impl_ab_c (3, &rpi->Arr, &rpi->Arr, &rpi->Arr_di, 
		      &rpc->Arr, &rpc->Arr, &rpc->Arr_di, cpi->kA4, cpi->kA5, dt);

/*
	d = rpi->Arr_di * rpi->Arr_di * cpi->kA4 * dt;
	rpi->Arr_di      = rpo->Arr_di      - 2 * d;
	rpi->Arr_tetra   = rpo->Arr_tetra   + d;

	d = rpi->Arr_tetra * cpi->kA5 * dt;
	rpi->Arr_tetra   = rpo->Arr_tetra   - d;
	rpi->Arr_di      = rpo->Arr_di      + 2 * d;
/* */

	rpd->Arr_di      = rpi->Arr_di;
	rpd->Arr_tetra   = rpi->Arr_tetra;
    	impl_ab_c (4, &rpi->Arr_di, &rpi->Arr_di, &rpi->Arr_tetra, 
		      &rpc->Arr_di, &rpc->Arr_di, &rpc->Arr_tetra, cpi->kA4, cpi->kA5, dt);

	for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn[n] * cpi->ktherm * dt;
	    rpi->Rn[n]    = rpo->Rn[n]    - d;
	    rpi->Ops      = rpo->Ops      + d;
	}
/*
	d = rpi->Ops * rpi->Gt * cpi->kOps * dt;
	rpi->Ops          = rpo->Ops       - d;
	rpi->Gt           = rpo->Gt        - d;
	rpi->Ops_Gt       = rpo->Ops_Gt    + d;

	d = rpi->Ops_Gt * cpi->kG2 * dt;
	rpi->Ops_Gt       = rpo->Ops_Gt    - d;
	rpi->Ops          = rpo->Ops       + d;
	rpi->Gt           = rpo->Gt        + d;
/* */
	rpd->Ops_Gt       = rpi->Ops_Gt;
	rpd->Ops          = rpi->Ops;
	rpd->Gt           = rpi->Gt;
    	impl_ab_c (5, &rpi->Ops, &rpi->Gt, &rpi->Ops_Gt, 
		      &rpc->Ops, &rpc->Gt, &rpc->Ops_Gt, cpi->kOps, cpi->kG2, dt);

/*
	d = rpi->Ops_Gt * cpi->kG3 * dt;
	rpi->Ops_Gt       = rpo->Ops_Gt    - d;
	rpi->Ops_G        = rpo->Ops_G     + d;
	d = rpi->Ops_G * cpi->kG4gdp * dt;
	rpi->Ops_G        = rpo->Ops_G     - d;
	rpi->Ops_Gt       = rpo->Ops_Gt    + d;
/* */
	rpd->Ops_G        = rpi->Ops_G;
	rpd->Ops_Gt       = rpi->Ops_Gt;
    	impl_a_b (6, &rpi->Ops_Gt, &rpi->Ops_G, 
		     &rpc->Ops_Gt, &rpc->Ops_G, cpi->kG3, cpi->kG4gdp, dt);

	d = rpi->Ops_G * cpi->kG5gtp * dt;
	rpi->Ops_G        = rpo->Ops_G      - d;
	rpi->Ops_Ggtp     = rpo->Ops_Ggtp   + d;
	
	d = rpi->Ops_Ggtp * cpi->kG6 * dt;
	rpi->Ops_Ggtp     = rpo->Ops_Ggtp   - d;
	rpi->Ops          = rpo->Ops        + d;
	rpi->Ggtp         = rpo->Ggtp       + d;

	d = rpi->Ops * cpi->kRrecyc * dt;
	rpi->Ops          = rpo->Ops        - d;
	rpi->Rh           = rpo->Rh         + d;
	if (rpi->Ops < 0) rpi->Ops = 0;
/*
  // fprintf (stderr,"R_Gt %g tc %g\n",rpi->R_Gt, tc);
 	   d = rpi->Gt * rpi->Rh * cpi->kG1_0 * 1.6 * dt;
	   rpi->Rh           = rpo->Rh        - d;
	   rpi->Gt           = rpo->Gt        - d;
	   rpi->R_Gt         = rpo->R_Gt      + d;
  // fprintf (stderr,"R_Gt %g dt %g\n",rpi->R_Gt,dt);

	   d = rpi->R_Gt * cpi->kG2 * 315 * dt;
	   rpi->R_Gt         = rpo->R_Gt      - d;
	   rpi->Rh           = rpo->Rh        + d;
	   rpi->Gt           = rpo->Gt        + d;
*/

      rpd->Rh      = rpi->Rh;
      rpd->Gt      = rpi->Gt;
      rpd->R_Gt    = rpi->R_Gt;
      impl_ab_c (7, &rpi->Rh, &rpi->Gt, &rpi->R_Gt, 
		    &rpc->Rh, &rpc->Gt, &rpc->R_Gt, cpi->kG1_0*1.6, cpi->kG2 * 315, dt);
    
   // fprintf (stderr,"R_Gt %g Rh %g \n",rpi->R_Gt, rpi->Rh);
   // fprintf (stderr,"\n",rpi->R_Gt, rpi->Rh);

/*
 	for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn[n] * rpi->Gt * cpi->kG1_0 * exp (-cpi->omGt*n) * dt;
	    rpi->Rn[n]     = rpo->Rn[n]     - d; 
	    rpi->Gt        = rpo->Gt        - d;
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  + d;
	}
	for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn_Gt[n] * cpi->kG2 * dt;
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  - d;
	    rpi->Rn[n]     = rpo->Rn[n]     + d; 
	    rpi->Gt        = rpo->Gt        + d;
	}
/* */

	for (n=0; n<NPHOS; n++) {
	    rpd->Rn_Gt[n]  = rpi->Rn_Gt[n];
	    rpd->Rn[n]     = rpi->Rn[n]; 
	    rpd->Gt        = rpi->Gt;
            impl_ab_c (8, &rpi->Rn[n], &rpi->Gt, &rpi->Rn_Gt[n], 
			  &rpc->Rn[n], &rpc->Gt, &rpc->Rn_Gt[n], cpi->kG1_0 * exp(-cpi->omGt*n), cpi->kG2, dt);
	}
/* */
/*
	for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn_Gt[n] * cpi->kG3 * dt;
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  - d;
	    rpi->Rn_G[n]   = rpo->Rn_G[n]   + d;
	 }
	 for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn_G[n] * cpi->kG4gdp * dt;
	    rpi->Rn_G[n]   = rpo->Rn_G[n]   - d;
	    rpi->Rn_Gt[n]  = rpo->Rn_Gt[n]  + d;
	}
/* */


 	for (n=0; n<NPHOS; n++) {
	    rpd->Rn_G[n]   = rpi->Rn_G[n];
	    rpd->Rn_Gt[n]  = rpi->Rn_Gt[n];
    	    impl_a_b (9, &rpi->Rn_Gt[n], &rpi->Rn_G[n], 
		         &rpc->Rn_Gt[n], &rpc->Rn_G[n], cpi->kG3, cpi->kG4gdp, dt);
	}
/* */

	for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn_G[n] * cpi->kG5gtp * dt;
	    rpi->Rn_G[n]    = rpo->Rn_G[n]    - d;
	    rpi->Rn_Ggtp[n] = rpo->Rn_Ggtp[n] + d;
	}

	for (n=0; n<NPHOS; n++) {
	    d = rpi->Rn_Ggtp[n] * cpi->kG6 * dt;
	    rpi->Rn_Ggtp[n] = rpo->Rn_Ggtp[n] - d;
	    rpi->Rn[n]      = rpo->Rn[n]      + d;
	    rpi->Ggtp       = rpo->Ggtp       + d;
	}

	d = rpi->Ggtp * cpi->kG7 * dt;
	rpi->Ggtp           = rpo->Ggtp       - d;
	rpi->Gagtp          = rpo->Gagtp      + d;
	rpi->Gbg            = rpo->Gbg        + d;

/*	
	d = rpi->PDE * rpi->Gagtp * cpi->kP1 * dt;
	rpi->PDE            = rpo->PDE        - d;
	rpi->Gagtp          = rpo->Gagtp      - d;
	rpi->PDE_Gagtp      = rpo->PDE_Gagtp  + d;
	d = rpi->PDE_Gagtp * cpi->kP1rev * dt;
	rpi->PDE_Gagtp      = rpo->PDE_Gagtp   - d;
	rpi->PDE            = rpo->PDE         + d;
	rpi->Gagtp          = rpo->Gagtp       + d;
/* */

	rpd->PDE            = rpi->PDE;
	rpd->Gagtp          = rpi->Gagtp;
	rpd->PDE_Gagtp      = rpi->PDE_Gagtp;
        impl_ab_c (10, 1, &rpi->PDE, &rpi->Gagtp, &rpi->PDE_Gagtp, 
		          &rpc->PDE, &rpc->Gagtp, &rpc->PDE_Gagtp, cpi->kP1, cpi->kP1rev, dt);
/* */

	// dark noise, check with tcomp18bd
	//  set to have s.d. of 1/5 pA, giving SNR of 5 for 1 pA SPR
	
	if (rpnti->dnois>0) {
                 if (rpnti->dstate) drand_setstate (rpnti->dstate);
		 rpi->PDE_Gagtp  += (drand() - 0.5) * 0.726 * rpnti->dnois * cpi->gain_spr_mul;
		 // rpi->PDE_Gagtp  += gasdev() * 4.55e-2 * rpnti->dnois * cpi->gain_spr_mul;
 		 restorstate();  	/* restore default random state */
	}

	d = rpi->PDE_Gagtp * cpi->kP2 * dt;
	rpi->PDE_Gagtp      = rpo->PDE_Gagtp    - d;
	rpi->PDEa_Gagtp     = rpo->PDEa_Gagtp   + d;

	d = rpi->PDEa_Gagtp * rpi->Gagtp * cpi->kP3 * dt;
	rpi->PDEa_Gagtp       = rpo->PDEa_Gagtp       - d;
	rpi->Gagtp            = rpo->Gagtp            - d;
	rpi->Gagtp_PDEa_Gagtp = rpo->Gagtp_PDEa_Gagtp + d;

	d = rpi->Gagtp_PDEa_Gagtp * cpi->kP4 * dt;
	rpi->Gagtp_PDEa_Gagtp  = rpo->Gagtp_PDEa_Gagtp  - d;
	rpi->Gagtp_aPDEa_Gagtp = rpo->Gagtp_aPDEa_Gagtp + d;

	d = rpi->RGS * rpi->Gagtp_aPDEa_Gagtp * cpi->kRGS1 * dt;	// 1st kRGS1
	rpi->RGS                   = rpo->RGS                  - d;
	rpi->Gagtp_aPDEa_Gagtp     = rpo->Gagtp_aPDEa_Gagtp    - d;
	rpi->RGS_Gagtp_aPDEa_Gagtp = rpo->RGS_Gagtp_aPDEa_Gagtp +d;

	d = rpi->RGS_Gagtp_aPDEa_Gagtp * cpi->kRGS2 * dt;
	rpi->RGS_Gagtp_aPDEa_Gagtp = rpo->RGS_Gagtp_aPDEa_Gagtp - d;
	rpi->PDEa_Gagtp            = rpo->PDEa_Gagtp            + d;
	rpi->RGS                   = rpo->RGS                   + d;
	rpi->Gagdp                 = rpo->Gagdp                 + d;

	d = rpi->RGS * rpi->PDEa_Gagtp * cpi->kRGS1 * dt;		// 2nd kRGS1
	rpi->RGS               = rpo->RGS               - d;
	rpi->PDEa_Gagtp        = rpo->PDEa_Gagtp        - d;
	rpi->RGS_PDEa_Gagtp    = rpo->RGS_PDEa_Gagtp    + d;


	d = rpi->RGS_PDEa_Gagtp * cpi->kRGS2 * dt;
	rpi->RGS_PDEa_Gagtp     = rpo->RGS_PDEa_Gagtp     - d;
	rpi->PDE                = rpo->PDE                + d;
	rpi->RGS                = rpo->RGS                + d;
	rpi->Gagdp              = rpo->Gagdp              + d;

	d = rpi->PDEa_Gagtp * cpi->kPDEshutoff * dt;
	rpi->PDEa_Gagtp         = rpo->PDEa_Gagtp         - d;
	rpi->PDE                = rpo->PDE                + d;
	rpi->Gagdp              = rpo->Gagdp              + d;

	d = rpi->Gagtp_aPDEa_Gagtp * cpi->kPDEshutoff * dt;
	rpi->Gagtp_aPDEa_Gagtp  = rpo->Gagtp_aPDEa_Gagtp  - d;
	rpi->PDEa_Gagtp         = rpo->PDEa_Gagtp         + d;
	rpi->Gagdp              = rpo->Gagdp              + d;

	d = rpi->Gagtp * cpi->kGshutoff * dt;
	rpi->Gagtp              = rpo->Gagtp              - d;
	rpi->Gagdp              = rpo->Gagdp              + d;

	d = rpi->Gagdp * rpi->Gbg * cpi->kGrecyc * dt;
	rpi->Gagdp              = rpo->Gagdp              - d;
	rpi->Gbg                = rpo->Gbg                - d;
	rpi->Gt                 = rpo->Gt                 + d;

/*
	d = rpi->Rect * rpi->Cafree * cpi->kRec1 * dt;
	rpi->Rect               = rpo->Rect               - d;
	rpi->Cafree             = rpo->Cafree             - d;
	rpi->Recr_Ca            = rpo->Recr_Ca            + d;
	d = rpi->Recr_Ca * cpi->kRec2 * dt;
	rpi->Recr_Ca            = rpo->Recr_Ca            - d;
	rpi->Rect               = rpo->Rect               + d;
	rpi->Cafree             = rpo->Cafree             + d;
/* */
	rpd->Recr_Ca            = rpi->Recr_Ca;
	rpd->Rect               = rpi->Rect;
	rpd->Cafree             = rpi->Cafree;
        impl_ab_c (11, &rpi->Rect, &rpi->Cafree, &rpi->Recr_Ca, 
	 	       &rpc->Rect, &rpc->Cafree, &rpc->Recr_Ca, cpi->kRec1, cpi->kRec2, dt);
/*
	d = rpi->Recr_Ca * rpi->RK * cpi->kRec3 * dt;
	rpi->Recr_Ca            = rpo->Recr_Ca            - d;
	rpi->RK                 = rpo->RK                 - d;
	rpi->Recr_Ca_RK         = rpo->Recr_Ca_RK         + d;
	d = rpi->Recr_Ca_RK * cpi->kRec4 * dt;
	rpi->Recr_Ca_RK         = rpo->Recr_Ca_RK         - d;
	rpi->Recr_Ca            = rpo->Recr_Ca            + d;
	rpi->RK                 = rpo->RK                 + d;
/* */

	rpd->Recr_Ca_RK         = rpi->Recr_Ca_RK;
	rpd->Recr_Ca            = rpi->Recr_Ca;
	rpd->RK                 = rpi->RK;
        impl_ab_c (12, &rpi->Recr_Ca, &rpi->RK, &rpi->Recr_Ca_RK, 
	 	       &rpc->Recr_Ca, &rpc->RK, &rpc->Recr_Ca_RK, cpi->kRec3, cpi->kRec4, dt);
/* */

/*
	d = rpi->Cafree * (cpi->eT - rpi->Cabuff) * cpi->k1 * dt;
	rpi->Cafree            = rpo->Cafree            - d;
	rpi->Cabuff            = rpo->Cabuff            + d;
	d = rpi->Cabuff * cpi->k2 * dt;
	rpi->Cabuff            = rpo->Cabuff            - d;
	rpi->Cafree            = rpo->Cafree            + d;
/* */


 	rpd->Cafree            = rpi->Cafree;
	rpd->Cabuff            = rpi->Cabuff;
    	impl_a_b (13, &rpi->Cafree, &rpi->Cabuff, 
		      &rpc->Cafree, &rpc->Cabuff, (cpi->eT - rpi->Cabuff) * cpi->k1, cpi->k2, dt);
/* */

	d = (rpi->Cafree - cpi->Ca_0) * cpi->gCa * dt;
	rpi->Cafree            = rpo->Cafree            - d;

	chan_cond = pow (rpi->cGMP / cpi->cGMPdark, cpi->ncg);

	d = 1e6 * cpi->fCa * cpi->Jdark / 
	    ((2 + cpi->fCa) * Fd * cpi->Vcyto) * chan_cond * dt;
	rpi->Cafree            = rpo->Cafree            + d;

	d = (cpi->Amax / (1 + pow(rpi->Cafree/cpi->Kc1, cpi->m1)) +
	     cpi->Amax / (1 + pow(rpi->Cafree/cpi->Kc2, cpi->m2))) * dt;
        rpi->cGMP               = rpo->cGMP               + d;

	d = (cpi->Bdark + cpi->Bsub * (rpi->PDEa_Gagtp + 2*rpi->Gagtp_aPDEa_Gagtp) + 
			                                   rpi->Gagtp_PDEa_Gagtp) * rpi->cGMP * dt;        
	rpi->cGMP              = rpo->cGMP              - d;

 	 // Ca2_frac = (rpi->Cafree - cpi->Ca_0) / (cpi->Ca_dark - cpi->Ca_0);
	 // Icgchan = 2/(2+cpi->fCa)* pow(rpi->cGMP/cpi->cGMPdark,cpi->ncg) * cpi->Jdark + 
	 //  					cpi->fCa/(cpi->fCa+2) * Ca2_frac * cpi->Jdark;
	 // cond = Icgchan * 1e-2 / 100e-3;
	 
	 cond = chan_cond / cpi->chan_eff;
      }
//  fprintf (stderr,"iflux %-7.3g Rh %-7.3g Rn0 %-7.3g R_Gt %-7.3g RK %-7.3g Arr %-7.3g Ggtp %-7.3g PDE %-7.3g Gagtp %-7.3g Cafree %-7.3g cGMP %-7.3g cond %-7.3g tc %g\n",
//	iflux, rpi->Rh, rpi->Rn[0], rpi->R_Gt, rpi->RK, rpi->Arr, rpi->Ggtp, rpi->PDE, rpi->Gagtp, rpi->Cafree, rpi->cGMP, cond, tc);
 // fprintf (stderr,"iflux %g Rn0 %g RK %g Arr %g Ggtp %g Cafree %g Cabuff %g cGMP %g Rh0 %g\n",
 // 	iflux, rpi->Rn[0], rpi->RK, rpi->Arr, rpi->Ggtp, rpi->Cafree,rpi->Cabuff,rpi->cGMP,rpi->Rh);
//  fprintf (stderr,"iflux %g Rn0 %g RK %g Arr %g Ggtp %g Cafree %g Cabuff %g cGMP %g cond %g\n",
// 	iflux, rpi->Rn[0], rpi->RK, rpi->Arr, rpi->Ggtp, rpi->Cafree,rpi->Cabuff,rpi->cGMP,cond);

	return cond;
}

/*------------------------------------*/
