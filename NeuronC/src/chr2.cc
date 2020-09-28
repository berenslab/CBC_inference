/* Module chr2 in Program nc */

/* sets up ChR2 channel-rhodopsin photoreceptor from Williams et al (2013)

Williams JC, Xu J, Lu Z, Klimas A, Chen X, Ambrosi CM, Cohen IS, Entcheva E.
(2013) Computational optogenetics: empirically-derived voltage- and
light-sensitive channelrhodopsin-2 model. PLoS Comput Biol. 2013;9(9):e1003220.
doi: 10.1371/journal.pcbi.1003220. https://www.ncbi.nlm.nih.gov/pubmed/24068903

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
int poisdev(double xm);
int binomdev(double p, int n);

/*------------------------------------*/

#define qtemp 22.0		/* base temperature for channel rate & cond */

double qcalc (double q10_val)
{
     return exp(log(q10_val) * (tempcel - qtemp) / 10.0);
}

/*------------------------------------*/

void chr2set(int pigm)

/* Set up constants for Williams (2013) channel-rhodopsin phototransuction */
/* Parameters taken from: https://www.ncbi.nlm.nih.gov/pubmed/24068903 */

/* Called by "initrec()" in ncmak.cc. */

{
    int n;
    double  qrate, sfactor;
    recparc *cpt;
    phvarsc *vpt;


  if (pigm>=NUMREC) pigm = NUMREC-1;
  if (pigm==ChR2) {
      cpt = (recparc*) rectypes[pigm];	/* set up cone cascade constants */
      cpt->stype = ChR2;	// Channel Rhodopsin 2 
  }
  else return;

  cpt = (recparc*)rectypes[pigm];  /* set up rodi cascade constants */
				/* last pigm type (pigm=NUMREC) is ksens */
  cpt->ctype = CHR;

  cpt->pigm = pigm;
  cpt->recslow = 1.0;		/* cone is slow by this factor */
  cpt->loopgain = 1;		/* gain for ca loop, sets stability */
  cpt->vrev  =  0.0;		// reversal pot for ChR channel

  cpt->qc = dqrec;		/* Q10 factor for conductance shape */
  cpt->qcond = dqcrec;		/* Q10 factor for conductance cond */
				/* not currently used */

  cpt->qt = 22.0;		/* base temperature for channel rate & cond */
  cpt->oldtempcel = tempcel;
  qrate = qcalc (cpt->qc);

  cpt->pathl = pigmlen[pigm];	/* path length through o.s. */

  cpt->q10_e1  = qcalc(1.46);
  cpt->q10_e2  = qcalc(2.77);
  cpt->q10_gd1 = qcalc(1.97);	/* from Williams et al (2013) */
  cpt->q10_gd2 = qcalc(1.77);
  cpt->q10_e12 = qcalc(1.1);	/* in dark */
  cpt->q10_e21 = qcalc(1.95);	/* in dark */
  cpt->q10_gr  = qcalc(2.56);

#define tchr2 1.3 		/* time const msec */
  cpt->tchr = tchr2;
  cpt->gamma = 0.1;

  cpt->gain_spr_mul = 1.0;	/* gain for single photon resp, sets peak ampl */
  cpt->chan_eff = 1;		

  vpt = &cpt->vars;

  vpt->SC1 = 1.0;
  vpt->SO1 = 0;
  vpt->SO2 = 0;
  vpt->SC2 = 0;
  vpt->p = 0;

}

#define runchr_euler runchr
// #define runchr_implicit runchr

/*------------------------------------*/

double runchr_euler(photrecc *rpntc, double iflux, double v)

/* channel-rhodopsin phototransduction from Williams et al (2013) */

/* Williams JC, Xu J, Lu Z, Klimas A, Chen X, Ambrosi CM, Cohen IS, Entcheva E.
   (2013) Computational optogenetics: empirically-derived voltage- and
   light-sensitive channelrhodopsin-2 model. PLoS Comput Biol. 2013;9(9):e1003220.
   doi: 10.1371/journal.pcbi.1003220. https://www.ncbi.nlm.nih.gov/pubmed/24068903
*/

/*  Channel-Rhodopsin 2 (ChR2) channel 

    Four states:

                e12
        O1     <--->     O2
                e21
      ^  |             ^ |
      k1 | Gd1        k2 | Gd2
           \/              \/
                Gr
        C1     <---     C2


        1     <--->     2

      ^ |             ^ |
        | \/            | \/

        0     <---      3
*/

{
	int i,j,k,n;
	double tc;
	double dt, vg, cond, nchans, so1, so2;
	double e12, e21, k1, k2, g, gd1, gd2, gr, s;
	double d1, d2, d3, d4, d5, d6, d7;
	double SC1old, SO1old, SO2old, SC2old, SO1_prev;
	double err, errmax;
	recparc *cpt;
	phvarsc *vpt;

     cpt = rpntc->chtyp;
     vpt = &rpntc->vars;	
     

     if      (iflux < 1) tc = 1;
     else if (iflux < 5) tc = 2;
     else if (iflux < 20) tc = 4;
     else if (iflux < 50) tc = 50;
     else if (iflux < 200) tc = 100;
     else if (iflux < 500) tc = 200;
     else if (iflux < 1000) tc = 500;
     else if (iflux < 2000) tc = 1000;
     else                   tc = 2000;

     tc = 1;

     dt = stiminc/(MSEC*tc);	// time increment 
     v *= VTOMV;		// volts calib in mV

/* See changes (below) to gd1, gd2, and gr that reduce gr rate and increase initial nose of response */
/*  for the forward-euler method now used instead of backwards euler */

#define e1 0.8535
#define e2 0.14

#define e12d 0.011
#define c1 0.005
#define c2 0.024

#define e21d 0.008
#define c3 0.004
#define c4 0.024

// #define gd2_val 0.05
#define gd2_val 0.001

     s = 0.5 * (1 + tanh(120 * (iflux*100.0 - 0.1)));
     // gd1 = (0.075 + 0.043 * tanh((v + 20) / -20)) * cpt->q10_gd1 * dt;
     gd1 = (0.025 + 0.043 * tanh((v + 20) / -20)) * cpt->q10_gd1 * dt;
     // gr  = (4.34587e5 * exp (-0.0211539274*v))    * cpt->q10_gr  * dt;
     gr  = (4.34587e2 * exp (-0.0211539274*v))    * cpt->q10_gr  * dt;
     e21 = (e21d + c3 * log(1.0 + iflux / c4))    * cpt->q10_e21 * dt;
     e12 = (e12d + c1 * log(1.0 + iflux / c2))    * cpt->q10_e12 * dt;
     gd2 = gd2_val * cpt->q10_gd2 * dt;

     vpt->p += (s - vpt->p) / cpt->tchr * dt*tc;

     k1  = e1 * iflux * 1e0 * vpt->p * cpt->q10_e1 * dt;
     k2  = e2 * iflux * 1e0 * vpt->p * cpt->q10_e2 * dt;

     for (i=0; i<tc; i++) {

	SC1old = vpt->SC1;
	SO1old = vpt->SO1;
	SO2old = vpt->SO2;
	SC2old = vpt->SC2;
	SO1_prev = vpt->SO1;

/*	// initial estimate

	vpt->SO1 += (k1 * vpt->SC1 + e21 * vpt->SO2 - (gd1 + e12) * vpt->SO1) * dt;
	vpt->SO2 += (k2 * vpt->SC2 + e12 * vpt->SO1 - (gd2 + e21) * vpt->SO2) * dt;
	vpt->SC2 += (gd2 * vpt->SO2 - k2 * vpt->SC2 - gr * vpt->SC2) * dt;
	// vpt->SC1 += gd1 * vpt->SO1 + gr * vpt->SC2 - k1 * vpt->SC1;
	vpt->SC1 = 1.0 - vpt->SO1 - vpt->SO2 - vpt->SC2;
	if (vpt->SC1 < 0) vpt->SC1 = 0;
	if (vpt->SC1 > 1) vpt->SC1 = 1;
*/

	   d1 = k1  * vpt->SC1;
	   d2 = gd1 * vpt->SO1;
	   d3 = gr  * vpt->SC2;
	   d4 = e21 * vpt->SO2;
	   d5 = e12 * vpt->SO1;
           d6 = k2  * vpt->SC2;
           d7 = gd2 * vpt->SO2;

	   // forward euler 

	   vpt->SO1 = (SO1old + (d1 + d4)) - d2 - d5;
	   if (vpt->SO1 > 1) vpt->SO1 = 1;
	   if (vpt->SO1 < 0) vpt->SO1 = 0;
	   vpt->SO2 = (SO2old + (d6 + d5)) - d7 - d4;
	   if (vpt->SO2 > 1) vpt->SO2 = 1;
	   if (vpt->SO2 < 0) vpt->SO2 = 0;
	   vpt->SC2 = (SC2old + d7) - d6 - d3;
	   if (vpt->SC2 > 1) vpt->SC2 = 1;
	   if (vpt->SC2 < 0) vpt->SC2 = 0;
	   vpt->SC1 = 1.0 - vpt->SO1 - vpt->SO2 - vpt->SC2;
	   // vpt->SC1 = SC1old + d2 + d3 - d1;
	   if (vpt->SC1 < 0) vpt->SC1 = 0;
	   if (vpt->SC1 > 1) vpt->SC1 = 1;

	  //  fprintf (stderr,"C1 %-10.5g p %-6.3g d1 %-10.6g d2d3 %-10.6g O1 %-10.5g O2 %-8.5g C2 %-8.5g err %-7.5g j %-2d tc %-3g\n",vpt->SC1,vpt->p,d1,d2+d3,vpt->SO1,vpt->SO2,vpt->SC2,err,j,tc);

//	   fprintf (stderr,"C1 %-10.5g d1 %-10.6g O1 %-10.5g O2 %-10.5g C2 %-10.5g err %-10.6g  e\n",vpt->SC1,d1,vpt->SO1,vpt->SO2,vpt->SC2,err);

     }
     // vg = v;					/* see this section in runrec() in ncsub.cc */
     // if (vg > -10) vg = -10;
     // g = (10.6408 - 14.6408 * exp(-vg/42.7671)) / vg;
     g = 1.0;
     cond = (vpt->SO1 + vpt->SO2 * cpt->gamma) * g;

// fprintf (stderr,"p k1 %g k2 %g gd1 %g gd2 %g e12 %g e21 %g gr %g iflux %g q %g cond %g C1 %g O1 %g O2 %g C2 %g\n",k1,k2,gd1,gd2*dt,e12,e21,gr,iflux,cpt->q10_e1*dt,cond,vpt->SC1,vpt->SO1,vpt->SO2,vpt->SC2);
	
	//    fprintf (stderr,"iflux %6.3g C1 %-10.5g p %-6.3g d1 %-10.6g d2d3 %-10.6g O1 %-10.5g O2 %-10.5g C2 %-10.5g err %-10.6g cond %-10.6g\n",iflux,vpt->SC1,vpt->p,d1,d2+d3,vpt->SO1,vpt->SO2,vpt->SC2,err, cond);

     return cond;
}

#undef gd2_val

/*------------------------------------*/

double runchr_implicit(photrecc *rpntc, double iflux, double v)

/* channel-rhodopsin phototransduction from Williams et al (2013) */

/* Williams JC, Xu J, Lu Z, Klimas A, Chen X, Ambrosi CM, Cohen IS, Entcheva E.
   (2013) Computational optogenetics: empirically-derived voltage- and
   light-sensitive channelrhodopsin-2 model. PLoS Comput Biol. 2013;9(9):e1003220.
   doi: 10.1371/journal.pcbi.1003220. https://www.ncbi.nlm.nih.gov/pubmed/24068903
*/

/*  Channel-Rhodopsin 2 (ChR2) channel 

    Four states:

                e12
        O1     <--->     O2
                e21
      ^  |             ^ |
      k1 | Gd1        k2 | Gd2
           \/              \/
                Gr
        C1     <---     C2


        1     <--->     2

      ^ |             ^ |
        | \/            | \/

        0     <---      3
*/

{
	int i,j,k,n;
	double tc;
	double dt, vg, cond, nchans, so1, so2;
	double e12, e21, k1, k2, g, gd1, gd2, gr, s;
	double d1, d2, d3, d4, d5, d6, d7;
	double SC1old, SO1old, SO2old, SC2old, SO1_prev;
	double err, errmax;
	recparc *cpt;
	phvarsc *vpt;

     cpt = rpntc->chtyp;
     vpt = &rpntc->vars;	
     

     if      (iflux < 1) tc = 1;
     else if (iflux < 5) tc = 2;
     else if (iflux < 20) tc = 4;
     else if (iflux < 50) tc = 50;
     else if (iflux < 200) tc = 100;
     else if (iflux < 500) tc = 200;
     else if (iflux < 1000) tc = 500;
     else if (iflux < 2000) tc = 1000;
     else                   tc = 2000;

     tc = 1;

     dt = stiminc/(MSEC*tc);	// time increment 
     v *= VTOMV;		// volts calib in mV

/* See changes (below) to gd1, gd2, and gr that reduce gr rate and increase initial nose of response */
/*  for the forward-euler method now used instead of backwards euler */

#define e1 0.8535
#define e2 0.14

#define e12d 0.011
#define c1 0.005
#define c2 0.024

#define e21d 0.008
#define c3 0.004
#define c4 0.024

// #define gd2_val 0.05
#define gd2_val 0.001

     s = 0.5 * (1 + tanh(120 * (iflux*100.0 - 0.1)));
     gd1 = (0.075 + 0.043 * tanh((v + 20) / -20)) * cpt->q10_gd1 * dt;
     // gd1 = (0.025 + 0.043 * tanh((v + 20) / -20)) * cpt->q10_gd1 * dt;
     gr  = (4.34587e5 * exp (-0.0211539274*v))    * cpt->q10_gr  * dt;
     // gr  = (4.34587e2 * exp (-0.0211539274*v))    * cpt->q10_gr  * dt;
     e21 = (e21d + c3 * log(1.0 + iflux / c4))    * cpt->q10_e21 * dt;
     e12 = (e12d + c1 * log(1.0 + iflux / c2))    * cpt->q10_e12 * dt;
     gd2 = gd2_val * cpt->q10_gd2 * dt;

     vpt->p += (s - vpt->p) / cpt->tchr * dt*tc;

     k1  = e1 * iflux * 1e0 * vpt->p * cpt->q10_e1 * dt;
     k2  = e2 * iflux * 1e0 * vpt->p * cpt->q10_e2 * dt;

     for (i=0; i<tc; i++) {

	SC1old = vpt->SC1;
	SO1old = vpt->SO1;
	SO2old = vpt->SO2;
	SC2old = vpt->SC2;
	SO1_prev = vpt->SO1;

/*	// initial estimate

	vpt->SO1 += (k1 * vpt->SC1 + e21 * vpt->SO2 - (gd1 + e12) * vpt->SO1) * dt;
	vpt->SO2 += (k2 * vpt->SC2 + e12 * vpt->SO1 - (gd2 + e21) * vpt->SO2) * dt;
	vpt->SC2 += (gd2 * vpt->SO2 - k2 * vpt->SC2 - gr * vpt->SC2) * dt;
	// vpt->SC1 += gd1 * vpt->SO1 + gr * vpt->SC2 - k1 * vpt->SC1;
	vpt->SC1 = 1.0 - vpt->SO1 - vpt->SO2 - vpt->SC2;
	if (vpt->SC1 < 0) vpt->SC1 = 0;
	if (vpt->SC1 > 1) vpt->SC1 = 1;
*/

	/* implicit calc of state variables */

	errmax = 1e-15;
	for (j=0,err=1; err > errmax && j<1000; j++) {

	   d1 = k1  * vpt->SC1;
	   d2 = gd1 * vpt->SO1;
	   d3 = gr  * vpt->SC2;
	   d4 = e21 * vpt->SO2;
	   d5 = e12 * vpt->SO1;
           d6 = k2  * vpt->SC2;
           d7 = gd2 * vpt->SO2;

 	   // backwards euler 

	   vpt->SO1 = (SO1old + (d1 + d4)) / (1 + d2 + d5);
	   vpt->SO2 = (SO2old + (d6 + d5)) / (1 + d7 + d4);
	   vpt->SC2 = (SC2old + d7) / (1 + d6 + d3);

	   // vpt->SC1 = (SC1old + (d2 + d3)) / (1 + d1);

	   vpt->SC1 = 1.0 - vpt->SO1 - vpt->SO2 - vpt->SC2;
	   if (vpt->SC1 < 0) vpt->SC1 = 0;


/*	   // not using d1 ... d7
 
	   vpt->SC1 = (SC1old + (gr * vpt->SC2 + gd1 * vpt->SO1) * dt) / (1 + k1 * vpt->SC1 * dt);
	   if (vpt->SC1 > 1) vpt->SC1 = 1;
	   vpt->SO1 = (SO1old + (k1 * vpt->SC1 + e21 * vpt->SO2) * dt) / (1 + (gd1 + e12) * vpt->SO1 * dt);
	   if (vpt->SO1 > 1) vpt->SO1 = 1;
	   vpt->SO2 = (SO2old + (k2 * vpt->SC2 + e12 * vpt->SO1) * dt) / (1 + (gd2 + e21) * vpt->SO2 * dt);
	   if (vpt->SO2 > 1) vpt->SO2 = 1;
	   vpt->SC2 = 1.0 - (vpt->SC1 + vpt->SO1 + vpt->SO2);
	   if (vpt->SC2 < 0) vpt->SC2 = 0;
*/
	   err = abs(vpt->SO1 - SO1_prev);
	   SO1_prev = vpt->SO1;
	  //  fprintf (stderr,"C1 %-10.5g p %-6.3g d1 %-10.6g d2d3 %-10.6g O1 %-10.5g O2 %-8.5g C2 %-8.5g err %-7.5g j %-2d tc %-3g\n",vpt->SC1,vpt->p,d1,d2+d3,vpt->SO1,vpt->SO2,vpt->SC2,err,j,tc);
	}

//	   fprintf (stderr,"C1 %-10.5g d1 %-10.6g O1 %-10.5g O2 %-10.5g C2 %-10.5g err %-10.6g  e\n",vpt->SC1,d1,vpt->SO1,vpt->SO2,vpt->SC2,err);

     }
     // vg = v;					/* see this section in runrec() in ncsub.cc */
     // if (vg > -10) vg = -10;
     // g = (10.6408 - 14.6408 * exp(-vg/42.7671)) / vg;
     g = 1.0;
     cond = (vpt->SO1 + vpt->SO2 * cpt->gamma) * g;

// fprintf (stderr,"p k1 %g k2 %g gd1 %g gd2 %g e12 %g e21 %g gr %g iflux %g q %g cond %g C1 %g O1 %g O2 %g C2 %g\n",k1,k2,gd1,gd2*dt,e12,e21,gr,iflux,cpt->q10_e1*dt,cond,vpt->SC1,vpt->SO1,vpt->SO2,vpt->SC2);
	
	//    fprintf (stderr,"iflux %6.3g C1 %-10.5g p %-6.3g d1 %-10.6g d2d3 %-10.6g O1 %-10.5g O2 %-10.5g C2 %-10.5g err %-10.6g cond %-10.6g\n",iflux,vpt->SC1,vpt->p,d1,d2+d3,vpt->SO1,vpt->SO2,vpt->SC2,err, cond);

     return cond;
}
#undef gd2_val

/*------------------------------------*/
