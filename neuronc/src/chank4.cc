/* segment chank4 in program nc */
/* Ih, the hyperpolarization-activated inward rectifier */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncomp.h"
#include "nconst.h"
#include "control.h"

#ifdef __cplusplus
extern "C" {
#endif

double log(double);
double exp(double);

#ifdef __cplusplus
}
#endif

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
				int nq, double bt);
void chanrate(chan *cpnt, chanparm *chp);

double al0 (chan *cpnt);
double bet0(chan *cpnt);

static double (*aln) (chan*) = al0;
static double (*betn)(chan*) = bet0;

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);

#ifdef __cplusplus
}
#endif

double ncabs(double x);
int set_perms (chantype *chtyp, iontab *ions, double svrev);

/*----------------------------------------*/

double calck4n (double v, int func)

/* Calculate K4 rate functions given voltage in mv.
/* for Ih, the hyperpolarization-activated inward rectifier */

/* Approximately matches behavior of channels in 

Hestrin, S. (1987)  J.  Physiol.  390, 319-333.
Maricq AV, and Korenbrot JI, (1990) J. Neurophysiol. 64:1917-28

*/

/* See tcomp31d for a test of this channel */

/* Calc rate in terms of 1/sec. 

   The "func" parameter defines:

    1	alpha n 
    2	beta  n 
*/

#define NUMKFUNC 2

{
   double val,x,y;

  if (func > NUMKFUNC) func = NUMKFUNC;
  switch (func) {				/* alpha functions */
				
  case 1:					/* alpha n */
    y = -0.3 * (-v-70.); 
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* singularity when x==0 */
          val =  y / x;
    else  val = 1.0;
    val *= 0.3;
    break;
						/* beta n */
  case 2:               
    val = 0.15 * exp ((-v-50) / -30.);

    break;
  }

  return val * HHRATE; 
}

/*----------------------------------------*/

chantype *makk4(void)

/* The inward rectifier in rod and cone inner segments: Ih */

{
  
 
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 3;                          /* 3 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(K,4,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dkihu;
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0.48;		/* permeability to Na+ ions */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions */
   ch->vrev = -0.03;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck4n;         /* default rate function */
   parm[0].funcname = (char *)"calck4n"; /* user rate function */
   if (dqn==0 && ((dqna+dqnb)>0))
        parm[0].dq[0] = (dqna+dqnb)*0.5;/* Q10 for n rate function */
   else parm[0].dq[0] = dqn;

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = aln;
   spnt[0].ratemul[0] = 1.0*n;
   spnt[0].rateo  [0] = 0;
                                     /*   0 <-> 1 <-> 2   */
   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = aln;
   spnt[1].ratemul[0] = 0.5*n;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = betn;
   spnt[1].ratemul[1] = 2.5*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 1;            /* state 2 = the open state */
   spnt[2].cond       = 1.0;
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = betn;
   spnt[2].ratemul[0] = 2.5*n;
   spnt[2].rateo  [0] = 1;

   return ch;
}

/*-----------------------------------------------------------------*/
   
/* Rate functions are taken from                             */
/* 	Altomare et al., 2001, J. Gen Phyhsiol 117: 519-532. */

/* For the purpose of normalizing kinetics to 22 deg C, */
/*  we assume Q10=2 */

/* #define K8A27RATE  (exp(log(2.0) * (BASETC- 27.0)/10.0); */
#define    K8A27RATE  0.70710678

/* #define K8G27RATE  (exp(log(2.0) * (BASETC- 27.0)/10.0); */
#define    K8G27RATE  0.70710678

#define    K8RTF  (RF*300*MVOLT)		// RT/F = 25.85 mV at 300 deg K

/* The "A" factor is change of one voltage sensor from the reluctant to willing state */

#define A 0.2

// #define F (1/sqrt(A))
#define F  2.236068
#define FF 5
#define F3 11.18034
#define F4 25

/* HCN1: */

#define a_hcn1  0.01406
#define b_hcn1  743.6
#define g_hcn1  2.289
#define d_hcn1  74.36
#define zb_hcn1 0.9369
#define zd_hcn1 1.018
#define za_hcn1 (-zb_hcn1)
#define zg_hcn1 (-zd_hcn1)

/* HCN2: */

#define a_hcn2  0.001308
#define b_hcn2  208.4 
#define g_hcn2  0.04214 
#define d_hcn2  86.66 
#define zb_hcn2 1.066 
#define zd_hcn2 1.021 
#define za_hcn2 (-zb_hcn2)
#define zg_hcn2 (-zd_hcn2)

/* HCN3: */

#define a_hcn3  a_hcn2 
#define b_hcn3  b_hcn2 
#define g_hcn3  g_hcn2 
#define d_hcn3  d_hcn2 
#define zb_hcn3 zb_hcn2
#define zd_hcn3 zd_hcn2
#define za_hcn3 za_hcn2
#define zg_hcn3 zg_hcn2

/* HCN4: */

#define a_hcn4  0.000496
#define b_hcn4  20.25 
#define g_hcn4  0.08444 
#define d_hcn4  6.196
#define zb_hcn4 0.8042
#define zd_hcn4 0.8934 
#define za_hcn4 (-zb_hcn4)
#define zg_hcn4 (-zd_hcn4)

/*--------------------------------------------*/

double calck8ab1 (double v, int func)

/* Calculate K type 8 (HCN) rate functions alpha, beta.
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 10) func = 10;
  switch (func) {

    case 1:  val = a_hcn1      * exp ( za_hcn1/K8RTF * v );  break; /* a */
    case 2:  val = a_hcn1 * F  * exp ( za_hcn1/K8RTF * v );  break; /* a1 */
    case 3:  val = a_hcn1 * FF * exp ( za_hcn1/K8RTF * v );  break; /* a2 */
    case 4:  val = a_hcn1 * F3 * exp ( za_hcn1/K8RTF * v );  break; /* a3 */
    case 5:  val = a_hcn1 * F4 * exp ( za_hcn1/K8RTF * v );  break; /* a4 */

    case 6:  val = b_hcn1      * exp ( zb_hcn1/K8RTF * v );  break; /* b */
    case 7:  val = b_hcn1 / F  * exp ( zb_hcn1/K8RTF * v );  break; /* b1 */
    case 8:  val = b_hcn1 / FF * exp ( zb_hcn1/K8RTF * v );  break; /* b2 */
    case 9:  val = b_hcn1 / F3 * exp ( zb_hcn1/K8RTF * v );  break; /* b3 */
    case 10: val = b_hcn1 / F4 * exp ( zb_hcn1/K8RTF * v );  break; /* b4 */

  }
  return val * K8A27RATE;
}

/*--------------------------------------------*/

double calck8gd1 (double v, int func)

/* Calculate K type 8 (HCN) rate functions gamma, delta. 
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 4) func = 4;
  switch (func) {

    case 1: val = g_hcn1     * exp (zg_hcn1/K8RTF * v); break; /* g */
    case 2: val = g_hcn1 * F * exp (zg_hcn1/K8RTF * v); break; /* gf */
    case 3: val = d_hcn1     * exp (zd_hcn1/K8RTF * v); break; /* d */
    case 4: val = d_hcn1 / F * exp (zd_hcn1/K8RTF * v); break; /* d/f */

  }
  return val * K8G27RATE;
}

/*--------------------------------------------*/

double calck8ab2 (double v, int func)

/* Calculate K type 8 (HCN) rate functions alpha, beta.
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 10) func = 10;
  switch (func) {

    case 1:  val = a_hcn2      * exp ( za_hcn2/K8RTF * v );  break; /* a */
    case 2:  val = a_hcn2 * F  * exp ( za_hcn2/K8RTF * v );  break; /* a1 */
    case 3:  val = a_hcn2 * FF * exp ( za_hcn2/K8RTF * v );  break; /* a2 */
    case 4:  val = a_hcn2 * F3 * exp ( za_hcn2/K8RTF * v );  break; /* a3 */
    case 5:  val = a_hcn2 * F4 * exp ( za_hcn2/K8RTF * v );  break; /* a4 */

    case 6:  val = b_hcn2      * exp ( zb_hcn2/K8RTF * v );  break; /* b */
    case 7:  val = b_hcn2 / F  * exp ( zb_hcn2/K8RTF * v );  break; /* b1 */
    case 8:  val = b_hcn2 / FF * exp ( zb_hcn2/K8RTF * v );  break; /* b2 */
    case 9:  val = b_hcn2 / F3 * exp ( zb_hcn2/K8RTF * v );  break; /* b3 */
    case 10: val = b_hcn2 / F4 * exp ( zb_hcn2/K8RTF * v );  break; /* b4 */

  }
  return val * K8A27RATE;
}

/*--------------------------------------------*/

double calck8gd2 (double v, int func)

/* Calculate K type 8 (HCN) rate functions gamma, delta. 
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 4) func = 4;
  switch (func) {

    case 1: val = g_hcn2     * exp (zg_hcn2/K8RTF * v); break; /* g */
    case 2: val = g_hcn2 * F * exp (zg_hcn2/K8RTF * v); break; /* gf */
    case 3: val = d_hcn2     * exp (zd_hcn2/K8RTF * v); break; /* d */
    case 4: val = d_hcn2 / F * exp (zd_hcn2/K8RTF * v); break; /* d/f */

  }
  return val * K8G27RATE;
}


/*--------------------------------------------*/

double calck8ab3 (double v, int func)

/* Calculate K type 8 (HCN) rate functions alpha, beta.
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 10) func = 10;
  switch (func) {

    case 1:  val = a_hcn3      * exp ( za_hcn3/K8RTF * v );  break; /* a */
    case 2:  val = a_hcn3 * F  * exp ( za_hcn3/K8RTF * v );  break; /* a1 */
    case 3:  val = a_hcn3 * FF * exp ( za_hcn3/K8RTF * v );  break; /* a2 */
    case 4:  val = a_hcn3 * F3 * exp ( za_hcn3/K8RTF * v );  break; /* a3 */
    case 5:  val = a_hcn3 * F4 * exp ( za_hcn3/K8RTF * v );  break; /* a4 */

    case 6:  val = b_hcn3      * exp ( zb_hcn3/K8RTF * v );  break; /* b */
    case 7:  val = b_hcn3 / F  * exp ( zb_hcn3/K8RTF * v );  break; /* b1 */
    case 8:  val = b_hcn3 / FF * exp ( zb_hcn3/K8RTF * v );  break; /* b2 */
    case 9:  val = b_hcn3 / F3 * exp ( zb_hcn3/K8RTF * v );  break; /* b3 */
    case 10: val = b_hcn3 / F4 * exp ( zb_hcn3/K8RTF * v );  break; /* b4 */

  }
  return val * K8A27RATE;
}

/*--------------------------------------------*/

double calck8gd3 (double v, int func)

/* Calculate K type 8 (HCN) rate functions gamma, delta. 
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 4) func = 4;
  switch (func) {

    case 1: val = g_hcn3     * exp (zg_hcn3/K8RTF * v); break; /* g */
    case 2: val = g_hcn3 * F * exp (zg_hcn3/K8RTF * v); break; /* gf */
    case 3: val = d_hcn3     * exp (zd_hcn3/K8RTF * v); break; /* d */
    case 4: val = d_hcn3 / F * exp (zd_hcn3/K8RTF * v); break; /* d/f */

  }
  return val * K8G27RATE;
}


/*--------------------------------------------*/

double calck8ab4 (double v, int func)

/* Calculate K type 8 (HCN) rate functions alpha, beta.
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 10) func = 10;
  switch (func) {

    case 1:  val = a_hcn4      * exp ( za_hcn4/K8RTF * v );  break; /* a */
    case 2:  val = a_hcn4 * F  * exp ( za_hcn4/K8RTF * v );  break; /* a1 */
    case 3:  val = a_hcn4 * FF * exp ( za_hcn4/K8RTF * v );  break; /* a2 */
    case 4:  val = a_hcn4 * F3 * exp ( za_hcn4/K8RTF * v );  break; /* a3 */
    case 5:  val = a_hcn4 * F4 * exp ( za_hcn4/K8RTF * v );  break; /* a4 */

    case 6:  val = b_hcn4      * exp ( zb_hcn4/K8RTF * v );  break; /* b */
    case 7:  val = b_hcn4 / F  * exp ( zb_hcn4/K8RTF * v );  break; /* b1 */
    case 8:  val = b_hcn4 / FF * exp ( zb_hcn4/K8RTF * v );  break; /* b2 */
    case 9:  val = b_hcn4 / F3 * exp ( zb_hcn4/K8RTF * v );  break; /* b3 */
    case 10: val = b_hcn4 / F4 * exp ( zb_hcn4/K8RTF * v );  break; /* b4 */

  }
  return val * K8A27RATE;
}

/*--------------------------------------------*/

double calck8gd4 (double v, int func)

/* Calculate K type 8 (HCN) rate functions gamma, delta. 
   Rate functions are taken from Altomare et al., 2001,
   J. Gen Phyhsiol 117: 519-532. */

/* Calibrated in mV, per sec */
/* Temperature 27 deg C  (300 deg K) */

{
   double val;

  if (func > 4) func = 4;
  switch (func) {

    case 1: val = g_hcn4     * exp (zg_hcn4/K8RTF * v); break; /* g */
    case 2: val = g_hcn4 * F * exp (zg_hcn4/K8RTF * v); break; /* gf */
    case 3: val = d_hcn4     * exp (zd_hcn4/K8RTF * v); break; /* d */
    case 4: val = d_hcn4 / F * exp (zd_hcn4/K8RTF * v); break; /* d/f */

  }
  return val * K8G27RATE;
}

/*--------------------------------------------*/

double k8a(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8a1(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8a2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8a3(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8a4(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[4];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8b(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[5];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8b1(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[6];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8b2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[7];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8b3(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[8];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8b4(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[9];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8g(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8gf(chan *cpnt)    /* g*f */

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8d(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double k8df(chan *cpnt)    /* d/f */

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*----------------------------------------*/

chantype *makHCN(int k_i, int hcn_i)

/* The HCN channel: Ih */

/*   Kinetic scheme and rate functions are taken from         */
/*       Altomare et al., 2001, J. Gen PHyhsiol 117: 519-532. */

/* Description of Altomare model using Neuron simulator:
   http://paynesnotebook.net/Research/NeuronSimulation/Notes */

{
  
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 10;                         /* 10 Markov states */
   nparm = 2;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(K,k_i,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dhcnu;
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0.48;		/* permeability to Na+ ions */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions */
   ch->vrev = -0.03;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 10;               /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   // begin joesterle
   if (hcn_i == 1) {
     //fprintf(stderr, "Case 1 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab1;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab1"; /* user rate function */
   } else if (hcn_i == 2) {
     //fprintf(stderr, "Case 2 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab2;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab2"; /* user rate function */
  } else if (hcn_i == 3) {
     //fprintf(stderr, "Case 3 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab3;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab3"; /* user rate function */
  } else if (hcn_i == 4) {
     //fprintf(stderr, "Case 4 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab4;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab4"; /* user rate function */
   }  else {
      fprintf(stderr, "HCN%d not implemented!\n", hcn_i);
   }
   parm[0].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[0].voff = 0;                    /* voltage offset from user */

   parm[1].nfval    = 4;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0;                /* don't use implicit vals (yet) */
   if (hcn_i == 1) {
      parm[1].chancalc =  calck8gd1;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd1"; /* user rate function */
    } else if (hcn_i == 2) {
      parm[1].chancalc =  calck8gd2;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd2"; /* user rate function */
    } else if (hcn_i == 3) {
      parm[1].chancalc =  calck8gd3;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd3"; /* user rate function */
    } else if (hcn_i == 4) {
      parm[1].chancalc =  calck8gd4;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd4"; /* user rate function */
    } else {
      fprintf(stderr, "HCN%d not implemented!\n", hcn_i); 
   }
   // end joesterle
   parm[1].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[1].voff = 0;                    /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = k8g;
   spnt[0].ratemul[0] = 4.0*n;
   spnt[0].rateo  [0] = 1;
   spnt[0].trans  [1] = 5;
   spnt[0].trate  [1] = k8a;
   spnt[0].ratemul[1] = 1.0*n;
   spnt[0].rateo  [1] = 0;

				     /*       4g      3g      2g      g                      */
                                     /*   0 <---> 1 <---> 2 <---> 3 <---> 4   closed states  */
				     /*       d       2d      3d      4d                     */
				     /*  /|      /|      /|      /|      /|                  */
				     /*   |       |       |       |       |                  */
				     /*  b|a    b1|a1   b2|a2   b3|a3   b4|a4                */
				     /*   |       |       |       |       |                  */
				     /*   |/      |/      |/      |/      |/                 */
				     /*      4gf     3gf      2gf     gf                     */
                                     /*   5 <---> 6 <---> 7  <---> 8 <---> 9   open states   */
				     /*      d/f    2d/f     3d/f    4d/f                    */
   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = k8g;
   spnt[1].ratemul[0] = 3.0*n;
   spnt[1].rateo  [0] = 1;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = k8d;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 6;
   spnt[1].trate  [2] = k8a1;
   spnt[1].ratemul[2] = 1.0*n;
   spnt[1].rateo  [2] = 0;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = k8g;
   spnt[2].ratemul[0] = 2.0*n;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = k8d;
   spnt[2].ratemul[1] = 2.0*n;
   spnt[2].rateo  [1] = 1;
   spnt[2].trans  [2] = 7;
   spnt[2].trate  [2] = k8a2;
   spnt[2].ratemul[2] = 1.0*n;
   spnt[2].rateo  [2] = 0;

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;	  	
   spnt[3].trate  [0] = k8g;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = k8d;
   spnt[3].ratemul[1] = 3.0*n;
   spnt[3].rateo  [1] = 1;
   spnt[3].trans  [2] = 8;
   spnt[3].trate  [2] = k8a3;
   spnt[3].ratemul[2] = 1.0*n;
   spnt[3].rateo  [2] = 0;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;	  	
   spnt[4].trate  [0] = k8d;
   spnt[4].ratemul[0] = 4.0*n;
   spnt[4].rateo  [0] = 1;
   spnt[4].trans  [1] = 9;
   spnt[4].trate  [1] = k8a4;
   spnt[4].ratemul[1] = 1.0*n;
   spnt[4].rateo  [1] = 0;
				     /*       4g      3g      2g      g                      */
                                     /*   0 <---> 1 <---> 2 <---> 3 <---> 4   closed states  */
				     /*       d       2d      3d      4d                     */
				     /*  /|      /|      /|      /|      /|                  */
				     /*   |       |       |       |       |                  */
				     /*  b|a    b1|a1   b2|a2   b3|a3   b4|a4                */
				     /*   |       |       |       |       |                  */
				     /*   |/      |/      |/      |/      |/                 */
				     /*      4gf     3gf      2gf     gf                     */
                                     /*   5 <---> 6 <---> 7  <---> 8 <---> 9   open states   */
				     /*      d/f    2d/f     3d/f    4d/f                    */

   spnt[5].numtrans   = 2;		/* open state */
   spnt[5].cond       = 1;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = k8gf;
   spnt[5].ratemul[0] = 4.0*n;
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 0;
   spnt[5].trate  [1] = k8b;
   spnt[5].ratemul[1] = 1.0*n;
   spnt[5].rateo  [1] = 0;

   spnt[6].numtrans   = 3;		/* open state */
   spnt[6].cond       = 1;
   spnt[6].trans  [0] = 7;	  	
   spnt[6].trate  [0] = k8gf;
   spnt[6].ratemul[0] = 3.0*n;
   spnt[6].rateo  [0] = 1;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = k8df;
   spnt[6].ratemul[1] = 1.0*n;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 1;
   spnt[6].trate  [2] = k8b1;
   spnt[6].ratemul[2] = 1.0*n;
   spnt[6].rateo  [2] = 0;

   spnt[7].numtrans   = 3;		/* open state */
   spnt[7].cond       = 1;
   spnt[7].trans  [0] = 8;	  	
   spnt[7].trate  [0] = k8gf;
   spnt[7].ratemul[0] = 2.0*n;
   spnt[7].rateo  [0] = 1;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = k8df;
   spnt[7].ratemul[1] = 2.0*n;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 2;
   spnt[7].trate  [2] = k8b2;
   spnt[7].ratemul[2] = 1.0*n;
   spnt[7].rateo  [2] = 0;

   spnt[8].numtrans   = 3;		/* open state */
   spnt[8].cond       = 1;
   spnt[8].trans  [0] = 9;	  	
   spnt[8].trate  [0] = k8gf;
   spnt[8].ratemul[0] = 1.0*n;
   spnt[8].rateo  [0] = 1;
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = k8df;
   spnt[8].ratemul[1] = 3.0*n;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 3;
   spnt[8].trate  [2] = k8b3;
   spnt[8].ratemul[2] = 1.0*n;
   spnt[8].rateo  [2] = 0;

   spnt[9].numtrans   = 2;		/* open state */
   spnt[9].cond       = 1;
   spnt[9].trans  [0] = 8;	  	
   spnt[9].trate  [0] = k8df;
   spnt[9].ratemul[0] = 4.0*n;
   spnt[9].rateo  [0] = 1;
   spnt[9].trans  [1] = 4;
   spnt[9].trate  [1] = k8b4;
   spnt[9].ratemul[1] = 1.0*n;
   spnt[9].rateo  [1] = 0;

   return ch;
}

/*----------------------------------------*/

chantype *makk9(void)

/* The HCN channel: Ih */

/*   Kinetic scheme and rate functions are taken from         */
/*       Altomare et al., 2001, J. Gen PHyhsiol 117: 519-532. */

/*   Modified for 6 states according to Michael van Rijssel's model */

/* Description of Altomare model using Neuron simulator:
   http://paynesnotebook.net/Research/NeuronSimulation/Notes */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 6;                          /*  Markov states */
   nparm = 2;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(K,9,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dkihu;
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0.48;		/* permeability to Na+ ions */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions */
   ch->vrev = -0.03;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 10;               /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calck8ab1;       /* default rate function */
   parm[0].funcname = (char *)"calck8ab1"; /* user rate function */
   parm[0].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[0].voff = 0;                    /* voltage offset from user */

   parm[1].nfval    = 4;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0;                /* don't use implicit vals (yet) */
   parm[1].chancalc =  calck8gd1;       /* default rate function */
   parm[1].funcname = (char *)"calck8gd1"; /* user rate function */
   parm[1].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[1].voff = 0;                    /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = k8g;
   spnt[0].ratemul[0] = 4.0*n;
   spnt[0].rateo  [0] = 1;

				     /*       4g      3g      2g      g                      */
                                     /*   0 <---> 1 <---> 2 <---> 3 <---> 4   closed states  */
				     /*       d       2d      3d      4d                     */
				     /*                                  /|                  */
				     /*                                   |                  */
				     /*                                 b4|a4                */
				     /*                                   |                  */
				     /*                                   |/                 */
				     /*                                                      */
                                     /*                                   9   open state     */
				     /*                                                      */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = k8g;
   spnt[1].ratemul[0] = 3.0*n;
   spnt[1].rateo  [0] = 1;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = k8d;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = k8g;
   spnt[2].ratemul[0] = 2.0*n;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = k8d;
   spnt[2].ratemul[1] = 2.0*n;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;	  	
   spnt[3].trate  [0] = k8g;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = k8d;
   spnt[3].ratemul[1] = 3.0*n;
   spnt[3].rateo  [1] = 1;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;	  	
   spnt[4].trate  [0] = k8d;
   spnt[4].ratemul[0] = 4.0*n;
   spnt[4].rateo  [0] = 1;
   spnt[4].trans  [1] = 5;
   spnt[4].trate  [1] = k8a4;
   spnt[4].ratemul[1] = 1.0*n;
   spnt[4].rateo  [1] = 0;
				     /*       4g      3g      2g      g                      */
                                     /*   0 <---> 1 <---> 2 <---> 3 <---> 4   closed states  */
				     /*       d       2d      3d      4d                     */
				     /*                                  /|                  */
				     /*                                   |                  */
				     /*                                 b4|a4                */
				     /*                                   |                  */
				     /*                                   |/                 */
				     /*                                                      */
                                     /*                                   9   open state     */
				     /*                                                      */

   spnt[5].numtrans   = 1;		/* open state */
   spnt[5].cond       = 1;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = k8b4;
   spnt[5].ratemul[0] = 1.0*n;
   spnt[5].rateo  [0] = 0;

   return ch;
}

/*----------------------------------------*/

chantype *makHCN_6ss(int k_i, int hcn_i)

/* The HCN channel: Ih */

/*   Kinetic scheme and rate functions are taken from         */
/*       Altomare et al., 2001, J. Gen PHyhsiol 117: 519-532. */

/*   Modified for 6 states according to Michael van Rijssel's model */

/* Description of Altomare model using Neuron simulator:
   http://paynesnotebook.net/Research/NeuronSimulation/Notes */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 6;                          /*  Markov states */
   nparm = 2;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(K,k_i,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dhcnu;
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0.48;		/* permeability to Na+ ions */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions */
   ch->vrev = -0.03;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 10;               /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   
   // begin joesterle
   if (hcn_i == 1) {
     //fprintf(stderr, "Case 1 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab1;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab1"; /* user rate function */
   } else if (hcn_i == 2) {
     //fprintf(stderr, "Case 2 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab2;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab2"; /* user rate function */
  } else if (hcn_i == 3) {
     //fprintf(stderr, "Case 3 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab3;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab3"; /* user rate function */
  } else if (hcn_i == 4) {
     //fprintf(stderr, "Case 4 <- %d\n", hcn_i); 
     parm[0].chancalc =  calck8ab4;       /* default rate function */
     parm[0].funcname = (char *)"calck8ab4"; /* user rate function */
   }  else {
      fprintf(stderr, "HCN%d not implemented!\n", hcn_i);
   }
   parm[0].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[0].voff = 0;                    /* voltage offset from user */

   parm[1].nfval    = 4;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0;                /* don't use implicit vals (yet) */
   if (hcn_i == 1) {
      parm[1].chancalc =  calck8gd1;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd1"; /* user rate function */
    } else if (hcn_i == 2) {
      parm[1].chancalc =  calck8gd2;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd2"; /* user rate function */
    } else if (hcn_i == 3) {
      parm[1].chancalc =  calck8gd3;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd3"; /* user rate function */
    } else if (hcn_i == 4) {
      parm[1].chancalc =  calck8gd4;       /* default rate function */
      parm[1].funcname = (char *)"calck8gd4"; /* user rate function */
    } else {
      fprintf(stderr, "HCN%d not implemented!\n", hcn_i); 
   }
   // end joesterle

   parm[1].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[1].voff = 0;                    /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = k8g;
   spnt[0].ratemul[0] = 4.0*n;
   spnt[0].rateo  [0] = 1;

				     /*       4g      3g      2g      g                      */
                                     /*   0 <---> 1 <---> 2 <---> 3 <---> 4   closed states  */
				     /*       d       2d      3d      4d                     */
				     /*                                  /|                  */
				     /*                                   |                  */
				     /*                                 b4|a4                */
				     /*                                   |                  */
				     /*                                   |/                 */
				     /*                                                      */
                                     /*                                   9   open state     */
				     /*                                                      */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = k8g;
   spnt[1].ratemul[0] = 3.0*n;
   spnt[1].rateo  [0] = 1;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = k8d;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = k8g;
   spnt[2].ratemul[0] = 2.0*n;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = k8d;
   spnt[2].ratemul[1] = 2.0*n;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;	  	
   spnt[3].trate  [0] = k8g;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = k8d;
   spnt[3].ratemul[1] = 3.0*n;
   spnt[3].rateo  [1] = 1;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;	  	
   spnt[4].trate  [0] = k8d;
   spnt[4].ratemul[0] = 4.0*n;
   spnt[4].rateo  [0] = 1;
   spnt[4].trans  [1] = 5;
   spnt[4].trate  [1] = k8a4;
   spnt[4].ratemul[1] = 1.0*n;
   spnt[4].rateo  [1] = 0;
				     /*       4g      3g      2g      g                      */
                                     /*   0 <---> 1 <---> 2 <---> 3 <---> 4   closed states  */
				     /*       d       2d      3d      4d                     */
				     /*                                  /|                  */
				     /*                                   |                  */
				     /*                                 b4|a4                */
				     /*                                   |                  */
				     /*                                   |/                 */
				     /*                                                      */
                                     /*                                   9   open state     */
				     /*                                                      */

   spnt[5].numtrans   = 1;		/* open state */
   spnt[5].cond       = 1;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = k8b4;
   spnt[5].ratemul[0] = 1.0*n;
   spnt[5].rateo  [0] = 0;

   return ch;
}


/*----------------------------------------*/

chantype *makk8(void)  { return makHCN_6ss(8, 1); }
chantype *makk11(void) { return makHCN_6ss(11, 1); }
chantype *makk12(void) { return makHCN_6ss(12, 2); }
chantype *makk13(void) { return makHCN_6ss(13, 3); }
chantype *makk14(void) { return makHCN_6ss(14, 4); }

/*----------------------------------------*/

#undef F
#undef FF
#undef F3
#undef F4

#undef A

#undef a1
#undef b1
#undef g1
#undef d1
#undef zb1
#undef zd1
#undef za1
#undef zg1

#undef a2
#undef b2
#undef g2
#undef d2
#undef zb2
#undef zd2
#undef za2
#undef zg2

#undef a4
#undef b4
#undef g4
#undef d4
#undef zb4
#undef zd4
#undef za4
#undef zg4

