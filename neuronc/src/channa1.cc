/* segment channa1 in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

chantype*makchantype(int ctype,int cnum,int nstates,int nparm,
			int nq, double dbasetc);
chantype *getchantype(int ctype, int cnum);
void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
                                double (*rrate)(chan *cpnt), double rate);
double rchanf(chan *cpnt);
double rchanr(chan *cpnt);
 
double al0 (chan *cpnt);
double bet0(chan *cpnt);
double al1 (chan *cpnt);
double bet1(chan *cpnt);

static double (*alm) (chan*) = al0;
static double (*betm)(chan*) = bet0;
static double (*alh) (chan*) = al1;
static double (*beth)(chan*) = bet1;

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

double ncabs(double x);

/*----------------------------------------*/

/*

The default Q10 for the unitary conductance of channels
is to 1.4, according to what has been measured:

Milburn T. Saint DA. Chung SH.   The temperature dependence of
conductance of the sodium channel: implications for mechanisms of
ion permeation.  Receptors & Channels. 3(3):201-211, 1995.

     Voltage-gated sodium channel currents were recorded in
cell-attached and inside-out membrane patches from rat
ventricular myocytes at temperatures ranging from 4 degrees C to
36 degrees C. The effects of changes in temperature on channel
conductance were accurately determined using a digital signal
processing technique based on hidden Markov models. We show that
the cardiac sodium channel has multiple conductance levels, with
the most frequently observed sublevel being either 11 pS or 22 pS
at room temperature in 280 mM Na+, depending on the dissociation
procedures adopted. The single channel current-voltage
relationship is ohmic at all of the temperatures studied. The
conductance increases steeply with temperature, with Q10 ranging
from 1.4 to 1.5. The proportional change in channel conductance
with increasing temperature was greater than the corresponding
change in bulk conductivity of electrolyte solutions, suggesting
that an ion traversing the channel needs to surmount a small
additional energy barrier. An activation energy deduced from a
plot of the logarithm of single channel conductance against the
inverse of temperature is about 28 kJ mole-1. We provide one
possible interpretation of the observed conductance-temperature
relationship in terms of the details of the microscopic
interactions operating between the protein wall, ions and water
molecules. 

*/


/*----------------------------------------*/

/* Vo is used here to define the resting potential
   that Hodgkin and Huxley (1952) used as a base 
   for their rate constant functions.  Vo is only
   used to normalize the modern definition of 
   membrane potential to their definition.
   "Vo" should not be modified as it refers only to
   the antiquated definition of membrane voltage: */

#define Vo (-65.)

/* The transformation from HH voltage to the standard
   used here is:

   new V =   - (oldV-Vo)
or
   new V =   (Vo - oldV)

*/

/*----------------------------------------*/

double calcna1m (double v, int func)

/* Calculate Na rate functions given voltage in mv.
   All rates are calculated exactly as in HH (1952) paper.
   Original rates were 1/msec, we multiply by 1000 here to 
   convert to 1/sec.

   The "func" parameter defines:

    1	alpha m 
    2	beta  m
*/

{
   double val,x,y;

  switch (func) {

  case 1:					/* alpha m */

/*    y = 0.1 * (v+25.); 			/* the old way */
    y = -0.1 * (v+40.);				/* the modern way */
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)			/* singularity when x==0 */
       val = y / x * MSSEC;
    else
       val = 1.0 * MSSEC;
    break;

  case 2:					/* beta m */
 /* val =  4 * exp (v / 18.);			/* the old way */
    val =  MSSEC * 4 * exp ((v-Vo) / -18.);		/* the modern way */
    break;
  }
  return val * dratehhm; 			/* normalize to dbasetc */
}						/*  i.e. 22 deg C */

/* - - - - - - - - - - - - - - - - - - -*/

double calcna1h (double v, int func)

/* Calculate Na rate functions given voltage in mv.
   All rates are calculated exactly as in HH (1952) paper.
   The "func" parameter defines:

    1	alpha h 
    2	beta  h
*/
 
{
   double val,x,y;

  switch (func) {

  case 1:					/* alpha h */
 /* val =        0.07 * exp (v / 20.);		/* old way */
    val =  MSSEC*0.07 * exp ((v-Vo) / -20.);	/* modern way */
    break;

  case 2:					/* beta h */
 /* y =  0.1 * (v+30.);				/* old way */
    y = -0.1 * (v+35.);				/* modern way */
    val = MSSEC * 1.0 / (exp (y) + 1.0);
    break;

   default: val = 0.0; break;
  }
  return val * dratehhh; 			/* normalize to dbasetc */
}						/*    i.e. 22 deg C */

/*--------------------------------------------*/

chantype *makna0(void)

/* Make channel with kinetics of Hodgkin-Huxley m3h description.  */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 2;				/* 2 Markov states */
   nparm = 2;				/* make 2 sets of params, m, h */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(NA,0,nstate,nparm,nq,dbasetc);/* make chan state const,parm */
   ch->hh = 1;				/* HH channel = uses m, h */
   ch->unitary = dnau;			/* conductance measured at dbasetc. */
					/* 32 pS @33 deg, Lipton & Tauk 1987 */
					/* = 22 pS @22 deg, Q10=1.4 */

   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / perm Na */
   ch->ions->ionp[PCA] = dpcana;		/* permeability to Ca++ ions/perm Na */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 2; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcna1m; 	/* default rate function */
   parm[0].funcname = (char *)"calcna0m"; /* user rate function */
   parm[0].dq[0]    = dqm; 		/* Q10 for m rate function */
   parm[0].voff = dnaoffsm; 		/* voltage offset from user */

   parm[1].nfval    = 2; 		/* number of func vals (ah, bh, etc.) */
   parm[1].nival    = 2; 		/* number of implicit vals (h1, h2) */
   parm[1].chancalc =  calcna1h;	/* default rate function */ 
   parm[1].funcname = (char *)"calcna0h"; /* user rate function */
   parm[1].dq[0]    = dqh; 		/* Q10 for h rate function */
   parm[1].voff = dnaoffsh; 		/* voltage offset from user */

   mak2state(spnt,rchanf,rchanr,CHANRATE / dnatauf * 6.6);
 
   return ch;
}

/*--------------------------------------------*/

chantype *makna1(void)

/* Sequential-state description of 
   Na+ channel that very closely mimics 
   kinetics of Hodgkin-Huxley m3h description.

   See Armstrong and Matteson, Curr. Topics Membranes Transport 22:331-352.

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans =   State to go to on this transition. 
    Set trate =   Function that returns basic rate for transition. 
    Set ratemul = Multiplier for rate function.
    Set rateo  =  Sets additional tau multiplier: 
                        0->am, 1->bm, 2->ah, 3->bh.

    dqm = Q10 for m = 2
    dqh = Q10 for h = 3
    dqc = Q10 for conductance = 1.4

    For info on Q10 of m and h and conductance, see:

    Frankenhaeuser, B., and Moore, L.E. (1963) J. Physiol 169: 431-467.

*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 8;				/* 8 Markov states */
   nparm = 2;				/* make 2 sets of params, m, h */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(NA,1,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dnau;			/* conductance measured at dbasetc. */
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = dpcana;		/* permeability to Ca++ ions/perm Na */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcna1m; 	/* default rate function */
   parm[0].funcname = (char *)"calcna1m"; /* user rate function */
   parm[0].dq[0]    = dqm; 		/* Q10 for m rate function */
   parm[0].voff = dnaoffsm; 		/* voltage offset from user */

   parm[1].nfval    = 2; 		/* number of func vals (ah, bh, etc.) */
   parm[1].nival    = 0; 		/* number of implicit vals (h1, h2) */
   parm[1].chancalc =  calcna1h;	/* default rate function */ 
   parm[1].funcname = (char *)"calcna1h"; 	/* user rate function */
   parm[1].dq[0]    = dqh; 		/* Q10 for m rate function */
   parm[1].voff = dnaoffsh; 		/* voltage offset from user */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alm;
   spnt[0].ratemul[0] = 3.0 * m;
   spnt[0].rateo  [0] = 0;	/* taum */
   spnt[0].trans  [1] = 4;
   spnt[0].trate  [1] = beth;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 3;	/* tauh */

                                        /*   0 <-> 1 <-> 2 <-> 3   */

                                        /*   |     |     |     |   */

                          		/*   |     |     |     |   */

                          		/*   4 <-> 5 <-> 6 <-> 7   */
   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = alm;
   spnt[1].ratemul[0] = 2.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = betm;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 5;
   spnt[1].trate  [2] = beth;
   spnt[1].ratemul[2] = 1.0 * h;
   spnt[1].rateo  [2] = 3;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = alm;
   spnt[2].ratemul[0] = 1.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betm;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = beth;
   spnt[2].ratemul[2] = 1.0 * h;
   spnt[2].rateo  [2] = 3;	/* tauh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 7;			/* state 3 = the open state */
   spnt[3].trate  [0] = beth;
   spnt[3].ratemul[0] = 1.0 * h;
   spnt[3].rateo  [0] = 3;	/* tauh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betm;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = alm;
   spnt[4].ratemul[0] = 3.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 0;
   spnt[4].trate  [1] = alh;
   spnt[4].ratemul[1] = 1.0 * h;
   spnt[4].rateo  [1] = 2;
                                        /*   0 <-> 1 <-> 2 <-> 3   */

                                        /*   |     |     |     |   */

                          		/*   |     |     |     |   */

                          		/*   4 <-> 5 <-> 6 <-> 7   */
   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = alm;
   spnt[5].ratemul[0] = 2.0 * m;
   spnt[5].rateo  [0] = 0;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = betm;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;
   spnt[5].trans  [2] = 1;
   spnt[5].trate  [2] = alh;
   spnt[5].ratemul[2] = 1.0;
   spnt[5].rateo  [2] = 2;

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 5;
   spnt[6].trate  [0] = betm;
   spnt[6].ratemul[0] = 2.0 * m;	
   spnt[6].rateo  [0] = 1;
   spnt[6].trans  [1] = 7;
   spnt[6].trate  [1] = alm;
   spnt[6].ratemul[1] = 1.0 * m;
   spnt[6].rateo  [1] = 0;
   spnt[6].trans  [2] = 2;
   spnt[6].trate  [2] = alh;
   spnt[6].ratemul[2] = 1.0 * h;
   spnt[6].rateo  [2] = 2;

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 3;			/* state 5 = inactivated state*/
   spnt[7].trate  [0] = alh;
   spnt[7].ratemul[0] = 1.0 * h;
   spnt[7].rateo  [0] = 2;	/* tauh */
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = betm;
   spnt[7].ratemul[1] = 3.0 * m;
   spnt[7].rateo  [1] = 1;	/* taum */

   return ch;
}
/*----------------------------------------*/

/* makna2() is in module "channa2.cc" */

/*----------------------------------------*/

/* makna3() is in module "channa3.cc" */

/*----------------------------------------*/

/* makna4() is in module "channa4.cc" */

/*----------------------------------------*/

/* makna6() is in module "channa6.cc" */

/*----------------------------------------*/

/* makif() is in module "chanif.cc", type 21 */


chantype *makna20(void)

/* 6-state model of Na+ channel that has similar 
   kinetics to Hodgkin-Huxley m3h description.

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans = state to go to on this transition. 
    Set trate = function that returns basic rate for transition. 
    Set ratemul = multiplier for rate function.
    Set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 6;				/* 6 Markov states */
   nparm = 2;				/* set up m, h */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(NA,20,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dnau;			/* conductance measured at dbasetc. */
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = dpcana;		/* permeability to Ca++ ions/perm Na */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcna1m; 	/* default rate function */
   parm[0].funcname = (char *)"calcna6m"; 	/* user function */
   parm[0].dq[0]    = dqm; 		/* Q10 for m rate function */

   parm[1].nfval    = 2; 		/* number of func vals (ah, bh, etc.) */
   parm[1].nival    = 0; 		/* number of implicit vals (h1, h2) */
   parm[1].chancalc =  calcna1h;	/* default rate function */ 
   parm[1].funcname = (char *)"calcna6h"; 	/* user function */
   parm[0].dq[0]    = dqh; 		/* Q10 for m rate function */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alm;
   spnt[0].ratemul[0] = 3.0 * m;
   spnt[0].rateo  [0] = 0;	/* taum */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  		/*   0 <-> 1 <-> 2 <-> 3   */
   spnt[1].trate  [0] = alm;
   spnt[1].ratemul[0] = 2.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;			/*               |     |   */
   spnt[1].trate  [1] = betm;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;			/*               4 <-> 5   */
   spnt[2].trate  [0] = alm;
   spnt[2].ratemul[0] = 1.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betm;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 4;
   spnt[2].trate  [2] = beth;
   spnt[2].ratemul[2] = 2.0 * h;
   spnt[2].rateo  [2] = 3;	/* taum */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 5;			/* state 3 = the open state */
   spnt[3].trate  [0] = beth;
   spnt[3].ratemul[0] = 2.0 * h;
   spnt[3].rateo  [0] = 3;	/* tauh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betm;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = alm;
   spnt[4].ratemul[0] = .010 * m; /*.001;  /* no return to deactivated state */
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 2;
   spnt[4].trate  [1] = alh;
   spnt[4].ratemul[1] = 1.0 * h;
   spnt[4].rateo  [1] = 2;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 3;			/* state 5 = inactivated state*/
   spnt[5].trate  [0] = alh;
   spnt[5].ratemul[0] = 1.0 * h;
   spnt[5].rateo  [0] = 2;	/* tauh */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = betm;
   spnt[5].ratemul[1] = 3.0 * m;
   spnt[5].rateo  [1] = 1;	/* taum */

   return ch;
}
/*----------------------------------------*/

