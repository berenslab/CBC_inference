/* segment chanca1 in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
				int nq, double bt);
chantype *getchantype(int ctype, int cnum);
void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
                                double (*rrate)(chan *cpnt), double rate);
double rchanf(chan *cpnt);
double rchanr(chan *cpnt);

double calcna2m(double v, int func);
double calcna2h(double v, int func);


double al0 (chan *cpnt);
double bet0(chan *cpnt);
double al1 (chan *cpnt);
double bet1(chan *cpnt);

static double (*alc) (chan*) = al0;
static double (*betc)(chan*) = bet0;
static double (*alch) (chan*) = al1;
static double (*betch)(chan*) = bet1;


#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

double ncabs(double x);

void chanrate(chan *cpnt, chanparm *chp);
double rt(chan *spnt);

/*----------------------------------------*/

double calcca0m(double v, int func)

/* Calculate alphac, betac per second as a function of mV.
   See "catype.n" for calibration script.
*/

{
   double val,x,y;

  switch (func) {
  case 1:					/* alpha c */
    y = -0.04 * (v+15.);
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* singularity when x==0 */
       val = 2.0 * y / x;
    else
       val = 2.0;
    break;

  case 2:					/* beta c */
    val = (5 * exp ((v+38) / -18.)); 
    break;
  }
  return val*MSSEC*HHRATE;
}

/*----------------------------------------*/

double calcca2m(double v, int func)

/* Calculate alphac, betac per second as a function of mV.

   Ca-T type channel from Kawai et al (1996) J. Gen. Physiol. 108:525-535.
   Corrected for voffset due to their use of 10 mM Ca++.
   See "catype.n" for calibration script.

*/

{
   double val,x,y;

  switch (func) {
  case 1:					/* alpha c */
    y = -0.1 * (v+49.);
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* check for singularity */
       val = y / x;
    else
       val = 1.0;
    val *= 40.0;
    break;

  case 2:					/* beta c */
    val = (60.0 * exp ((v+69) / -25.)); 
    break;
  }
  return val * HHRATE;
}

/*----------------------------------------*/

double calcca2h(double v, int func)

/* Calculate alphach, betach as a function of mv.

   Ca-T type channel from Kawai et al (1996) J. Gen. Physiol.  108:525-535.
   Corrected for voffset due to their use of 10 mM Ca++.
   See "catype.n" for calibration script.
*/

{
   double val,x,y;

  switch (func) {
  case 1:					/* alpha c */
    val = (12 * exp ((v+73.) / -9.)); 
    break;

  case 2:					/* beta c */
    y = -0.05 * (v+70.);
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* check for singularity */
       val = y / x;
    else
       val = 1.0;
    val *= 24.; 
    val *= 1- 1 / (1 + ( exp ((v+50.) / -25.)));

    break;
  }
  return val * HHRATE;
}

/*--------------------------------------------*/

chantype *makca0(void)

/* Make L-type Ca channel with kinetics of c3 like HH.  */

/* Unitary current 20 pS @ 35 deg C from:

     Karschin and Lipton, (1989) J. Physiol.  418: 379-396.
*/
{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 2;				/* 2 Markov states */
   nparm = 1;				/* make 1 set of params, c */
   nq = 1;				/* number of Q10 values */
   ch = makchantype(CA, 0, nstate, nparm, nq, dbasetc);/* make chan const,parm */
   ch->hh = 1;				/* HH chan = uses "m" */
   ch->unitary = dcalu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* permeability to K ions / perm Ca */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 2; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca0m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca0m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */
   parm[0].voff = dcaoffs; 		/* voltage offset from user */

   mak2state(spnt,rchanf,rchanr,CHANRATE);

   return ch;
}

/*--------------------------------------------*/

chantype *makca1(void)

/* Sequential-state description of L-type Ca++ channel, c3 like m3h.

   See Armstrong and Matteson, Curr. Topics Membranes Transport 22:331-352.

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans =   State to go to on this transition. 
    Set trate =   Function that returns basic rate for transition. 
    Set ratemul = Multiplier for rate function.
    Set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.


                  0 <-> 1 <-> 2 <-> 3 
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 4;				/* 4 Markov states */
   nparm = 1;				/* make 1 set of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch = makchantype(CA, 1, nstate, nparm, nq, dbasetc);/* make chan const,parm */
   ch->unitary = dcalu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* permeability to K ions / perm Ca */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca0m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca1m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alc;
   spnt[0].ratemul[0] = 3.0 * m;
   spnt[0].rateo  [0] = 0;	/* taum */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = alc;
   spnt[1].ratemul[0] = 2.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = betc;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = alc;
   spnt[2].ratemul[0] = 1.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betc;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = betc;
   spnt[3].ratemul[0] = 3.0 * m;
   spnt[3].rateo  [0] = 1;	/* taum */

   return ch;
}

/*--------------------------------------------*/

chantype *makca2(void)

/* Make T-type channel with kinetics of c2h like HH.  */

/* Unitary current 8 pS @ 35 deg C from:

     Karschin and Lipton, (1989) J. Physiol.  418: 379-396.

     qfactor = exp (log(Q10=1.4) * (35-6.3)/10) = 2.626

     conductance @ 6.3 deg C = 8 pS / 2.626 = 3 pS @ 6.3 deg
*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 2;				/* 2 Markov states */
   nparm = 2;				/* make 2 set of params, c */
   nq = 1;				/* number of Q10 values */
   ch = makchantype(CA, 2, nstate, nparm, nq, dbasetc);/* make chan const,parm */
   ch->hh = 1;				/* HH chan = uses "m", "h" */
   ch->unitary = dcatu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* permeability to K ions / perm Ca */
   ch->cabnd  = 1;			/* needs Ca binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 2; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca2m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca2m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   parm[1].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[1].nival    = 2; 		/* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calcca2h; 	/* default rate function */
   parm[1].funcname = (char *)"calcca2h"; 	/* user rate function */
   parm[1].dq[0] = dqca; 		/* Q10 for c rate function */

   mak2state(spnt,rchanf,rchanr,CHANRATE);

   return ch;
}

/*--------------------------------------------*/

chantype *makca3(void)

/* Sequential-state description of T-type Ca++ channel, c2 like m2h.

   See Armstrong and Matteson, Curr. Topics Membranes Transport 22:331-352.

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans =   State to go to on this transition. 
    Set trate =   Function that returns basic rate for transition. 
    Set ratemul = Multiplier for rate function.
    Set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.


                  0 <-> 1 <-> 2 

                  |     |     |

                  3 <-> 4 <-> 5 
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 6;				/* 6 Markov states */
   nparm = 2;				/* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch = makchantype(CA, 3, nstate, nparm, nq, dbasetc);/* make chan const,parm */
   ch->unitary = dcatu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* rel permeability to K ions */
   ch->cabnd  = 1;			/* needs Ca binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca2m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca3m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   parm[1].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calcca2h; 	/* default rate function */
   parm[1].funcname = (char *)"calcca3h"; 	/* user rate function */
   parm[1].dq[0] = dqca; 		/* Q10 for c rate function */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alc;
   spnt[0].ratemul[0] = 2.0 * m;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 3;
   spnt[0].trate  [1] = betch;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 3;	/* tauh */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = alc;
   spnt[1].ratemul[0] = 1.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = betc;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 4;	
   spnt[1].trate  [2] = betch;
   spnt[1].ratemul[2] = 1.0 * h;
   spnt[1].rateo  [2] = 3;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 1.0;
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = betc;
   spnt[2].ratemul[0] = 2.0 * m;
   spnt[2].rateo  [0] = 1;	/* taum */
   spnt[2].trans  [1] = 5;
   spnt[2].trate  [1] = betch;
   spnt[2].ratemul[1] = 1.0 * h;
   spnt[2].rateo  [1] = 3;	/* tauh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = alc;
   spnt[3].ratemul[0] = 2.0 * m;
   spnt[3].rateo  [0] = 0;	/* taum */
   spnt[3].trans  [1] = 0;
   spnt[3].trate  [1] = alch;
   spnt[3].ratemul[1] = 1.0 * h;
   spnt[3].rateo  [1] = 2;	/* tauh */

   spnt[4].numtrans   = 3;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;	  	
   spnt[4].trate  [0] = alc;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;	
   spnt[4].trate  [1] = betc;
   spnt[4].ratemul[1] = 1.0 * m;
   spnt[4].rateo  [1] = 1;
   spnt[4].trans  [2] = 1;	
   spnt[4].trate  [2] = alch;
   spnt[4].ratemul[2] = 1.0 * h;
   spnt[4].rateo  [2] = 2;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = betc;
   spnt[5].ratemul[0] = 2.0 * m;
   spnt[5].rateo  [0] = 1;	/* taum */
   spnt[5].trans  [1] = 2;
   spnt[5].trate  [1] = alch;
   spnt[5].ratemul[1] = 1.0 * h;
   spnt[5].rateo  [1] = 2;	/* tauh */


   return ch;
}

/*--------------------------------------------*/

chantype *makca4(void)

/* Make T-type channel with kinetics of c3h like HH.  */

/* Unitary current 8 pS @ 35 deg C from:

     Karschin and Lipton, (1989) J. Physiol.  418: 379-396.

     qfactor = exp (log(Q10=1.4) * (35-6.3)/10) = 2.626

     conductance @ 6.3 deg C = 8 pS / 2.626 = 3 pS @ 6.3 deg
*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 2;				/* 2 Markov states */
   nparm = 2;				/* make 2 set of params, c */
   nq = 1;				/* number of Q10 values */
   ch = makchantype(CA, 4, nstate, nparm, nq, dbasetc);/* make chan const,parm */
   ch->hh = 1;				/* HH chan = uses "m", "h" */
   ch->unitary = dcatu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* permeability to K ions / perm Ca */
   ch->cabnd  = 1;			/* needs Ca binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 2; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca2m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca4m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   parm[1].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[1].nival    = 2; 		/* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calcca2h; 	/* default rate function */
   parm[1].funcname = (char *)"calcca4h"; 	/* user rate function */
   parm[1].dq[0] = dqca; 		/* Q10 for c rate function */

   mak2state(spnt,rchanf,rchanr,CHANRATE);

   return ch;
}

/*--------------------------------------------*/

chantype *makca5(void)

/* Sequential-state description of T-type Ca++ channel, c3 like m3h.

   See Armstrong and Matteson, Curr. Topics Membranes Transport 22:331-352.

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans =   State to go to on this transition. 
    Set trate =   Function that returns basic rate for transition. 
    Set ratemul = Multiplier for rate function.
    Set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.


                  0 <-> 1 <-> 2 <-> 3

                  |     |     |     |

                  4 <-> 5 <-> 6 <-> 7 
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 8;				/* 8 Markov states */
   nparm = 2;				/* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch = makchantype(CA, 5, nstate, nparm, nq, dbasetc);/* make chan const,parm*/
   ch->unitary = dcatu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* rel permeability to K ions */
   ch->cabnd  = 1;			/* needs Ca binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca2m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca5m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   parm[1].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calcca2h; 	/* default rate function */
   parm[1].funcname = (char *)"calcca5h"; 	/* user rate function */
   parm[1].dq[0] = dqca; 		/* Q10 for c rate function */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alc;
   spnt[0].ratemul[0] = 3.0 * m;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 4;
   spnt[0].trate  [1] = betch;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 3;	/* tauh */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = alc;
   spnt[1].ratemul[0] = 2.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = betc;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 5;	
   spnt[1].trate  [2] = betch;
   spnt[1].ratemul[2] = 1.0 * h;
   spnt[1].rateo  [2] = 3;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = alc;
   spnt[2].ratemul[0] = 1.0 * m;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betc;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = betch;
   spnt[2].ratemul[2] = 1.0 * h;
   spnt[2].rateo  [2] = 3;	/* tauh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = betc;
   spnt[3].ratemul[0] = 3.0 * m;
   spnt[3].rateo  [0] = 1;	/* taum */
   spnt[3].trans  [1] = 7;
   spnt[3].trate  [1] = betch;
   spnt[3].ratemul[1] = 1.0 * h;
   spnt[3].rateo  [1] = 3;	/* tauh */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = alc;
   spnt[4].ratemul[0] = 3.0 * m;
   spnt[4].rateo  [0] = 0;	/* taum */
   spnt[4].trans  [1] = 0;
   spnt[4].trate  [1] = alch;
   spnt[4].ratemul[1] = 1.0 * h;
   spnt[4].rateo  [1] = 2;	/* tauh */

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;	  	
   spnt[5].trate  [0] = alc;
   spnt[5].ratemul[0] = 2.0 * m;
   spnt[5].rateo  [0] = 0;
   spnt[5].trans  [1] = 4;	
   spnt[5].trate  [1] = betc;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;
   spnt[5].trans  [2] = 1;	
   spnt[5].trate  [2] = alch;
   spnt[5].ratemul[2] = 1.0 * h;
   spnt[5].rateo  [2] = 2;

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;	  	
   spnt[6].trate  [0] = alc;
   spnt[6].ratemul[0] = 1.0 * m;
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;	
   spnt[6].trate  [1] = betc;
   spnt[6].ratemul[1] = 2.0 * m;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 2;	
   spnt[6].trate  [2] = alch;
   spnt[6].ratemul[2] = 1.0 * h;
   spnt[6].rateo  [2] = 2;

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 6;
   spnt[7].trate  [0] = betc;
   spnt[7].ratemul[0] = 3.0 * m;
   spnt[7].rateo  [0] = 1;	/* taum */
   spnt[7].trans  [1] = 3;
   spnt[7].trate  [1] = alch;
   spnt[7].ratemul[1] = 1.0 * h;
   spnt[7].rateo  [1] = 2;	/* tauh */


   return ch;
}

/*--------------------------------------------*/

/* Rates for Ca type 6 */

/*             sec-1                         */

#define kvo   2500
#define k_vo     8
#define Vkv     25
#define Vk_v   -18
#define Vk_o   -34

#define ko    4000
#define k_o     25 
#define ki      70
#define k_i      1.4

double calcca6m(double v, int func)

/* Calculate kinetic rates per second per mV.
   From Serrano et al (1999), J Gen Physiol 114:185-201

*/

{
   double val,x,y;

  if (func > 4) func = 4;
  switch (func) {
  case 1: val = kvo  *  exp(v/Vkv);     break;
  case 2: val = k_vo *  exp(v/Vk_v);    break;
  case 3: val = ko;                     break;
  case 4: val = k_o  *  exp(v/Vk_o);    break;
  }
  return val;
}

/*--------------------------------------------*/

double calcca6h(double v, int func)

/* Calculate kinetic rates per second per mV.
*/

{
   double val,x,y;

  switch (func) {
  case 1: val = ki;        break;
  case 2: val = k_i;       break;
  }
  return val;
}

/*--------------------------------------------*/

double ca6kv(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double ca6k_v(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double ca6ko(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double ca6k_o(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double ca6ki(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double ca6k_i(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}


/*--------------------------------------------*/

chantype *makca6(void)

/* Markov sequential-state description of T-type Ca++ channel.
   Note that inactivation is not voltage-dependent.
   From Serrano JR, Perez-Reyes E, and Jones SW (1999), J Gen Physiol 114:185-201

            4kv          3kv           2kv         kv            ko
         0 <->        1  <->       2  <->      3   <->      4    <->   5     Open
            k_v          2k_v          3k_v        4k_v          k_o

     /\           /\           /\           /\           /\         /\
     k_i | kif3   k_i | kif2   k_i | kif   k_i | ki     k_i | ki   k_i | ki
     ---   \/     ---   \/     ---   \/          \/           \/         \/      
      h3           h2           h

             4kv/f         3kv/f        2kv/f       kv          ko
          6  <->       7  <->      8   <->     9   <->      10  <->    11    Inact
             k_vh         2k_vh         3k_vh       4k_v        k_o 
*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double dbcat;
     double m,h;
     int nstate, nparm, nq;

#define f      0.2  
#define h1     0.5

#define f2     (f*f)
#define f3     (f*f*f)

#define h2     (h1*h1)
#define h3     (h1*h1*h1)

   nstate = 12;				/* 12 Markov states */
   nparm = 2;				/* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   dbcat = 20.0;			/* temperature for orig measurements */
   ch = makchantype(CA, 6, nstate, nparm, nq, dbcat);/* make chan const,parm*/
   ch->unitary = dcatu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* rel permeability to K ions */
   ch->cabnd  = 1;			/* needs Ca binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 4; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca6m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca6m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   parm[1].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calcca6h; 	/* default rate function */
   parm[1].funcname = (char *)"calcca6h"; 	/* user rate function */
   parm[1].dq[0] = dqca; 		/* Q10 for c rate function */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = ca6kv;
   spnt[0].ratemul[0] = 4.0 * m;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 6;
   spnt[0].trate  [1] = ca6ki;
   spnt[0].ratemul[1] = 1.0 * f3 * h;
   spnt[0].rateo  [1] = 3;	/* tauh */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = ca6kv;
   spnt[1].ratemul[0] = 3.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = ca6k_v;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 7;	
   spnt[1].trate  [2] = ca6ki;
   spnt[1].ratemul[2] = 1.0 * f2 * h;
   spnt[1].rateo  [2] = 3;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = ca6kv;
   spnt[2].ratemul[0] = 2.0 * m;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = ca6k_v;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 8;
   spnt[2].trate  [2] = ca6ki;
   spnt[2].ratemul[2] = 1.0 * f * h;
   spnt[2].rateo  [2] = 3;	/* tauh */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;	  	
   spnt[3].trate  [0] = ca6kv;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = ca6k_v;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 9;
   spnt[3].trate  [2] = ca6ki;
   spnt[3].ratemul[2] = 1.0 * h;
   spnt[3].rateo  [2] = 3;	/* tauh */

   spnt[4].numtrans   = 3;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;	  	
   spnt[4].trate  [0] = ca6ko;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = ca6k_v;
   spnt[4].ratemul[1] = 4.0 * m;
   spnt[4].rateo  [1] = 1;	/* taum */
   spnt[4].trans  [2] = 10;
   spnt[4].trate  [2] = ca6ki;
   spnt[4].ratemul[2] = 1.0 * h;
   spnt[4].rateo  [2] = 3;	/* tauh */

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = ca6k_o;
   spnt[5].ratemul[0] = 1.0 * m;
   spnt[5].rateo  [0] = 1;	/* taum */
   spnt[5].trans  [1] = 11;
   spnt[5].trate  [1] = ca6ki;
   spnt[5].ratemul[1] = 1.0 * h;
   spnt[5].rateo  [1] = 3;	/* tauh */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = ca6kv;
   spnt[6].ratemul[0] = 4.0 / f * m;
   spnt[6].rateo  [0] = 0;	/* taum */
   spnt[6].trans  [1] = 0;
   spnt[6].trate  [1] = ca6k_i;
   spnt[6].ratemul[1] = 1.0 / h3 * h;
   spnt[6].rateo  [1] = 2;	/* tauh */

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;	  	
   spnt[7].trate  [0] = ca6kv;
   spnt[7].ratemul[0] = 3.0 / f * m;
   spnt[7].rateo  [0] = 0;
   spnt[7].trans  [1] = 6;	
   spnt[7].trate  [1] = ca6k_v;
   spnt[7].ratemul[1] = 1.0 * h1 * m;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 1;	
   spnt[7].trate  [2] = ca6k_i;
   spnt[7].ratemul[2] = 1.0 / h2 * h;
   spnt[7].rateo  [2] = 2;

   spnt[8].numtrans   = 3;
   spnt[8].cond       = 0;
   spnt[8].trans  [0] = 9;	  	
   spnt[8].trate  [0] = ca6kv;
   spnt[8].ratemul[0] = 2.0 / f * m;
   spnt[8].rateo  [0] = 0;
   spnt[8].trans  [1] = 7;	
   spnt[8].trate  [1] = ca6k_v;
   spnt[8].ratemul[1] = 2.0 * h1 * m;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 2;	
   spnt[8].trate  [2] = ca6k_i;
   spnt[8].ratemul[2] = 1.0 / h1 * h;
   spnt[8].rateo  [2] = 2;

   spnt[9].numtrans   = 3;
   spnt[9].cond       = 0;
   spnt[9].trans  [0] = 10;	  	
   spnt[9].trate  [0] = ca6kv;
   spnt[9].ratemul[0] = 1.0 * m;
   spnt[9].rateo  [0] = 0;
   spnt[9].trans  [1] = 8;	
   spnt[9].trate  [1] = ca6k_v;
   spnt[9].ratemul[1] = 3.0 * h1 * m;
   spnt[9].rateo  [1] = 1;
   spnt[9].trans  [2] = 3;	
   spnt[9].trate  [2] = ca6k_i;
   spnt[9].ratemul[2] = 1.0 * h;
   spnt[9].rateo  [2] = 2;

   spnt[10].numtrans   = 3;
   spnt[10].cond       = 0;
   spnt[10].trans  [0] = 11;	  	
   spnt[10].trate  [0] = ca6ko;
   spnt[10].ratemul[0] = 1.0 * m;
   spnt[10].rateo  [0] = 0;
   spnt[10].trans  [1] = 9;	
   spnt[10].trate  [1] = ca6k_v;
   spnt[10].ratemul[1] = 4.0 * m;
   spnt[10].rateo  [1] = 1;
   spnt[10].trans  [2] = 4;	
   spnt[10].trate  [2] = ca6k_i;
   spnt[10].ratemul[2] = 1.0 * h;
   spnt[10].rateo  [2] = 2;

   spnt[11].numtrans   = 2;
   spnt[11].cond       = 0;
   spnt[11].trans  [0] = 10;
   spnt[11].trate  [0] = ca6k_o;
   spnt[11].ratemul[0] = 1.0 * m;
   spnt[11].rateo  [0] = 1;	/* taum */
   spnt[11].trans  [1] = 5;
   spnt[11].trate  [1] = ca6k_i;
   spnt[11].ratemul[1] = 1.0 * h;
   spnt[11].rateo  [1] = 2;	/* tauh */


   return ch;
}

#undef kvo
#undef k_vo
#undef Vkv
#undef Vk_v

#undef ko
#undef k_o
#undef ki
#undef k_i

#undef f2
#undef f3

#undef h2
#undef h3

#undef f
#undef h1

/*--------------------------------------------*/

/*             msec-1                         */

#define kvo      3.7
#define k_vo    0.02
#define Vkv     25.5
#define Vk_v   -15.3

#define ko        20
#define k_o     1.65 
#define ki     0.062
#define k_i     0.12e-3

double calcca7m(double v, int func)

/* Calculate kinetic rates per ms per mV.

   See "catype.n" for calibration script.
*/

{
   double val,x,y;

  if (func > 4) func = 4;
  switch (func) {
  case 1: val = kvo  *  exp(v/Vkv);     break;
  case 2: val = k_vo *  exp(v/Vk_v);    break;
  case 3: val = ko;                     break;
  case 4: val = k_o;                    break;
  }
  return val * MSSEC;
}

/*--------------------------------------------*/

double calcca7h(double v, int func)

/* Calculate alphach, betach per ms a function of mV.
   See "catype.n" for calibration script.
*/

{
   double val,x,y;

  switch (func) {
  case 1: val = ki;        break;
  case 2: val = k_i;       break;
  }
  return val * MSSEC;
}

/*--------------------------------------------*/

chantype *makca7(void)

/* Markov sequential-state description of T-type Ca++ channel
   from Lee SC, Hayashida Y, and Ishida AT (2003) J Neurophysiol 90:3888-3901.
   Note that inactivation is not voltage-dependent.

   Like Serrano JR, Perez-Reyes E, and Jones SW (1999), J Gen Physiol 114:185-201
   expt k_o is not voltage dependent.

            4kv          3kv           2kv         kv            ko
         0 <->        1  <->       2  <->      3   <->      4    <->   5     Open
            k_v          2k_v          3k_v        4k_v          k_o

     /\           /\           /\           /\           /\         /\
     k_i | kif3   k_i | kif2   k_i | kif   k_i | ki     k_i | ki   k_i | ki
     ---   \/     ---   \/     ---   \/          \/           \/         \/      
      h3           h2           h

             4kv/f         3kv/f        2kv/f       kv          ko
          6  <->       7  <->      8   <->     9   <->      10  <->    11    Inact
             k_vh         2k_vh         3k_vh       4k_v        k_o 
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double dbcat;
     double m,h;
     int nstate, nparm, nq;

/* sec-1    */

#define f      0.305  
#define h1     0.226

#define f2     (f*f)
#define f3     (f*f*f)

#define h2     (h1*h1)
#define h3     (h1*h1*h1)

   nstate = 12;				/* 12 Markov states */
   nparm = 2;				/* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   dbcat = 23.0;			/* temperature for orig measurements */
   ch = makchantype(CA, 7, nstate, nparm, nq, dbcat);/* make chan const,parm*/
   ch->unitary = dcatu;
   ch->vrev = R2F * ktemp * log(dcao/(dcai+dpki));   /* GHK equation */
   ch->ions->ionp[PCA] = 1.0;			/* rel permeability to Ca ions */
   ch->ions->ionp[PNA] = dpnaca;		/* rel permeability to Na ions */
   ch->ions->ionp[PK] = dpkca;		/* rel permeability to K ions */
   ch->cabnd  = 1;			/* needs Ca binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 4; 		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calcca7m; 	/* default rate function */
   parm[0].funcname = (char *)"calcca7m"; 	/* user rate function */
   parm[0].dq[0] = dqca; 		/* Q10 for c rate function */

   parm[1].nfval    = 2; 		/* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0; 		/* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calcca7h; 	/* default rate function */
   parm[1].funcname = (char *)"calcca7h"; 	/* user rate function */
   parm[1].dq[0] = dqca; 		/* Q10 for c rate function */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = ca6kv;
   spnt[0].ratemul[0] = 4.0 * m;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 6;
   spnt[0].trate  [1] = ca6ki;
   spnt[0].ratemul[1] = 1.0 * f3 * h;
   spnt[0].rateo  [1] = 3;	/* tauh */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = ca6kv;
   spnt[1].ratemul[0] = 3.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = ca6k_v;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 7;	
   spnt[1].trate  [2] = ca6ki;
   spnt[1].ratemul[2] = 1.0 * f2 * h;
   spnt[1].rateo  [2] = 3;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = ca6kv;
   spnt[2].ratemul[0] = 2.0 * m;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = ca6k_v;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 8;
   spnt[2].trate  [2] = ca6ki;
   spnt[2].ratemul[2] = 1.0 * f * h;
   spnt[2].rateo  [2] = 3;	/* tauh */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;	  	
   spnt[3].trate  [0] = ca6kv;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = ca6k_v;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 9;
   spnt[3].trate  [2] = ca6ki;
   spnt[3].ratemul[2] = 1.0 * h;
   spnt[3].rateo  [2] = 3;	/* tauh */

   spnt[4].numtrans   = 3;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;	  	
   spnt[4].trate  [0] = ca6ko;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = ca6k_v;
   spnt[4].ratemul[1] = 4.0 * m;
   spnt[4].rateo  [1] = 1;	/* taum */
   spnt[4].trans  [2] = 10;
   spnt[4].trate  [2] = ca6ki;
   spnt[4].ratemul[2] = 1.0 * h;
   spnt[4].rateo  [2] = 3;	/* tauh */

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = ca6k_o;
   spnt[5].ratemul[0] = 1.0 * m;
   spnt[5].rateo  [0] = 1;	/* taum */
   spnt[5].trans  [1] = 11;
   spnt[5].trate  [1] = ca6ki;
   spnt[5].ratemul[1] = 1.0 * h;
   spnt[5].rateo  [1] = 3;	/* tauh */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = ca6kv;
   spnt[6].ratemul[0] = 4.0 / f * m;
   spnt[6].rateo  [0] = 0;	/* taum */
   spnt[6].trans  [1] = 0;
   spnt[6].trate  [1] = ca6k_i;
   spnt[6].ratemul[1] = 1.0 / h3 * h;
   spnt[6].rateo  [1] = 2;	/* tauh */

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;	  	
   spnt[7].trate  [0] = ca6kv;
   spnt[7].ratemul[0] = 3.0 / f * m;
   spnt[7].rateo  [0] = 0;
   spnt[7].trans  [1] = 6;	
   spnt[7].trate  [1] = ca6k_v;
   spnt[7].ratemul[1] = 1.0 * h1 * m;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 1;	
   spnt[7].trate  [2] = ca6k_i;
   spnt[7].ratemul[2] = 1.0 / h2 * h;
   spnt[7].rateo  [2] = 2;

   spnt[8].numtrans   = 3;
   spnt[8].cond       = 0;
   spnt[8].trans  [0] = 9;	  	
   spnt[8].trate  [0] = ca6kv;
   spnt[8].ratemul[0] = 2.0 / f * m;
   spnt[8].rateo  [0] = 0;
   spnt[8].trans  [1] = 7;	
   spnt[8].trate  [1] = ca6k_v;
   spnt[8].ratemul[1] = 2.0 * h1 * m;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 2;	
   spnt[8].trate  [2] = ca6k_i;
   spnt[8].ratemul[2] = 1.0 / h1 * h;
   spnt[8].rateo  [2] = 2;

   spnt[9].numtrans   = 3;
   spnt[9].cond       = 0;
   spnt[9].trans  [0] = 10;	  	
   spnt[9].trate  [0] = ca6kv;
   spnt[9].ratemul[0] = 1.0 * m;
   spnt[9].rateo  [0] = 0;
   spnt[9].trans  [1] = 8;	
   spnt[9].trate  [1] = ca6k_v;
   spnt[9].ratemul[1] = 3.0 * h1 * m;
   spnt[9].rateo  [1] = 1;
   spnt[9].trans  [2] = 3;	
   spnt[9].trate  [2] = ca6k_i;
   spnt[9].ratemul[2] = 1.0 * h;
   spnt[9].rateo  [2] = 2;

   spnt[10].numtrans   = 3;
   spnt[10].cond       = 0;
   spnt[10].trans  [0] = 11;	  	
   spnt[10].trate  [0] = ca6ko;
   spnt[10].ratemul[0] = 1.0 * m;
   spnt[10].rateo  [0] = 0;
   spnt[10].trans  [1] = 9;	
   spnt[10].trate  [1] = ca6k_v;
   spnt[10].ratemul[1] = 4.0 * m;
   spnt[10].rateo  [1] = 1;
   spnt[10].trans  [2] = 4;	
   spnt[10].trate  [2] = ca6k_i;
   spnt[10].ratemul[2] = 1.0 * h;
   spnt[10].rateo  [2] = 2;

   spnt[11].numtrans   = 2;
   spnt[11].cond       = 0;
   spnt[11].trans  [0] = 10;
   spnt[11].trate  [0] = ca6k_o;
   spnt[11].ratemul[0] = 1.0 * m;
   spnt[11].rateo  [0] = 1;	/* taum */
   spnt[11].trans  [1] = 5;
   spnt[11].trate  [1] = ca6k_i;
   spnt[11].ratemul[1] = 1.0 * h;
   spnt[11].rateo  [1] = 2;	/* tauh */


   return ch;
}

#undef kvo
#undef k_vo
#undef Vkv
#undef Vk_v

#undef ko
#undef k_o
#undef ki
#undef k_i

#undef f2
#undef f3

#undef h2
#undef h3

#undef f
#undef h1

/*--------------------------------------------*/

/* check:

   Ueda Y, Kaneko, A., and Kaneda, M. (1992) Voltage -dependent
ionic currents in solitary horizontal cells isolated from cat
retina.  J. Neurophysiol. 68: 1143-1150.

   Usui, S., Kamiyama, Y, Ishii, H, and Ikeno, H (1996a)
Reconstruction of retinal horizontal cell responses by the ionic
current model.  Vis. Research 36: 1711-1719.

   Usui, S., Kamiyama, Y, Ishii, H, and Ikeno, H (1996b) Ionic
current model of bipolar cells in the lower vertebrate retina.
Vis. Research 36: 4069-4076.

   Sullivan, JM, and Lasater, EM (1990) Sustained and transient
potassium currents of cultured horizontal cells isolated from
white bass retina.  J. Neurophysiol. 64: 1758-1766.

Satoh H. Aoki K. Watanabe SI. Kaneko A. L-type calcium
channels in the axon terminal of mouse bipolar cells.
Neuroreport. 9:2161-5, 1998.

de la Villa P. Vaquero CF. Kaneko A. Two types of calcium
currents of the mouse bipolar cells recorded in the retinal slice
preparation. European Journal of Neuroscience. 10(1):317-23, 1998

Kaneda M. Mochizuki M. Aoki K. Kaneko A. Modulation of GABAC
response by Ca2+ and other divalent cations in horizontal cells
of the catfish retina.  Journal of General Physiology.
110(6):741-7, 1997. 

0. Sasaki T. Kaneko A. L-Glutamate-induced responses in OFF-type
bipolar cells of the cat retina. Vision Research. 36(6):787-95,
1996 

Downing JE. Kaneko A. Cat retinal ganglion cells show transient
responses to acetylcholine and sustained responses to
L-glutamate. Neuroscience Letters.  137(1):114-8, 1992 

Kaneda M. Kaneko A. Voltage-gated calcium currents in isolated
retinal ganglion cells of the cat. Japanese Journal of
Physiology. 41(1):35-48, 1991.

Kaneko A. Tachibana M. Pinto LH. Transient calcium current of
retinal bipolar cells of the mouse. [Review] [12 refs]
Neuroscience Research - Supplement.  10:S67-76, 1989

Kaneko A. Pinto LH. Tachibana M. Transient calcium current of
retinal bipolar cells of the mouse. Journal of Physiology.
410:613-29, 1989 Mar.

Tachibana M. Kaneko A. Retinal bipolar cells receive negative
feedback input from GABAergic amacrine cells. [Review] [37 refs]
Visual Neuroscience.  1(3):297-305, 1988.  

Kaneko A. The functional role of retinal horizontal cells.
[Review] [86 refs] Japanese Journal of Physiology. 37(3):341-58,
1987.   Complete Reference 

Lohrke S. Hofmann HD. Voltage-gated currents of rabbit A- and
B-type horizontal cells in retinal monolayer cultures. Visual
Neuroscience. 11(2):369-78, 1994 

Fenwick, E.M, Marty, A., and Neher, E. 1982b Sodium and calcium
channels in bovine chromaffin cells.  J. Physiol. 331::599-635.

Nowycky M.C., Fox, A.P. and Tsien, R.W. 1985b Long opening mode
of neuronal calcium channels and its promotion by the
dihydropyridine calcium agonist Bay K 8644. PNAS 82:2178-2182.

Serrano et al (1999), J Gen Physiol 114:185-201.
   
Lee SC, Hayashida Y, and Ishida AT (2003) J Neurophysiol 90:3888-3901.

*/

