/* segment chank3 in program nc */

/* type KA inactivating K channel */

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

double al0 (chan *cpnt);
double bet0(chan *cpnt);
double al1 (chan *cpnt);
double bet1(chan *cpnt);

static double (*aln) (chan*) = al0;
static double (*betn)(chan*) = bet0;
static double (*ald) (chan*) = al1;
static double (*betd)(chan*) = bet1;

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

double calck3n (double v, int func)

/* Calculate K3 rate functions given voltage in mv.

   Original was in terms of per msec.
   Calc in terms of sec.

   The "func" parameter defines:

    1	alpha n 
    2	beta  m
*/

{
   double val,x,y;

  switch (func) {

  case 1:					/* alpha n type KA chan */
    y = -0.1 * (v+90.);                         /* the modern way */
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* singularity when x==0 */
       val = 0.0107 * y / x;     /* taken from Fohlmeister + Miller 1998 */
    else
       val = 0.0107;
    break;

  case 2:
    val =  0.0178 * exp ((v+30) / -10.);        /* type KA chan */
    break;              /* taken from Fohlmeister and Miller 1998 */
  }
  return val*MSSEC*HHRATE; 			/* normalize to dbasetc */
}						/* i.e. 22 deg C */

/* - - - - - - - - - - - - - - - - - - - -*/

double calck3d (double v, int func)

/* Calculate K3 inactivation rate functions given voltage in mv.

   Original was in terms of per msec.
   Calc in terms of sec.

   The "func" parameter defines:

    1	alpha d 
    2	beta  d
*/

{
   double val,y;

  switch (func) {
  case 1:					/* alpha d (K inact) */
    val = (0.0071 * exp (-0.05*(v+70)));
    break;

  case 2:					/* beta d */
    y = -0.1 * (v+40.);
    val = (0.107 / (exp (y) + 1.0));
    break;
  }
  return val*MSSEC*HHRATE; 			/* normalize to dbasetc */
}

/*--------------------------------------------*/

chantype *makk2(void)

/* HHdescription of KA channel. Based on kinetics of Hodgkin-Huxley 
   m3h description.
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 2;                          /* 2 Markov states for noise */
   nparm = 2;                           /* make 2 sets of params, n, d */
   nq = 2;				/* make 1 Q10 value */
   ch=makchantype(K,2,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->hh=1;				/* HH chan - uses "m", "h" */
   ch->unitary = dkau;			/* 22 pS = 30 pS @30 deg, Lipton&Tauck*/
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions / perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 2;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck3n;         /* default rate function */
   parm[0].funcname = (char *)"calck2n"; /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for n rate function */
   parm[0].dq[1] = dqnb;                /* Q10 for n rate function */
   parm[0].voff = dkoffsn;              /* voltage offset set by user */

   parm[1].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 2;                /* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calck3d;         /* default rate function */
   parm[1].funcname = (char *)"calck2d"; /* user rate function */
   parm[1].dq[0] = dqd;                 /* Q10 for n rate function */
   parm[1].dq[1] = dqd;                 /* Q10 for n rate function */
   parm[1].voff = dkoffsh;              /* voltage offset set by user */

   mak2state(spnt,rchanf,rchanr,CHANRATE);

   return ch;
}

/*--------------------------------------------*/

chantype *makk3(void)

/* Sequential-state description of 
   KA channel. Based on kinetics of Hodgkin-Huxley 
   m3h description.

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

*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n,h;
     int nstate, nparm, nq;

   nstate = 8;                          /* 8 Markov states */
   nparm = 2;                           /* make 2 sets of params, n, d */
   nq = 2;				/* make 1 Q10 value */
   ch=makchantype(K,3,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dkau;			/* 22 pS@22 = 30pS @30, Lipton&Tauck */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions / perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
  
   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck3n;         /* default rate function */
   parm[0].funcname = (char *)"calck3n";  /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for alphan rate function */
   parm[0].dq[1] = dqnb;                /* Q10 for betan rate function */
   parm[0].voff = dkoffsn;              /* voltage offset set by user */

   parm[1].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[1].chancalc =  calck3d;         /* default rate function */
   parm[1].funcname = (char *)"calck3d";  /* user rate function */
   parm[1].dq[0] = dqd;                 /* Q10 for n rate function */
   parm[1].dq[1] = dqd;                 /* Q10 for n rate function */
   parm[1].voff = dkoffsh;              /* voltage offset set by user */


   n = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = aln;
   spnt[0].ratemul[0] = 3.0 * n;
   spnt[0].rateo  [0] = 0;	/* taum */
   spnt[0].trans  [1] = 4;
   spnt[0].trate  [1] = betd;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 3;	/* tauh */

                                        /*   0 <-> 1 <-> 2 <-> 3   */

                                        /*   |     |     |     |   */

                          		/*   |     |     |     |   */

                          		/*   4 <-> 5 <-> 6 <-> 7   */
   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = aln;
   spnt[1].ratemul[0] = 2.0 * n;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = betn;
   spnt[1].ratemul[1] = 1.0 * n;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 5;
   spnt[1].trate  [2] = betd;
   spnt[1].ratemul[2] = 1.0 * h;
   spnt[1].rateo  [2] = 3;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = aln;
   spnt[2].ratemul[0] = 1.0 * n;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betn;
   spnt[2].ratemul[1] = 2.0 * n;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = betd;
   spnt[2].ratemul[2] = 1.0 * h;
   spnt[2].rateo  [2] = 3;	/* tauh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 7;			/* state 3 = the open state */
   spnt[3].trate  [0] = betd;
   spnt[3].ratemul[0] = 1.0 * h;
   spnt[3].rateo  [0] = 3;	/* tauh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betn;
   spnt[3].ratemul[1] = 3.0 * n;
   spnt[3].rateo  [1] = 1;	/* taum */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = aln;
   spnt[4].ratemul[0] = 3.0 * n;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 0;
   spnt[4].trate  [1] = ald;
   spnt[4].ratemul[1] = 1.0 * h;
   spnt[4].rateo  [1] = 2;
                                        /*   0 <-> 1 <-> 2 <-> 3   */

                                        /*   |     |     |     |   */

                          		/*   |     |     |     |   */

                          		/*   4 <-> 5 <-> 6 <-> 7   */
   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = aln;
   spnt[5].ratemul[0] = 2.0 * n;
   spnt[5].rateo  [0] = 0;	/* taum, taua */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = betn;
   spnt[5].ratemul[1] = 1.0 * n;
   spnt[5].rateo  [1] = 1;	/* taum, taub */
   spnt[5].trans  [2] = 1;
   spnt[5].trate  [2] = ald;
   spnt[5].ratemul[2] = 1.0;
   spnt[5].rateo  [2] = 2;	/* taum, tauc */

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 5;
   spnt[6].trate  [0] = betn;
   spnt[6].ratemul[0] = 2.0 * n;	
   spnt[6].rateo  [0] = 1;	/* taum, taub */
   spnt[6].trans  [1] = 7;
   spnt[6].trate  [1] = aln;
   spnt[6].ratemul[1] = 1.0 * n;
   spnt[6].rateo  [1] = 0;	/* taum, taua */
   spnt[6].trans  [2] = 2;
   spnt[6].trate  [2] = ald;
   spnt[6].ratemul[2] = 1.0 * h;
   spnt[6].rateo  [2] = 2;	/* tauh */

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 3;			/* state 5 = inactivated state*/
   spnt[7].trate  [0] = ald;
   spnt[7].ratemul[0] = 1.0 * h;
   spnt[7].rateo  [0] = 2;	/* tauh */
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = betn;
   spnt[7].ratemul[1] = 3.0 * n;
   spnt[7].rateo  [1] = 1;	/* taum */

   return ch;
}

/*----------------------------------------*/
