/* segment changly in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
					int nq, double bt);

double rnt(chan *spnt);
double rt(chan *spnt);

/*--------------------------------------------*/

double calcgly1(double v, int func) {}

/*--------------------------------------------*/

chantype *makgly1(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

/* Taken from:

		Busch and Sakmann (1990) 
			[cited by Destexhe, Mainen and Sejnowski (1998),
                              Kinetic Models of Synaptic Transmission,
                              in Methods in Neuronal Modeling, 
                              Ed. by Koch and Segev ]

   in which a 5-state GABA-A Markov diagram is provided:

                  gaba      gaba
                0  <->  1   <->   2

                        |         |

                        3         4     Open states

*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     int nstate, nparm, nq;


   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(GLY,1,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PCL] = 1.0;			/* permeability to Cl- ions */
   ch->ions->ionp[PK]  = 0.0;			/* permeability to K+ ions */
   ch->unitary = dglyu;
   ch->trconc = 100e-6;			/* default max trans conc */
   ch->vrev = vcl;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcgly1;       /* default rate function */
   parm[0].funcname = (char *)"calcgly1";   /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GLY-GLU] = 1.0;		/* response to GLY */
   respamp[STRY-GLU] = -1.0;		/* response to strychnine */
   // respamp[PTX-GLU] = -1.0;		/* response to picrotoxin */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 20e6;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 10e6;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 4.6e3;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 3;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = 3.3e3;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 4;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 10.6e3;
   spnt[2].rateo  [0] = 2;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 9.2e3;
   spnt[2].rateo  [1] = 1;	

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 1;			/* state 3 = the open state */
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 9.8e3;
   spnt[3].rateo  [0] = 2;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 2;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 410;	
   spnt[4].rateo  [0] = 2;

   return ch;

}

/*--------------------------------------------*/

