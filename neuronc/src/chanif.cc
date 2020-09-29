/* segment initchan in program nc */

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
double ccavoff(chan *ch);
double vext(conn *ch);

/*--------------------------------------------*/

#define ITAU 2e3

static int vi = 0;

double aif(chan *cpnt)

/* alpha for integrate-and-fire channel */

{
   double v,v1,v2,retval;

   v1 = cpnt->conc[0].cest; 
   v2 = cpnt->conc[1].cest; 
   v = cpnt->comp1->v - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsm;
   if (v > 0 && ! v2) v1 += timinc * ITAU * v;
   if (v1 > 1.0) {
   	 v2 = 1; 
	 retval = .3;
   }
   else  retval = 0;
   cpnt->conc[0].cest = v1; 
   cpnt->conc[1].cest = v2; 
  return retval;
}

/*----------------------------------------*/

double bif(chan *cpnt)

/* beta for integrate-and-fire channel */

{
   double v1, v2, retval;

   v1 = cpnt->conc[0].cest; 
   v2 = cpnt->conc[1].cest; 
   if (v2) {
         v1 -= timinc * 1e4; 
   }
   if (v1 < 0.1) {
       v2 = 0; 
       retval = 0.3;
   }
   else  retval= 0;
   cpnt->conc[0].cest = v1; 
   cpnt->conc[1].cest = v2; 
   return retval;
}

/*--------------------------------------------*/

chantype *makif(void)

/* Sequential-state description of 
   integrate-and-fire channel.

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans =   State to go to on this transition. 
    Set trate =   Function that returns basic rate for transition. 
    Set ratemul = Multiplier for rate function.
    Set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 2;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(NA,21,nstate,nparm,nq,dbasetc); /* make chan const, parm */
   ch->unitary = 20e-12;
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqm;                    /* Q10 for m rate function */

                                /*   0 <-> 1 */
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = aif;
   spnt[0].ratemul[0] = 1.0;
   spnt[0].rateo  [0] = 0;	/* offsm, taum */

   spnt[1].numtrans   = 1;
   spnt[1].cond       = 1;
   spnt[1].trans  [0] = 0;
   spnt[1].trate  [0] = bif;
   spnt[1].ratemul[0] = 1.0;
   spnt[1].rateo  [0] = 1;	/* offsm, taum */

   return ch;
}

/*----------------------------------------*/
