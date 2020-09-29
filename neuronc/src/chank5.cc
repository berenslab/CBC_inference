/* segment chank5 in program nc */
/* The K channel in horizontal cells */
/* Like Ih, activated by hyperpolarization, but has a K reversal potential */

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

/*----------------------------------------*/

double calck5n (double v, int func)

/* Calculate K 5 rate functions given voltage in mv.
/* for K ir channel in fish horizontal cells */

/* Approximately matches behavior of Ih channels in 
   Dong C-J, Werblin FS (1995)

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
    y = -0.1 * (-v-60.); 
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* singularity when x==0 */
          val =  y / x;
    else  val = 1.0;
    val *= 0.6;
    break;
						/* beta n */
  case 2:               
    val = 0.5 * exp ((-v-40) / -30.);
    break;
  }
  return val * 30 * HHRATE; 
}

/*----------------------------------------*/

chantype *makk5(void)

/* The K inward rectifier in horizontal cells: Ihz */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 3;                          /* 3 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(K,5,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dkihu;
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0.0001; 		/* permeability to Na+ ions */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck5n;         /* default rate function */
   parm[0].funcname = (char *)"calck5n"; /* user rate function */
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

/*----------------------------------------*/
