/* segment chansyn in program nc */

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
void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
                                double (*rrate)(chan *cpnt), double rate);
double rchanf(chan *cpnt);
double rchanr(chan *cpnt);

double qrate(chanparm *chp);

/*--------------------------------------------*/

double calcsyn (double v, int func)

/* Calculate rate function for synaptic channel,
   given voltage in mv.  Return rate calib in 1/sec;

   The "func" parameter is not used here.
*/

{
  return 1.0;
}

/*--------------------------------------------*/

double rsynf(chan *cpnt)

/* Forward binding rate constant for channel opening.
   Return conductance timinc and scaling factor for use 
   in dochan2().

   Not used unless trconc (transmitter conc.) is required
   (i.e. because state transition functions require it).

 */

{
   double nt;
   synap *spnt;

 if (!cpnt) return (0.0);
 if (spnt=(synap*)cpnt->comp2) {
   nt=spnt->conduct * spnt->trconc;     	/* find postsyn conductance */
 }
 else nt = 0;
 return (nt*qrate(cpnt->chtyp->parm));		/* rate per second */
}

/*--------------------------------------------*/

double rsynr(chan *cpnt)

/* Reverse binding rate constant for channel opening.
   Return conductance timinc and scaling factor for use 
   in dochan2().

   Not used unless trconc (transmitter conc.) is required
   (i.e. because state transition functions require it).

 */

{
   double nt;
   synap *spnt;

 if (!cpnt) return (0.0);
 if (spnt=(synap*)cpnt->comp2) {
   nt=(1.0 - spnt->conduct) * spnt->trconc;     /* find postsyn conductance */
 }
 else nt = 0;
 return (nt*qrate(cpnt->chtyp->parm));		/* rate per second */
}

/*--------------------------------------------*/

chantype *maksyn1(void)

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

  {
     double rsyn(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double a,b,c,d;
     int nstate, nparm, nq;

   nstate = 2;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(SYN2,1,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PNA] = 2.0/3;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 0.5;		/* permeability to K+ ions */
					/*  set to give vrev=0 */
   ch->ions->ionp[PCA] = dpcasyn2;	/* permeability to Ca++ ions */

   ch->unitary = 15e-12;
   ch->cabnd   = 1;			/* needs Ca binding */
   ch->cgbnd   = 1;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcsyn;         /* default rate function */
   parm[0].funcname = (char *)"calcsyn";        /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   mak2state(spnt,rchanf,rchanr,CHANRATE);
   return ch;
}

