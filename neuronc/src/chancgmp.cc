/* segment chancgmp in program nc */

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

double al0 (chan *cpnt);
double bet0 (chan *cpnt);
				/* Opening rate constant for CNG channel */
				/*  Dependent only on Vm.  */
static double (*betcgmp)(chan*)  = al0; 

static double (*cgmp2_k3)(chan*)  = al0; 
static double (*cgmp2_k4)(chan*)  = bet0; 

double rnt(chan *cp);     	/* rate functions */
double rt(chan *cp);    	/*  function of time only */ 
double rca(chan *cpnt);		/* Ca gates transition func */

chantype *makchantype(int ctype,int cnum,int nstates,int nparm, 
					int nq, double bt);
double getntresp (chan *ch, comp *cpnt);
int     isntresp (chan *ch, comp *cpnt);

double qrate(chanparm *chp);

/*--------------------------------------------*/

double calccgmp (double v, int func)

/* Calculate "beta" rate function for cGMP-gated channel,
   given voltage in mv.  Return rate calib in 1/sec;

   The "func" parameter is not used here because there is only
   one function for the cGMP channel.

   z = .23

*/

{
    double rate;

  rate = exp (0.23 * 0.001 * v * FR / ktemp);
  return rate;

}

/*--------------------------------------------*/

double calccgmp2 (double v, int func)

/* Calculate opening and closing rate function for cGMP-gated channel,
   given voltage in mv.  Return rate calib in 1/sec;

   The "func" parameter is used here because there are 
   2 voltage-gated functions for the cGMP channel.

   z = .23

     rate = exp (0.5 * 0.23 * 0.001 * v * FR / ktemp);
     rate = exp (1.15e-4            * v * FR / ktemp);
*/

{
    double rate;

  switch (func) {

   case 1:
     rate = exp ( (0.5 * 0.23 * 0.001 * FR) * v / ktemp);
     break;

   case 2:
     rate = exp (-(0.5 * 0.23 * 0.001 * FR) * v / ktemp);
     break;

  }
  return rate;
}

/*--------------------------------------------*/

chantype *makcgmp1(void)

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
                Taylor, W.R., and Baylor, D.A. (1995) 
                J. Physiol. 483: 567-582.

   in which a 7-state cGMP-gated markov diagram is provided.
Only alpha and beta for the transition between the open and
closed state were provided.  Rates for cGMP binding were inferred
(RGS) from probability of the activated state (includes both open
and closed state) at different ligand concentrations and from
kinetics (Baylor, 1988).  Inactivated states were arranged
according to frequency and duration reported (RGS 10/98).

For info on conductance, see: 

  Ruiz, M. and Karpen, J.W. (1997), Nature 389: 389-92.

  Karpen has recent manuscript (now in review) on cGMP binding rates.
       (RGS 11/98)

                                           7
                                           
                                           |
                                            
                  0   <->  1  <->  2  <->  3  <->  6
                                           
                                          
                                   |       | 
                                    
                                   
                                   4  <->  5

State 0: Closed non-bound
      1: Closed 1 cGMP bound
      2: Closed 2 cGMP bound
      3: Closed 3 cGMP bound
      4: partially Open  (25%)
      5: fully Open  (100%)
      6: short Inactivated
      7: long Inactivated
*/

  {
     double rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 8;                          /* Number of Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,1,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  5e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 1;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   a = 250e6;
   b = 8.0e2;
   c = 10;
   d = 0.2;

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 1 * a;
   spnt[0].rateo  [0] = 0;
                                     /*                        7            */
   spnt[1].numtrans   = 2;           /*                                     */
   spnt[1].cond       = 0;           /*                        |            */ 
   spnt[1].trans  [0] = 2;           /*                                     */
   spnt[1].trate  [0] = rnt;         /* 0 <->  1  <->  2  <->  3  <->  6    */
   spnt[1].ratemul[0] = 1 * a;       /*                                     */
   spnt[1].rateo  [0] = 0;	     /*                                     */
   spnt[1].trans  [1] = 0;	     /*                |       |            */
   spnt[1].trate  [1] = rt;          /*                                     */
   spnt[1].ratemul[1] = 5 * b;       /*                                     */
   spnt[1].rateo  [1] = 0;	     /*                4  <->  5            */

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 4 * a;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 5 * b;
   spnt[2].rateo  [1] = 0;
   spnt[2].trans  [2] = 4;
   spnt[2].trate  [2] = betcgmp;	/* open rate const from calccgmp() */
   spnt[2].ratemul[2] = 5e4;		//1.51e4;
   spnt[2].rateo  [2] = 1;

   spnt[3].numtrans   = 4;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 6;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 4 * c;		/* transition to short I */
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 1 * b;
   spnt[3].rateo  [1] = 0;
   spnt[3].trans  [2] = 5;
   spnt[3].trate  [2] = betcgmp;	/* opening rate constant */
   spnt[3].ratemul[2] = 5e4;		//1.51e4;
   spnt[3].rateo  [2] = 1;
   spnt[3].trans  [3] = 7;
   spnt[3].trate  [3] = rt;
   spnt[3].ratemul[3] = 4 * d;		/* transition to long I */
   spnt[3].rateo  [3] = 1;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0.25;		/* partially open state */
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 4 * a;	
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 2;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 2.1e4;		/* closing rate constant */
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;		/* the open state */
   spnt[5].trans  [0] = 3;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 2.1e4;		/* closing rate constant */
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1 * b;
   spnt[5].rateo  [1] = 0;

   spnt[6].numtrans   = 1;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 3;		/* state 6 = short Inactivated state */
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = 30 * c;
   spnt[6].rateo  [0] = 1;

   spnt[7].numtrans   = 1;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 3;		/* state 7 = long Inactivated state */
   spnt[7].trate  [0] = rt;
   spnt[7].ratemul[0] = 20 * d;
   spnt[7].rateo  [0] = 1;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp2(void)

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
		Benndorf, K., Koopmann, R., Eismann, E., and Kaupp, U.B. (1999)
		J. Gen. Physiol. 114: 477-489.

                He, Y., Ruiz, M.L., and Karpen J.W. (2000) 97: 895-900.

   in which allosteric, sequential 5-state, and sequential
6-state (G3O and G4O) models of the rod cGMP-gated channel made
from its alpha subunit are described.  Since the channel is now
known to have 4 binding sites, and the G3O model was "as good" as
the allosteric model, and better than the G4O model, they state
that the allosteric model is the best (since it has 4 cGMP
binding transitions.)  We implement here the G3O model because it
has fewer states and will run faster.  It probably doesn't
correctly represent the noise correctly but that can wait.  Note
that the wild-type channel has a configuration a-b-a-b (He et al,
2000), implying that the channel implemented here is probably a
bit different than the wild type.

 The 5-state cGMP-gated markov diagram (G30):

                                            open
                 k1                      k3
            0   <->  1  <->  2  <->  3  <->  4
                 k2                      k4    

State 0: Closed non-bound
      1: Closed 1 cGMP bound
      2: Closed 2 cGMP bound
      3: Closed 3 cGMP bound
      4: Open 

*/  

  {
     double rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 5;                          /* Number of Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,2,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 2;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp2;       /* default rate function */
   parm[0].funcname = (char *)"calccgmp2";      /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e7;
   k2 = 1.5e3;
   k3 = 1.39e3;
   k4 = 1.67e2;

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;
                                 
   spnt[1].numtrans   = 2;      
   spnt[1].cond       = 0;     
   spnt[1].trans  [0] = 2;    
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 0;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 0;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = cgmp2_k3;
   spnt[3].ratemul[0] = 1 * k3;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 3 * k2;
   spnt[3].rateo  [1] = 0;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 1.0;		/* open state */
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = cgmp2_k4;
   spnt[4].ratemul[0] = 1 * k4;	
   spnt[4].rateo  [0] = 1;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp3(void)

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

   cGMP type 2 (above), but simplified, and with 1 Ca-bound state added.


               cG      cG      cG 
            0  <->  1  <->  2  <->  3
 
                                    |  Ca
                                       
                                    4
 
State 0: Closed non-bound
      1: Closed 1 cGMP bound
      2: Closed 2 cGMP bound
      3: Open, 3 cGMP 
      5: Closed, 1 Ca bound

*/  

  {
     double rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 5;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,3,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e8;
   k2 = 1.5e3;
   k3 = 5e6;		/* [Ca] * Ca binding rate */
   k4 = 1e0;		/* [Ca] * Ca unbinding rate */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 3 * k2;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 4;
   spnt[3].trate  [1] = rca;	
   spnt[3].ratemul[1] = k3;
   spnt[3].rateo  [1] = 2;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 0;		/* open state, 1 Ca bound */
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = k4;		/* closing rate constant */
   spnt[4].rateo  [0] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp4(void)

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

   cGMP type 2 (above), but simplified, and with 2 Ca-bound states added.


               cG      cG      cG 
            0  <->  1  <->  2  <->  3
 
                                    |  Ca
                                       
                                    4
 
                                    |  Ca
                                       
                                    5
State 0: Closed non-bound
      1: Closed 1 cGMP bound
      2: Closed 2 cGMP bound
      3: Open,  3 cGMP bound
      4: partially Closed, 1 Ca bound
      6: Inactivated, 2 Ca bound

*/  

  {
     double rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,c,d,e;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 6;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,4,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e7;
   k2 = 1.5e3;
   c = 10;
   d = 0.2;
   e = 5e5;		/* [Ca] * Ca binding rate */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 3 * k2;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 4;
   spnt[3].trate  [1] = rca;	
   spnt[3].ratemul[1] = e;
   spnt[3].rateo  [1] = 2;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0.5;		/* open state, 1 Ca bound */
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 1e0;		/* closing rate constant */
   spnt[4].rateo  [0] = 3;
   spnt[4].trans  [1] = 5;
   spnt[4].trate  [1] = rca;
   spnt[4].ratemul[1] = e;
   spnt[4].rateo  [1] = 2;

   spnt[5].numtrans   = 1;
   spnt[5].cond       = 0;		/* closed, 2 Ca bound */
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 1e0;
   spnt[5].rateo  [0] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp5(void)

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

   cGMP type 4 (above), but added states with 2 cGMP and 1, 2 Ca bound.


               cG      cG      cG 
            0  <->  1  <->  2  <->  3
 
                            |       |  Ca
                                       
                            4  <->  5
 
State 0: Closed, non-bound
      1: Closed, 1 cGMP, bound
      2: Closed, 2 cGMP, bound
      3: Open,   3 cGMP, bound
      4: Closed, 2 cGMP, 1 Ca bound
      5: Closed, 3 cGMP, 1 Ca bound

*/  

  {
     double rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 6;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,5,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * .3;     /* cGMP binding rate */
   k2 = 1.5e4 * .3;     /* cGMP unbinding rate */
   k3 =   5e6;		/* Ca binding rate */
   k4 =   .5e0;         /* Ca unbinding rate */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 1;
   spnt[2].trans  [2] = 4;
   spnt[2].trate  [2] = rca;
   spnt[2].ratemul[2] = 1 * k3;
   spnt[2].rateo  [2] = 2;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 3 * k2;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 5;
   spnt[3].trate  [1] = rca;	
   spnt[3].ratemul[1] = 1 * k3;
   spnt[3].rateo  [1] = 2;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 1 * k1;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 2;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;		/* 1 Ca bound */
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 3 * k2;
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 3;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[5].rateo  [1] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp6(void)

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

   cGMP type 4 (above), but added states with 2, 3 cGMP and 1 Ca bound.


               cG      cG      cG 
            0  <->  1  <->  2  <->  3
 
            |       |       |       |  Ca
                                       
            4  <->  5  <->  6  <->  7

State 0: Closed, non-bound
      1: Closed, 1 cGMP, bound
      2: Closed, 2 cGMP, bound
      3: Open,   3 cGMP, bound
      4: Closed, 0 cGMP, 1 Ca bound
      5: Closed, 1 cGMP, 1 Ca bound
      6: Closed, 2 cGMP, 1 Ca bound
      7: Closed, 3 cGMP, 1 Ca bound, Inactivated
*/  

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 8;                          /* 8 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,6,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * .3;
   k2 = 1.5e4 * .3;
   k3 =   5e6;			/* [Ca]  Ca binding rate */
   k4 =   0.5e0;		/* [Ca]  Ca unbinding rate */

   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 4;
   spnt[0].trate  [1] = rca;
   spnt[0].ratemul[1] = k3;
   spnt[0].rateo  [1] = 2;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 5;
   spnt[1].trate  [2] = rca;
   spnt[1].ratemul[2] = k3;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 1;
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = rca;
   spnt[2].ratemul[2] = k3;
   spnt[2].rateo  [2] = 2;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 3 * k2;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 7;
   spnt[3].trate  [1] = rca;	
   spnt[3].ratemul[1] = k3;
   spnt[3].rateo  [1] = 2;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 3 * k1;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 0;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = k4;
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = rnt;
   spnt[5].ratemul[0] = 2 * k1;	
   spnt[5].rateo  [0] = 0;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1 * k2;
   spnt[5].rateo  [1] = 1;
   spnt[5].trans  [2] = 1;
   spnt[5].trate  [2] = rt;
   spnt[5].ratemul[2] = k4;
   spnt[5].rateo  [2] = 3;

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = rnt;
   spnt[6].ratemul[0] = 1 * k1;
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = 2 * k2;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 2;
   spnt[6].trate  [2] = rt;
   spnt[6].ratemul[2] = k4;
   spnt[6].rateo  [2] = 3;

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 6;
   spnt[7].trate  [0] = rt;
   spnt[7].ratemul[0] = 3 * k2;
   spnt[7].rateo  [0] = 1;
   spnt[7].trans  [1] = 3;
   spnt[7].trate  [1] = rt;
   spnt[7].ratemul[1] = k4;
   spnt[7].rateo  [1] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp7(void)

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

   cGMP type 4 (above), but added states with 2, 3 cGMP and 1 Ca bound.


               cG            cG      cG 
            0  <->     1    <->  2  <->  3
 
        *A3 | /A3  *A2 | /A2  *A | /A    |  Ca

               *A           *A       *A    
            4  <->     5    <->  6  <->  7
               /A           /A       /A

State 0: Closed, non-bound
      1: Closed, 1 cGMP, bound
      2: Closed, 2 cGMP, bound
      3: Open,   3 cGMP, bound
      4: Closed, 0 cGMP, 1 Ca bound
      5: Closed, 1 cGMP, 1 Ca bound
      6: Closed, 2 cGMP, 1 Ca bound
      7: Closed, 3 cGMP, 1 Ca bound, Inactivated

   Same as type 6 above but factor A added in Ca recovery.

*/  

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4,a,a2,a3;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 8;                          /* 8 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,7,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;                 /* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * .3;
   k2 = 1.5e4 * .3;
   k3 =   5e6;			/* [Ca]  Ca binding rate */
   k4 =   0.5e0;		/* [Ca]  Ca unbinding rate */
   a  = 1.5;
   a2 = a*a;
   a3 = a*a*a;

   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 4;
   spnt[0].trate  [1] = rca;
   spnt[0].ratemul[1] = k3 / a3;
   spnt[0].rateo  [1] = 2;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 5;
   spnt[1].trate  [2] = rca;
   spnt[1].ratemul[2] = k3/a2;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 1;
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = rca;
   spnt[2].ratemul[2] = k3/a;
   spnt[2].rateo  [2] = 2;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 3 * k2;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 7;
   spnt[3].trate  [1] = rca;	
   spnt[3].ratemul[1] = k3;
   spnt[3].rateo  [1] = 2;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 3 * k1 * a;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 0;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = k4*a3;
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = rnt;
   spnt[5].ratemul[0] = 2 * k1 * a;	
   spnt[5].rateo  [0] = 0;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1 * k2 / a;
   spnt[5].rateo  [1] = 1;
   spnt[5].trans  [2] = 1;
   spnt[5].trate  [2] = rt;
   spnt[5].ratemul[2] = k4*a2;
   spnt[5].rateo  [2] = 3;

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = rnt;
   spnt[6].ratemul[0] = 1 * k1 * a;
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = 2 * k2 / a;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 2;
   spnt[6].trate  [2] = rt;
   spnt[6].ratemul[2] = k4*a;
   spnt[6].rateo  [2] = 3;

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 6;
   spnt[7].trate  [0] = rt;
   spnt[7].ratemul[0] = 3 * k2 / a;
   spnt[7].rateo  [0] = 1;
   spnt[7].trans  [1] = 3;
   spnt[7].trate  [1] = rt;
   spnt[7].ratemul[1] = k4;
   spnt[7].rateo  [1] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp8(void)

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

   cGMP type 4 (above), but added states with 2, 3 cGMP and 1 Ca bound.


               cG      cG      cG 
            0  <->  1  <->  2  <->  3
 
                    |               |  Ca
                                       
                    4  <->  5  <->  6

State 0: Closed, non-bound
      1: Closed, 1 cGMP, bound
      2: Closed, 2 cGMP, bound
      3: Open,   3 cGMP, bound
      4: Closed, 1 cGMP, 1 Ca bound
      5: Closed, 2 cGMP, 1 Ca bound
      6: Closed, 3 cGMP, 1 Ca bound, Inactivated
*/  

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pna, pca;
     int nstate, nparm, nq;

   nstate = 7;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,8,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;			/* permeability to Ca++ ions */
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;			/* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * .3;
   k2 = 1.5e4 * .3;
   k3 =   5e6;			/* [Ca]  Ca binding rate */
   k4 =   0.5e0;		/* [Ca]  Ca unbinding rate */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * k1;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 4;
   spnt[1].trate  [2] = rca;
   spnt[1].ratemul[2] = 1 * k3;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1 * k1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2 * k2;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 3 * k2;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 5;
   spnt[3].trate  [1] = rca;	
   spnt[3].ratemul[1] = 1 * k3;
   spnt[3].rateo  [1] = 2;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 2 * k1;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 1;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;		/* 1 Ca bound */
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = rnt;
   spnt[5].ratemul[0] = 1 * k1;
   spnt[5].rateo  [0] = 0;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 2 * k2;		/* closing rate constant */
   spnt[5].rateo  [1] = 1;

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;		/* 1 Ca bound */
   spnt[6].trans  [0] = 3;
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = 1 * k4;
   spnt[6].rateo  [0] = 3;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = 3 * k2;		/* closing rate constant */
   spnt[6].rateo  [1] = 1;


   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp9(void)

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

   cGMP type 4 (above), but added states with 2, 3 cGMP and 1 Ca bound.


               cG      cG 
            0  <->  1  <->  2 
 
                    |       |  Ca
                                       
                    3  <->  4  

State 0: Closed, non-bound
      1: Closed, 1 cGMP  bound
      2: Open,   2 cGMP  bound
      3: Closed, 1 cGMP, 1 Ca bound
      4: Closed, 2 cGMP, 1 Ca bound
*/  

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 5;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,9,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  5e-6;			/* max conc of nt */
   pca = dpcacgmp;			/* permeability to Ca++ ions */
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;			/* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * .3;
   k2 = 1.5e4 * .3;
   k3 =   5e6;			/* [Ca]  Ca binding rate */
   k4 =   0.5e0;		/* [Ca]  Ca unbinding rate */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 2 * k1;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 1 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 3;
   spnt[1].trate  [2] = rca;
   spnt[1].ratemul[2] = 1 * k3;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 1.0;		/* the open state */
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 2 * k2;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 4;
   spnt[2].trate  [1] = rca;	
   spnt[2].ratemul[1] = 1 * k3;
   spnt[2].rateo  [1] = 2;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rnt;
   spnt[3].ratemul[0] = 1 * k1;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 1;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[3].rateo  [1] = 3;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;		/* 1 Ca bound */
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 2 * k2;
   spnt[4].rateo  [0] = 1;
   spnt[4].trans  [1] = 2;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[4].rateo  [1] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp10(void)

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

   cGMP type 4 (above), but added states with 2, 3 cGMP and 1 Ca bound.


               cG      cG 
            0  <->  1  <->  2 
 
            |               |  Ca
                                       
            3  <->  4  <->  5  

State 0: Closed, non-bound
      1: Closed, 1 cGMP  bound
      2: Open,   2 cGMP  bound
      3: Closed, 0 cGMP, 1 Ca bound
      4: Closed, 1 cGMP, 1 Ca bound
      5: Closed, 2 cGMP, 1 Ca bound
*/  

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 6;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,10,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;			/* permeability to Ca++ ions */
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;			/* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * .3;
   k2 = 1.5e4 * .3;
   k3 =   5e6;			/* [Ca]  Ca binding rate */
   k4 =   0.5e0;		/* [Ca]  Ca unbinding rate */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 2 * k1;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 3;
   spnt[0].trate  [1] = rca;
   spnt[0].ratemul[1] = 1 * k3;
   spnt[0].rateo  [1] = 2;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 1 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 1.0;		/* the open state */
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 2 * k2;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 5;
   spnt[2].trate  [1] = rca;	
   spnt[2].ratemul[1] = 1 * k3;
   spnt[2].rateo  [1] = 2;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rnt;
   spnt[3].ratemul[0] = 1 * k1;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 0;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[3].rateo  [1] = 3;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;		/* 1 Ca bound */
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 1 * k1;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 1 * k2;		/* closing rate constant */
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;		/* 1 Ca bound */
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 1 * k2;
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 2;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[5].rateo  [1] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makcgmp11(void)

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

   cGMP type 4 (above), but added states with 2, 3 cGMP and 1 Ca bound.


               cG        cG 
            0  <->    1  <->    2 
 
       *A^2 | / A^2             |  Ca

                *A        *A                
            3  <->    4  <->    5  
               /A         /A

State 0: Closed, non-bound
      1: Closed, 1 cGMP  bound
      2: Open,   2 cGMP  bound
      3: Closed, 0 cGMP, 1 Ca bound
      4: Closed, 1 cGMP, 1 Ca bound
      5: Closed, 2 cGMP, 1 Ca bound

    Same as type 10 above but factor A added in Ca recovery.
*/  

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k1,k2,k3,k4,A;
     double pca, pna;
     int nstate, nparm, nq;

   nstate = 6;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(CGMP,11,nstate,nparm,nq,BASETSYN);/* make chan const, parm */
   ch->unitary = 25e-12;
   ch->trconc =  50e-6;			/* max conc of nt */
   pca = dpcacgmp;			/* permeability to Ca++ ions */
   pna = 0.5089;
   ch->ions->ionp[PNA] = pna;                 /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0 - pna;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = pca;			/* permeability to Ca++ ions */
					/*  See Yau & Baylor, 1989 */
   ch->vrev = 0;			/* default reversal potential */
   ch->cabnd  = 0;			/* needs Ca binding */
   ch->cgbnd  = 0;			/* needs cGMP binding */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 1;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calccgmp;        /* default rate function */
   parm[0].funcname = (char *)"calccgmp";       /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */

   respamp[CGMP-GLU] = 1.0;              /* response to cGMP */
   respamp[CAMP-GLU] = 1e-2;             /* response to cAMP */

   k1 =   3e9 * 0.3;
   k2 = 1.5e4 * 0.3;
   k3 =   5e6 * 0.25;		/* [Ca]  Ca binding rate */
   k4 =   0.5e0 * 0.25;		/* [Ca]  Ca unbinding rate */
   A = 4;			/* recovery factor */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 2 * k1;
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 3;
   spnt[0].trate  [1] = rca;
   spnt[0].ratemul[1] = 1 * k3 /A/A;
   spnt[0].rateo  [1] = 2;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 1 * k1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1 * k2;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 1.0;		/* the open state */
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 2 * k2;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 5;
   spnt[2].trate  [1] = rca;	
   spnt[2].ratemul[1] = 1 * k3;
   spnt[2].rateo  [1] = 2;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rnt;
   spnt[3].ratemul[0] = 1 * k1 * A;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 0;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 1 * k4 * A*A;	/* closing rate constant */
   spnt[3].rateo  [1] = 3;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;		/* 1 Ca bound */
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 1 * k1 * A;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 1 * k2 / A;	/* closing rate constant */
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;		/* 1 Ca bound */
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 1 * k2 / A;
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 2;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1 * k4;		/* closing rate constant */
   spnt[5].rateo  [1] = 3;

   return ch;
}

