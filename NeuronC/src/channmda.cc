/* segment channmd in program nc */

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
void chanrate(chan *cpnt, chanparm *chp);
double qrate(chanparm *chp);
double rnt(chan *spnt);
double rt(chan *spnt);
int set_perms (chantype *chtyp, iontab *ions, double svrev);

/*--------------------------------------------*/

double calcnmda1(double v, int func) 

{
   double val, kplus, kminus, vm, mg;

switch (func) {

   case 1:
    val = 6.1e5 * exp (-v / 17.0);		/* kplus */
    break;

   case 2:
    val = 5.4e3 * exp (v  / 47.0);		/* kminus */
    break;
  }
  return val;

}

/*--------------------------------------------*/

double almg(chan *cpnt)

/* Binding rate contants taken from 
	Ascher and Nowak, (1988), J. Physiol 399: 247.
	v in volts, mg in M, kplus is in /M/Sec.

  Multiply output rate by qrate for use in dochan2().  

*/

{
   double kplus, mg;
   chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  mg = ((ntchan *)cpnt)->mg;
  kplus = chp->fval[0] * mg;
  return (kplus);
}

/*--------------------------------------------*/

double betmg(chan *cpnt)

{
   double kminus, vm;
   chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  kminus = chp->fval[1];
  return (kminus);
}

/*--------------------------------------------*/

chantype *maknmda1(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.

    Markov state diagram:  0 <-> 1 <-> 2 <-> 3 <-> 4
                                       |
                                       |
                                       5

    Here we use Clements and Westbrook (1991) for all the
rates except 3<->4 which are from Ascher and Nowak (1988).

    state 0 = unbound receptor
          1 = receptor bound with 1 ligand
          2 = receptor bound with 2 ligands
          3 = open state
          4 = blocked by Mg
          5 = inactivated
          
*/
  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     int nstate, nparm, nq;
     double *respamp;
     double rb,ru;


   nstate = 6;                          /* 6 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(NMDA,1,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->unitary = 20e-12;
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = dpcanmda;	/* Ca perm, Jahr & Stevens 1993 */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   ch->trconc = 100e-6;
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 2;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcnmda1;       /* default rate function */
   parm[0].funcname = (char *)"calcnmda1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */
					/* Q10 for NMDA receptor cond & open */
					/* time is "in excess of 2.5" */
					/*  McLarnon JG. Curry K. (1990) */
					/*  Exp Br Res 82:82-88 */

   respamp[GLU-GLU] = 1.0;              /* response to glutamate */
   respamp[AMPA-GLU] = 0.01;            /* response to AMPA */
   respamp[NMDA-GLU] = 0.2;             /* response to NMDA */
   respamp[CNQX-GLU] = 0.0;             /* response to CNQX */

   rb = 5e6;				/* 4.9 /uM/sec */
   ru = 5e0;				/* set so that Kd is 1 uM */
 
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 2*rb;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = rb;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ru;
   spnt[1].rateo  [1] = 0;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 100;
   spnt[2].rateo  [0] = 3;	/* taud sets flicker rate */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2*ru;
   spnt[2].rateo  [1] = 0;
   spnt[2].trans  [2] = 5;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = 4;
   spnt[2].rateo  [2] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;	/* the conducting state */
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = almg;
   spnt[3].ratemul[0] = 1;
   spnt[3].rateo  [0] = 2;	/* tauc sets Mg flicker rate */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 300;
   spnt[3].rateo  [1] = 3;	/* taud sets flicker rate */

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = betmg;
   spnt[4].ratemul[0] = 1;
   spnt[4].rateo  [0] = 2; 	/* tauc sets Mg flicker rate */

   spnt[5].numtrans   = 1;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 2;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 0.3;
   spnt[5].rateo  [0] = 1;

   return ch;
}

/*--------------------------------------------*/

chantype *maknmda2(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.

                                                       10    11 
                                                       |  / mg
                                                       | /
                                                       4   <->  5
                                                    / 
                             nmda  nmda          / 
    Markov state diagram:  0 <-> 1 <->  2 <->  3  <->  6   <->  7
                                               | \     |  \
                                               |    \  10  11 mg
                                               |      \
                                               |       8   <->  9 
                                               |    /    \
                                               | /         \ mg
                                               10           11 
    State 0 is unbound 
          1,2 are bound
          3 is bound but closed
          4, 6 and 8 are open
	  5, 7, 9 are flicker closed
          10 is inactivated
          11 is Mg blocked

    taum (taua) = activation rate multiplier   (0)
    tauh (taub) = deactivation rate multiplier (1)
    tauc        = Mg flicker rate              (2)
    taud        = Open-closed flicker rate     (3)

   Kleckner and Pallotta have three "open, burst" arms to set
burst kinetics. Clements and Westbrook and Lin and Stevens show
the desensitized state but don't have multiple open states, so we
divide rates for the several paths out of the desensitized state.
We use Lin and Stevens' rates to and from the desensitized state.
We use Asher and Nowak's rates for mg block, dividing the return
rates in thirds because there are 3 open states.

   In the scheme used here, 0->1 and 1->2 are set to be a ligand
binding transition, whereas 2->3 is set to be a non ligand
binding transition.  Clements and Westbrook set these to 2*rb and
1*rb respectively whereas Jackson sets them to 2.56 and 14.7.
Assuming a 1 uM concentration of ligand (used by Klekner and
Pallotta and implicitly by Jackson), the rates are rb =
4.9e6/M/sec and 0->1 uses 0.5*rb and 1->2 uses 3*rb.
Alternately, the 0->1 transition could be assigned a fixed rate
and the 2->3 transition could be a ligand-binding rate. 

   We assume a constant glycine concentration of 10 uM. This
implies a glycine binding rate of 60-110/sec (Lester et al. 1993).

    See: Jackson MB, (1997) Biophys J. 73: 1382-1394.
         Kleckner NW and Pallotta BS, (1995) J. Physiol. 486:411-426. 
         Jahr and Stevens CF, (1990a) J Neurosci. 10: 1830-1837.
         Jahr and Stevens CF, (1990b) J Neurosci. 10: 3178-3182.
         Jahr and Stevens CF, (1993) PNAS 90: 11573-11577.
         Lin and Stevens CF, (1994) J. Neurosci. 14: 2153-2160. 	
	 Clements JD and Westbrook DL (1991) Neuron 7: 605-613.
	 Ascher P and Nowak L (1988) J. Physiol 399: 247-266. 
	 Lester RAJ, Tong G, Jahr CE (1993) 13:1088-1096. 
*/
  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double rb,ru;
     int nstate, nparm, nq;

   rb = 4.9e6;				/* 4.9 /uM/sec */
   ru = 4.9;				/* set so that Kd is 1 uM */
 
   nstate = 12;                         /* 12 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(NMDA,2,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->unitary = 20e-12;
   ch->ions->ionp[PNA] = 1;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1;		/* permeability to K+ ions */
   ch->ions->ionp[PCA] = dpcanmda;	/* Ca perm, Jahr & Stevens 1993 */
   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   ch->trconc = 100e-6;
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 2;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcnmda1;       /* default rate function */
   parm[0].funcname = (char *)"calcnmda1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */
					/* Q10 is "in excess of 2.5" */
					/*  McLarnon JG. Curry K. (1990) */
					/*  Exp Br Res 82:82-88 */
   respamp[GLU-GLU] = 1.0;              /* response to glutamate */
   respamp[AMPA-GLU] = 0.01;            /* response to AMPA */
   respamp[NMDA-GLU] = 0.2;             /* response to NMDA */
   respamp[CNQX-GLU] = 0.0;             /* response to CNQX */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 0.5*rb;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 3*rb;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 3.14;
   spnt[1].rateo  [1] = 0;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 98.3;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 47.6;
   spnt[2].rateo  [1] = 0;

   spnt[3].numtrans   = 5;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 170;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 6;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 600;
   spnt[3].rateo  [1] = 0;
   spnt[3].trans  [2] = 8;
   spnt[3].trate  [2] = rt;
   spnt[3].ratemul[2] = 230;
   spnt[3].rateo  [2] = 0;
   spnt[3].trans  [3] = 2;
   spnt[3].trate  [3] = rt;
   spnt[3].ratemul[3] = 983;
   spnt[3].rateo  [3] = 0;
   spnt[3].trans  [4] = 10;
   spnt[3].trate  [4] = rt;
   spnt[3].ratemul[4] = 1.5;	
   spnt[3].rateo  [4] = 1;	/* offseth, tauh sets deactivation rate */

   spnt[4].numtrans   = 4;
   spnt[4].cond       = 1.0;	/* conducting state */
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 6250;
   spnt[4].rateo  [0] = 3;	/* taud sets flicker rate */
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 6250;
   spnt[4].rateo  [1] = 3;
   spnt[4].trans  [2] = 10;
   spnt[4].trate  [2] = rt;
   spnt[4].ratemul[2] = 48;	/* total with 6->10, 8->10 to = 47.3/sec */
   spnt[4].rateo  [2] = 1;	/*  Lin and Stevens (1994) */
   spnt[4].trans  [3] = 11;
   spnt[4].trate  [3] = almg;
   spnt[4].ratemul[3] = 1;
   spnt[4].rateo  [3] = 2;	/* tauc sets Mg flicker rate */

   spnt[5].numtrans   = 1;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 16670;
   spnt[5].rateo  [0] = 3;	/* taud sets flicker rate */

   spnt[6].numtrans   = 4;
   spnt[6].cond       = 1.0;	/* conducting state */
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = 210;
   spnt[6].rateo  [0] = 3;	/* taud sets flicker rate */
   spnt[6].trans  [1] = 3;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = 210;
   spnt[6].rateo  [1] = 3;
   spnt[6].trans  [2] = 10;
   spnt[6].trate  [2] = rt;
   spnt[6].ratemul[2] = 48;	/* total with 4->10, 8->10 to = 47.3/sec */
   spnt[6].rateo  [2] = 1;	/*  Lin and Stevens (1994) */
   spnt[6].trans  [3] = 11;
   spnt[6].trate  [3] = almg;
   spnt[6].ratemul[3] = 1;
   spnt[6].rateo  [3] = 2;	/* tauc sets Mg flicker rate */

   spnt[7].numtrans   = 1;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 6;
   spnt[7].trate  [0] = rt;
   spnt[7].ratemul[0] = 16670;
   spnt[7].rateo  [0] = 3;

   spnt[8].numtrans   = 4;
   spnt[8].cond       = 1.0;	/* conducting state */
   spnt[8].trans  [0] = 9;
   spnt[8].trate  [0] = rt;
   spnt[8].ratemul[0] = 110;
   spnt[8].rateo  [0] = 3;	/* taud sets flicker rate */
   spnt[8].trans  [1] = 3;
   spnt[8].trate  [1] = rt;
   spnt[8].ratemul[1] = 110;
   spnt[8].rateo  [1] = 3;
   spnt[8].trans  [2] = 10;
   spnt[8].trate  [2] = rt;
   spnt[8].ratemul[2] = 48;	/* total with 4->10, 6->10 = 47.3/sec */
   spnt[8].rateo  [2] = 1;
   spnt[8].trans  [3] = 11;
   spnt[8].trate  [3] = almg;
   spnt[8].ratemul[3] = 1;	
   spnt[8].rateo  [3] = 2;	/* tauc sets Mg flicker rate */

   spnt[9].numtrans   = 1;
   spnt[9].cond       = 0;
   spnt[9].trans  [0] = 8;
   spnt[9].trate  [0] = rt;
   spnt[9].ratemul[0] = 16670;
   spnt[9].rateo  [0] = 3;	/* taud sets flicker rate */

   spnt[10].numtrans   = 4;
   spnt[10].cond       = 0;
   spnt[10].trans  [0] = 3;
   spnt[10].trate  [0] = rt;
   spnt[10].ratemul[0] = .2;	/* total with 10->4,6,8 to = .5/sec */
   spnt[10].rateo  [0] = 1;	/*  Lin and Stevens, (1994) */ 
   spnt[10].trans  [1] = 4;
   spnt[10].trate  [1] = rt;
   spnt[10].ratemul[1] = 0.1;
   spnt[10].rateo  [1] = 1;	/* tauh sets deactivation rate */
   spnt[10].trans  [2] = 6;
   spnt[10].trate  [2] = rt;
   spnt[10].ratemul[2] = 0.1;
   spnt[10].rateo  [3] = 1;	/* tauh sets deactivation rate */
   spnt[10].trans  [3] = 8;
   spnt[10].trate  [3] = rt;
   spnt[10].ratemul[3] = 0.1;
   spnt[10].rateo  [3] = 1;	/* tauh sets deactivation rate */

   spnt[11].numtrans   = 3;
   spnt[11].cond       = 0;
   spnt[11].trans  [0] = 4;
   spnt[11].trate  [0] = betmg;
   spnt[11].ratemul[0] = 0.333;
   spnt[11].rateo  [0] = 0;
   spnt[11].trans  [1] = 6;
   spnt[11].trate  [1] = betmg;
   spnt[11].ratemul[1] = 0.333; 
   spnt[11].rateo  [1] = 0;
   spnt[11].trans  [2] = 8;
   spnt[11].trate  [2] = betmg;
   spnt[11].ratemul[2] = 0.333; 
   spnt[11].rateo  [2] = 0;

   return ch;
}

/*--------------------------------------------*/


/* see
 
  Popescu G, Auerbach A. (2004) The NMDA receptor gating machine: lessons from single channels. Neuroscientist. 2004 Jun;10(3):192-198.

*/

