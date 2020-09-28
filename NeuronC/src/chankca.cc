/* segment chankca in program nc */

/* sets up KCa channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "control.h"
#include "ncio.h"

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
					int nq, double bt);
void makchanimpl (chanparm *chp);
void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
                                double (*rrate)(chan *cpnt), double rate);
double rchanf(chan *cpnt);
double rchanr(chan *cpnt);

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

extern int interp;		/* running interpreter */

double ncabs(double x);
double qrate(chanparm *chp);
double rca(chan *cpnt);
double rca2(chan *cpnt);
double rt(chan *cpnt);
void varcopyu(void);
attrib *getattrib(int elnum);
chantype *getchantype(int ctype, int stype);
double calcchaninf(double v, chanparm *chp, double arate, double brate);
double calcchantau(double v, chanparm *chp, double arate, double brate);
double ccavoff (chan *ch);
double vext(conn *ch);

/*----------------------------------------*/

double calckca1 (double v, int func)

{
  return 1.0;
}

/*----------------------------------------*/

double calckca2 (double v, int func)

{
  return 1.0;
}

/*----------------------------------------*/

#define KCA_RATE 50 

double akcacalc(double v, double ca, double rate, double d1, double k1)

/* Alpha for calcium-activated potassium channel,
      from Hines, 1989 and Moczydlowski and Latorre (1983).

   Given Vm in volts, return rate.
   (Original given units of mV.)

   Does not include effect of timestep and temperature.
*/


{
   double xalpha,r,vm;

    vm = v * MVOLT;
    if (ca==0.0) ca = 1e-9;
    r = rate;
    xalpha = r / (1.0 + k1/ca * exp (d1* -vm/10.));
//ncfprintf (stderr,"%g %g %g %g %g %g\n",v,ca,rate,d1,k1,xalpha);

    return xalpha * HHRATE;
}

/*----------------------------------------*/

double bkcacalc(double v, double ca, double rate, double d2, double k2)

/* Beta for calcium-activated potassium channel,
      from Hines, 1989 and Moczydlowski and Latorre (1983).

   Given Vm in volts, return rate.
   (Original given units of mV.)

   Does not include effect of timestep and temperature.
*/

{
   double xbeta,r,vm;

    vm = v * MVOLT;
    if (ca==0.0) ca = 1e-9;
    r = rate;
    xbeta = r / (1.0 + ca / (k2 * exp (d2* -vm/10.)));
    return xbeta * HHRATE;
}

double bkcacalc_new(double v, double ca, double rate, double d2, double k2)

/* Beta for KCa channels, returns constant rate. This gives
    decreasing activation tau with increasing activation. */ 

#define KCA_BETA 1.0 

{
    return rate * KCA_BETA * HHRATE;
}

/*----------------------------------------*/

double bkcainf(double v, double ca, mattrib *mpnt, double arate, double brate)

/* Calculate ninf for KCa channels, assuming that [Ca]i = dcai.  */
/* Used by kcainf() below. */

{
  double alpha, beta, inf;
  chattrib *kcapnt;

  kcapnt = (chattrib *)mpnt;
  alpha = akcacalc(v, ca, arate, kcapnt->d1, kcapnt->k1);
  beta  = bkcacalc(v, ca, brate, kcapnt->d2, kcapnt->k2);
  inf = alpha / (alpha + beta); 
  return inf; 
}

/*----------------------------------------*/

double bkcatau(double v, double ca, mattrib *mpnt, double arate, double brate, 
					chanparm *chp)

/* Calculate ntau for KCa channels, assuming that [Ca]i = dcai.  */
/* Used by kcatau() below. */

{
  double alpha, beta, tau, qkca;
  chattrib *kcapnt;

  kcapnt = (chattrib *)mpnt;
  qkca = qrate(chp);				/* calculate rate factor */
  alpha = qkca * akcacalc(v, ca, arate, kcapnt->d1, kcapnt->k1) * KCA_RATE;
  beta  = qkca * bkcacalc(v, ca, brate, kcapnt->d2, kcapnt->k2) * KCA_RATE;
  tau = timinc / (alpha + beta); 
  return tau; 
}

/*--------------------------------------------*/

datum kcainf (datum &v, datum &ca, datum &elname) 

/* Return the equilibrium activation for KCa chans
   at a given voltage v and [Ca]. */

{
  attrib *apnt;
  mattrib *mpnt;
  chantype *cht;
  chanparm *chp;
  double voff, arate, brate, val;
  int elnum;

  if (interp) varcopyu();
  elnum = (int)elname.val;
  if (!(apnt=getattrib(elnum))) {
    ncfprintf (stderr,"pinf: can't find element %d\n",elnum);
    v.val = 0; return v;
  }
  if (!(cht=getchantype(apnt->ctype, apnt->stype))) {
     ncfprintf (stderr,"pinf: incorrect channel type for element %d\n",elnum);
     v.val = 0; return v;
  }
  chp = &cht->parm[0];
  switch (apnt->ctype) {
   case CA:
   case K:
   case KCa:
   case NA:
	mpnt = (mattrib *)apnt; break;
   default: v.val=0; return v; break;
  }
  if (chp->pn) { 		/* h */
    voff =  mpnt->voffsh;
    arate = mpnt->tauc;
    brate = mpnt->taud;
  }
  else {			/* m */
    voff =  mpnt->voffsm;
    arate = mpnt->taua;
    brate = mpnt->taub;
  } 
  if (voff==NULLVAL)  voff = 0;
  if (arate==NULLVAL) arate = 1;
  if (brate==NULLVAL) brate = 1;
  arate = 1/arate;		/* rate is inverse tau */
  brate =  1/brate;
  switch (apnt->ctype) {	/* special case for Kca chans */

   case KCa: val = bkcainf(v.val-voff, ca.val, mpnt, arate, brate); break;

   default: val = calcchaninf(v.val-voff,chp,arate,brate);
	    break;
  }
  v.val = val;
  return v;
}

/*--------------------------------------------*/

datum kcatau (datum &v, datum &ca, datum &elname) 

/* Return the time constant of activation for KCa chans
   at a given voltage v and [Ca]. */

{
  attrib *apnt;
  mattrib *mpnt;
  chantype *cht;
  chanparm *chp;
  double voff, arate, brate, val;
  int elnum;

  if (interp) varcopyu();
  elnum = (int)elname.val;
  if (!(apnt=getattrib(elnum))) {
    ncfprintf (stderr,"pinf: can't find element %d\n",elnum);
    v.val = 0; return v;
  }
  if (!(cht=getchantype(apnt->ctype, apnt->stype))) {
     ncfprintf (stderr,"pinf: incorrect channel type for element %d\n",elnum);
     v.val = 0; return v;
  }
  chp = &cht->parm[0];
  switch (apnt->ctype) {
   case CA:
   case K:
   case KCa:
   case NA:
	mpnt = (mattrib *)apnt; break;
   default: v.val=0; return v; break;
  }
  if (chp->pn) { 		/* h */
    voff =  mpnt->voffsh;
    arate = mpnt->tauc;
    brate = mpnt->taud;
  }
  else {			/* m */
    voff =  mpnt->voffsm;
    arate = mpnt->taua;
    brate = mpnt->taub;
  } 
  if (voff==NULLVAL)  voff = 0;
  if (arate==NULLVAL) arate = 1;
  if (brate==NULLVAL) brate = 1;
  arate = 1/arate;		/* rate is inverse tau */
  brate =  1/brate;
  switch (apnt->ctype) {	/* special case for Kca chans */

   case KCa: val = bkcatau(v.val-voff, ca.val, mpnt, arate, brate, chp); break;

   default: val = calcchantau(v.val-voff,chp,arate,brate);
	    break;
  }
  v.val = val;
  return v;
}

/*----------------------------------------*/

void kcatab(double v, double ca, int typ, double nrate, 
		double d1, double d2, double k1, double k2, chanparm *chp)
{
     double a,b,ab,ab2,tval;
     double kca1, kca2, qkca;

    if (! chp->ival) {
	makchanimpl(chp);			/* make dummy tables */
    }
    if (nrate==0) nrate = 1e-4;
    qkca = qrate(chp);				/* calculate rate factor */
    a = qkca * akcacalc(v,ca,nrate,d1,k1) * KCA_RATE;
    b = qkca * bkcacalc(v,ca,nrate,d2,k2) * KCA_RATE;
    ab  = a+b;
    ab2 = ab*0.5;

    tval = (1.0 - ab2) /
    	   (1.0 + ab2); 

    if (tval>0) {
      kca1 = tval;				/* C-N implicit */
      kca2 = a / (1.0+ab2); 
    }
    else {					/* purely implicit */
      kca1 = 1.0 / (1.0+ab);
      kca2 = a   / (1.0+ab); 
    }
    chp->ival[0] = kca1;                        /* save values for impl calc */
    chp->ival[1] = kca2;
}

/*----------------------------------------*/

double alnca(chan *cpnt)
{
     double v;
     kcachan *kcapnts;
     cacomp *capnt;
     double ca,qkca,rval;

  if (!cpnt) return 0.0;
  v = cpnt->comp1->v - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsm;
  kcapnts = (kcachan*)cpnt;
  if (!(capnt=cpnt->comp1->capnt)) {
        ca = dcai;
  }
  else ca = capnt->cais[0];
  qkca = qrate(cpnt->chtyp->parm);		/* calculate rate factor */
  rval = (qkca*akcacalc(v,ca,1.0,kcapnts->d1,kcapnts->k1)); 
  return (rval);
}

/*----------------------------------------*/

double betnca(chan *cpnt)
{
     double v;
     kcachan *kcapnts;
     cacomp *capnt;
     double ca,qkca;

  if (!cpnt) return 0.0;
  v = cpnt->comp1->v - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsm;
  kcapnts = (kcachan*)cpnt;
  if (!(capnt=cpnt->comp1->capnt)) {
         ca = dcai;
  }
  else ca = capnt->cais[0];
  qkca = qrate(cpnt->chtyp->parm);		/* calculate rate factor */
  return (qkca*bkcacalc(v,ca,1.0,kcapnts->d2,kcapnts->k2));
}

/*----------------------------------------*/

double alnkca6(chan *cpnt)

/* Like alnca, but with only voltage sensitivity, no Ca sensitivity */

/*  Use for rate func for g (to open state). */
/*  Normalized correctly for 22 deg C */

/* According to Horrigan et al (1999), alpha and beta are same,
     but d and g differ by factor of 2.
*/

{
     double v;
     kcachan *kcapnts;
     cacomp *capnt;
     double qkca;
     double x,y,xalpha;

//#define GF 0.7
#define GF 1.0

  if (!cpnt) return 0.0;
  v = cpnt->comp1->v - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsm;
  kcapnts = (kcachan*)cpnt;
  qkca = qrate(cpnt->chtyp->parm);		/* calculate rate factor */
  y =  MVOLT * 0.05 * -(v-.015) * kcapnts->d1;
  x = exp (y) - 1;
  if (ncabs(x) > 1e-8) xalpha =  y/x;
  else                 xalpha =  1.0;
  return qkca * xalpha;
}

/*----------------------------------------*/

double betnkca6(chan *cpnt)

/* Like betnca, but with only voltage sensitivity, no Ca sensitivity */

/*  Use for rate func for d (to closed state). */

{
     double v;
     kcachan *kcapnts;
     cacomp *capnt;
     double qkca;

//#define DF 0.4
#define DF 1.0

  if (!cpnt) return 0.0;
  v = cpnt->comp1->v - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsm;
  kcapnts = (kcachan*)cpnt;
  qkca = qrate(cpnt->chtyp->parm);		/* calculate rate factor */
  return qkca * exp (-(v-.015) * MVOLT * 0.05 * kcapnts->d2);
}

/*----------------------------------------*/

chantype *makkca0(void)

/* Kca channel, type 0 */
/* by default set to be Ca sensitive only, no V sensitivity */

/*  SK type channel, "apamin-sensitive" */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 2;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,0,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->hh = 1;				/* HH chan - uses "m" */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->unitary = dkcasu;		/* 14.2 pS @22 deg makes 22 ps @ 5 deg */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca1;	/* user function */
   parm[0].funcname = (char *)"calckca1"; /* default rate func */
   parm[0].nfval    = 0; 
   parm[0].nival    = 2;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqkca;               /* Q10 for rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   mak2state(spnt,rchanf,rchanr,1000);	/* sets noise properties only */
					/*  (but adds slight delay, too) */

					/* Rate set to approximate openings */
					/*  of sKCa channel in Wang, Robinson */
					/*  and Chalupa, (1998). */
   return ch;
}

/*----------------------------------------*/

chantype *makkca1(void)

/* Sequential-state sKCa channel, like K type 0 */
/* by default set to be Ca sensitive only, no V sensitivity */

/*  SK type channel, "apamin-sensitive" */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define ca1f1 500.0
#define ca1r1 500.0
#define ca1f2 400.0
#define ca1r2 100.0

   nstate = 3;                          /* 3 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,1,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->unitary = dkcasu;		/* 14.2 pS @22 makes 22 ps @ 35 deg */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca1;        /* dummy rate func */
   parm[0].funcname = (char *)"calckca1";       /* dummy user func */
   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqkca;               /* Q10 for rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alnca;
   spnt[0].ratemul[0] = ca1f1;
   spnt[0].rateo  [0] = 0;      /* arate */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rt;
   spnt[1].ratemul[0] = ca1f2;
   spnt[1].rateo  [0] = 3;	/* drate */
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = betnca;
   spnt[1].ratemul[1] = ca1r1;
   spnt[1].rateo  [1] = 1;	/* brate */

   spnt[2].numtrans   = 1;
   spnt[2].cond       = 1;
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = rt;             /* function of ca conc, time */
   spnt[2].ratemul[0] = ca1r2;
   spnt[2].rateo  [0] = 3;       /* drate */

/*  The old way: */ 
/*   mak2state(spnt,alnca,betnca,10);     /* Sets both macro and micro */
					/*  channel properties */
   return ch;
}

/*--------------------------------------------*/

chantype *makkca2(void)

/* bKCa channel, type 2 */ 
/*  set to be Ca- and V-sensitive. */

/*   Like type 0, except BK type channel, voltage-sensitive */

/* Q10 for BKCa cond =~ 1.4, Barret, Magleby and Pallotta, 1982,
                             J. Physiol. 331: 211-230
*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 2;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,2,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->hh = 1;				/* HH chan - uses "m" */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->unitary = dkcabu;		/* 74 pS @22 makes 115 ps @ 22 deg */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca2;	/* user function */
   parm[0].funcname = (char *)"calckca2";	/* default rate func */
   parm[0].nfval    = 0; 
   parm[0].nival    = 2;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqkca;               /* Q10 for rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   mak2state(spnt,rchanf,rchanr,1000);	/* sets noise properties only */
					/*  (but adds slight delay, too) */
   return ch;
}

/*----------------------------------------*/

chantype *makkca3(void)

/* Sequential-state bKCa channel, like KCa type 2 */
/* by default set to be Ca- and V-sensitive. */

/*  BK type channel */

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define ca3f1 KCA_RATE
#define ca3r1 KCA_RATE

   nstate = 2;                          /* 3 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,3,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dkcabu;		/* 74 pS @22 makes 115 ps @ 35 deg */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca2;        /* dummy rate func */
   parm[0].funcname = (char *)"calckca2";       /* dummy user func */
   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqkca;               /* Q10 for rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alnca;
   spnt[0].ratemul[0] = ca3f1;
   spnt[0].rateo  [0] = 0;      /* arate */

   spnt[1].numtrans   = 1;
   spnt[1].cond       = 1;
   spnt[1].trans  [0] = 0;
   spnt[1].trate  [0] = betnca;
   spnt[1].ratemul[0] = ca3r1;
   spnt[1].rateo  [0] = 1;	/* brate */

   return ch;
}

/*--------------------------------------------*/

chantype *makkca4(void)

/* Sequential-state sKCa channel, 

   from:  Hirschberg et al., (1998) J. Gen. Physiol. 111: 565-581

  in which a 6-state markov diagram is provided to describe the
cloned rSK2 channel. Two sets of rate constants were provided, one for
normal high-probability gating, and another for simulating
"low-probability" gating behavior (not used here).  This channel
has no gating voltage sensitivity.

        f1      f2        f3
    0  <->  1  <->   2   <->    3
        r1      r2        r3

                       ^           ^
                 f4  | r4   f5  |  r5
                 \/         \/           

                     4          5

State 0 is unbound
      1 is bound with 1 Ca++
      2 is bound with 2 Ca++
      3 is bound with 3 Ca++
      4 is open, short duration
      5 is open, long duration

From Hirschberg et al (1998): Single channel currents were
recorded from Xenopus oocytes expressing the apamin-sensitive
clone rSK2. Channel activity was detectable in 0.2 uM [Ca]i and
was maximal above 2 uM [Ca]i.  Analysis of stationary currents
revealed two open times and three closed times, with only the
longest closed time being ca dependent, decreasing with
increasing [Ca]i concentrations.  In addition, elevated [Ca]i
concentrations resulted in a larger percentage of long openings
and short closures.  Membrane voltage did not have significant
effects on either open or closed times.  The open probability was
~0.6 in 1 uM free [Ca].

To set greater [Ca] sensitivity, increase "arate" by setting 
to a value less than 1, i.e. "taua = 0.1" will increase the [Ca] 
sensitivity by a factor of 10.

*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define ca4f1 200e6
#define ca4r1 80.
#define ca4f2 160e6
#define ca4r2 80.0
#define ca4f3 80e6
#define ca4r3 200.0
#define ca4f4 160.0
#define ca4r4 1000.0
#define ca4f5 1200.0
#define ca4r5 100.0

   nstate = 6;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,4,nstate,nparm,nq,dbasetc);  /* make chan state info */
   ch->unitary = dkcasu;		/* 14.2 pS @22 makes 22 ps @ 35 deg */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca2;        /* dummy rate func */
   parm[0].funcname = (char *)"calckca2";       /* dummy user func */
   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0] = dqkca;               /* Q10 for rca rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rca;
   spnt[0].ratemul[0] = ca4f1;
   spnt[0].rateo  [0] = 0;      /* arate */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rca;
   spnt[1].ratemul[0] = ca4f2;
   spnt[1].rateo  [0] = 0;	/* arate */
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ca4r1;
   spnt[1].rateo  [1] = 1;	/* brate */

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rca;             /* function of ca conc, time */
   spnt[2].ratemul[0] = ca4f3;
   spnt[2].rateo  [0] = 0;      /* arate */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;             /* function of time */
   spnt[2].ratemul[1] = ca4r2;
   spnt[2].rateo  [1] = 1;      /* brate */
   spnt[2].trans  [2] = 4;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = ca4f4;
   spnt[2].rateo  [2] = 2;      /* crate */
 
   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 5;
   spnt[3].trate  [0] = rt;             /* function of time */
   spnt[3].ratemul[0] = ca4f5;
   spnt[3].rateo  [0] = 2;      /* crate */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;             /* function of time */
   spnt[3].ratemul[1] = ca4r3;
   spnt[3].rateo  [1] = 1;      /* brate */
 
   spnt[4].numtrans   = 1;
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 2;
   spnt[4].trate  [0] = rt;             /* function of time */
   spnt[4].ratemul[0] = ca4r4;
   spnt[4].rateo  [0] = 3;      /* drate */
 
   spnt[5].numtrans   = 1;
   spnt[5].cond       = 1.0;
   spnt[5].trans  [0] = 3;
   spnt[5].trate  [0] = rt;             /* function of time */
   spnt[5].ratemul[0] = ca4r5;
   spnt[5].rateo  [0] = 3;      /* drate */

   return ch;
}

/*--------------------------------------------*/

chantype *makkca5(void)

/* Sequential-state sKCa channel, 

   from:  Sah and Clements (1999) J. Neurosci. 19: 3657-3664 

  in which a 6-state markov diagram is provided.  This is the
channel that underlies the "slow afterhyperpolarization" in
hippocampal pyramidal neurons. Its identity is unknown at
present, but it is apamin-insensitive and has no gating voltage
sensitivity.

        4rb     3rb       2rb       rb      ro 
    0  <->  1  <->   2   <->    3  <->  4  <->  5
        ru      2ru       3ru       4ru     rc

State 0 is unbound
      1 is bound with 1 Ca++
      2 is bound with 2 Ca++
      3 is bound with 3 Ca++
      4 is bound with 4 Ca++
      5 is open, flickering

From Sah and Clements, (1999): "The rates were chosen based on
the following three assumptions concerning the Ca-activated K
channels: 1) the steady-state dose-response curve has a steep
activation above the resting [Ca]i of 50 nM and has an EC50 of
150 nM, so that it is efficiently activated by small increases in
[Ca]i;  2) when [Ca]i falls rapidly, the decay of sIAHP is
limited by the channel closing and Ca dissociation rates to give
a time constant of ~1.5 sec; and 3) the peak open probability of
the channel is ~0.6, and its mean open time is 2.5 msec based on
estimates from noise analysis of sIAHP (Sah and Isaacson, 1995).

*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define ca5rb 10e6 
#define ca5ru 0.5
#define ca5ro 600.0
#define ca5rc 400.0

   nstate = 6;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,5,nstate,nparm,nq,dbasetc);  /* make chan state info */
   ch->unitary = dkcasu;		/* 14.2 pS @22 makes 22 ps @ 35 deg */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca2;	/* dummy rate func */
   parm[0].funcname = (char *)"calckca2";	/* dummy user func */
   parm[0].nfval    = 0;		/* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* number of implicit vals (m1, m2) */
   parm[0].dq[0]    = dqkca;            /* Q10 for rca rate function */
   parm[0].voff     = dkoffsn;		/* voltage offset from user */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rca;
   spnt[0].ratemul[0] = ca5rb * 4.0;
   spnt[0].rateo  [0] = 0;      /* arate */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rca;
   spnt[1].ratemul[0] = ca5rb * 3.0;
   spnt[1].rateo  [0] = 0;	/* arate */
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ca5ru;
   spnt[1].rateo  [1] = 1;	/* brate */

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rca;
   spnt[2].ratemul[0] = ca5rb * 2.0;
   spnt[2].rateo  [0] = 0;	/* arate */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = ca5ru * 2.0;
   spnt[2].rateo  [1] = 1;	/* brate */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rca;
   spnt[3].ratemul[0] = ca5rb;
   spnt[3].rateo  [0] = 0;	/* arate */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = ca5ru * 3.0;
   spnt[3].rateo  [1] = 1;	/* brate */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = ca5ro;
   spnt[4].rateo  [0] = 3;	/* arate */
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = ca5ru * 4.0;
   spnt[4].rateo  [1] = 3;	/* drate */

   spnt[5].numtrans   = 1;
   spnt[5].cond       = 1.0;	/* the open state */
   spnt[5].trans  [0] = 4;
   spnt[5].trate  [0] = rt;             /* function of time */
   spnt[5].ratemul[0] = ca5rc;
   spnt[5].rateo  [0] = 3;      /* drate */

   return ch;
}

/*----------------------------------------*/

chantype *makkca6(void)

/* Sequential-state bKCa channel, 

   from:  Cox, DH, Cui, J. and Aldrich, RW (1997)
   J. Gen. Physiol. 110: 257-281. 

  in which a 10-state markov diagram is provided.

        4a      3a        2a        a
    0  <->   1  <->   2   <->  3   <->  4 
        b       2b        3b        4b

      ^        ^        ^        ^        ^
 go | do  g1 | d1  g2 | d2  g3 | d3  g4 | d4
 \/       \/       \/       \/       \/  

       4af      3af      2af       af 
    5  <->   6  <->   7  <->   8  <->   9
       b/f      2b/f     3b/f     4b/f


State 0 is unbound
      1 is bound with 1 Ca++
      2 is bound with 2 Ca++
      3 is bound with 3 Ca++
      4 is bound with 4 Ca++
      5 is open, short duration
      6 is open, short duration
      7 is open, short duration
      8 is open, long duration
      9 is open, long duration

Original Cox rate constants (23 deg C):

Steady state parameters

L(0) = 1647
Kc   = 11.0 uM
Ko   =  1.1 uM

- - - - - - - - - - - - -- 
Kinetic parameters

Vertical rate constants

      sec-1

a  = 1e9 / [Ca] * unbound sites
bc = 1e9 * Kc * bound sites
bo = 1e9 * Ko * bound sites
Q  = 1.40e
 
c05 = 2.39
c16 = 7.0
c27 = 40
c38 = 295
c49 = 557
charge C->O = .73e

c50 = 3936
c61 = 1152
c72 = 659
c83 = 486
c94 = 92
charge O->C = -.67e

Horizontal rate constants

Ca++ on rates = 1e9/M/S for each binding site in both open and
closed channel conformations and were adjusted appropriately for
the number of available sites.

Ca++ off rates equaled (1e9 Kc) for each binding site in the
closed conformation and (1e9 Ko) for each binding site in the
open conformation, and were adjusted appropriately for the number
of bound Ca++.

*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define kcaKc  11.0e-6 
#define kcaKo  1.1e-6 

#define kca6f  1e9		/* rate of Ca per binding site */
#define kca6a  kca6f*1.0	/* binding rate of Ca per unbound site */
#define kca6bc (kca6f * kcaKc)	/* unbinding of Ca per bound closed site */
#define kca6bo (kca6f * kcaKo)	/* unbinding of Ca per bound open site */

#define kca6c05 2.39
#define kca6c16 7.0
#define kca6c27 40
#define kca6c38 295
#define kca6c49 557

#define kca6c50 3936
#define kca6c61 1152
#define kca6c72 659
#define kca6c83 486
#define kca6c94 92

#define kf 1.0
#define kr 1.0
#define kvf 1.0
#define kvr 1.0

#define BASETKCA6 23.0

   nstate = 10;                         /* 10 Markov states */
   nparm = 2;                           /* make 2 sets of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,6,nstate,nparm,nq,dbasetc);  /* make chan state info */
   ch->unitary = dkcabu;		/* 74 pS  @22 makes 115 ps @ 35 deg */
   					/* 44 pS @6.3 makes 115 ps @ 35 deg */
   					/* 68 pS @6.3 makes 115 ps @ 22 deg */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
					/* second param, for Ca binding */
   parm[0].chancalc =  calckca2;        /* dummy rate func */
   parm[0].funcname = (char *)"calckca2";       /* dummy user func */
   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0]    = dqkca;            /* Q10 for rca rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

					/* second param, for Ca binding */
   parm[1].chancalc =  calckca2;        /* dummy rate func */
   parm[1].funcname = (char *)"calckca2";       /* dummy user func */
   parm[1].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[1].dq[0] = dqkca;               /* Q10 for rca2 rate function */
   parm[1].voff = dkoffsn;              /* voltage offset from user */

   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rca2;
   spnt[0].ratemul[0] = kca6a*4.0*kf;	/* rate func multiplier */
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 5;
   spnt[0].trate  [1] = alnkca6;
   spnt[0].ratemul[1] = kca6c05*kvf;
   spnt[0].rateo  [1] = 2;      /* crate */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rca2;
   spnt[1].ratemul[0] = kca6a*3.0*kf;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = kca6bc;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 6;
   spnt[1].trate  [2] = alnkca6;
   spnt[1].ratemul[2] = kca6c16*kvf;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rca2;               /* function of ca */
   spnt[2].ratemul[0] = kca6a*2.0*kf;
   spnt[2].rateo  [0] = 0;      /* arate */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = kca6bc*2.0;
   spnt[2].rateo  [1] = 1;      /* brate */
   spnt[2].trans  [2] = 7;
   spnt[2].trate  [2] = alnkca6;	    /* function of v */
   spnt[2].ratemul[2] = kca6c27*kvf;
   spnt[2].rateo  [2] = 2;      /* crate */
 
   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = rca2;                /* function of ca */
   spnt[3].ratemul[0] = kca6a*kf;
   spnt[3].rateo  [0] = 0;      /* arate */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = kca6bc*3.0;
   spnt[3].rateo  [1] = 1;      /* brate */
   spnt[3].trans  [2] = 8;
   spnt[3].trate  [2] = alnkca6;             /* function of v */
   spnt[3].ratemul[2] = kca6c38*kvf;
   spnt[3].rateo  [2] = 2;      /* crate */
 
   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = kca6bc*4.0;
   spnt[4].rateo  [0] = 1;      /* brate */
   spnt[4].trans  [1] = 9;
   spnt[4].trate  [1] = alnkca6;             /* function of v */
   spnt[4].ratemul[1] = kca6c49*kvf;
   spnt[4].rateo  [1] = 2;      /* crate */
 
   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = rca2;
   spnt[5].ratemul[0] = kca6a*4.0*kf;
   spnt[5].rateo  [0] = 0;      /* arate */
   spnt[5].trans  [1] = 0;
   spnt[5].trate  [1] = betnkca6;
   spnt[5].ratemul[1] = kca6c50*kvr;
   spnt[5].rateo  [1] = 3;      /* drate */

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 1.0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = rca2;
   spnt[6].ratemul[0] = kca6a*3.0*kf;
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = kca6bo;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 1;
   spnt[6].trate  [2] = betnkca6;
   spnt[6].ratemul[2] = kca6c61*kvr;
   spnt[6].rateo  [2] = 3;

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 1.0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = rca2;             /* function of ca */
   spnt[7].ratemul[0] = kca6a*2.0*kf;
   spnt[7].rateo  [0] = 0;      /* arate */
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = rt;
   spnt[7].ratemul[1] = kca6bo*2.0;
   spnt[7].rateo  [1] = 1;      /* brate */
   spnt[7].trans  [2] = 2;
   spnt[7].trate  [2] = betnkca6;	    /* function of v */
   spnt[7].ratemul[2] = kca6c72*kvr;
   spnt[7].rateo  [2] = 3;      /* drate */
 
   spnt[8].numtrans   = 3;
   spnt[8].cond       = 1.0;
   spnt[8].trans  [0] = 9;
   spnt[8].trate  [0] = rca2;                /* function of ca */
   spnt[8].ratemul[0] = kca6a*kf;
   spnt[8].rateo  [0] = 0;      /* arate */
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = rt;
   spnt[8].ratemul[1] = kca6bo*3.0;
   spnt[8].rateo  [1] = 1;      /* brate */
   spnt[8].trans  [2] = 3;
   spnt[8].trate  [2] = betnkca6;             /* function of v */
   spnt[8].ratemul[2] = kca6c83*kvr;
   spnt[8].rateo  [2] = 3;      /* drate */
 
   spnt[9].numtrans   = 2;
   spnt[9].cond       = 1.0;
   spnt[9].trans  [0] = 8;
   spnt[9].trate  [0] = rt;
   spnt[9].ratemul[0] = kca6bo*4.0;
   spnt[9].rateo  [0] = 1;      /* brate */
   spnt[9].trans  [1] = 4;
   spnt[9].trate  [1] = betnkca6;             /* function of v */
   spnt[9].ratemul[1] = kca6c94*kvr;
   spnt[9].rateo  [1] = 3;      /* drate */
 
   return ch;
}

#undef kf 
#undef kr 

/*--------------------------------------------*/

chantype *makkca8(void)

/* Sequential-state bKCa channel, 

   from:  Horrigan, FT, Cui, J, and Aldrich, RW (1999)
   J. Gen. Physiol. 114: 277-304. 

  in which a 10-state markov diagram is provided.

        4a      3a        2a        a
    0  <->   1  <->   2   <->  3   <->  4 
        b       2b        3b        4b

      ^        ^        ^        ^        ^
 go | do  g1 | d1  g2 | d2  g3 | d3  g4 | d4
 \/       \/       \/       \/       \/  

       4af      3af      2af       af 
    5  <->   6  <->   7  <->   8  <->   9
       b/f      2b/f     3b/f     4b/f


State 0 is unbound
      1 is bound with 1 Ca++
      2 is bound with 2 Ca++
      3 is bound with 3 Ca++
      4 is bound with 4 Ca++
      5 is open, short duration
      6 is open, short duration
      7 is open, short duration
      8 is open, long duration
      9 is open, long duration

Original Horrigan rate constants (20 deg C):

     sec-1

a  = 1276  / [Ca] (uM)
b  = 35370 / [Ca] (um)
d0 =  0.0074
d1 =  0.126
d2 =  2.14
d3 = 25.7
d4 = 49.3
g0 = 3700
g1 = 3700
g2 = 3700
g3 = 2612
g4 =  295

za =  0.275
zb = -0.275
zd =  0.262
zg = -0.138
f  = sqrt(17) = 4.123
d0/g0 = 2e-6

We modify the Horrigan model by adding Ca sensitivity
from states 0 -> 4, taken from the 8-state model of McManus and
Magleby (1991).

The rate functions for voltage sensitivity are slightly modified
from the "alnca()" and "betnca()" functions above. The idea is
that the voltage sensitivity of the vertical transitions is less
than for the [Ca]-dependent transitions. 

This channel is untested. The states are set up, and the rate
functions might work OK, but the rate multipliers are probably
not correct.  

The reason is that the rate constants from Horrigan et al. (1999)
are for voltage sensitivity, and the rate constants from McManus
and Magleby are for [Ca] sensitivity.  The Markov states from
both studies look similar but they are not the same and cannot
easily be incorporated into one model. The system is very complex
and the kinetics are not well understood (as of 4/2000).  

For example, Cui et al (1997) tried (as an exercise) to add
voltage sensitivity to the vertical transitions of McManus and
Magleby's model, and found that they could get it to match at
high and low [Ca] but at medium [Ca] it didn't match.

It appears that more states will be necessary to fully represent
the effect of both [Ca] and voltage.  It is already known that
[Ca] has a bigger effect on kinetics, and but to incorporate the
voltage sensitivity, one needs more than the 8 - 10 state model
that works for [Ca] alone.

*/

{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

#define kca8a  1276.0
#define kca8b  35370.0
#define kca8d0 0.0074
#define kca8d1 0.126
#define kca8d2 2.14
#define kca8d3 25.7
#define kca8d4 49.3
#define kca8g0 3700.0
#define kca8g1 3700.0
#define kca8g2 3700.0
#define kca8g3 2612.0
#define kca8g4 295.0
#define kca8f  4.123
#define kf    1.0
#define kr    1.0

   nstate = 10;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(KCa,8,nstate,nparm,nq,dbasetc);  /* make chan state info */
   ch->unitary = dkcabu;		/* 74 pS @22 makes 115 ps @ 35 deg */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak*.5;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak*.5;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].chancalc =  calckca2;        /* dummy rate func */
   parm[0].funcname = (char *)"calckca2";       /* dummy user func */
   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].dq[0]    = dqkca;            /* Q10 for rca rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = alnca;
   spnt[0].ratemul[0] = kca8a*4.0*kf;	/* rate func multiplier */
   spnt[0].rateo  [0] = 0;
   spnt[0].trans  [1] = 5;
   spnt[0].trate  [1] = alnkca6;
   spnt[0].ratemul[1] = kca8g0;
   spnt[0].rateo  [1] = 2;      /* crate */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = alnca;
   spnt[1].ratemul[0] = kca8a*3.0*kf;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = betnca;
   spnt[1].ratemul[1] = kca8b;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 6;
   spnt[1].trate  [2] = alnkca6;
   spnt[1].ratemul[2] = kca8g1;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = alnca;             /* function of ca, v */
   spnt[2].ratemul[0] = kca8a*2.0*kf;
   spnt[2].rateo  [0] = 0;      /* arate */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betnca;             /* function of v */
   spnt[2].ratemul[1] = kca8b*2.0;
   spnt[2].rateo  [1] = 1;      /* brate */
   spnt[2].trans  [2] = 7;
   spnt[2].trate  [2] = alnkca6;
   spnt[2].ratemul[2] = kca8g2;
   spnt[2].rateo  [2] = 2;      /* drate */
 
   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = alnca;             /* function of ca, v */
   spnt[3].ratemul[0] = kca8a*kf;
   spnt[3].rateo  [0] = 0;      /* arate */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betnca;             /* function of time */
   spnt[3].ratemul[1] = kca8b*3.0;
   spnt[3].rateo  [1] = 1;      /* brate */
   spnt[3].trans  [2] = 8;
   spnt[3].trate  [2] = alnkca6;             /* function of v */
   spnt[3].ratemul[2] = kca8g3;
   spnt[3].rateo  [2] = 2;      /* drate */
 
   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = betnca;             /* function of ca, v */
   spnt[4].ratemul[0] = kca8b*4.0;
   spnt[4].rateo  [0] = 1;      /* brate */
   spnt[4].trans  [1] = 9;
   spnt[4].trate  [1] = alnkca6;             /* function of v */
   spnt[4].ratemul[1] = kca8g4;
   spnt[4].rateo  [1] = 2;      /* drate */
 
   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = alnca;
   spnt[5].ratemul[0] = kca8a*kca8f*4.0*kf;
   spnt[5].rateo  [0] = 0;      /* arate */
   spnt[5].trans  [1] = 0;
   spnt[5].trate  [1] = betnkca6;
   spnt[5].ratemul[1] = kca8d0;
   spnt[5].rateo  [1] = 3;      /* drate */

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 1.0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = alnca;
   spnt[6].ratemul[0] = kca8a*kca8f*3.0*kf;
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = betnca;
   spnt[6].ratemul[1] = kca8b/kca8f;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 1;
   spnt[6].trate  [2] = betnkca6;
   spnt[6].ratemul[2] = kca8d1;
   spnt[6].rateo  [2] = 3;

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 1.0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = alnca;             /* function of ca, v */
   spnt[7].ratemul[0] = kca8a*kca8f*2.0*kf;
   spnt[7].rateo  [0] = 0;      /* offsm */
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = betnca;             /* function of v */
   spnt[7].ratemul[1] = kca8b/kca8f/2.0;
   spnt[7].rateo  [1] = 1;      /* offsm */
   spnt[7].trans  [2] = 2;
   spnt[7].trate  [2] = betnkca6;
   spnt[7].ratemul[2] = kca8d2;
   spnt[7].rateo  [2] = 3;      /* drate */
 
   spnt[8].numtrans   = 3;
   spnt[8].cond       = 1.0;
   spnt[8].trans  [0] = 9;
   spnt[8].trate  [0] = alnca;             /* function of ca, v */
   spnt[8].ratemul[0] = kca8a*kca8f;
   spnt[8].rateo  [0] = 0;      /* offsm */
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = betnca;             /* function of v */
   spnt[8].ratemul[1] = kca8b/kca8f*3.0;
   spnt[8].rateo  [1] = 1;      /* offsm */
   spnt[8].trans  [2] = 3;
   spnt[8].trate  [2] = betnkca6;             /* function of v */
   spnt[8].ratemul[2] = kca8d3;
   spnt[8].rateo  [2] = 3;      /* drate */
 
   spnt[9].numtrans   = 2;
   spnt[9].cond       = 1.0;
   spnt[9].trans  [0] = 8;
   spnt[9].trate  [0] = betnca;             /* function of v */
   spnt[9].ratemul[0] = kca8b/kca8f*4.0;
   spnt[9].rateo  [0] = 1;      /* offsm */
   spnt[9].trans  [1] = 4;
   spnt[9].trate  [1] = betnkca6;             /* function of v */
   spnt[9].ratemul[1] = kca8d4;
   spnt[9].rateo  [1] = 3;      /* drate */
 
   return ch;
}

/*--------------------------------------------*/
