/* Segment ncomp in Program nc */

/* Simulates neuronal circuits */
/*  by numerically integrating difference equations */
/*  using iterative relaxation technique */
/*  to implement the "implicit" integration method */ 

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncelem.h" 
#include "ncomp.h"
#include "control.h"

#define DEBUG
#define NOSKIP 128

#define SQRT10 3.16227766
#define SSQRT10 1.77827942

extern comp *compnt,*compend;

double exp();
double sqrt();

extern double alpham,betam;
extern double alphah,betah;
extern double alphan,betan;

static double maxerr,oldmaxerr;
static int ncomp,tcomp;

/*------------------------------------*/

double dochan(chpnt)
  chan *chpnt;
{
   int i,nstate,nsites;
   static stconc *conc;
   double p,bnldev(),dochan2();
   static double conduct;

  nstate = chpnt->chtyp->numstate;
  conc = chpnt->conc;
  for (i=0; i<nstate; i++,conc++) {  		/* copy old conc estimate */ 
    conc->cval = conc->cest;	   
  }
  conduct=dochan2(chpnt);
  if (nsites=chpnt->nchan) {			  /* */
    p = conduct / (nsites * chpnt->cq);           /* find open prob */
    conduct = bnldev(p,nsites) * chpnt->cq;       /* conduct w/noise */
  } 

 /* fprintf (stderr,"typ %d c %10.4g cond %10.4g\n",
		chpnt->ctype,conduct,chpnt->conduct); /* */
 return conduct;
}
/*------------------------------------*/

double dochan2(chpnt)
  chan *chpnt;
{
   int i,j,nstate,numtr;
   static chanstate *stpnt;
   static stconc *conc;
   static double (**f_rate)(), trrate, conduct;
   static float *trmul;
   static char *ratev;
   double trv;			/* voltage value for channel transition */
   char *trans;
   double cerr,cmaxerr,critc;
   static stconc oconc[NUMSTATE];
   static stconc concold[NUMSTATE];

  critc = 1e-8;
  nstate = chpnt->chtyp->numstate;
  conc = chpnt->conc;
  for (i=0; i<nstate; i++,conc++) {  /* zero deltas, set conc val = prev est */
    concold[i].dcon = conc->dcon;
    conc->dcon = 0;
    oconc[i].dcon = 0;
  }	
	/* calculate delta concentrations from rate functions: */

  for (cmaxerr=1.0; cmaxerr>critc; ) {
    cmaxerr = 0.0;
    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++) {  /* find change in conc */
       trans = stpnt->trans;	     /* get parameters */
       f_rate  = stpnt->trate;
       trmul = stpnt->ratemul;
       numtr = stpnt->numtrans;
       ratev = stpnt->ratev;
       for (j=0; j<numtr; j++,trans++,f_rate++,trmul++,ratev++) {
         trv = (*ratev ? chpnt->comp2->v : chpnt->comp1->v);
         trrate = conc->cest * (*f_rate)(trv) * *trmul;	      /* per timinc */
	 if (trrate > conc->cest) trrate = conc->cest;
	 if (trrate + chpnt->conc[*trans].cest > 1.0)
		 trrate = 1.0 - chpnt->conc[*trans].cest;
    /*     trrate = trrate * 2.0 / (trrate + 2.0); */
 /*if (chpnt->ctype==NA)
  fprintf (stderr,"i %d j %d trrate %10.4g cest %10.4g f_rate %10.4g mul %g\n",
	i,j,trrate,conc->cest,(*f_rate)(trv),*trmul); /* */
         chpnt->conc[*trans].dcon += trrate;
         conc->dcon -= trrate;
       }
    }
	/* after calculating all deltas, find new conc estimates: */

    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++) {
/*       if      (conc->dcon< -conc->cest) conc->dcon = -conc->cest;
       else if (conc->dcon+conc->cval>1.0) conc->dcon = 1.0 - conc->cval; /* */
       conc->cest = conc->cval +
		 0.3*conc->dcon + 0.3*oconc[i].dcon+0.4*concold[i].dcon;
/*       if      (conc->cest<0.0) conc->cest = 0.0;
       else if (conc->cest>1.0) conc->cest = 1.0;   /* */
       cerr = conc->dcon - oconc[i].dcon; 
/* fprintf (stderr,"typ %d conc %10.4g dconc %10.4g err %10.4g cond %g  i %d\n",
		chpnt->ctype,conc->cest,conc->dcon,cerr,stpnt->cond,i); /* */
       if (cerr < 0.0) cerr = -cerr; 	/* absolute value */
       if (cerr > cmaxerr) cmaxerr = cerr;
       oconc[i].dcon = conc->dcon;       /* save old value of change */
       conc->dcon = 0;			/* erase delta for next iteration */
    }
/* fprintf (stderr,"cmax %10.4g\n",cmaxerr);   /* */
  }    /* for (cmaxerr;;) */

  conduct = 0.0;
  conc = chpnt->conc;
  stpnt = chpnt->chtyp->state;
  for (i=0; i<nstate; i++,stpnt++,conc++) { /* integrate and total cond */
    conduct += conc->cest * stpnt->cond;
/*    if (chpnt->ctype==NA) fprintf (stderr,"%10.4g ",conc->cest); /* */
  }
/* if (chpnt->ctype==NA) 
	 fprintf (stderr,"\n"); /* */

  return (conduct);
}

/*------------------------------------*/

runcomp()

/* Advance the network of compartments one step in time,
using a second-order implicit numerical integration.
Iterate the approximation until the error is
reduced below criterion.
*/

{
  static int left,niter,titer,less,altstep=0;
  static double nodcur;                    /* total node current */
  static double nodcr;                     /* total unvarying node curr */
  static comp *pnt,*ocpnt;
  static conn *conpnt;
  static double delcrit,tdelcrit;
  static double relmult,critx,critxx,critk,crita;
  static double tcond;			   /* total conductance in comp */
  static double tcondn;			   /* tot cond not varying w/vest */
  static double conduct;		   /* conductance of channel */
  static double vest;		   	   /* conductance of channel */
  static conlst *lpnt;
  static recpar *cpt;                      /* receptor constants */
  static hhchan *hhpnt;
  static chan *chpnt;
  static double m,n,h,dm,dn,dh;

#ifdef DEBUG
  if (debug & 1) fprintf (stderr,"runcomp\n");
#endif

 
 if (timinc > 1.0) runsyn(0);			/* if static, calc synapses */


 if (euler) implicit = 0;

 for (pnt=compnt; pnt; pnt=pnt->cnext)                  /* First estimate */
  {                                                     /* for all comps */
   nodcur = nodcr = 0.0;
   tcond = tcondn = 0.0;
   for (lpnt=pnt->clst; lpnt; lpnt=lpnt->lnext) 
    {                                           /* check all connections */
     conpnt = lpnt->conpnt;
     if (! conpnt) break;
     if (conpnt->comp1 == pnt) ocpnt = conpnt->comp2; /* get other compartment*/
     else                      ocpnt = conpnt->comp1;
     switch (conpnt->ctype)
      {                         /* add connection conductances which change */

       case GJ:                                 /* may vary gj in future */
                nodcur += conpnt->conduct * ocpnt->v;   /* current thru g j */
                tcond += conpnt->conduct;              /* total cond */
               break;
       case SYNAPSE:					 /* cur thru syn chan */
		if (pnt==((synap *)conpnt)->comp2) {
                  nodcr += ((synap *)conpnt)->conduct *((synap *)conpnt)->vrev;
                  tcondn += ((synap *)conpnt)->conduct;   /* total cond */
  	/* fprintf (stderr,"vrev %10.3g conduct %10.3g\n",
		((synap *)conpnt)->vrev,((synap *)conpnt)->conduct); /* */
		}
               break;

      case ROD:
      case CONE:                                /* find receptor instantan */
                                                /* conductance */
          cpt = ((recep *)conpnt)->chtyp;
          nodcr += ((recep *)conpnt)->conduct * cpt->vrev;/* cur thru rec chan*/
          tcondn  += ((recep *)conpnt)->conduct;               /* total cond */
/*   fprintf (stderr,"vrev %10.3g conduct %10.3g\n",
		cpt->vrev,((recep *)conpnt)->conduct); /* */
         break;

      case AXIALRES:
      case RESISTOR:
                nodcur += conpnt->conduct * ocpnt->v;  /* cur to next comp */
                break;              /* conductance is collected in pnt->tcond */
      case BATT:
     		if (ocpnt == conpnt->comp1)
                  nodcur += (ocpnt->v - conpnt->conduct) * BATTCOND;
		else 
                  nodcur += (ocpnt->v + conpnt->conduct) * BATTCOND;
                break;		/* cur to next comp through batt */

      case BUF:	if (pnt == ((dbuf *)conpnt)->comp2) {
		   if (((dbuf *)conpnt)->delay) {
		     pnt->extv = *((dbuf *)conpnt)->delpnt * DELCONV;
		   }
		   else pnt->extv = ((dbuf *)conpnt)->comp1->v;
		  }
		break;

      case CAP: nodcr += (pnt->oldv - ocpnt->oldv) * conpnt->conduct/timinc;
                break;			/* -(v3 - v2) * c / t */

      case NA:  hhpnt = (hhchan *)conpnt;

		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

      		narate(pnt->v);			/* set Na rate chtyp */ 
		m = hhpnt->mest;		/* use previous best est */
  		h = hhpnt->hest;
		hhpnt->m = m;
		hhpnt->h = h;
  		dm = alpham - (alpham + betam) * m;	/* per timinc */
  		dh = alphah - (alphah + betah) * h;
     /* fprintf (stderr,"Na: alpham %10.4g betam %g v %10.4g\n",
					alpham,betam,pnt->v); /* */
  		hhpnt->dm = dm;
  		hhpnt->dh = dh;
  		m += dm; 
  		h += dh; 
		hhpnt->mest = m;
		hhpnt->hest = h;
  		conduct = m*m*m*h * hhpnt->maxcond;
    /* fprintf (stderr,"Na: m %10.4g h %10.4g cond %10.4g max %10.4g\n",
	m,h,conduct,hhpnt->maxcond); /* */
                nodcur += conduct * hhpnt->vrev;
                tcond += conduct;    			/* total cond */
		break;

		case 1:				/* Sequential state */
		case 2:

	 	chpnt = (chan *)conpnt;
  		conduct = dochan((chan *)conpnt) * chpnt->maxcond;
                nodcur += conduct * chpnt->vrev;
                tcond += conduct;    			/* total cond */
		break;

		}	/* switch (stype) */
                break;

      case K:   
	 	hhpnt = (hhchan *)conpnt;

		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

  		krate(pnt->v);			/* set K rate chtyp */
		n = hhpnt->mest;
		hhpnt->m = n;
  		dn = alphan - (alphan + betan) * n;
  		hhpnt->dm = dn;
  		n += dn; 			/* per timinc */
		hhpnt->mest = n;
  		conduct = n*n*n*n * hhpnt->maxcond;
   /*  fprintf (stderr,"K: n %10.4g cond %10.4g max %10.4g\n",
	n,conduct,hhpnt->maxcond); /* */
                nodcur += conduct * hhpnt->vrev;
                tcond += conduct;    			/* total cond */
                break;

		case 1:				/* Sequential state */
		case 2:
	 	chpnt = (chan *)conpnt;
  		conduct = dochan((chan *)conpnt) * chpnt->maxcond;
                nodcur += conduct * chpnt->vrev;
                tcond += conduct;    		/* total cond */
		break;

		}	/* switch (stype) */
		break;

      default:  break;

      } /* switch */
    }      /* for (lpnt=; ;) */

   tcond += tcondn + pnt->tcond;                 /* total conductance */
   nodcur -= pnt->v * tcond;                     /* current due to comp volts */
   nodcr  += pnt->vrev * pnt->rm;                /* membrane conductance */

/*  fprintf (stderr,"pnt->vrev %10.3g\n",pnt->vrev); /* */
/*  fprintf (stderr,"nodcr %10.3g nodcur %10.3g pntcond %g tcond %g\n",
			nodcr,nodcur,pnt->tcond,tcond); /* */

   if (pnt->miscfl & IEXT)                    /* comp has a current src */
      nodcr += pnt->exti;                     /* external current */
   nodcur += nodcr;			      /* total comp current */

 if (euler)	{     			      /* forward-Euler */
     vest = pnt->v + nodcur * pnt->k;	      /* simple first-order est */
     if (pnt->miscfl & VEXT) {                /* comp is a battery */
       pnt->extvi = (vest - pnt->extv) / (pnt->k);
       vest = pnt->extv;		      /*  battery voltage  */
     }
     else {
       pnt->extvi = 0.0;
     }
     pnt->vest = vest;				/* just save the estimate */
 }						/* (not quite) end of Euler */
						/*   (see below) */

 else {						/* Implicit modes */

   pnt->tcondn = tcondn;			/* cond not varying w/vest */
   pnt->implf = 1. / (1. + tcond * pnt->k);   /* implicit factor   */

   if (implicit) {			/* Backward Euler, fully implicit */
     pnt->nodc = nodcr;				/* curr for estimate */
     nodcur *=  pnt->k;                      	/* first-ord only */
     pnt->vest = pnt->v + nodcur;             /* simple first-order est */
     if (pnt->miscfl & VEXT) {                   /* comp is a battery */
          pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k);
          pnt->vest = pnt->extv;
       }
     else {
       pnt->extvi = 0.0;
     }
     pnt->verr = 1.;
   }						/* end of implicit */

   else	{     /* Crank-Nicolson */
     pnt->nodc = nodcur + nodcr;	      /* curr for estimate */
     nodcur *=  pnt->k * 2.;                  /* first-ord only */
     pnt->vest = pnt->v + nodcur;             /* simple first-order est */
     if (pnt->miscfl & VEXT) {                /* comp is a battery */
       pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k*2.);
       pnt->vest = pnt->extv;		      /*  battery voltage  */
     }
     else {
       pnt->extvi = 0.0;
     }
     pnt->verr = 1.;
   }    /* crank-nicholson */
  }  /* CN or implicit */
 }   /* for (pnt=;;) */

 if (euler) {					/* after everything else, */
   for (pnt=compnt; pnt; pnt=pnt->cnext) {
	pnt->v = pnt->vest;			/* save as new voltage */
   }
   return;
 }


 if (!altstep) {
     crita = 1e-1;
     crita = crit;
 }
 else {
     crita = crit;
 }

 altstep = !altstep;
 relmult = relincr * timinc / 1e-4;
 oldmaxerr=maxerr = 1.0;                         /* set maxerr so loop runs */ 
 ncomp = tcomp = 0; 
 delcrit =   10.0;		/* SQRT10; /* */
 left = 1;
for (titer=0,critx=1.000001e-2;
		 critx >= crita;
		 critx = maxerr,
		  critx /= delcrit, titer+=niter, left++) /* */
 {
 oldmaxerr = 1.0;
 critxx = critx;
 for (niter=0; maxerr>critx; niter++)	/* Iterate until convergence */
  {
   if (timinc > 1.0) runsyn(0);		/* if static, calc synapses */

/* fprintf (stderr,"maxerr %g critx %g\n",maxerr,critx); /* */


   maxerr = 0.0;
   if (left&1) 
     for (ncomp=0,pnt=compend; pnt; pnt=pnt->clast)     /* left estimate */
      {
#ifdef DEBUG
        if (!(debug&NOSKIP))
#endif
          if (pnt->verr<critxx) continue; /* */
        docomp(pnt);
      }
   else 
     for (ncomp=0,pnt=compnt; pnt; pnt=pnt->cnext)      /* right estimate */
      {
#ifdef DEBUG
        if (!(debug&NOSKIP))
#endif
          if (pnt->verr<critxx) continue; /* */
        docomp(pnt);
      }           /* for (ncomp=0,pnt=compnt;) */

#ifdef DEBUG
  if (debug & 4 && debugz & 1)
	 fprintf (stderr,"ncomp  %4d maxerr %8.3g critxx %8.3g\n",
					ncomp,maxerr,critxx);
#endif

   /* Reduce the criterion if error increases, and also at least
       once every 10 iterations.  This causes compartments
       outside the region of iteration (i.e. skipped over)
       to be recomputed and reduces the total number of iterations
       necessary for complete convergence.
   */
       
   if (maxerr>oldmaxerr || ((niter>0) && (niter%10==0))) {
      critxx /= sqrt(delcrit);
#ifdef DEBUG
  if (debug & 4 && debugz & 1)
	 fprintf (stderr,"maxerr %8.3g, reducing critxx to %8.3g\n",
					maxerr,critxx);
#endif
   }
   oldmaxerr = maxerr; 
   if (critxx < 1e-30) {
      fprintf (stderr,"Ncomp: panic, no convergence, continuing...");
      return;
   }
  
 }            /* for (niter=0; maxerr>critx; niter++)  */


 if (niter<2) {
 /*     delcrit *= delcrit;
      if (delcrit>10.00001) delcrit=10.0;
#ifdef DEBUG
      if (debug & 4 && debugz & 1)
	 fprintf (stderr,"Increasing delcrit to %8.3g\n",delcrit);
#endif
*/
    }

 if (niter>10) {
    tdelcrit = sqrt(delcrit);
    if (tdelcrit>1.002) {
      delcrit=tdelcrit;			/* */
#ifdef DEBUG
      if (debug & 4 && debugz & 1)
	 fprintf (stderr,"Reducing delcrit to %8.3g\n",delcrit);
#endif
    }
 }

/* Reset error after each level of critx */
/*  This has the effect of spreading the region of iteration,
    and causes compartments just outside to be recomputed.
*/

 critk = critx / (delcrit * sqrt(delcrit))  + crit;
 for (pnt=compnt; pnt; pnt=pnt->cnext)
  {
   pnt->verr += critk;
  }			/* */


#ifdef DEBUG
  if (debug & 2 && debugz & 1)
   fprintf (stderr,"niter %2d                    critx  %8.3g delcrit %8.3g\n",
				niter,critx,delcrit);
#endif

  }         /* for (critx ; ; ) */

#ifdef DEBUG
  if (debug & 2 && debugz & 1)
	fprintf (stderr,"*tcomp %4d titer %d\n",tcomp,titer);
#endif

/* Just save it and calculate relax for next time step. */

 for (pnt=compnt; pnt; pnt=pnt->cnext)
  {
   pnt->oldv = pnt->v = pnt->vest;		/* save it */
  }

/* This code runs only once per time step, */
/* so it adds negligible time to computataion */
/* if model is large and requires many iterations */
/* per time step. */

if (relincr && titer > 10)			/* if lots of iterations */
 for (pnt=compnt; pnt; pnt=pnt->cnext) {
   if (!(pnt->miscfl & VEXT)) {			/* if comp is not a battery */
      less = pnt->t - pnt->g;			/* then compute oscillations */
      m = ((pnt->g < less) ? less : pnt->g);	/* max of less, greater */
      if (m==pnt->t) pnt->relax += relmult;	/* adjust comp relax */
      else           pnt->relax -= (1.0 - ((double) m)/pnt->t) * 10.0 * relmult;
      if      (pnt->relax > 1.0) pnt->relax = 1.0;
      else if (pnt->relax < 0.0) pnt->relax = 0.0;
    }

#ifdef DEBUG
  if (debug & 8 && debugz & 1)
     fprintf (stderr,"c %-4d relax %-8.4g g %-4d t %-4d misc %d\n",
			   pnt->num,pnt->relax,pnt->g,pnt->t,pnt->miscfl); /* */
#endif
   pnt->g = pnt->t = 0;				/* reset greater and total */
  }
}

/*------------------------------------*/

docomp(pnt)
   comp *pnt;

{
  static double nodcur;                         /* total comp current */
  static double tcondv;                         /* total comp conductance */
  static double conduct;                        /* channel conductance */
  static double implf;				/* new implicit factor */
  static comp *ocpnt;
  static conn *conpnt;
  static double err;
  static double icomp;
  static conlst *lpnt;
  static hhchan *hhpnt;
  static chan *chpnt;
  static double m,n,h,dm,dn,dh;
  int i;

#ifdef DEBUG
     tcomp++;
     ncomp++;
#endif
     nodcur = pnt->nodc;                        /* start with old node cur */
     tcondv = 0.0;				/* non-varying conductance */
     for (lpnt=pnt->clst; lpnt; lpnt=lpnt->lnext) 
      {
       conpnt = lpnt->conpnt;                   /* check connection */
       if (! conpnt) break;
       if (conpnt->comp1 == pnt) ocpnt = conpnt->comp2; /* get other comp */
       else                      ocpnt = conpnt->comp1;
       switch (conpnt->ctype)
        {
         case AXIALRES:
         case RESISTOR:
         case GJ:
                nodcur += ocpnt->vest * conpnt->conduct;  /* curr thru resist */
                break;
         case BATT:
     		if (ocpnt == conpnt->comp1) 
                  nodcur += (ocpnt->vest - conpnt->conduct) * BATTCOND;
		else
                  nodcur += (ocpnt->vest + conpnt->conduct) * BATTCOND;
                break;		     /* curr thru batt small resistor */

         case BUF: if (pnt == ((dbuf *)conpnt)->comp2) {
		     if (((dbuf *)conpnt)->delay) {
		       pnt->extv = *((dbuf *)conpnt)->delpnt * DELCONV;
		     }
		     else pnt->extv = ((dbuf *)conpnt)->comp1->vest;
		   }
		break;
         case CAP:
		nodcur += ocpnt->vest * conpnt->conduct / timinc;	
                break; 
         case NA:
  		hhpnt = (hhchan *)conpnt;
		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

		narate(pnt->v);			/* set Na rate chtyp */ 
	/* 	for (i=0; i<2; i++) {
		  m = hhpnt->m;
  		  h = hhpnt->h;
  		  dm = alpham - (alpham + betam) * hhpnt->mest;
  		  dh = alphah - (alphah + betah) * hhpnt->hest;
  		  m += (.5*dm + .5*hhpnt->dm);
  		  h += (.5*dh + .5*hhpnt->dh);
		  hhpnt->mest = m;
		  hhpnt->hest = h;
		}
	*/
		m=(hhpnt->m + (alpham+hhpnt->dm)*.5) / (1.0+(alpham+betam)*.5);
		h=(hhpnt->h + (alphah+hhpnt->dh)*.5) / (1.0+(alphah+betah)*.5);

/*
		m=(hhpnt->m + alpham) / (1.0+alpham+betam);
		h=(hhpnt->h + alphah) / (1.0+alphah+betah);
*/
		hhpnt->mest = m;
		hhpnt->hest = h;
  		conduct = m*m*m*h * hhpnt->maxcond;
                nodcur += conduct * hhpnt->vrev;
		tcondv += conduct;
/*     fprintf (stderr,"Na: m %10.4g h %10.4g cond %10.4g max %10.4g\n",
	m,h,conduct,hhpnt->maxcond); /* */
		break;

		case 1:				/* Sequential state */
		case 2:

	 	chpnt = (chan *)conpnt;
  		conduct = dochan2((chan *)conpnt) * chpnt->maxcond;
                nodcur += conduct * chpnt->vrev;
		tcondv += conduct;
		break;

		}  /* switch (stype) */
		break;

      	case K:   
  		hhpnt = (hhchan *)conpnt;
		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

		krate(pnt->v);			/* set K rate chtyp */ 
	/* 	for (i=0; i<2; i++) {
		  n = hhpnt->m;
  		  dn = alphan - (alphan + betan) * hhpnt->mest;
  		  n += (.5*dn + .5*hhpnt->dm)
		  hhpnt->mest = n;
		}
	*/
		n=(hhpnt->m + (alphan+hhpnt->dm)*.5) / (1.0+(alphan+betan)*.5);

/*		n=(hhpnt->m + alphan) / (1.0+alphan+betan);
*/
		hhpnt->mest = n;
  		conduct = n*n*n*n * hhpnt->maxcond;
                nodcur += conduct * hhpnt->vrev;
		tcondv += conduct;
		break;

		case 1:				/* Sequential state */
		case 2:

	 	chpnt = (chan *)conpnt;
  		conduct = dochan2((chan *)conpnt) * chpnt->maxcond;
                nodcur += conduct * chpnt->vrev;
		tcondv += conduct;
		break;

		}  /* switch (stype) */
		break;

         case SYNAPSE:                            /* currents don't need est */
         case ROD:
         case CONE:
         default:
                break;
        }  /* switch */
       }  /* for (lpnt= ; ; ) */

     err = pnt->vest;
     if (tcondv) {				/* if any new conductances */
        implf = 1. / (1. + (tcondv+pnt->tcondn+pnt->tcond) * pnt->k);
        pnt->vest = (pnt->v + nodcur * pnt->k) * implf;
     }
     else  {
 	implf = pnt->implf;
        pnt->vest = (pnt->v + nodcur * pnt->k) * implf;
     }
     if (pnt->miscfl & VEXT) {			/* comp is a battery */
	icomp = (pnt->vest - pnt->extv) / (pnt->k * implf);
	if (implicit) {
           pnt->extvi = icomp;
        }
	else {    			/* crank-nicholson */
           pnt->extvi = icomp * 0.5;
	}   
       pnt->vest = pnt->extv;
       pnt->verr = 1.0;			/* set verr so vclamp always runs */
     }
     else {
        err -= pnt->vest;                         /* find convergence deriv */
        pnt->vest -= err * pnt->relax;
        pnt->t++;
        if (err>0.0) pnt->g++;
	else err = -err;
        pnt->verr = err;
        if (err > maxerr) maxerr = err;
     }
}

/*------------------------------------*/

