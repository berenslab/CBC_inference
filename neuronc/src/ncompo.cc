/* Segment ncomp in Program nc */

/* Models neuronal circuits */


#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncelem.h" 
#include "ncomp.h"
#include "control.h"

#define DEBUG

extern comp *compnt;

double exp();
double sqrt();

extern double alpham,betam;
extern double alphah,betah;
extern double alphan,betan;

/*------------------------------------*/

double dochan(chpnt)
{
   int i,nstate,numtr,nsites;
   static chanstate *stpnt;
   static stconc tconc[NUMSTATE], *conc;
   static double (**rate)(), trrate;
   static float *trmul;
   double binomdev();

  nstate = chpnt->chtyp->numstate;
  stpnt = chpnt->chtyp->state;
  conc = chpnt->conc;
  for (i=0; i<nstate; i++,conc++) {  /* clear old estimate */
    conc->cest = 0;
    conc->dcon = 0;
  }	
  conc = chpnt->conc;
  for (i=0; i<nstate; i++,conc++) {  /* find change in conc */
     trans = stpnt->trans;
     rate  = stpnt->trate;
     trmul = stpnt->ratemul;
     numtr = stpnt->numtrans;
     for (j=0; j<numtr; j++,trans++,rate++,trmul++) {
       trrate = conc->cval * (*rate)() * *trmul * timinc;
       chpnt->conc[*trans].dcon += trrate;
       conc->dcon -= trrate;
     }
  }
  conduct = 0.0;
  conc = chpnt->conc;
  for (i=0; i<nstate; i++,conc++) { /* integrate and total cond */
    conc->cest = conc->cval + conc->dcon;
    conduct += conc->cest;			/* conduct <= 1.0 */
  }
  if (nsites=chpnt->nchan) {
    p = conduct / (nsites * chpnt->cq);           /* find binding prob */
    conduct = binomdev(p,nsites) * chpnt->cq;       /* bound  w/noise */
  }
     switch (conpnt->ctype)
 /* fprintf (stderr,"typ %d c %10.4g cond %10.4g\n",
		chpnt->ctype,conduct,chpnt->conduct); /* */
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
  int i,j,ncomp,tcomp,niter,nstate,ns,numtr;
  char *trans;
  static double nodcur;                         /* total node current */
  static double nodcr;                          /* total unvarying node curr */
  static comp *pnt,*cmpnt;
  static conn *conpnt;
  static double critx,maxerr,err,vest,delcrit;
  static double rpar,nodc,icomp,xxx;
  static conlst *lpnt;
  static recpar *cpt;                           /* receptor constants */
  static hhchan *hhpnt;
  static chanstate *stpnt;
  static chan *chpnt;
  static double (**rate)(), trrate;
  static float *trmul;
  static double m,n,h,dm,dn,dh;
  static double conduct;
  static stconc tconc[NUMSTATE], *conc;

#ifdef DEBUG
  if (debug & 1) fprintf (stderr,"runcomp\n");
#endif

 for (pnt=compnt; pnt; pnt=pnt->cnext)                  /* First estimate */
  {                                                     /* for all comps */
   nodcur = nodcr = 0.0;
   rpar = 0.0;
   for (lpnt=pnt->clst; lpnt; lpnt=lpnt->lnext) 
    {                                           /* check all connections */
     conpnt = lpnt->conpnt;
     if (! conpnt) break;
     if (conpnt->comp1 == pnt) cmpnt = conpnt->comp2; /* get other compartment*/
     else                      cmpnt = conpnt->comp1;
     switch (conpnt->ctype)
      {                         /* add connection conductances which change */

       case GJ:                                 /* may vary gj in future */
                nodcur += conpnt->conduct * cmpnt->v;   /* current thru g j */
                rpar   += conpnt->conduct;              /* total cond */
               break;
       case SYNAPSE:                                    /* cur thru syn chan */
                nodcr += ((synap *)conpnt)->conduct * ((synap *)conpnt)->vrev;
                rpar  += ((synap *)conpnt)->conduct;    /* total cond */
               break;

      case ROD:
      case CONE:                                /* find receptor instantan */
                                                /* conductance */
          cpt = ((recep *)conpnt)->chtyp;
          nodcr += ((recep *)conpnt)->conduct * cpt->vrev;/* cur thru rec chan*/
          rpar  += ((recep *)conpnt)->conduct;               /* total cond */
         break;

      case RESISTOR:
                nodcur += conpnt->conduct * cmpnt->v;  /* cur to next comp */
                break;                  /* conductance is collected in trt */
      case BATT:
     		if (cmpnt == conpnt->comp1)
                  nodcur += (cmpnt->v - conpnt->conduct) * BATTCOND;
		else 
                  nodcur += (cmpnt->v + conpnt->conduct) * BATTCOND;
                break;		/* cur to next comp through batt */

      case CAP: nodcr += (pnt->oldv - cmpnt->oldv) * conpnt->conduct/timinc;
                break;			/* -(v3 - v2) * c / t */

      case NA:  narate(pnt->v);			/* set Na rate chtyp */ 
  		hhpnt = (hhchan *)conpnt;

		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

		m = hhpnt->mest;		/* use previous best est */
  		h = hhpnt->hest;
		hhpnt->m = m;
		hhpnt->h = h;
  		dm = alpham - (alpham + betam) * m;
  		dh = alphah - (alphah + betah) * h;
  		hhpnt->dm = dm;
  		hhpnt->dh = dh;
  		m += dm * timinc; 
  		h += dh * timinc; 
		hhpnt->mest = m;
		hhpnt->hest = h;
  		hhpnt->conduct = m*m*m*h * hhpnt->maxcond;
  /* fprintf (stderr,"m %10.4g h %10.4g cond %10.4g max %10.4g\n",
	m,h,hhpnt->conduct,hhpnt->maxcond); /* */
                nodcur += hhpnt->conduct * hhpnt->vrev;
                rpar   += hhpnt->conduct;    		/* total cond */
		break;

		case 1:				/* Sequential state */
		case 2:

  		chpnt->conduct = dochan(conpnt) * chpnt->maxcond;
                nodcur += chpnt->conduct * chpnt->vrev;
                rpar   += chpnt->conduct;    		/* total cond */
		break;

		}	/* switch (stype) */
                break;

      case K:   krate(pnt->v);			/* set K rate chtyp */ 
  		hhpnt = (hhchan *)conpnt;

		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

		n = hhpnt->mest;
		hhpnt->m = n;
  		dn = alphan - (alphan + betan) * n;
  		hhpnt->dm = dn;
  		n += dn * timinc; 
		hhpnt->mest = n;
 /* fprintf (stderr,"n %10.4g cond %10.4g\n",n,hhpnt->conduct); /* */
  		hhpnt->conduct = n*n*n*n * hhpnt->maxcond;
                nodcur += hhpnt->conduct * hhpnt->vrev;
                rpar   += hhpnt->conduct;    		/* total cond */
                break;

		case 1:				/* Sequential state */
		case 2:

  		chpnt->conduct = dochan(conpnt) * chpnt->maxcond;
                nodcur += chpnt->conduct * chpnt->vrev;
                rpar   += chpnt->conduct;    		/* total cond */
		break;

		}	/* switch (stype) */
		break;

      default:  break;

      } /* switch */
    }      /* for (lpnt=; ;) */

   rpar += pnt->trt;                            /* total conductance */
   pnt->implf = 1. / (1. + rpar * pnt->k);      /* implicit factor   */
   nodcur -= pnt->v * rpar;                     /* current due to comp volts */
   nodcr  += pnt->vrev * pnt->rm;               /* membrane conductance */

/* fprintf (stderr,"nodcr %10.3g nodcur %10.3g\n",nodcr,nodcur); /* */

   if (pnt->extyp & IEXT)                       /* comp is a current src */
      nodcr += pnt->exti;                       /* external current */
   nodcur += nodcr;

   if (implicit) {				/* fully implicit */
     pnt->nodc = nodcr;				/* curr for estimate */
     nodcur *=  pnt->k;                      	/* first-ord only */
     pnt->vest = pnt->v + nodcur;            	/* simple first-order est */
     if (pnt->extyp & VEXT) {                   /* comp is a battery */
       pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k);
       pnt->vest = pnt->extv;
     }
     else {
       pnt->extvi = 0.0;
     }
   }						/* end of implicit */

   else	{     /* Crank-Nicholson or Euler */
     pnt->nodc = nodcur + nodcr;		/* curr for estimate */
     nodcur *=  pnt->k * 2.;                  /* first-ord only */
     pnt->vest = pnt->v + nodcur;             /* simple first-order est */
     if (pnt->extyp & VEXT) {                   /* comp is a battery */
       pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k*2.);
       pnt->vest = pnt->extv;			/*  battery voltage  */
     }
     else {
       pnt->extvi = 0.0;
     }
   }    /* crank-nicholson */

   pnt->verr = 1.;
 }

/* just save it:  only first order estimate */

if (euler) {
 for (pnt=compnt; pnt; pnt=pnt->cnext)
  {
   pnt->v = pnt->vest;
  }
   return;                              /* if return here, only first-order */
}


 maxerr = 1.;                           /* set maxerr so loop runs */ 
 tcomp = 0; 
 delcrit = sqrt(10.0);

 for (critx=1.000001e-3; critx >= crit; critx /= delcrit)
 {
 for (niter=0; maxerr > critx; niter++)         /* Iterate until convergence */
  {
   maxerr = 0.0;
   for (ncomp=0,pnt=compnt; pnt; pnt=pnt->cnext)        /* new estimate */
    {
    if (pnt->verr < critx) continue;  		/* */
#ifdef DEBUG
     tcomp++; ncomp++;
#endif
     nodcur = pnt->nodc;                        /* start with old node cur */
     for (lpnt=pnt->clst; lpnt; lpnt=lpnt->lnext) 
      {
       conpnt = lpnt->conpnt;                   /* check connection */
       if (! conpnt) break;
       if (conpnt->comp1 == pnt) cmpnt = conpnt->comp2;/* get other compartmnt*/
       else                      cmpnt = conpnt->comp1;
       switch (conpnt->ctype)
        {
         case RESISTOR:
         case GJ:
                nodcur += cmpnt->vest * conpnt->conduct;  /* curr thru resist */
                break;
         case BATT:
     		if (cmpnt == conpnt->comp1) 
                  nodcur += (cmpnt->vest - conpnt->conduct) * BATTCOND;
		else
                  nodcur += (cmpnt->vest + conpnt->conduct) * BATTCOND;
                break;		     /* curr thru batt small resistor */
         case CAP:
		nodcur += cmpnt->vest * conpnt->conduct / timinc;	
                break; 
         case NA:
		narate(pnt->vest);		/* set Na rate chtyp */ 
  		hhpnt = (hhchan *)conpnt;
		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

		m = hhpnt->m;
  		h = hhpnt->h;
  		dm = alpham - (alpham + betam) * hhpnt->mest; /* use new m,h */
  		dh = alphah - (alphah + betah) * hhpnt->hest;
  		m += (dm + hhpnt->dm) * 0.5 * timinc; 
  		h += (dh + hhpnt->dh) * 0.5 * timinc; 
		hhpnt->mest = m;
		hhpnt->hest = h;
  		hhpnt->conduct = m*m*m*h * hhpnt->maxcond;
                nodcur += hhpnt->conduct * hhpnt->vrev;
		break;

		case 1:				/* Sequential state */
		case 2:

  		chpnt = (chan *)conpnt;
		nstate = chpnt->chtyp->numstate;
		stpnt = chpnt->chtyp->state;
		for (i=0; i<nstate; i++) {	/* clear temp concs */
		  tconc[i].cest = 0;
		  tconc[i].dcon = 0;
		}
		conc = chpnt->conc;
		for (i=0; i<nstate; i++,conc++) { /* find change in conc */
		  trans = stpnt->trans;		  /* but use new conc est */
		  rate  = stpnt->trate;
		  trmul = stpnt->ratemul;
		  numtr = stpnt->numtrans;
		  for (j=0; j<numtr; j++,trans++,rate++,trmul++) {
		    trrate = conc->cest * (*rate)() * *trmul * timinc;
		    tconc[*trans].dcon += trrate;
		    tconc[i].dcon -= trrate;
		  }
		}
		conduct = 0.0;
		conc = chpnt->conc;
		for (i=0; i<nstate; i++,conc++) { /* integrate and total cond */
		   conc->cest = conc->cval + (conc->dcon + tconc[i].dcon) * 0.5;
		   conduct += conc->cest;
		}
  		chpnt->conduct = conduct * chpnt->maxcond;
 /* fprintf (stderr,"c %10.4g cond %10.4g\n",conduct,chpnt->conduct); /* */
                nodcur += chpnt->conduct * chpnt->vrev;
                rpar   += chpnt->conduct;    		/* total cond */
		break;

		}  /* switch (stype) */
		break;

      	case K:   
		krate(pnt->v);			/* set K rate chtyp */ 
  		hhpnt = (hhchan *)conpnt;
		switch (hhpnt->stype) {

		case 0:				/* Hodgkin-Huxley */

		n = hhpnt->m;
  		dn = alphan - (alphan + betan) * hhpnt->mest; /* use new n */
  		n += (dn + hhpnt->dm) * 0.5 * timinc; 
		hhpnt->mest = n;
  		hhpnt->conduct = n*n*n*n * hhpnt->maxcond;
                nodcur += hhpnt->conduct * hhpnt->vrev;
		break;

		case 1:				/* Sequential state */
		case 2:

  		chpnt = (chan *)conpnt;
		nstate = chpnt->chtyp->numstate;
		stpnt = chpnt->chtyp->state;
		for (i=0; i<nstate; i++) {	/* clear temp concs */
		  tconc[i].cest = 0;
		  tconc[i].dcon = 0;
		}
		conc = chpnt->conc;
		for (i=0; i<nstate; i++,conc++) { /* find change in conc */
		  trans = stpnt->trans;		  /* but use new conc est */
		  rate  = stpnt->trate;
		  trmul = stpnt->ratemul;
		  numtr = stpnt->numtrans;
		  for (j=0; j<numtr; j++,trans++,rate++,trmul++) {
		    trrate = conc->cest * (*rate)() * *trmul * timinc;
		    tconc[*trans].dcon += trrate;
		    tconc[i].dcon -= trrate;
		  }
		}
		conduct = 0.0;
		conc = chpnt->conc;
		for (i=0; i<nstate; i++,conc++) { /* integrate and total cond */
		   conc->cest = conc->cval + (conc->dcon + tconc[i].dcon) * 0.5;
		   conduct += conc->cest;
		}
  		chpnt->conduct = conduct * chpnt->maxcond;
 /* fprintf (stderr,"c %10.4g cond %10.4g\n",conduct,chpnt->conduct); /* */
                nodcur += chpnt->conduct * chpnt->vrev;
                rpar   += chpnt->conduct;    		/* total cond */
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
     pnt->vest = (pnt->v + nodcur * pnt->k) * pnt->implf;
     if (pnt->extyp & VEXT) {			/* comp is a battery */
	icomp = (pnt->vest - pnt->extv) / (pnt->k * pnt->implf);
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
        pnt->vest -= err * relax;
        err = abs (err);
        pnt->verr = err;
        if (err > maxerr) maxerr = err;
     }
    }           /* for (ncomp=0,pnt=compnt;) */

#ifdef DEBUG
  if (debug & 4) fprintf (stderr,"ncomp %d maxerr %8.3g\n",ncomp,maxerr);
#endif
  
   }            /* for ( ; maxerr>critx; )  */

/* reset error each level of critx */

/* for (pnt=compnt; pnt; pnt=pnt->cnext)
  {
   pnt->verr = 1.0; 
  }
*/

#ifdef DEBUG
  if (debug & 2) fprintf (stderr,"niter %d critx %8.3g\n",niter,critx);
#endif

  }         /* for (critx ; ; ) */

#ifdef DEBUG
  if (debug & 2) fprintf (stderr,"tcomp %d\n",tcomp);
#endif
/* just save it */

 for (pnt=compnt; pnt; pnt=pnt->cnext)
  {
   pnt->oldv = pnt->v = pnt->vest;
  }
}

/*------------------------------------*/

