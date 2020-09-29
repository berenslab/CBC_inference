/* Module ncconvert in Program NC */

/* Translates neuron cable description into
   the equivalent compartmental model. */

/* Assume a dendrite is made up of small isopotential
compartments which are modeled as electrical nodes.
Each node is connected to its neighbors through axial 
resistances, and has conduction pathways to ground
through membrane resistance and capacitance.
Additional conduction pathways are allowed for any node. 

Unless otherwise noted, all length values are calibrated
in meters, all voltages are in volts, all currents are
in amperes, and capacitance is in farads.

	Oct 87			R.G. Smith

*/

extern "C" {

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#endif
#include <math.h>
#include <string.h>

}

#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "ncelem.h"
#include "control.h"
#include "ncio.h"

#define CMMSCL  1e-2			/* scale from cm to meters */
#define CMMSCL2 1e-4			/* scale from cm2 to meters2 */
#define UMMSCAL 1e-6			/* scale from um to meters */
#define CMUMSCL  1e4			/* scale from cm to microns */
#define CMUMSCL2 1e8			/* scale from cm2 to microns2 */

#define DEBUG

#ifdef DEBUG
#include "ncdebug.h"
#endif

extern int cumelem,cumnode;
extern int cumcomp;			/* total number of compartments */
extern int cumconn;			/* total number of connections */
extern int cumchan;			/* total number of channels */
extern int cumcacomp;			/* total number of Ca comps */
extern int cumcycacomp;			/* total number of cycA comps */
extern int cumcycgcomp;			/* total number of cycG comps */
extern elem *elempnt;			/* pointer to element list */
extern elem *elemend;			/* pointer to end of element list */
extern comp *compnt;			/* pointer to compartment list */
extern comp *compend;			/* pointer to end of compartment list */
extern conn *connpnt;			/* pointer to connection list */
extern conn *connend;			/* pointer to end of connection list */
extern node *nodepnt;			/* pointer to node list */
extern synap *synpnt;			/* pointer to synapse list */

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double sqrt(double);

#ifdef __cplusplus
}
#endif

comp *makcomp(elem *epnt, double rm, double cap, double vrest, double vrev);
char *findsym(int);
char *nodstr(int a, int b, int c, int d);
chan *makchan(attrib *apnt, comp *comp1, comp *comp2, int saved);
chan *makca(elem *epnt, cattrib *apnt, comp *comp1, double area, int saved);
chan *makcacomp(elem *epnt, cattrib *apnt, comp *comp1, double area);
void redocacomp(comp *ocpnt, double oarea);
synap *maksynap(synapse *epnt, comp *comp1, comp *comp2);
photrec *makphotrec(photorec *epnt, comp *cpnt);
load *makload(loadelem *epnt, comp *comp1);
conn *makconn(comp *comp1, comp *comp2, double s, int type);
double dist3d(node *n1, node *n2);

void execerror (const char*, const char*);
void calcnod (double Ri, double Rm, double Cap, double d, double clen, \
		double *rri, double *rrm, double *rcap);
void maklst(conlst **head, conn *cpnt);
void chandens(elem *epnt, comp *cpnt, double area);
void modload(loadelem *lepnt);
void modconn(elem *epnt, double conduct);
synap *modsynap(synapse *sepnt);
photrec *modphotrec(photorec *pepnt);
void checkcomp(elem *epnt, comp **pnt1, comp **pnt2);
dbuf *makdelbuf(elem *epnt, comp *comp1, comp *comp2);
ndbuf *makndbuf(elem *epnt, comp *comp1, comp *comp2);
ntcomp *makntcomp (comp *cpnt, int type);

void freelst (conlst *clst);
void initsyn(void);
void setsynv(synap *spnt, double v);
chan *addchan(comp *pnt, chan *chpnt);
void prelem(elem *epnt);
comp *othercomp (comp *pnt,conn *cpnt);
void delelem (elem *epnt);
void setnodlst(node *nodp, elem *epnt);
void unsetnodlst(node *nodp, elem *epnt);
char *prnode(int n1, int n2, int n3, int n4);
void delhashcomp (comp *cpnt);
gj *makgj(gapjunc *cepnt, comp *comp1, comp *comp2, int type);
gj *modgj(gapjunc *cepnt);
pnx *makpnx(pannex *pepnt, comp *comp1, comp *comp2, int type);
chan *modchan(elem *cepnt, attrib *apnt, double lcap);
ntcomp *makaddntcomp (comp *cpnt, int type, double val, int n);
void delntcomp (ntcomp *npnt);
void condense_comp_chan(comp *pnt);
elem *savelem (elem *epnt);
void efree (void *ptr);
void delcacomp(cacomp *capnt);

/*------------------------------------*/

double calc_complen(double Ri,double Rm,double dia, double cplam)

{ 
    double ri, rm, lambda, complen;

  ri = Ri * CMUMSCL / (MPI * dia * dia / 4);	/* calc ri (ohm/m) */
  rm = Rm * CMUMSCL2 / (MPI * dia);		/* calc rm (ohm*m) */
  lambda = sqrt (rm / ri); 
  complen = lambda * cplam;			/* compartment length */
  if (complen==0) complen = .001;

    if (debug & NCCONVERT && debugz & 4)
	ncfprintf (stderr,"calc_complen Ri %g Rm %g lam %g cplam %g complen %g\n",
		Ri,Rm, lambda,cplam,complen);

  return (complen);
}

/*------------------------------------*/

void initcomp(void)

/* Add new compartments to the model */
/* or modify existing compartments. */
/* If element is to be saved for future reference, */
/* store its low-level structure in epnt->lpnt. */
/* If element is actually a modification of an existing */
/* element, then don't remake it, but call the */
/* modification routine. */
   
{
    int i,measrcomp,numcomp,ncomp,typ;
    static double ri,rm,cri,crm,lambda,ccap,lcap,complen;
    static double len,dia,dia2;
    static double Ri,Rm,cplam,ldia,len_ratio,cumlen;
    double taper_factor, maxd, mind;
    elem *epnt,*tepnt;
    attrib *apnt;
    node  *npnt;
    comp  *cpnt,*cpntn,*cpntold,*cpntlast;
    synap *tsynap;
    photrec *trecep;
    conn *tconn;

    cable *cepnt;
    sphere *sepnt;
    loadelem *lepnt;
    capac *cappnt;
    electrode *elecpnt;
    gapjunc *gjepnt;
    pannex *pannexpnt;

if (debug & NCCONVERT) ncfprintf (stderr,"initcomp\n");

			/* reduce cable length by sphere radius */

 for (epnt=elempnt; epnt; epnt=epnt->next) {		/* for all elements */
        conlst *lpnt;
        double mindia1, mindia2;
        int found1,found2;

   if (epnt->lptr) continue;	/* skip if element has already been processed */

   if (epnt->ctype==CABLE) {

      found1 = 0;
      if (npnt=epnt->nodp1) {		/* check for spheres at node 1 */
        mindia1 = 1e10;
        for (lpnt=npnt->elemlst; lpnt; lpnt=lpnt->next) {
           if (sepnt=(sphere *)lpnt->conpnt) {
	     if (sepnt->ctype==SPHERE) {
               found1 = 1;
               if (mindia1 > sepnt->dia) mindia1 = sepnt->dia;
             } 
           }
        }  /* for (lpnt;;) */
	if (!found1) mindia1 = 0.0;
      }  /* nodp1 */

      found2 = 0;
      if (npnt=epnt->nodp2) {		/* check for spheres at node 2 */
        mindia2 = 1e10;
        for (lpnt=npnt->elemlst; lpnt; lpnt=lpnt->next) {
           if (sepnt=(sphere *)lpnt->conpnt) {
	     if (sepnt->ctype==SPHERE) {
               found2 = 1;
               if (mindia2 > sepnt->dia) mindia2 = sepnt->dia;
             } 
           }
        }  /* for (lpnt;;) */
	if (!found2) mindia2 = 0.0;
      }  /* nodp2 */

      cepnt=(cable*)epnt;
      if ((len=cepnt->length)==NULLVAL) {             /* if length not set */
          if (!epnt->nodp1 || !epnt->nodp2) len = 10.0; /* default length */
          else len = dist3d(epnt->nodp1,epnt->nodp2);   /* calc length */
	  cepnt->length = len - (mindia1 + mindia2) * 0.5;
	  if (cepnt->length < 0.0) cepnt->length = 1e-3;
      }
   }  /* if CABLE */
 }  /* for (epnt;;) */

 for (epnt=elempnt; epnt; ) {				/* for all elements */

   if (!epnt->newelem && !epnt->modif) {	/* skip if old and unchanged */
	epnt = epnt->next;
	continue;
   }
   typ = epnt->ctype; 

/*  if (epnt->modif) 
	ncfprintf (stderr,"modify elem: num %d typ %d '%s'\n",
		epnt->modif,epnt->ctype,findsym(epnt->ctype)); /* */

   if ((debug & NCCONVERT) && (debugz & 4)) {

        prelem(epnt);
	ncfprintf (stderr,"elem: num %d typ %d '%s' lptr %x\n",
		epnt->elnum,epnt->ctype,findsym(epnt->ctype),epnt->lptr); /* */
  }

  if (epnt->newelem || epnt->modif) { /* if element hasn't been processed yet */

  switch (typ) {

  case CABLE:			
  {
     double cvrev, cvrest;
				/* initial calculations */
  cepnt = (cable *)epnt;
				/* Don't allow modification of cable */
  if (epnt->modif) { }		/* because that would change number of */	
				/* compartments. Actually, it's do-able... */
				/* Just find old voltages, erase the old comps */
				/*  and make new ones. But no time now...*/
  else {
   if (epnt->lptr==NULL) { 		/* if not made (not saved) yet */
     if ((len=cepnt->length)==NULLVAL) {    /* and if length not set */
	if (!epnt->nodp1 && !epnt->nodp2) break;
	len = dist3d(epnt->nodp1,epnt->nodp2); /* calc length from locs */
     }
      if (len < 1e-6) { 
         char buf1[40], buf2[40];

      len = 0.1;
  /*
      strcpy (buf1,prnode(epnt->node1a,epnt->node1b,
			  epnt->node1c,epnt->node1d));
      strcpy (buf2,prnode(epnt->node2a,epnt->node2b,
		          epnt->node2c,epnt->node2d));
      ncfprintf (stderr,
	"Cable %d from node %s to %s: endpoints are at same location.\n",
            epnt->elnum,buf1,buf2);
	  execerror ("Initcomp: ","Cable needs location or length");
	  break;
*/
     }
    if ((Ri=cepnt->Ri)<=0) { 				  /* use default values if nec. */
      Ri = dri / exp(log(dqri) * (tempcel-dbasetc)/10.0); /* add temp coeff for Ri, Rm */
    }
    if ((Rm=cepnt->Rm)<=0) { 
      Rm = drm / exp(log(dqrm) * (tempcel-dbasetc)/10.0); /* def val = 1, no effect */
    }
    if ((cvrev=cepnt->vrev)   == NULLVAL) cvrev  = vcl;
    if ((cvrest=cepnt->vrest) == NULLVAL) cvrest = vcl;
    if ((cplam = cepnt->cplam)==NULLVAL) cplam = complam;


    dia  = cepnt->dia;
    dia2 = cepnt->dia2;

    if (dia  != NULLVAL && dia  < 1e-2) dia = 1e-2;		/* check for zero diameter */
    if (dia2 != NULLVAL && dia2 < 1e-2) dia2 = 1e-2;

    if (dia==NULLVAL) {
	   char buf1[40],buf2[40];

      strcpy (buf1,prnode(epnt->node1a,epnt->node1b,
			  epnt->node1c,epnt->node1d));
      strcpy (buf2,prnode(epnt->node2a,epnt->node2b,
		          epnt->node2c,epnt->node2d));
	  ncfprintf (stderr,
	"# Diameter unspecified for cable %d from node %s to %s, using 1 um.\n", 
		epnt->elnum,buf1,buf2);
	  dia = 1.0;
	}
    if (dia2==NULLVAL) dia2 = dia;

    if (dia > dia2) { maxd = dia;  mind = dia2; }
    else            { maxd = dia2; mind = dia; }

    complen = calc_complen(Ri,Rm,maxd,cplam);	/* complen for near end */

    if (dia != dia2) {		/* normalize by integral of sqrt(d) */
				/*  Much thanks to MCW van Rossum */
      taper_factor = 2/3.0 * (maxd*sqrt(maxd)-mind*sqrt(mind)) / (maxd-mind);
      complen *= taper_factor;
    }

    measrcomp = (int)(len / complen) + 1;	/* number of comps in cable */
    numcomp = measrcomp + 1;
    len_ratio = len / measrcomp / complen;
    complen = len / measrcomp;

    if (debug & NCCONVERT && debugz & 4)
	ncfprintf (stderr,"cable num  %d Ri %g Rm %g lam %g cplam %g len %g complen %g numcomp %d\n",
		epnt->elnum,cepnt->Ri,cepnt->Rm,
		lambda,cplam,len, complen,numcomp);

    if ((lcap=cepnt->Cm)==NULLVAL) lcap = dcm;
    if (lcap==0) lcap = 1e-6;

    cpntold=(comp*)NULL;
    ncomp = numcomp-1;
    /* cumlen = 0; */
    for (i=0; i<numcomp; i++) {		/* make comps for one cable seg */

					/* calc complen for this comp only*/; 
      ldia = ((ncomp-i)*dia + i*dia2) / ncomp;
      if (ldia <= 0) ldia = 1e-20;	/* if taper to zero (um), make non-zero */

      complen = calc_complen(Ri,Rm,ldia,cplam) * len_ratio;

	/* if ((i==0) || (i==numcomp-1)) cumlen += complen*0.5;/
	else cumlen += complen;
	ncfprintf (stderr,"i %d len %g cumlen %g\n",i,len,cumlen);
	*/

      calcnod (Ri,Rm,lcap,ldia,complen,&cri,&crm,&ccap);

      if (i==0) {			/* first end (beginning) */
	  				/* first, find start node */
	 if (! (npnt=epnt->nodp1)) {
	   ncfprintf(stderr,"initcomp: elem %d: can't find node 1: %s\n",
			epnt->elnum,nodstr(epnt->node1a,epnt->node1b,
				           epnt->node1c,epnt->node1d));
 	    execerror("Error","");
		break;
	  }
	  cpntn = npnt->comptr;		/* check for compartment at start */
	  if (cpntn) {			/* node already has compartment   */
	    if (cpntn->ctype != SPHERE) {
	      cpntn->vrev = cvrev;      /* cable controls node reversal v */
	      cpntn->v = cvrest;        /* cable controls node initial v  */
	    }
	    cpntn->rm  += 0.5/crm;	/* add to node conductance        */
            cpntn->cap += ccap/2;	/* membrane capacitance           */
	    cpntlast = cpntn;		/* save node comp pntr for connection */
	  }
	  else {			/* this is first comp for node */
					/* make the new compartment */
	    cpnt = makcomp(epnt,0.5/crm,ccap/2,cvrest,cvrev);
	    cpnt->cplam = cplam;	/* save complam for condense() below */
	    npnt->comptr = cpnt;	/* save comp pntr for node */ 
	    npnt->compnum = cpnt->num;	/* save comp number for node */
            maklst(&cpnt->nodlst,(conn*)npnt);	/* save node pntr for comp */
            cpntlast = cpnt;		/* last pointer for next connection */
	  }
	 if (epnt->attpnt)
	      chandens(epnt,cpntlast,ccap*0.5/lcap); /* volt-sens chan */
      }
      else if(i==(numcomp-1)) {		/* second end */
	 				/* find far end node */
	 if (! (npnt=epnt->nodp2)) {
	   ncfprintf(stderr,"initcomp: elem %d: can't find node 2: %s\n",
			epnt->elnum,nodstr(epnt->node2a,epnt->node2b,
					   epnt->node2c,epnt->node2d));
 	     execerror("Error","");
  	     break;
	 }
	 cpntn = npnt->comptr;		/* check for compartment at end */
	 if (cpntn) {			/* node already has compartment */
	    tconn = makconn(cpntn,cpntold,
				1/cri,AXIALRES);/* axial resis to last comp */
	    if (cpntn->ctype != SPHERE) {
	      cpntn->vrev = cvrev;      /*cable controls node rev */
	      cpntn->v = cvrest;        /*cable controls node init */
	    }
	    cpntn->rm  += 0.5/crm;	/* add conductance */
            cpntn->cap += ccap/2;	/* membrane capacitance         */
	    cpnt = cpntn;
	  }
	  else {			/* this is first comp for node */
					/* make the new compartment */
	    cpnt = makcomp(epnt,0.5/crm,ccap/2,cvrest,cvrev);
	    cpnt->cplam = cplam;	/* save complam for condense() below */
	    tconn = makconn(cpnt,cpntold,
			1/cri,AXIALRES); /* axial resis to last comp */
	    npnt->comptr = cpnt;	/* save comp pntr for node */ 
	    npnt->compnum = cpnt->num;	/* save comp number for node */
            maklst(&cpnt->nodlst,(conn*)npnt);/* save node pntr for comp */
	  }
	 if (epnt->attpnt)
	   chandens(epnt,cpnt,ccap*0.5/lcap); /* if volt-sens chan */

					/* even if there are only 2 comps, */
    	if ( numcomp==2) 		/* save the conn to remember cable */ 
	     epnt->lptr = (elem *)tconn; 	     /* save low-lev ptr */
		
      }
      else {	/* i>0 && i<numcomp-1 */   /* this is middle of branch */
					   /* make the new compartment */
	cpnt = makcomp(epnt,1.0/crm,ccap,cvrest,cvrev);
	cpnt->cplam = cplam;		/* save complam for condense() below */
	tconn = makconn(cpnt,cpntold,
			1/cri,AXIALRES); /* axial resis to last comp */
	if (epnt->attpnt)
		 chandens(epnt,cpnt,ccap/lcap); /* voltage-sens chans */
	cpntlast = cpnt;
					/* if there are more than 2 comps, */
    	if (numcomp>2) 			/* save first conn to remember cable */
    	 	  if (i==1) epnt->lptr = (elem *)tconn;   /* save low-lev ptr */
      }
      cpntold = cpntlast;		/* save the current pointer */

      }   /* for ( ; i<numcomp; ) */
    }	/* if (lptr==NULL) */
   }   /* else if not (modif) */
  }    /* case CABLE: */
   break;

     case SPHERE:
	   {
              double svrev, svrest;

	   sepnt = (sphere *) epnt;
    	   if ((svrev=sepnt->vrev)  == NULLVAL)  svrev  = vcl;
    	   if ((svrest=sepnt->vrest) == NULLVAL) svrest = vcl;
	   if ((Rm=sepnt->Rm) <= 0) Rm = drm;
  	   dia = sepnt->dia * UMMSCAL;
 	   rm = Rm * CMMSCL2 / (MPI * dia * dia);
    	   if ((lcap=sepnt->Cm)==NULLVAL) lcap = dcm;	/* lcap is Cm */
 	   ccap = lcap * MPI * dia * dia / CMMSCL2;  

	 if (! (npnt=epnt->nodp1)) {
	  ncfprintf(stderr,"initcomp: sphere %d: can't find node %s\n",
			epnt->elnum,nodstr(epnt->node1a,epnt->node1b,
				    	   epnt->node1c,epnt->node1d));
 	    execerror("Error","");
		break;
	  }
	  if (epnt->lptr==NULL) {
	   cpnt = npnt->comptr;		/* check for existing compartment */
	   if (cpnt) {			/* node already has compartment */
	     if (typ == SPHERE) {
	      cpnt->vrev  = svrev;      /* sphere controls node reversal v */
	      cpnt->v = svrest;         /* sphere controls node initial v */
	     }
	    cpnt->rm  += 1/rm;		/* add conductance 		*/
            cpnt->cap += ccap;		/* membrane capacitance         */
	   }
	   else {
	    cpnt = makcomp(epnt,1/rm,ccap,svrest,svrev);
	    npnt->comptr = cpnt;	/* save comp pntr for node 	*/ 
	    npnt->compnum = cpnt->num;	/* save comp number for node */
            maklst(&cpnt->nodlst,(conn *)npnt);	/* save node pntr for comp */
	   }
	   if (epnt->attpnt)
		 chandens(epnt,cpnt,ccap/lcap); /* voltage-sens chans */
	   epnt->lptr = (elem *)cpnt; /* save comp ptr */
	   } /* if (epnt->lptr==NULL) */
	  }
	  break;

      case LOAD:
           {
	       double lvrev, lvrest;

	   lepnt = (loadelem*)epnt;
	   ccap = SMALLCAP;
	   if ((lvrev=lepnt->vrev)   == NULLVAL) lvrev = 0;
	   if ((lvrest=lepnt->vrest) == NULLVAL) lvrest = 0;
  	   if (epnt->modif) {
	       modload(lepnt);			/* modify old load */
	   }
	   else {
		load *tload;

	  if (epnt->lptr==NULL) { 			/* check lowlev ptr */
	    if (! (npnt=epnt->nodp1)) {
	     ncfprintf(stderr,"initcomp: load %d: can't find node %s\n",
			epnt->elnum,nodstr(epnt->node1a,epnt->node1b,
				    epnt->node1c,epnt->node1d));
 	        execerror("Error","");
		break;
	    }
	    cpnt = npnt->comptr;	/* check for existing compartment */
	    if (!cpnt) {		/* node doesn't have compartment */
	        cpnt = makcomp(epnt,0.0,ccap,lvrest,lvrev);
	        npnt->comptr = cpnt;	       /* save comp pntr for node  */ 
	        npnt->compnum = cpnt->num;     /* save comp number for node */
                maklst(&cpnt->nodlst,(conn *)npnt); /* save node pntr for comp*/
	    }
	    tload = makload((loadelem*)epnt,cpnt); 
	    epnt->lptr = (elem *)tload; /* save synap ptr */
	   }
	  }
	 }
	   break;

      case GNDCAP:
	  {
		double capvrest;

	   cappnt = (capac *)epnt;
	   rm = LARGERES;
	   ccap = cappnt->c;
	   if ((capvrest=cappnt->vrest) == NULLVAL) capvrest = 0;
           if (! (npnt=epnt->nodp1)) {
	   ncfprintf(stderr,"initcomp: sphere %d: can't find node %s\n",
			epnt->elnum,nodstr(epnt->node1a,epnt->node1b,
				           epnt->node1c,epnt->node1d));
 	     execerror("Error","");
		break;
	   }
	  cpnt = npnt->comptr;		/* check for existing compartment */
	  if (cpnt) {			/* node already has compartment */
	   if (cappnt->vrest != NULLVAL) /* gndcap controls node initial v */
	          cpnt->v = capvrest;    
	    cpnt->rm  += 1/rm;		/* add conductance 		*/
            cpnt->cap += ccap;		/* membrane capacitance         */
	  }
	  else {
	    cpnt = makcomp(epnt,1/rm,ccap,capvrest,(double)NULLVAL);
	    npnt->comptr = cpnt;	/* save comp pntr for node 	*/ 
	    npnt->compnum = cpnt->num;  /* save comp number for node    */
            maklst(&cpnt->nodlst,(conn *)npnt);	/* save node pntr for comp */
	  }
	  if (epnt->attpnt)
		 chandens(epnt,cpnt,cpnt->cap/lcap); /* voltage-sens chans */
          }
	  epnt->lptr = (elem *)cpnt; /* save ptr */
	 break;

      case RESISTOR:

	 ri = ((resistor *)epnt)->r;
	 if (ri == 0.0) ri = 1.0;
	 if (epnt->modif) {
		modconn(epnt,1/ri);			/* modify old resist */
	 }
	 else {
	  if (epnt->lptr==NULL) { 			/* check lowlev ptr */
	   checkcomp(epnt,&cpntn,&cpntlast);
	   tconn=makconn(cpntn,cpntlast,1/ri,typ); /* axial resis to comp */
	   epnt->lptr = (elem *)tconn; /* save ptr */
	  }
	 }
 	 break;

      case ELECTRODE:

	 elecpnt = (electrode *)epnt;
	 if ((ri=elecpnt->r)==NULLVAL) ri = drs;
	 if ((ccap=elecpnt->c)==NULLVAL) ccap = deleccap;
	 if (ri == 0.0) ri = 1.0;
	 if (epnt->modif) {
	    //	modconn(epnt,1/ri);			/* modify old resist */
	 }
	 else {
	    double vrest, ccomp;
	  if (epnt->lptr==NULL) { 			/* check lowlev ptr */
	   checkcomp(epnt,&cpntn,&cpntlast);
	   tconn=makconn(cpntn,cpntlast,1/ri,RESISTOR); /* axial resis to comp */
	   epnt->lptr = (elem *)tconn; /* save ptr */
	   if ((ccomp=elecpnt->ccomp)!=NULLVAL) {	
 	     if (cpntlast) {
	 	cpntlast->cap += ccomp; 	/* add cap to comp at elec tip, can be neg */
		if (cpntlast->cap <= 0) cpntlast->cap = 1e-13; // SMALLCAP;
	     }
	   }

	   cpnt = cpntn;		/* get compartment */
	   if (cpnt) {			/* node already has compartment */
	     if ((vrest=elecpnt->vrest) != NULLVAL) /* gndcap controls node initial v */
	          cpnt->v = vrest;    
              cpnt->cap += ccap;	/* electrode capacitance         */
	      cpnt->rm  += 1/(LARGERES*0.1); /* add conductance 		*/
	   }
	  }
	 }
 	 break;

      case DIODE:

	 ri = ((diode *)epnt)->r;
	 if (ri == 0.0) ri = 1.0;
	 if (epnt->modif) {
		modconn(epnt,1/ri);			/* modify old resist */
	 }
	 else {
	  if (epnt->lptr==NULL) { 			/* check lowlev ptr */
	   checkcomp(epnt,&cpntn,&cpntlast);
	   tconn=makconn(cpntn,cpntlast,1/ri,typ); /* axial resis to comp */
	   epnt->lptr = (elem *)tconn; /* save ptr */
	  }
	 }
 	 break;

      case CAP:
	 cappnt = (capac *)epnt;	 
	 checkcomp(epnt,&cpntn,&cpntlast);
	 makconn(cpntn,cpntlast,cappnt->c,typ); /* series capac */
	 epnt->lptr = (elem *)cpntn; /* save ptr */
 	 break;

      case BATT:

	 checkcomp(epnt,&cpntn,&cpntlast);
	 makconn(cpntn,cpntlast,((batt *)epnt)->v,typ); /* series batt */
	 epnt->lptr = (elem *)cpntn; /* save ptr */
 	 break;

      case BUF:			/* voltage buffer (voltage follower) */
	{
         comp  *comp1,*comp2;
	 dbuf *tdpnt;
	 int nvclamps;

	 checkcomp(epnt,&comp1,&comp2);
	 tdpnt = makdelbuf(epnt,comp1,comp2);
	 nvclamps = comp2->miscfl | VEXT;
	 nvclamps++;
	 if (nvclamps > VEXT) nvclamps=VEXT;
	 comp2->miscfl |= (VEXT & nvclamps)|VBAT|VBUF;/* set flags*/
	 epnt->lptr = (elem *)comp1; /* save ptr */
 	 break;
	}

      case NBUF:		/* nt buffer (voltage sets ntrans) */
	{
         comp  *comp1,*comp2;
	 ndbuf *tnpnt;

	 checkcomp(epnt,&comp1,&comp2);
	 tnpnt = makndbuf(epnt,comp1,comp2);
	 makntcomp(comp2,((nbuf*)epnt)->ntrans);
	 epnt->lptr = (elem *)comp1; /* save ptr */
 	 break;
	}

      case GJ:				/* gap junction */

	 gjepnt = (gapjunc *)epnt;
	 if (epnt->modif) {
		modgj(gjepnt);		/* modify old gj */
	 }
	 else {
	  if (epnt->lptr==NULL) { 	/* check lowlev ptr */
                 gj *tgj;

	     checkcomp(epnt,&cpntn,&cpntlast);
	     tgj=makgj(gjepnt,cpntn,cpntlast,typ); /* make new gj */
	     epnt->lptr = (elem *)tgj; /* save gj ptr */
	  }
	 }
 	 break;

      case PNX:				/* Pannexin */
	 pannexpnt = (pannex *)epnt;
	 if (epnt->modif) {
	 	modgj(gjepnt);		/* modify old gj */
	 }
	 else {
	  if (epnt->lptr==NULL) { 	/* check lowlev ptr */
                 pnx *tpnx;

	     checkcomp(epnt,&cpntn,&cpntlast);
	     tpnx=makpnx(pannexpnt,cpntn,cpntlast,typ); /* make new pnx */
	     epnt->lptr = (elem *)tpnx; /* save pnx ptr */
	  }
	 }
 	 break;

      case CACOMP:
      case CHAN:		/* chan or ca comp with compartment */	
	 if (epnt->attpnt) apnt = epnt->attpnt;
         else { 
		ncfprintf (stderr,"initcomp chan: can't find chan attrib\n");
 	        execerror("Error","");
		break;
	 }
	 switch (apnt->ctype) {
	    case AMPA:
	    case KAINATE:
	    case GABA:
	    case GLY:
	    case CHRC:
	    case CGMP:
	    case NA:
	    case K:
	    case ClCa:
	    case KCa:
	    case CA: break;

	    case VESNOISE:
	    case CCHNOISE:
		   ncfprintf (stderr,"initcomp chan: Illegal noise.\n");	  
		  break;
	 }
	 if (epnt->modif) {
	   modchan(epnt,apnt,lcap); /* volt-sens chan*/
	 }
	 else {
	  if (epnt->lptr==NULL) { 	/* check lowlev ptr */
	   if (! (npnt=epnt->nodp1)) {
	    ncfprintf (stderr,"initcomp chan: can't find node %s\n",
			nodstr(epnt->node1a,epnt->node1b,
			       epnt->node1c,epnt->node1d));
 	        execerror("Error","");
		break;
	   }				/* Check back at node, to */
	   cpnt = npnt->comptr;		/* look for existing compartment */
	   if (!cpnt) {
	     ccap = 0.0;
    	     cpnt = makcomp(epnt,1.0/LARGERES,SMALLCAP,NULLVAL,
			((chattrib *)apnt)->vrev);
	     npnt->comptr = cpnt;	/* save comp pntr for node 	*/ 
	     npnt->compnum = cpnt->num;  /* save comp number for node    */
             maklst(&cpnt->nodlst,(conn *)npnt);  /* save node pntr for comp */
	     cpnt->ctype = apnt->ctype;	/* set chan type only if new */
	   }
	   else {
	     ccap = cpnt->cap;
	   }

           for ( ; apnt; apnt=apnt->attpnt) {
	     chan *chpnt;
	     int chtyp;
             chattrib *chapnt;
	     double area;

	   chapnt = (chattrib *)apnt;
    	   switch (chapnt->ctype) {

	   case VESNOISE:
	   case CCHNOISE: break;	/* do nothing here if noise attrib, etc. */

	     default:
	     case  NA:
	     case   K:
	     case ClCa:
	     case KCa:
	     case  CA:
	    case CGMP: 
				/* When "density" specified, we need "area". */

	      if ((ccap<=SMALLCAP) && 
                  (chapnt->maxcond==NULLVAL) && 
		  ((chapnt->ndensity!=NULLVAL) ||
		   (chapnt->density !=NULLVAL))) {
                   double r;

                r = CACOMPRAD * 1e-4;		/* calib in cm2 */
               area = 4.0 * MPI * r * r;
     ncfprintf (stderr,"# nc: initcomp: compartment for channel has not been defined.\n");
     ncfprintf (stderr,"# Using diameter of %g um for compartment.",CACOMPRAD*2); 
             }
             else
	       area = ccap/lcap;     /* Area used when comp already specified */
				     /*  calibrated in cm2 */

				     /* convert density to N */

	     if (chapnt->ndensity!=NULLVAL) {
	          chapnt->n = chapnt->ndensity * area * UM2CM2;
	     } 
	     if (chapnt->density != NULLVAL) {
	        chapnt->maxcond = chapnt->density * area;
	     }

             chtyp = chapnt->ctype;
	     switch (chtyp) {

	      case NA:	
	      case K:
	      case ClCa:
	      case KCa:
	      case KAINATE:
	      case GABA:
	      case CHRC:
	      case GLY:
	             chpnt=makchan(apnt,cpnt,(comp*)NULL, epnt->nocondens);
		     break;

	      case CACOMP:
		    chpnt=(chan*)makcacomp(epnt,(cattrib*)apnt,cpnt,area); 
		    break;

	      case AMPA:
	      case CGMP:
	      case CA:
		    chpnt=makca(epnt,(cattrib*)apnt,cpnt,area, epnt->nocondens); 
		    break;

	      }  /* switch (chapnt->ctype) */

	      epnt->lptr = (elem *)chpnt; /* save chan ptr */
	      apnt->lptr = (elem *)chpnt; /* save chan ptr */
	     } /* switch (chapnt->ctype) */
           }  /* for (;apnt;) */
	  }  /* if (lptr==NULL) */
	 }  /* else */
 	 break;

      case SYNAPSE:

	 if (epnt->modif) {
	    modsynap((synapse*)epnt);		/* modify old synapse */
	 }
	 else {
	  if (epnt->lptr==NULL) { 		/* check lowlev ptr */
	    checkcomp(epnt,&cpntn,&cpntlast); 
	    tsynap = maksynap((synapse*)epnt,cpntn,cpntlast); 
	    epnt->lptr = (elem *)tsynap; /* save synap ptr */
	  }
	 }
	 break;

      case ROD:
      case CONE:
      case CHR:
      case VTRANSDUCER:
      case ITRANSDUCER:

	 if (epnt->modif) {
	   modphotrec((photorec*)epnt);		/* modify old photrec */
	 }
	 else {
	  if (epnt->lptr==NULL) { 			/* check lowlev ptr */
	     if (! (npnt=epnt->nodp1)) {
	 	ncfprintf (stderr,"check photrec node: can't find node %s\n",
		  nodstr(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d));
 	        execerror("Error","");
		break;
	     }
  	    cpnt = npnt->comptr;		/* check for comp at node */
  	    if (!cpnt) {			/* new node for comp */
    	      cpnt = makcomp(epnt,1/LARGERES,SMALLCAP,NULLVAL,NULLVAL);
    	      npnt->comptr = cpnt;		/* save comp pntr for node */ 
	      npnt->compnum = cpnt->num;  /* save comp number for node    */
              maklst(&cpnt->nodlst,(conn *)npnt);/* save node pntr for comp */
  	    }
	    trecep = makphotrec((photorec*)epnt,cpnt);
	    epnt->lptr = (elem *)trecep; /* save photrec ptr */
	    if (npnt->xloc>=LARGNOD) npnt->xloc = trecep->xloc;
	    if (npnt->yloc>=LARGNOD) npnt->yloc = trecep->yloc;
	    // if (npnt->zloc>=LARGNOD) npnt->zloc = trecep->zloc;
	  }
	 }
	 break;
 
      case GNDBATT:

	 rm = 1e6;			/* small leakage current */
	 ccap = 1.0;			/* one farad for battery cap */
	  if (! (npnt=epnt->nodp1)) {
		ncfprintf (stderr,"initcomp battery: can't find node %s\n",
		  nodstr(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d));
 	        execerror("Error","");
		break;
	  }
	  cpnt = npnt->comptr;		/* check for existing compartment */
	  if (cpnt) {			/* node already has compartment */
	   cpnt->rm += 1/rm;		/* add conductance 		*/
	   cpnt->cap += ccap;
	  }
	  else {
	    cpnt = makcomp(epnt,1/rm,ccap,NULLVAL,NULLVAL);
	    npnt->comptr = cpnt;	/* save comp pntr for node 	*/ 
	    npnt->compnum = cpnt->num;  /* save comp number for node    */
            maklst(&cpnt->nodlst,(conn*)npnt);	/* save node pntr for comp */
	  }
	  cpnt->ctype = epnt->ctype;	
	  cpnt->extv  = ((batt *)epnt)->v;/* set up compartment as battery */
	  cpnt->v     = ((batt *)epnt)->v;
	  {
	     int nvclamps;
	  nvclamps = cpnt->miscfl | VEXT;
	  nvclamps++;
	  if (nvclamps > VEXT) nvclamps=VEXT;
	  cpnt->miscfl |= (VEXT & nvclamps)|VBAT;/* set vclamp and batt flgs*/
  	  epnt->lptr = (elem *)cpnt; /* save ptr */
	  }
	 break;

      case ELEMENT:
	 break;

   }  /* switch */ 
  }  /* if (!epnt->newelem || epnt->modif) */	/* if element hasn't been processed yet */

  if (epnt->modif) { 		/* always erase modify elems */
    tepnt = epnt;
    epnt = epnt->next;
    if (tepnt->nodp1)unsetnodlst(tepnt->nodp1,tepnt); /* free nodpnt to elem */
    if (tepnt->nodp2)unsetnodlst(tepnt->nodp2,tepnt); /* free nodpnt to elem */
    delelem(tepnt);                     	      /* erase the element */
  }
  else {
    epnt->newelem = 0;		/* don't run this one again unless modified */
    epnt = epnt->next;
  }
 } /* for (epnt;; ) */

if (debug & NCCONVERT) ncfprintf (stderr,"initcomp end\n");
}

/*-------------------------------------*/

void chandens(elem *epnt, comp *cpnt, double area)
              
/* Add ion channels to a compartment.
   If channel "maxcond" is not specified, calculate
   channel max conductance by multiplying channel
   density by compartment membrane area.  Surface area
   is calibrated in cm2.
*/

{
    chan *chpnt;
    attrib *apnt;
    chattrib *chapnt;
    double ddensity;

  for (apnt=epnt->attpnt; apnt; apnt=apnt->attpnt) {

    chapnt = (chattrib *)apnt;
    switch (chapnt->ctype) {

   case VESNOISE:
    case CCHNOISE: break;	/* do nothing here if noise attrib, etc. */

     default:
      case  NA:
      case   K:
      case ClCa:
      case KCa:
      case  CA:
     case CGMP: 
			/* if "ndensity" has been specified, use it */
    if (chapnt->ndensity!=NULLVAL) { 
	chapnt->n = chapnt->ndensity*area*UM2CM2;
    }

    else {		/* look for "density" */
      switch (apnt->ctype) {
	case NA: ddensity = dnadens; break;
	case ClCa:
	case KCa:
	case K: ddensity = dkdens; break;
	case CGMP: 
        case CA: ddensity = dcadens; break;
	default: break;		/* do nothing for noise or other params */
      }	
      if (chapnt->density == NULLVAL) {
	 if (chapnt->n==NULLVAL) {	/* if N specified, use it */
	   chapnt->density = ddensity;	/* default density is last resort */
           chapnt->maxcond = chapnt->density * area;
	}
      }
     else chapnt->maxcond=chapnt->density * area; /* use density if specified */
    }
    switch (apnt->ctype) {
	case NA: 
          chpnt = makchan(apnt,cpnt,(comp*)NULL, epnt->nocondens); 
	  break;
	case K: 
	case ClCa: 
	case KCa: 
          chpnt = makchan(apnt,cpnt,(comp*)NULL, epnt->nocondens); 
	  break;
	case CGMP: 
        case CA:
          chpnt = makca(epnt,(cattrib*)apnt,cpnt,area, epnt->nocondens); 
	  break;
	default: break;		/* do nothing for noise or other params */

      }	    /* switch (ctype) */
    if (epnt->nocondens) apnt->lptr = (elem *)chpnt; /* save chan ptr */

   }  /* switch (chapnt->ctype) */

  }  /* for (apnt= ; ; ) */
} 

/*------------------------------------*/

void checkcomp(elem *epnt, comp **pnt1, comp **pnt2)
               
/* find the compartments that the synapse connects to,
   and make new ones if necessary.
*/

{
    node  *npnt;
    comp  *cpnt,*cpntn,*cpntl;

if (debug & NCCONVERT && debugz & 4) ncfprintf (stderr,"checkcomp \n");

 				/* find first node */
  if (! (npnt=epnt->nodp1)) {
	ncfprintf (stderr,"Checkcomp: elem '%s', ",findsym(epnt->ctype));
	ncfprintf (stderr,"can't find node %s\n",
	  nodstr(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d));
        execerror("Error: incorrect node spec.","");
	return;
  }
  cpntn = npnt->comptr;		/* check for compartment at end */
  if (cpntn) {			/* node already has compartment */
           cpntl = cpntn;
  }
  else {			/* this is first comp for node */
    cpnt = makcomp(epnt,1/LARGERES,SMALLCAP,NULLVAL,NULLVAL);
    npnt->comptr = cpnt;	/* save comp pntr for node */ 
    npnt->compnum = cpnt->num;  /* save comp number for node    */
    maklst(&cpnt->nodlst,(conn*)npnt);	/* save node pntr for comp */
    cpntl = cpnt;
  }

 				/* next, find second node */
  if (! (npnt=epnt->nodp2)) {
	ncfprintf (stderr,"checkcomp: elem '%s',\n",findsym(epnt->ctype));
	ncfprintf (stderr,"      Found first node: %s,\n",
	  nodstr(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d));
	ncfprintf (stderr,"can't find second node: %s\n",
	  nodstr(epnt->node2a,epnt->node2b,epnt->node2c,epnt->node2d));
        execerror("Error: incorrect node spec.","");
	return;
  }
  cpntn = npnt->comptr;		 /* check for compartment at end */
  if (!cpntn) {			 /* this is first comp for node */
    cpnt = makcomp(epnt,1/LARGERES,SMALLCAP,NULLVAL,NULLVAL);
    npnt->comptr = cpnt;	 /* save comp pntr for node */ 
    npnt->compnum = cpnt->num;  /* save comp number for node    */
    maklst(&cpnt->nodlst,(conn*)npnt);	/* save node pntr for comp */
    cpntn = cpnt;
  }
  *pnt1 = cpntl;		 /* return the compartment pointers */
  *pnt2 = cpntn;

  if (debug & NCCONVERT && debugz & 4) ncfprintf (stderr,"checkcomp end\n");
}

/*------------------------------------*/

char *nodstr (int a, int b, int c, int d)
             

/* make string for printing node in readable format */

{
   static char nodst[48] = {0};

 if (a == NULLVAL) a = NULND;
 if (b == NULLVAL) b = NULND;
 if (c == NULLVAL) c = NULND;
 if (d == NULLVAL) d = NULND;
 if (a>=0) {
    if (b>=0) {
       if (c>=0) {
          if (d>=0)
             sprintf (nodst,"%3d %3d %3d %3d",a,b,c,d);
           else
             sprintf (nodst,"%3d %3d %3d",a,b,c);
       }
       else sprintf (nodst,"%3d %3d",a,b);
    }
    else sprintf (nodst,"%3d",a);
 }
 else sprintf (nodst,"%3d %3d %3d %3d",a,b,c,d);
 return (nodst);
}

/*------------------------------------*/

void compcon(void)

/* calculate the total unvarying conductance for all compartments. */

{
    static double tcond, tcondm;
    extern double timinc;
    comp *cpnt;
    conn *conptr;
    conlst *lpnt;

if (debug & NCCONVERT) ncfprintf (stderr,"compcon\n");

 for (cpnt=compnt; cpnt; cpnt=cpnt->next)
  {
   cpnt->tcond = 0;			/* total constant conductance */
   cpnt->tcondm = 0;			/* total membrane conductance */
   for (lpnt=cpnt->clst; lpnt; lpnt=lpnt->next) 
    {					/* calc total conductance through */
     conptr = lpnt->conpnt;		/*   all connections       */
     if (! conptr) break;
     tcond=tcondm=0.0;
     switch (conptr->ctype)
      {
	case AXIALRES:
     	  tcond = conptr->conduct;		/* neighbor's axial resist */
		    break;
 	//case RESISTOR:
     	//  tcond = conptr->conduct;		/* neighbor's axial resist */
	//  tcondm = tcond;
		    break;
	//case LOAD:				/* load */
     	//  tcond  = ((load *)conptr)->conduct;	/* conductance of load */
	//  tcondm = tcond;
		    break;

	case BATT:				/* series battery */
     	  tcond  = BATTCOND;			/* battery conductance */
     	  tcondm = tcond;			/* battery conductance */
		    break;

	case CAP:
     	  tcond += conptr->conduct/timinc;	/* capacitor reactance */
		    break;
	case GJ:
	case ROD:
	case CONE:
	case CHR:
	case VTRANSDUCER:
	case ITRANSDUCER:
	case SYNAPSE:
	case DIODE:
	default:
		    break;
      }
     cpnt->tcond += tcond;			/* total constant conductance */
     cpnt->tcondm += tcondm;			/* total membrane conductance */
    }
   cpnt->tcond += cpnt->rm;			/* total constant conductance */
   cpnt->tcondm += cpnt->rm;			/* total membrane conductance */
   if (cpnt->miscfl & VBAT) 			/* if comp is battery */
     cpnt->tcond += BATTCOND;			/* battery conductance */
   cpnt->k = timinc / cpnt->cap;		/* constant for node */
   cpnt->relax = relax;				/* set initial over-relax */
   cpnt->vest = cpnt->v;			/* set initial volt est */
   if  (djnoise) {
	if (cpnt->jnoise) cpnt->jnoise *= djnoise;
	else              cpnt->jnoise = djnoise;
   }
  }
  initsyn();					/* initialize all synapses */

if (debug & NCCONVERT) ncfprintf (stderr,"compcon end\n");
}

/*------------------------------------*/

void initsyn(void)

/* initialize all synaptic connections' */
/*  internal activity according to the presynaptic */
/*  voltage */

{
   synap *spnt;

 for (spnt=synpnt; spnt; spnt= (synap*)spnt->next)
   setsynv(spnt, spnt->comp1->v-spnt->comp1->vext);
}

/*------------------------------------*/

void calcnod (double Ri, double Rm, double Cap, double d, double clen, \
		double *rri, double *rrm, double *rcap)
                                           

/* Calculate the conductance and capacitance for an incremental
node.  Ri, Rm and Cap (on input) are scaled to centimeters.
clen is length of compartments.
*/

{
  double a;
  double ri;			/* normalized to ohms / meter */
  double rm;			/* normalized to ohms * meter */
  double cap;			/* normalized to uF / meter */

 clen *= UMMSCAL;		/* normalize comp len to meters */
 a = d * UMMSCAL / 2;
 ri = Ri * CMMSCL / (MPI * a * a);
 rm = Rm * CMMSCL2 /  (2 * MPI * a);
 cap = Cap * 2 * MPI * a  / CMMSCL2;
 *rrm = rm / clen;
 *rri = ri * clen;
 *rcap = cap * clen;

/*
if (1)
 {
  ncfprintf (stderr,"a %g Ri %g Rm %g Cap %g clen %g\n",a,Ri,Rm,Cap,clen);
  ncfprintf (stderr,"a %g ri %g rm %g cap %g\n",a,ri,rm,cap);
  ncfprintf (stderr,"rri %g rrm %g rcap %g\n",*rri,*rrm,*rcap);
 }
*/

}

/*------------------------------------*/

void condense(void)

/*   This algorithm condenses closely coupled compartments
   with their neighbors.  This function is useful because
   the algorithm of "initcomp()" (which generates compartments
   from the element list) leaves a minimum of 2 compartments
   per cable segment, even if the 2 compartments are closer
   (smaller ri/rm) than specified by "complam".

     The algorithm: For all compartments, check whether 
   (for each of their axial resistive connections)
   sqrt(rm/ri) > lamcrit.  If it isn't (i.e compartment
   is more closely coupled than criterion), then condense the
   compartment into its neighbor.  Add gm and cm, and add
   ri in proportion to the gm's of the two compartments.
   Move any attributes and synapses to the neighbor, and
   erase the connection between the two compartments.  Also,
   the node that points to the compartment is redirected to
   the neighbor as well.  This operation, of course, may
   leave the neighbor compartment closely coupled to another
   neighbor so that it also needs to be condensed.  Since
   the neighbor may already have been passed over, the whole
   process needs to be repeated.  Therefore, do this condensation
   operation on all compartments, and go back to the start of
   the compartment list to iterate until no more changes.

   Since all compartments in a cable (except ends) are always
   made with proper size rm/ri, we only need to check those 
   compartments associated with nodes.  Therefore, the main 
   loop of the algorithm is done for nodes, not compartments.
   Thus it is fairly easy to move the nodes' compartment pointers
   from the old erased compartment to the new compartment.
*/

{
    int change,done,n,oldcomp;
    comp *cpnt,*ocpnt,*tcpnt;
    conn *conpnt;
    node *npnt,*tnpnt;
    cacomp *capnt,*ocapnt;
    static conlst *lpnt,*nlpnt,*olpnt,*tlpnt;
    double ri,rm,orm,lam,cond,crit,varycrit;
    double ocond,tcond,minlam,ominlam,maxcplam;

#ifdef DEBUG
 if (debug & NCCONVERT && debugz && 8) ncfprintf (stderr,"condense\n");
#endif

maxcplam= 0;
minlam=1e10;

if (lamcrit==0) return;			/* no condense on zero lamcrit */
if (lamcrit<0) lamcrit = - lamcrit;
					/* condense lowest lambda comps first */
					/* by varying crit from .2 to 1 */
for (varycrit=.2; varycrit<=1; varycrit+= 0.2) {

 crit = lamcrit * varycrit;		/* variable criterion for lambda */
#ifdef DEBUG
    if ((debug & NCCONVERT) && (debugz & 8))
     ncfprintf(stderr,"condense: lambda crit %g...\n", crit);
#endif
 for (change=1; change; ) {	/* iterate last crit while changes made */
  change = 0;
  ominlam = minlam;
  minlam=1e10;
  for (npnt=nodepnt; npnt; npnt=npnt->next) {	/* for all nodes */
    cpnt = npnt->comptr;			/* get node's compartment */
    if (!cpnt) {
  ncfprintf (stderr,"condense: node %s has no compartment. continuing...\n",
	   nodstr(npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4));
        continue;			/* skip null node */
    }
    oldcomp = 0;
    for (done=0,lpnt=cpnt->clst; lpnt && !done; lpnt=lpnt->next) {
       if (!(conpnt=lpnt->conpnt)) {
	  ncfprintf (stderr,"condense: missing connection in comp %d\n",
						cpnt->num); /* */
	  continue;
       }
       if (conpnt->ctype==AXIALRES) {	/* check all axial connections */

     if (conpnt->comp1==cpnt) ocpnt = conpnt->comp2;  /* find other comp */
     else                     ocpnt = conpnt->comp1;

    rm  =  cpnt->rm;			/* save comp rm */
    orm = ocpnt->rm;			/* save neighbor's comp rm */
    ri = conpnt->conduct;
    if (ri==0) ri = 1e-30;
    tcond = orm/(rm+orm) * 1/ri;	/* this compartment's share of ri */
    ocond =  rm/(rm+orm) * 1/ri;	/* other compartment's share of ri */
    lam = sqrt (rm/ri);
    ominlam = minlam; 
    if (lam < minlam) minlam = lam;	/* remember minimum lambda */

    if (cpnt->cplam > maxcplam) maxcplam = cpnt->cplam;	/* maximum complam */

#ifdef DEBUG
    if ((debug & NCCONVERT) && (debugz & 8))
    ncfprintf(stderr,"node %s, conn from comp %d to %d lam %-4.2g cplam %g\n",
			nodstr(npnt->nodenm1,npnt->nodenm2,
			       npnt->nodenm3,npnt->nodenm4),
			cpnt->num,ocpnt->num,lam,cpnt->cplam); /* */
#endif
//if (0)
    if (cpnt==ocpnt) {
#ifdef DEBUG
//    if ((debug & NCCONVERT) && (debugz & 8))
    ncfprintf(stderr,
"node %s, conn from comp %d to %d, type %d, cond %g is in tight loop ****\n",
			nodstr(npnt->nodenm1,npnt->nodenm2,
			       npnt->nodenm3,npnt->nodenm4),
			cpnt->num,ocpnt->num,conpnt->ctype,conpnt->conduct);
			 /* */
#endif
        lpnt->conpnt = (conn *)NULL;   /* remember we're not erasing conlst */
        efree (conpnt);		 /* erase connection from comp to itself */
	continue;
    }
    if (lam < cpnt->cplam * crit) {		/* if this comp is too close */
						/* proceed with condensation */
#ifdef DEBUG
    if ((debug & NCCONVERT) && (debugz & 8))
     ncfprintf(stderr,"condensing comp %d lam %-4.2g into comp %d\n",
			cpnt->num,lam,ocpnt->num); /* */
#endif
       if (varycrit > .99) change = 1;		/* only repeat last crit */
       ocpnt->rm = ocpnt->rm + rm;		/* condense rm("gm"), cap */
       ocpnt->cap = ocpnt->cap + cpnt->cap;
       minlam = ominlam;			/* set minlam to prev. value */
						/* because comp is condensed */

	/* At this point, we've decided to condense a compartment */
	/*  into its neighbor. First, we move all the compartment's */
	/*  connections (except the one already connecting it to */
	/*  the neighbor we're condensing) to connect to the compartment */
	/*  we're condensing into.  Once we've moved the connections, */
	/*  we delete the old connection pointers if new ones have been *
	/*  made (i.e. for channels). Otherwise all the old connection */
	/*  pointers get moved to the neighbor's conlst. */

	/* Special calcium stuff: */

	if (capnt=cpnt->capnt) {  	    /* if ca comp */
	   if (ocapnt=ocpnt->capnt) {       /* use "other" ca comp */
	      if (capnt->cashell > ocapnt->cashell) {
	         /* if this cacomp's number of shells is greater, use */	
		  ocapnt->cashell = capnt->cashell;
		  ocapnt->caoshell = capnt->caoshell;
	 	  efree (ocapnt->cais);	
		  ocapnt->cais = capnt->cais;	
		  capnt->cais = (double *)NULL;
	 	  efree (ocapnt->caos);	
		  ocapnt->caos = capnt->caos;	
		  capnt->caos = (double *)NULL;
		  if (capnt->cab) {
		     efree (ocapnt->cab);	
		     ocapnt->cab = capnt->cab;	
		     capnt->cab = (double *)NULL;
		     ocapnt->cabt = capnt->cabt;
		     ocapnt->cabf = capnt->cabf;
		     ocapnt->cabr = capnt->cabr;
		  }
		  if (capnt->vm2>0) {			/* using CICR */
        	     ocapnt->cas = capnt->cas;
        	     ocapnt->vm2 = capnt->vm2;
        	     ocapnt->vm3 = capnt->vm3;
        	     ocapnt->ncicr = capnt->ncicr;
        	     ocapnt->mcicr = capnt->mcicr;
        	     ocapnt->pcicr = capnt->pcicr;
        	     ocapnt->kacicr = capnt->kacicr;
        	     ocapnt->kfcicr = capnt->kfcicr;
        	     ocapnt->krcicr = capnt->krcicr;
        	     ocapnt->k1cicr = capnt->k1cicr;
        	     ocapnt->k2cicr = capnt->k2cicr;
        	     ocapnt->c1cicr = capnt->c1cicr;
		  }
		  if (capnt->bip3>0) {			/* using IP3 release */
		     ocapnt->cas2 = capnt->cas2;
		     ocapnt->ip3i = capnt->ip3i;
		     ocapnt->bip3 = capnt->bip3;
		     ocapnt->vip3 = capnt->vip3;
		     ocapnt->b2ip3 = capnt->b2ip3;
		     ocapnt->a2ip3 = capnt->a2ip3;
		     ocapnt->a3ip3 = capnt->a3ip3;
		     ocapnt->oip3 = capnt->oip3;
		     ocapnt->v2ip3 = capnt->v2ip3;
		     ocapnt->v3ip3 = capnt->v3ip3;
		     ocapnt->v4ip3 = capnt->v4ip3;
		     ocapnt->k3ip3 = capnt->k3ip3;
		     ocapnt->d1ip3 = capnt->d1ip3;
		     ocapnt->d2ip3 = capnt->d2ip3;
		     ocapnt->d3ip3 = capnt->d3ip3;
		     ocapnt->d4ip3 = capnt->d4ip3;
		     ocapnt->mip3 = capnt->mip3;
		     ocapnt->hip3 = capnt->hip3;
		     ocapnt->mtypeip3 = capnt->mtypeip3;
		  }
	      }
	      redocacomp(ocpnt, capnt->area);	/* redo cacomp in ocpnt */
	      delcacomp (capnt);   /* erase ca comp */
	      cpnt->capnt = (cacomp*)NULL; 
           }
	   else {			 	/* use "this" ca comp */
	     ocpnt->capnt=ocapnt=capnt;  	/* move cacomp pointer */
	     ocapnt->comp1 = ocpnt;		/* ca comp pntr to v comp */
	   }
	}

	for (n=0,olpnt=tlpnt=(conlst*)NULL,nlpnt=cpnt->clst;
		 nlpnt;
		 nlpnt=nlpnt->next, (tlpnt?efree(tlpnt),1:0),
					tlpnt=(conlst*)NULL) {
	   if (!(conpnt=nlpnt->conpnt)) continue;
           switch (conpnt->ctype) {

             case AXIALRES:		
		if (nlpnt==lpnt) {	/* don't want to move orig conn */
	          if (olpnt) olpnt->next = nlpnt->next;  
		  else cpnt->clst = nlpnt->next;
		  tlpnt = nlpnt;	/* remember this one to free later */
		}
		/* It is unlikely that there might be 2 axial conns */
		/* between the 2 compartments being condensed,  */
		/* but it is possible. */
		/* In this case, there should no big problem because */
		/* this loop will move the second connection's */
		/* pointer so the connection points from the comp */
		/* to the same comp. This will, of course, forget */
		/* about the second conn. But if such a conn is */
		/* important, then one should not condense at all. */ 

		else {				/* move other axial res */
     	          if (conpnt->comp1==cpnt) conpnt->comp1 = ocpnt;
     	          else                     conpnt->comp2 = ocpnt;
		  n++;				/* count up other axial conns */
		}
	       break;

	     case GJ:
	     case RESISTOR:
     	       if (conpnt->comp1==cpnt) conpnt->comp1 = ocpnt;
     	       else                     conpnt->comp2 = ocpnt;
	       n++;			/* count up resistive conns */
	       break;
 
	     case SYN2: 	/* move channel to other comp */
	     case GABA:		/*  then erase conn ptr after done */
	     case GLY:
	     case NA:
	     case ClCa:
	     case KCa:
	     case K:
	     case CGMP:	
	     case AMPA:
	     case KAINATE:
	     case NMDA:
	     case CA:

	       if (olpnt) olpnt->next = nlpnt->next; 
	       else cpnt->clst = nlpnt->next;
	       tlpnt = nlpnt;		/* free conn pointer later (above) */
	       if (!addchan (ocpnt,(chan*)conpnt)){ /* add chan to other comp*/
		  ncfprintf (stderr,"can't add chan to other comp.\n");
	        }
		break;

	     case DIODE:
	     case SYNAPSE:
     	       if (((synap *)conpnt)->comp1==cpnt)
				 ((synap *)conpnt)->comp1 = ocpnt;
     	       else              ((synap *)conpnt)->comp2 = ocpnt;
		break;

	     case ROD:
	     case CONE:
	     case CHR:
	     case VTRANSDUCER:
	     case ITRANSDUCER:
     	       if (((photrec *)conpnt)->comp1==cpnt)
				 ((photrec *)conpnt)->comp1 = ocpnt;
		break;

	     case LOAD:
     	       if (((load *)conpnt)->comp1==cpnt)
				 ((load *)conpnt)->comp1 = ocpnt;
		break;

             case BATT:
	     case CAP:
     	       if (conpnt->comp1==cpnt) conpnt->comp1 = ocpnt;
     	       else                     conpnt->comp2 = ocpnt;
	       break;

	     case BUF:
     	       if (((dbuf *)conpnt)->comp1==cpnt)
				 ((dbuf *)conpnt)->comp1 = ocpnt;
	       else		 ((dbuf *)conpnt)->comp2 = ocpnt;
	       break;

	     default: {
                        conlst *xlpnt;
		        node *xnpnt;

	       ncfprintf (stderr,"condense: unknown connection type %d comp %d\n",
						conpnt->ctype,cpnt->num);
 		ncfprintf (stderr,
	"This error probably originates from a loop that was condensed.\n");
 		ncfprintf (stderr,
	"Please check your model for small loops between cables.\n");
 		ncfprintf (stderr,
	"The offending node is probably one of the following:\n");

        for (xlpnt=cpnt->nodlst; xlpnt; xlpnt=xlpnt->next) {
           ncfprintf (stdout,"node "); 
           if (xlpnt->conpnt) xnpnt = (node *)xlpnt->conpnt;
           if (xnpnt->nodenm4 > NULLVAL)
                 ncfprintf (stdout,"%-2d %-2d %-2d %-2d",
                  xnpnt->nodenm1,xnpnt->nodenm2,xnpnt->nodenm3,xnpnt->nodenm4);
           else if (xnpnt->nodenm3 > NULLVAL)
                 ncfprintf (stdout,"%-2d %-2d %-2d",
                        xnpnt->nodenm1,xnpnt->nodenm2,xnpnt->nodenm3);
           else if (xnpnt->nodenm2 > NULLVAL)
                 ncfprintf (stdout,"%-4d %-3d",xnpnt->nodenm1,xnpnt->nodenm2);
           else
                 ncfprintf (stdout,"%-4d    ",xnpnt->nodenm1);
           ncfprintf (stdout,"\n"); 
         }

        }  /* default */

	break;
	   }   /* switch (conpnt->ctype) */

	   if (nlpnt!=tlpnt)            /* if not erased */
	      olpnt = nlpnt;		/* remember previous one in list */

	}  /* for (nlpnt=cpnt->clst;;) */
	
		   /* Special nt, cyca/g comp stuff: */

       if (cpnt->ntpnt) {	/* if we need to condense a nt comp. */
	      ntcomp *npt;				

	  for (npt=cpnt->ntpnt; npt; npt=npt->ntpnt) { 
              makaddntcomp(ocpnt,npt->ctype,npt->val,npt->n);
	  }
	  delntcomp(cpnt->ntpnt);	/* delete */
       }

		/* now add comp's share of ri to its other axial connections */
		/* partition ri evenly between other conns */

       if (n)			/* if there are any other axial conns */
	for (olpnt=(conlst*)NULL,nlpnt=cpnt->clst; nlpnt; nlpnt=nlpnt->next) {
	   if (!(conpnt=nlpnt->conpnt)) continue;
           switch (conpnt->ctype) {
             case AXIALRES:			/* add ri to other axials */
             case GJ:				/* and even to gj's */
		if (ri==0.0) break;
		if ((cond=conpnt->conduct)==0.0) break;
		conpnt->conduct = 1 / (1/cond + tcond/n);
	   }  /* switch */
	}  /* for */
			/* erase the connection, count resistors  */ 
			/* for "other" compartment (ocpnt) */

       for (n=0,tlpnt=olpnt=(conlst*)NULL,nlpnt=ocpnt->clst;
			 nlpnt;
			 nlpnt=nlpnt->next, (tlpnt?efree(tlpnt),1:0),
					tlpnt=(conlst*)NULL) {
	   if (!(conpnt=nlpnt->conpnt)) continue;
	   if (conpnt->comp1==ocpnt) tcpnt = conpnt->comp2;
	   else                      tcpnt = conpnt->comp1;
           switch (conpnt->ctype) {
             case AXIALRES:
	       if (tcpnt==cpnt) {     		/* if this is right one */

		/* patch the connection list pointers */
		/*  and erase the connection */

      		  if (conpnt->last) conpnt->last->next = conpnt->next;
      		  else connpnt = conpnt->next;
      		  if (conpnt->next) conpnt->next->last = conpnt->last;
      		  else connend = conpnt->last;
	          efree (conpnt);		/* erase the connection */
		  cumconn--;

	          if (olpnt) olpnt->next = nlpnt->next;
		  else ocpnt->clst = nlpnt->next;  /* erase conn pointer */
		  tlpnt = nlpnt;
	       }
	       else n++;			/* count other axial conn's */
	       break;

	     case RESISTOR:
	     case GJ:
		n++;
		break;
	    }
	    if (nlpnt!=tlpnt)                   /* if not erased */
	        olpnt = nlpnt;			/* remember previous one */

 	} /* for (nlpnt=ocpnt->clst;;) */

	  /* now add "other" comp's share of ri to its other axial conns */
	  /* partition ri evenly between other conns */

	if (n) 				/* if neighbor has any axial conns */
	 for (nlpnt=ocpnt->clst; nlpnt; nlpnt=nlpnt->next) {
	   if (!(conpnt=nlpnt->conpnt)) continue;
           switch (conpnt->ctype) {
             case AXIALRES:		/* add extra ri to neighbor's conns */
             case GJ:			/*  and even to gap junctions */
			if (ri==0.0) break;
			if ((cond=conpnt->conduct)==0.0) break;
			conpnt->conduct = 1 / (1/cond + ocond/n);
			break;
	   } /* switch */
	 }  /* for (nlpnt) */

			/* find the end of the neighbor's conn list again */

	for (olpnt=(conlst*)NULL,nlpnt=ocpnt->clst; nlpnt; nlpnt=nlpnt->next) 
		olpnt = nlpnt;

	/* now move the remaining conn list from the erased compartment */
	/*  to the new (i.e. neighbor) compartment */

	if (olpnt) olpnt->next = cpnt->clst;
	else ocpnt->clst = cpnt->clst;
	
	/* move all node compartment pointers from old comp to new */

      for (nlpnt=cpnt->nodlst; nlpnt; nlpnt=nlpnt->next) {
	  tnpnt = (node *)nlpnt->conpnt;
	  if (tnpnt->comptr == cpnt) {
		tnpnt->comptr = ocpnt;
		tnpnt->compnum = ocpnt->num;
	  }
      }
		/* find the end of the neighbor's node list */

	for (olpnt=(conlst*)NULL,nlpnt=ocpnt->nodlst; nlpnt; nlpnt=nlpnt->next) 
		olpnt = nlpnt;

		/* now move the node list from the erased compartment */
		/*  to the new (i.e. neighbor) compartment */

	if (olpnt) olpnt->next = cpnt->nodlst;
	else ocpnt->nodlst = cpnt->nodlst;

		/* finally, patch the compartment list pointers */
		/*  and erase the compartment */

      if (cpnt->last) cpnt->last->next = cpnt->next;  /* patch pointers */
      else compnt = cpnt->next;
      if (cpnt->next) cpnt->next->last = cpnt->last;
      else compend = cpnt->last;

      delhashcomp(cpnt);		/* delete from the hash list */ 
      efree (cpnt);			/* erase the compartment */
      cumcomp--;
      done = 1;		/* stop when a comp and its conns have been erased */

      }	      /* if (lam;;) */
     }	  /* if (conpnt->ctype==AXIALRES) */
    if (done) break;	/* go on to check next compartment */

    oldcomp = cpnt->num;
   }   /*  for (done=0,lpnt;; ) */
  }	  /* for (npnt=nodepnt;;) */
#ifdef DEBUG
    if ((debug & NCCONVERT) && (debugz & 8))
     if (change) ncfprintf(stderr,"condense: starting over...\n");
#endif
 }      /* for (change=0;; ) */

 if (minlam > (maxcplam*lamcrit)) break;	/* stop when obviously done */
}   /* for (crit=0.2;;) */

if (lamcrit > 1e-6) {
  for (cpnt=compnt; cpnt; cpnt=cpnt->next) {	/* compress comps' channels */
    condense_comp_chan(cpnt);
  }
}

#ifdef DEBUG
  if (debug & NCCONVERT && debugz & 8) ncfprintf (stderr,"condense end\n");
#endif

}     /* condense */


/*------------------------------------*/

comp *othercomp (comp *pnt,conn *cpnt)

/* Return pointer to other comp that connection points to. 
   For reliability, we should modify code in "condense()" above 
    to use this.
*/

{
  if (!pnt || !cpnt) {
      ncfprintf (stderr,"Othercomp: invalid pointer\n");
      return (comp*)NULL;
  }
  if (cpnt->comp1==pnt) return cpnt->comp2;
  else return cpnt->comp1; 
}

