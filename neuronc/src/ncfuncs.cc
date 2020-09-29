
/* functions for C++ NC script language */

#include "ncsub.h"
#include "ncelem.h"
#include "control.h"
#include "nc.h"
#include "y.tab.h"
#include "ncomp.h"
#include "ncplot.h"
#include "gprim.h"
#include "stim.h"


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <setjmp.h>
#include <sys/times.h>
#include <sys/types.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
//#include "stdplt.h"

#ifdef __cplusplus
}
#endif

#include "ncio.h"

extern plotfr plotnod[PLOTNODSIZ];  /* holds nodes to be recorded */
extern int plotnum;                /* number of plots to be displayed */
extern int runyet;
extern char *rootframe;             /* name of root plot frame */
extern int debug;
extern int debugz;
extern int setdebug;
extern int setdebgz;
extern int setdisp;
extern int *istat;
extern int makestim;            /* = 1 -> says we're running "stim", not "nc" */
extern int prcomps;		/* =1 -> says create node file in addition to display nodes */
extern char *progname;
extern char *infile;		/* input file name, used by "plotmod" */
extern char *einfile;
extern elem *elempnt;
extern double dsize;

extern FILE *stimout;
extern FILE *compout;


#define SECPERMIN 60.0          /* number of seconds per minute */

/* #define TIMECALIB 60.0       /* times() system call gives 1/60  sec */
/* #define TIMECALIB 100.0      /* times() system call gives 1/100 sec */
/* #define TIMECALIB 1024.0     /* times() system call gives 1/1024 sec */

#ifdef DEC_ALPHA
#define TIMECALIB 1024.0        /* times() system call gives .9765625 msec */
#else
#define TIMECALIB 100.0         /* times() system call gives 1/100 sec */
#endif

extern struct tms timebuf;
extern struct tms *timepnt;
extern double runmin;
extern double timecalib;
extern double startclk,startu,totmin;

recstim *makvstim(double start, int nodnm1, int nodnm2, int nodnm3, int nodnm4,
                        double value, const char *action);
recstim *makrstim(double start, int recnm1, int recnm2, int recnm3, int recnm4,
                double inten, double wavel, double mask, const char *action);

extern elem *elempnt;           /* pointer to current element */

node *maknod(nodeint nodea, nodeint nodeb, nodeint nodec,
        nodeint noded);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, nodeint node1d, const char *s);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, const char *s);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c);
int findphotrec (int nodenm1, int nodenm2, int nodenm3, photrec **rpnt, const char *str);
int findphotrec (int nodenm1, int nodenm2, int nodenm3, int nodenm4, photrec **rpnt, const char *str);

elem *makelem(int etype, elem *epnt);
void eninstall (elem *nelem);
void checkelemnode(elem *elpnt);
void plotfunc(double (*func)(double v, double t));
void plotfunc(double (*func)(double v, double v2, double t));
void plotfunc(double (*func)(double v, double v2, double v3, double t));
void plotfunc(double (*func)(double v, double v2, double v3, double v4, double t));
void plotfunc(double (*func)(double v, double v2, double v3, double v4, double v5, double t));
void plotfunc(double (*func)(double v, double v2, double v3, double v4, double v5, double v6, double t));
void plotfunc(double (*func)(double v, double v2, double v3, double v4, double v5, double v6, double v7, double t));
void plotarr(double *var, int maxindex);
void plotvar(double *var);
void plotval (double plotvl, int plotnum);
void plotval (double plotvl, double plotvl2, int plotnum);
void plotval (double plotvl, double plotvl2, double plotvl3, int plotnum);
void plotval (double plotvl, double plotvl2, double plotvl3, double plotvl4, int plotnum);
void plotval (double plotvl, double plotvl2, double plotvl3, double plotvl4, double plotvl5, int plotnum);
void plotval (double plotvl, double plotvl2, double plotvl3, double plotvl4, double plotvl5, 
	      double plotvl6, int plotnum);
void plotval (double plotvl, double plotvl2, double plotvl3, double plotvl4, double plotvl5, 
	      double plotvl6, double plotvl7, int plotnum);
void onintr(int i);
void fpecatch(int i);
char *print_version(double version);
void ncleanup();
void initpl(int charind);
void init();
void initchan();
void findconnect(int donefl);
void runsim(double time, int run);
void prtime(FILE *stimout);
void sortstimfile();
void ncmdlin(int argc, char **argv);
void bexit();
char *prnode(int n1, int n2, int n3, int n4);
char *prnode(node *npnt);
void execerror (const char *s, const char *t);
elem *findelem (int num);
void setnodlst(node *nodp, elem *epnt);
void unsetnodlst(node *nodp, elem *epnt);
double dist3d(node *n1, node *n2);
chattrib *chanattr(elem *elpnt, short int ctype);
void copyattr(attrib *src, attrib *dest);
double getnt (comp *cptr,int nt);
void vinit(void);
char *emalloc(unsigned int n);
int gausnn(double mean, double stdev, double density, double ms, int grseed,
        double framex, double framey, double xcent, double ycent,
        int numcells, double **xarr, double **yarr, int first_center, int filout, int textfl,
        int info);
void ncdispe (int elemnum, double zrange1, double zrange2,
        int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
        double vmax, double vmin, double dscale, int hide, int cmap);
void ncdispe (elem *epnt, int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
        double vmax, double vmin, double dscale, int hide, int cmap);
void dispstim(double stime, double dscale, int cmap, double flux_max, double flux_min);
void dispstim(double sttime, double stime, double dscale, int cmap, double flux_max, double flux_min);
void setptrs(void);
node *erasenode(node *npnt);
void eraseelem(elem *epnt);
double recchan(int cnum, int nf, int pval);
double recchan(elem *elpnt, int nf, int pval);
double rec_node_chan(node *npnt, int nf, int pval);
double recsynapse(int snum, int nf);
double recsynap(synap *spnt, int nf, int snum);
double reccable(int elemnum, double fracdist);
void efree(void *p);
double system_speed(void);
chantype *getchantype(int ctype, int cnum);
double rechan(conn *cpnt, int nf, int pval, int cnum);
elem *getepnt (nodeint nodenm1, nodeint nodenm2);

#ifdef __cplusplus
extern "C" {
#endif

char *gfrname(void);

#ifdef __cplusplus
}
#endif


/*---------------------------------------------------------------------*/

node *nd(int na, int nb, int nc, int nd)

{
    int found;
    node *npnt;

  npnt=findnode(na,nb,nc,nd,NULL); 		/* look for existing node */
  if (npnt==NULL) {				/* if not found, then make it */
    npnt=maknod(na,nb,nc,nd);	
  }
  return npnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *nd(int na, int nb, int nc)

{
    int found;
    node *npnt;

  npnt=findnode(na,nb,nc);	 		/* look for existing node */
  if (npnt==NULL) {				/* if not found, then make it */
    npnt=maknod(na,nb,nc,NULLVAL);	
  }
  return npnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *nd(int na)

{
  return (nd (na,NULLVAL,NULLVAL,NULLVAL));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *nd(int na, int nb)

{
  return (nd (na,nb,NULLVAL,NULLVAL));
}


/*---------------------------------------------------------------------*/

node *ndn(int na, int nb, int nc, int nd)
{
    int found;
    node *npnt;
    char sbuf[50];

  npnt=findnode(na,nb,nc,nd,"ndn"); 		/* look for existing node */
  if (npnt==NULL) {				/* if not found, don't make it */
     sprintf (sbuf,"Missing node %d %d %d; ",na,nb,nc);
     execerror (sbuf,"stopping...");
  }
  return npnt;
}

/*---------------------------------------------------------------------*/

node *ndn(int na, int nb, int nc)
{
    int found;
    node *npnt;
    char sbuf[50];

  npnt=findnode(na,nb,nc); 			/* look for existing node */
  if (npnt==NULL) {				/* if not found, don't make it */
     sprintf (sbuf,"Missing node %d %d %d; ",na,nb,nc);
     execerror (sbuf,"stopping...");
  }
  return npnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *ndn(int na)

{
  return (ndn (na,NULLVAL,NULLVAL));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *ndn(int na, int nb)

{
  return (ndn (na,nb,NULLVAL));
}

/*---------------------------------------------------------------------*/

node *nde(int na, int nb, int nc, int nd)
{
  return(findnode(na,nb,nc,nd,NULL)); 		/* look for existing node, */
}						/* no error messages */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *nde(int na, int nb, int nc)
{
  return(findnode(na,nb,nc));			/* look for existing node */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *nde(int na, int nb)
{
  return(findnode(na,nb,NULLVAL,NULL));		/* look for existing node */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *nde(int na)
{
  return(findnode(na,NULLVAL,NULLVAL,NULL)); 	/* look for existing node */
}

/*---------------------------------------------------------------------*/

node *ndt(int na, int nb, int nc, int nd)

/* Make node numbers but don't find or make node. */
/* Must call efree() after */

{
    node *npnt;

 if ((npnt=(node *)emalloc(sizeof(node))) == NULL) {
     ncfprintf (stderr,"no space left for temp node.\n");
     return ((node *)NULL);
  }
  npnt->ctype = NNODE;		/* not to be used for valid node */
  npnt->nodenm1 = na;
  npnt->nodenm2 = nb;
  npnt->nodenm3 = nc;
  npnt->nodenm4 = nd;
  return npnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *ndt(int na)

{
  return (ndt (na,NULLVAL,NULLVAL,NULLVAL));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *ndt(int na, int nb)

{
  return (ndt (na,nb,NULLVAL,NULLVAL));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *ndt(int na, int nb, int nc)
{
  return (ndt (na,nb,nc,NULLVAL));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

node *ndt(void)
{
  return (ndt (NULLVAL,NULLVAL,NULLVAL,NULLVAL));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void ndt_free(node *nd)

{ 
  if (nd!=NULL) {
    if (nd->ctype==NNODE) efree (nd);
  }
}

/*---------------------------------------------------------------------*/

node *loc(node *npnt,double xloc,double yloc,double zloc)
{
  npnt->xloc = xloc;
  npnt->yloc = yloc;
  npnt->zloc = zloc;
  return (npnt);
}

/*---------------------------------------------------------------------*/

node *loc(node *npnt,double xloc,double yloc,double zloc, int region, int dendn)
{
  npnt->xloc = xloc;
  npnt->yloc = yloc;
  npnt->zloc = zloc;
  npnt->region = region;
  npnt->dendn = dendn;
  return (npnt);
}

/*---------------------------------------------------------------------*/

void label(node *npnt, int color, const char *text)
{
  if (npnt!=NULL) {
      npnt->label = color;
      npnt->labeltext = text;
  }
}

/*---------------------------------------------------------------------*/

void label(node *npnt, int color)
{
  label(npnt, color, NULL);
}

/*---------------------------------------------------------------------*/

void label(synapse *spnt, int color, const char *text)
{
  if (spnt!=NULL) {
      if (spnt->nodp1!=NULL) {
	  spnt->nodp1->label = color;
          spnt->nodp1->labeltext = text;
     }
  }
}
/*---------------------------------------------------------------------*/

void label(synapse *spnt, int color)
{
  label(spnt, color, NULL);
}

/*---------------------------------------------------------------------*/

void erase (node *npnt)

{
  erasenode(npnt);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void erase (elem *epnt)

{
 eraseelem(epnt);
}

/*---------------------------------------------------------------------*/

elem *at (node *node1, int elemtype )

{
   elem *epnt;

  epnt = makelem (elemtype,NULL);
  if (node1!=NULL) {
    epnt->node1a = node1->nodenm1;        /* set node number */
    epnt->node1b = node1->nodenm2;
    epnt->node1c = node1->nodenm3;
    epnt->node1d = node1->nodenm4;
  }
  epnt->node2a = NULLVAL;
  epnt->node2b = NULLVAL;
  epnt->node2c = NULLVAL;
  epnt->node2d = NULLVAL;
  checkelemnode(epnt);
  eninstall (epnt);            /* install new element in hash node table (req node1a,1b) */
  return epnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *at (int node1a, int node1b, int node1c, int node1d, int elemtype )
{
   elem *epnt;

  epnt = makelem (elemtype,NULL);
  epnt->node1a = node1a;        /* set node number */
  epnt->node1b = node1b; 
  epnt->node1c = node1c;
  epnt->node1d = node1d;
  epnt->node2a = NULLVAL;
  epnt->node2b = NULLVAL;
  epnt->node2c = NULLVAL;
  epnt->node2d = NULLVAL;
  checkelemnode(epnt);
  eninstall (epnt);            /* install new element in hash node table (req node1a,1b) */
  return epnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *at (int node1, int elemtype)

{
  return at(node1, NULLVAL, NULLVAL, NULLVAL, elemtype);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *at (int node1a, int node1b, int elemtype)

{
  return at(node1a, node1b, NULLVAL, NULLVAL, elemtype);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *at (int node1a, int node1b, int node1c, int elemtype )

{
  return at(node1a, node1b, node1c, NULLVAL, elemtype);
}

/*---------------------------------------------------------------------*/

void at(elem *epnt, node *pnode, double distfrac, node *nnode)

/* "at(ndn(node), elnum offset=frac nd(node));" */

/* make a new node along a cable element */

{
  int found, elnum;
  nodeint nod1a,nod1b,nod1c,nod1d;
  node *npnt,*npnt2,*npnt3;
  elem *elpnt,*epnt2;
  attrib *apnt,*napnt;
  conlst *lpnt;
  double xdist,ydist,zdist,len,ndia;


  nod1a=nnode->nodenm1;
  nod1b=nnode->nodenm2;
  nod1c=nnode->nodenm3;
  nod1d=nnode->nodenm4;

  elpnt=at (pnode, ELEMENT); /* this gets parent node number, */
			     /*  makes new temp. element, and connects parent node */
			     /*  at one end but does not set type of element, */
			     /*  nor its loc */

  if (!(npnt=findnode(elpnt->node1a,elpnt->node1b,
			elpnt->node1c,elpnt->node1d,"Atmake"))) {
     ncfprintf (stderr,"conn1m: can't find node %s\n",
			 prnode(elpnt->node1a,elpnt->node1b,
				elpnt->node1c,elpnt->node1d));
    execerror ("Missing parent node: ","stopping...");
  }
  				/* find previously made cable element */
  if (epnt==NULL) {
      ncfprintf (stderr,"Edist: can't find element %d\n",epnt->elnum);
      return;
  }
  else checkelemnode(epnt);
  elnum = epnt->elnum;

  for (found=0,lpnt=npnt->elemlst; lpnt && !found; lpnt=lpnt->next) {
    if (!(epnt=(elem *)lpnt->conpnt)) continue;
    if (epnt->elnum==elnum) found=1;
  }

  if (!found) {
     ncfprintf (stderr,"At node: can't find elem %d connected to node %s\n",
		elnum,prnode(elpnt->node1a,elpnt->node1b,
			     elpnt->node1c,elpnt->node1d));
     execerror ("Missing element: ","stopping...");
   }
				/* find the other end of element */

/* ncfprintf (stderr,"conn1m elem %d elnum %d  n1 %d n2 %d\n",
		epnt->elnum, elnum,epnt->nodp1,epnt->nodp2); /* */

  if ((npnt2=epnt->nodp1)==npnt) npnt2 = epnt->nodp2;
  if (!npnt2) {
      ncfprintf (stderr,"Atmake: can't find second node for elem %d\n", elnum);
      execerror("Missing node: ","stopping..."); 
  }
				/* stop if element is not a cable */
  if (epnt->ctype != CABLE) {
      ncfprintf (stderr,"At node: element %d is not a cable\n", elnum);
      execerror("Can't make node:","stopping..."); 
  } 

  elpnt = (cable *)makelem(CABLE,elpnt);         /* copy temp elem into cable */

			/* Next, find or make the new node,  */
			/* connect it properly to the new cable, */
			/*  and reconnect old cable to the new node. */ 
			/* Use old node if exists: */
if (!(npnt3=findnode(nod1a,nod1b,nod1c,nod1d,(char *)NULL))){ 
    setnodlst((npnt3=maknod(nod1a,nod1b,nod1c,nod1d)),elpnt);/* new node, pntr*/ 
  }
  elpnt->nodp2 = npnt3;		/* set pointer from new cable to new node */
  elpnt->node2a = nod1a;
  elpnt->node2b = nod1b;
  elpnt->node2c = nod1c;
  elpnt->node2d = nod1d;
  elpnt->nodp1 = npnt;		/* set pointer from new cable to old node */
  eninstall (elpnt);            /* install new element in hash node table (req node1a,1b) */
  unsetnodlst(npnt,epnt);	/* del pointer from parent node to old cable */
  setnodlst(npnt,elpnt);	/* set pointer from parent node to new cable */
  setnodlst(npnt3,epnt);	/* set pointer from new node to old cable */
  
  if (epnt->nodp1==npnt) {		/* set node pointer for old cable */
	epnt->nodp1 = npnt3;		/*  to point to "new" node (npnt3) */
	epnt->node1a = nod1a;
	epnt->node1b = nod1b;
	epnt->node1c = nod1c;
	epnt->node1d = nod1d;
  }
  else {
	epnt->nodp2 = npnt3;
	epnt->node2a = nod1a;
	epnt->node2b = nod1b;
	epnt->node2c = nod1c;
	epnt->node2d = nod1d;
  }

  if (distfrac > 1.0) distfrac = 1.0;
  if (distfrac < 0.0) distfrac = 0.0;

  if ((len=((cable*)epnt)->length)==NULLVAL)
     len = dist3d(npnt,npnt2);
  ((cable *)elpnt)->length = len * distfrac;	  /*set length of new cable */
  ((cable *)epnt)->length  = len * (1 - distfrac); /*set new length of old cabl*/

						/* find new dia */
  ndia = ((cable *)epnt)->dia * (1 - distfrac) + ((cable *)epnt)->dia2*distfrac;

  ((cable *)elpnt)->dia  = ((cable *)epnt)->dia;    /* set new dia first */
  ((cable *)epnt)->dia = ndia;	    		    /* set old dia */
  ((cable *)elpnt)->dia2   = ndia;                  /* set new dia2 */
  ((cable *)elpnt)->cplam = ((cable *)epnt)->cplam; /* make same cplam */
  ((cable *)elpnt)->Rm    = ((cable *)epnt)->Rm;
  ((cable *)elpnt)->Ri    = ((cable *)epnt)->Ri;
  ((cable *)elpnt)->Cm    = ((cable *)epnt)->Cm;
  ((cable *)elpnt)->vrev  = ((cable *)epnt)->vrev;
  ((cable *)elpnt)->vrest = ((cable *)epnt)->vrest;
  ((cable *)elpnt)->elabl = ((cable *)epnt)->elabl;
  ((cable *)elpnt)->region = ((cable *)epnt)->region;

 				/* copy attributes from old cable to new */ 
for (apnt=epnt->attpnt; apnt; apnt=apnt->attpnt) {
   napnt = chanattr(elpnt,apnt->ctype);		/* alloc new attribute */
   copyattr(apnt,napnt);
} 
  				/* Finally, calculate location of new node */
  xdist = npnt2->xloc - npnt->xloc;
  ydist = npnt2->yloc - npnt->yloc;
  zdist = npnt2->zloc - npnt->zloc;

  npnt3->xloc = npnt->xloc + xdist * distfrac;
  npnt3->yloc = npnt->yloc + ydist * distfrac;
  npnt3->zloc = npnt->zloc + zdist * distfrac;

}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, int ncells, double density, int grseed, double reg,
		double xcent, double ycent, int ginfo)
{
  int filout, textfl, n, first_center;
  double mean, stdev;
  double xsize, ysize;

  if (density > 0 && ncells > 0) {
     xsize = ysize = sqrt (ncells/density);
  }
  else xsize = ysize = 0;
  ncells = gausnn(mean=0, stdev=0, density, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells, xarr, yarr, first_center=0, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, double density, int grseed, double reg,
		double xsize, double ysize, double xcent, double ycent, int ginfo)
{
  int filout, textfl, ncells, first_center;
  double mean, stdev;

  ncells = gausnn(mean=0, stdev=0, density, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells=0, xarr, yarr, first_center=0, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, int ncells, int grseed, double reg,
		double xsize, double ysize, double xcent, double ycent, int ginfo)
{
  int filout, textfl, n, first_center;
  double mean, stdev, density;

  ncells = gausnn(mean=0, stdev=0, density=0, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells, xarr, yarr, first_center=0, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, double nnd, double reg, int grseed, 
		double xsize, double ysize, double xcent, double ycent, int ginfo)
{
  int filout, textfl, ncells, first_center;
  double mean, stdev, density;

  ncells = gausnn(mean=nnd, stdev=0, density=0, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells=0, xarr, yarr, first_center=0, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

/* These versions can make first cell in center */

int gausnn (double **xarr, double **yarr, int ncells, double density, int grseed, double reg,
		double xcent, double ycent, int first_center, int ginfo)
{
  int filout, textfl, n;
  double mean, stdev;
  double xsize, ysize;

  if (density > 0 && ncells > 0) {
     xsize = ysize = sqrt (ncells/density);
  }
  else xsize = ysize = 0;
  ncells = gausnn(mean=0, stdev=0, density, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells, xarr, yarr, first_center, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, double density, int grseed, double reg,
		double xsize, double ysize, double xcent, double ycent, int first_center, int ginfo)
{
  int filout, textfl, ncells;
  double mean, stdev;

  ncells = gausnn(mean=0, stdev=0, density, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells=0, xarr, yarr, first_center, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, int ncells, int grseed, double reg,
		double xsize, double ysize, double xcent, double ycent, int first_center, int ginfo)
{
  int filout, textfl, n;
  double mean, stdev, density;

  ncells = gausnn(mean=0, stdev=0, density=0, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells, xarr, yarr, first_center, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

int gausnn (double **xarr, double **yarr, double nnd, double reg, int grseed, 
		double xsize, double ysize, double xcent, double ycent, int first_center, int ginfo)
{
  int filout, textfl, ncells;
  double mean, stdev, density;

  ncells = gausnn(mean=nnd, stdev=0, density=0, reg, grseed, xsize, ysize, 
	   xcent, ycent, ncells=0, xarr, yarr, first_center, filout=0, textfl=0, ginfo);
  return ncells;
}

/*---------------------------------------------------------------------*/

sphere *make_sphere (node *node1, double dia, double Rm, double Cm)

{
   sphere *s;

   s = (sphere *)at(node1, SPHERE);
   s->dia = dia;
   s->Rm = Rm;
   s->Cm = Cm;
   return s;
}

/*---------------------------------------------------------------------*/

sphere *make_sphere (node *node1, double dia, double Rm)

{
   sphere *s;

   s = (sphere *)at(node1, SPHERE);
   s->dia = dia;
   s->Rm = Rm;
   return s;
}

/*---------------------------------------------------------------------*/

sphere *make_sphere (node *node1, double dia)

{
   sphere *s;

   s = (sphere *)at(node1, SPHERE);
   s->dia = dia;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_photorec (int type, node *node1, double xpos, double ypos, int stimchan)

/* make a photoreceptor and add it to the stimulus list (for stimfile) */

{
   photorec *s;
   recnod *maksnod(),*rpnt;

   s = (photorec *)at(node1, type);
   s->xpos = xpos;		    /* allows photorec to be at different loc than node */
   s->ypos = ypos;

   if (makestim) {
      rpnt = maksnod();
      rpnt->recnm1 = s->node1a;     /* node numbers are defined by conn() */
      rpnt->recnm2 = s->node1b;
      rpnt->recnm3 = s->node1c;
      rpnt->recnm4 = s->node1d;
      if (xpos!=LARGENODE) rpnt->xpos = xpos;
      else                 rpnt->xpos = 0;
      if (ypos!=LARGENODE) rpnt->ypos = ypos;
      else                 rpnt->ypos = 0;
   }
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_photorec (int type, node *node1, double xpos, double ypos)

{
     int stimchan;

   return make_photorec (type, node1, xpos, ypos, stimchan=0);
}

/*---------------------------------------------------------------------*/

photorec *make_photorec (int type, node *node1, int stimchan)

{
   return make_photorec (type, node1, node1->xloc, node1->yloc, stimchan);
}

/*---------------------------------------------------------------------*/

photorec *make_photorec (int type, node *node1)

{
     int stimchan;

   return make_photorec (type, node1, node1->xloc, node1->yloc, stimchan=0);
}

/*---------------------------------------------------------------------*/

photorec *make_cone (node *node1)

{
   photorec *s;

   s = (photorec *)at(node1, CONE);
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_cone (node *node1, int pigm)

{
   photorec *s;

   s = (photorec *)at(node1, CONE);
   s->pigm = pigm;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_cone (node *node1, double dia)

{
   photorec *s;

   s = (photorec *)at(node1, CONE);
   s->dia = dia;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_cone (node *node1, int pigm, double dia)

{
   photorec *s;
   s = (photorec *)at(node1, CONE);
   s->pigm = pigm;
   s->dia = dia;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_cone (node *node1, double xloc, double yloc, int stimchan)

{
  photorec *s;
  s =  make_photorec(CONE,node1,xloc,yloc,stimchan);
  return s;
}

/*---------------------------------------------------------------------*/

photorec *make_cone (node *node1, double xloc, double yloc)

{
  photorec *s;
  s =  make_photorec(CONE,node1,xloc,yloc);
  return s;
}

/*---------------------------------------------------------------------*/

photorec *make_chr (node *node1)

{
   photorec *s;

   s = (photorec *)at(node1, CHR);
   s->pigm = ChR2;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_chr (node *node1, double xloc, double yloc)

{
   photorec *s;
   s =  make_photorec(CHR,node1,xloc,yloc);
   s->pigm = ChR2;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_rod (node *node1)

{
   photorec *s;

   s = (photorec *)at(node1, ROD);
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_rod (node *node1, double dia)

{
   photorec *s;

   s = (photorec *)at(node1, ROD);
   s->dia = dia;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_rod (node *node1, int pigm, double dia)

{
   photorec *s;

   s = (photorec *)at(node1, ROD);
   s->pigm = pigm;
   return s;
}

/*---------------------------------------------------------------------*/

photorec *make_rod (node *node1, double xloc, double yloc, int stimchan)

{
   return make_photorec(ROD,node1,xloc,yloc,stimchan);
}

/*---------------------------------------------------------------------*/

photorec *make_rod (node *node1, double xloc, double yloc)

{
   return make_photorec(ROD,node1,xloc,yloc);
}

/*---------------------------------------------------------------------*/

photorec *make_transducer (node *node1)

{

   return make_photorec(VTRANSDUCER,node1);
}

/*---------------------------------------------------------------------*/

photorec *make_transducer (node *node1, double xpos, double ypos)

{
   return make_photorec(VTRANSDUCER,node1,xpos,ypos);
}

/*---------------------------------------------------------------------*/

photorec *make_itransducer (node *node1, double xpos, double ypos)

{
   return make_photorec(ITRANSDUCER,node1,xpos,ypos);
}

/*---------------------------------------------------------------------*/

photorec *make_itransducer (node *node1)

{

   return make_photorec(ITRANSDUCER,node1);
}

/*---------------------------------------------------------------------*/

elem *conn (node *node1, node *node2, int elemtype )
{
   elem *epnt;

  epnt = makelem (elemtype,NULL);
  if (node1!=NULL) {
    epnt->node1a = node1->nodenm1;        /* set node number */
    epnt->node1b = node1->nodenm2; 
    epnt->node1c = node1->nodenm3;
    epnt->node1d = node1->nodenm4;
  }
  if (node2!=NULL) {
    epnt->node2a = node2->nodenm1;
    epnt->node2b = node2->nodenm2;
    epnt->node2c = node2->nodenm3;
    epnt->node2d = node2->nodenm4;
  }
  checkelemnode(epnt);
  eninstall (epnt);                      /* install new element in hash node table (req nodenm1,2) */
  return epnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *conn (int node1a, int node1b, int node1c, int node1d,
		 int node2a, int node2b, int node2c, int node2d, int elemtype )
{
   elem *epnt;

  epnt = makelem (elemtype,NULL);
  epnt->node1a = node1a;        /* set node number */
  epnt->node1b = node1b; 
  epnt->node1c = node1c;
  epnt->node1d = node1d;
  epnt->node2a = node2a;
  epnt->node2b = node2b;
  epnt->node2c = node2c;
  epnt->node2d = node2d;
  checkelemnode(epnt);
  eninstall (epnt);                      /* install new element in hash node table (req nodenm1,2) */
  return epnt;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *conn (int node1, int node2, int elemtype)

{
  return conn(node1, NULLVAL, NULLVAL, NULLVAL,
		   node2, NULLVAL, NULLVAL, NULLVAL, elemtype);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *conn (int node1a, int node1b, int node2a, int node2b, int elemtype)

{
  return conn(node1a, node1b, NULLVAL, NULLVAL,
		   node2a, node2b, NULLVAL, NULLVAL, elemtype);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

elem *conn (int node1a, int node1b, int node1c, 
		 int node2a, int node2b, int node2c, int elemtype )

{
  return conn(node1a, node1b, node1c, NULLVAL,
		   node2a, node2b, node2c, NULLVAL, elemtype);
}

/*---------------------------------------------------------------------*/

cable *make_cable (node *node1, node *node2)

/* Make raw cable. */

{
   elem *epnt; 

  epnt = conn (node1, node2, CABLE);
  return (cable *)epnt;
}

/*---------------------------------------------------------------------*/

cable *make_cable (node *node1, node *node2, double dia)

/* make untapered cable */

{
   elem *epnt; 

  epnt = conn (node1, node2, CABLE);
  ((cable *)epnt)->dia  = dia;
  ((cable *)epnt)->dia2 = dia;
  return (cable *)epnt;
}

/*---------------------------------------------------------------------*/

cable *make_cable (node *node1, node *node2, double dia, double dia2)

/* make tapered cable */

{
   elem *epnt; 

  epnt = conn (node1, node2, CABLE);
  ((cable *)epnt)->dia  = dia;
  ((cable *)epnt)->dia2 = dia2;
  return (cable *)epnt;
}

/*---------------------------------------------------------------------*/

electrode *make_electrode (node *node1, node *node2)

{
   electrode *s;

   s = (electrode *)conn (node1, node2, ELECTRODE);
   return s;
}

/*---------------------------------------------------------------------*/

electrode *make_electrode (node *node1, node *node2, double rs, double cap)

/* make electrode with series resistance and grounded capacitance */

{
   electrode *s;

   s = (electrode *)conn (node1, node2, ELECTRODE);
   s->r  = rs;
   s->c  = cap;
   return s;
}
/*---------------------------------------------------------------------*/

electrode *make_electrode (node *node1, node *node2, double rs, double cap, double ccomp)

/* make electrode with series resistance and grounded capacitance */

{
   electrode *s;

   s = (electrode *)conn (node1, node2, ELECTRODE);
   s->r  = rs;
   s->c  = cap;
   s->ccomp  = ccomp;       // cap added to compartment at tip, can be negative
   return s;
}
/*---------------------------------------------------------------------*/

loadelem *make_shunt (node *node1, double rs)

/* make shunt load with series resistance, grounded to vrev */

{
   loadelem *s;

   s = (loadelem *)at (node1, LOAD);
   s->r  = rs;
   s->vrev  = 0;
   return s;
}
/*---------------------------------------------------------------------*/

loadelem *make_shunt (node *node1, double rs, double vrev)

/* make shunt load with series resistance, grounded to vrev */

{
   loadelem *s;

   s = (loadelem *)at (node1, LOAD);
   s->r  = rs;
   s->vrev  = vrev;
   return s;
}
/*---------------------------------------------------------------------*/

void ename (elem *epnt, int *num)

{
  *num  = epnt->elnum;
  epnt->nocondens = 1; 
}

/*---------------------------------------------------------------------*/

void ename (elem *epnt, double *num)

{
  *num  = epnt->elnum;
  epnt->nocondens = 1; 
}

/*---------------------------------------------------------------------*/

void plot (int param, node *npnt)

{
    int pmod, pval;

 if (param !=FUNCTION && param !=VAR && param !=S && param !=CABLE) {

   if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

   if (npnt!=NULL) {
     plotnod[plotnum].cnod1 = npnt->nodenm1;
     plotnod[plotnum].cnod2 = npnt->nodenm2;
     plotnod[plotnum].cnod3 = npnt->nodenm3;
     plotnod[plotnum].cnod4 = npnt->nodenm4;
   }
 }
 pval = pmod = 0;

  switch (param) {
   case V:
        pmod  = VREC;
        break;
   case I:
        pmod  = IREC;
        break;
   case IM:
        pmod  = MREC;
        break;
   case L:
        pmod  = LREC;
        break;
   case VM:
        pmod  = WREC;
        break;
   case FA0: case FA1: case FA2: case FA3: case FA4: case FA8: case FA9:
   case FH0: case FH1: case FH2: case FH3: case FH4:
   case FB0: case FB1: case FB2: case FB3: case FB4:
   case FC0: case FC1: case FC2: case FC3: case FC4: case FC9:
        pmod  = NRECA0 + (param - FA0);
        break;
   case G: 
        pmod  = NRECG;
        break;
   case G0: case G1: case G2: case G3:
        pmod  = NRECG0 + (param - G0);
        break;
   //case S:
   //     plotnod[plotnum].spnt  = (char *)d1.sym;
   //     pmod  = SREC;
   //     break;
   //case CABLE:
   //     pmod  = CABLE;
   //     plotnod[plotnum].cnod2 = (int)(d2.val * CABLFRAC);
   //     break;
   case GLU:
   case AMPA:
   case KAINATE:
   case NMDA:
   case CNQX:
   case GABA:
   case PTX:
   case BIC:
   case GLY:
   case STRY:
   case CGMP:
   case CAMP:
   case PH:
   case ATP:
	pmod  = param - GLU + CREC_GLU; 
	      break;

 }
 plotnod[plotnum].pval  = pval;
 plotnod[plotnum].pmod  = pmod;

 if (strcmp(gfrname(),rootframe))                /* if diff than root frame */
   strcpy(plotnod[plotnum].plframe,gfrname());  /* copy current frame name */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int ename)

{
    int pmod, pval;
    elem *epnt=NULL;

 if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

 plotnod[plotnum].cnod1 = ename;
 plotnod[plotnum].cnod2 = NULLVAL;
 plotnod[plotnum].cnod3 = NULLVAL;
 plotnod[plotnum].cnod4 = NULLVAL;
 epnt=findelem(ename);
 
 pval = pmod = 0;

  switch (param) {
   case V:
        pmod  = VREC;
        break;
   case I:
	if (epnt) pmod = NRECGI;
	else      pmod = IREC;
        break;
   case IM:
        pmod  = MREC;
        break;
   case L:
        pmod  = LREC;
        break;
   case VM:
        pmod  = WREC;
        break;
   case FA0: case FA1: case FA2: case FA3: case FA4: case FA8: case FA9:
   case FH0: case FH1: case FH2: case FH3: case FH4:
   case FB0: case FB1: case FB2: case FB3: case FB4:
   case FC0: case FC1: case FC2: case FC3: case FC4: case FC9:
        pmod  = NRECA0 + (param - FA0);
        break;
   case G: 
        pmod  = NRECG;
        break;
   case G0: case G1: case G2: case G3:
        pmod  = NRECG0 + (param - G0);
        break;
   //case S:
   //     plotnod[plotnum].spnt  = (char *)d1.sym;
   //     pmod  = SREC;
   //     break;
   //case CABLE:
   //     pmod  = CABLE;
   //     plotnod[plotnum].cnod2 = (int)(d2.val * CABLFRAC);
   //     break;
   case GLU:
   case AMPA:
   case KAINATE:
   case NMDA:
   case CNQX:
   case GABA:
   case PTX:
   case BIC:
   case GLY:
   case STRY:
   case CGMP:
   case CAMP: pmod  = param - GLU + CREC_GLU; 
	      break;

 }
 plotnod[plotnum].pval  = pval;
 plotnod[plotnum].pmod  = pmod;

 if (strcmp(gfrname(),rootframe))                /* if diff than root frame */
   strcpy(plotnod[plotnum].plframe,gfrname());  /* copy current frame name */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int ename, double frac)

{
    int pmod, pval;
   
 if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

 plotnod[plotnum].cnod1 = ename;
 plotnod[plotnum].cnod2 = NULLVAL;
 plotnod[plotnum].cnod3 = NULLVAL;
 plotnod[plotnum].cnod4 = NULLVAL;
 
 pval = pmod = 0;
  switch (param) {
    case CABLE:
        pmod  = CABLE;
        plotnod[plotnum].cnod2 = (int)(frac * CABLFRAC);
	break;
    default:
        pmod  = 0;
	break;
  }
 plotnod[plotnum].pval  = pval;
 plotnod[plotnum].pmod  = pmod;

 if (strcmp(gfrname(),rootframe))                /* if diff than root frame */
   strcpy(plotnod[plotnum].plframe,gfrname());  /* copy current frame name */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double t), double plval, double maxx, double minx, int pmod, int pmod2)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plotnum);
  plotnod[plotnum].pmod  = pmod;
  plotnod[plotnum].pmod2 = pmod2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double t), double plval, double maxx, double minx)

{
   plot_func (func, plval, maxx, minx, FREC,0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_funcg (double(*func)(double val, double t), double plval, double maxx, double minx)

{
   plot_func (func, plval, maxx, minx, FREC, NRECG0);
}
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_funci (double(*func)(double val, double t), double plval, double maxx, double minx)

{
   plot_func (func, plval, maxx, minx, FREC, IREC);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double t), double plval, double plval2, double maxx, double minx, int pmod, int pmod2)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plval2,plotnum);
  plotnod[plotnum].pmod  = pmod;
  plotnod[plotnum].pmod2 = pmod2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double t), double plval, double plval2, double maxx, double minx) 
{
   plot_func (func, plval, plval2, maxx, minx, FREC,0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double t), double plval, double plval2, double plval3, double maxx, double minx, int pmod, int pmod2)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plval2,plval3,plotnum);
  plotnod[plotnum].pmod  = pmod;
  plotnod[plotnum].pmod2 = pmod2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double t), double plval, double plval2, double plval3, double maxx, double minx)

{
   plot_func (func, plval, plval2, plval3, maxx, minx, FREC,0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double t), double plval, double plval2, double plval3, double plval4, double maxx, double minx, int pmod, int pmod2)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plval2,plval3,plval4,plotnum);
  plotnod[plotnum].pmod  = pmod;
  plotnod[plotnum].pmod2 = pmod2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double t), double plval, double plval2, double plval3, double plval4, double maxx, double minx)

{
    plot_func (func, plval, plval2, plval3, plval4, maxx, minx, FREC,0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double t), double plval, double plval2, double plval3, double plval4, double plval5, double maxx, double minx, int pmod, int pmod2)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plval2,plval3,plval4,plval5,plotnum);
  plotnod[plotnum].pmod  = pmod;
  plotnod[plotnum].pmod2 = pmod2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double t), double plval, double plval2, double plval3, double plval4, double plval5, double maxx, double minx)

{
   plot_func (func, plval, plval2, plval3, plval4, plval5, maxx, minx, FREC, 0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double maxx, double minx, int pmod)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plval2,plval3,plval4,plval5,plval6,plotnum);
  plotnod[plotnum].pmod  = pmod;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double maxx, double minx)

{
   plot_func (func, plval, plval2, plval3, plval4, plval5, plval6, maxx, minx, FREC);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double val7, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double plval7, double maxx, double minx, int pmod)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotfunc(func);
  plotval(plval,plval2,plval3,plval4,plval5,plval6,plval7,plotnum);
  plotnod[plotnum].pmod  = pmod;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double val7, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double plval7, double maxx, double minx)

{
   plot_func (func, plval, plval2, plval3, plval4, plval5, plval6, plval7, maxx, minx, FREC);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_var (double *var, double plval, double maxx, double minx, int pmod)
{
   double t;
 
  if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

  plotnod[plotnum].pymax  = maxx;
  plotnod[plotnum].pymin  = minx;
  if (!makestim) plotvar(var);
  plotval(plval,plotnum);
  plotnod[plotnum].pmod  = pmod;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_var (double *var, double plval, double maxx, double minx)

{
   plot_var (var, plval, maxx, minx, SREC);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot_arr (double *arr, int maxindex)

{
  if (!makestim) plotarr(arr, maxindex);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, node *npnt, double maxx, double minx)

{
   double t;
 
 plot (param, npnt);
 if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
 plotnod[plotnum].pymax  = maxx;
 plotnod[plotnum].pymin  = minx;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int cagparm, int ename)

{
   int pmod;
   int pval;

 if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;

 plotnod[plotnum].cnod1 = ename;
 plotnod[plotnum].cnod2 = NULLVAL;
 plotnod[plotnum].cnod3 = NULLVAL;
 plotnod[plotnum].cnod4 = NULLVAL;
 pval = 0;
 switch(param) {
  case G:
        if (cagparm > 0 && cagparm <30) {
               pval = cagparm;
               pmod = NRECG;                /* record state fraction */
	}else
        switch(cagparm) {

          case 0:   pmod = NRECG0; break;        /* record conductance */
          case VVREV:pmod = NRECGV; break;
          case I:   pmod = NRECGI; break;
          case M:   pmod = NRECGM; break;
          case H:   pmod = NRECGH; break;
          case CA:  pmod = NRECGC; break;
        }
        break;

   default: break;
 }
 plotnod[plotnum].pval  = pval;
 plotnod[plotnum].pmod  = pmod;

 if (strcmp(gfrname(),rootframe))                /* if diff than root frame */
   strcpy(plotnod[plotnum].plframe,gfrname());  /* copy current frame name */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int cagparm, node *npnt)

{
   int pmod;
   int pval;

  if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
  if (npnt!=NULL) {
    plotnod[plotnum].cnod1 = npnt->nodenm1;
    plotnod[plotnum].cnod2 = npnt->nodenm2;
    plotnod[plotnum].cnod3 = npnt->nodenm3;
    plotnod[plotnum].cnod4 = npnt->nodenm4;
  }
  pval = 0;
  switch(param) {

   case CA:
        if (cagparm > -1000 && cagparm <=1000) {
               pval = cagparm;
               pmod = CREC_CACONC;           /* record state fraction */
	}else
        switch(cagparm) {

          case 0:  pmod = CREC_CACONC;
                   pval = 1; break;          /* record cai ([Ca] inside) */
          case VVREV:pmod = CREC_CAV; break;
          case I:   pmod = CREC_CAIT; break;
          case IP:  pmod = CREC_CAIP; break;
          case IE:  pmod = CREC_CAIE; break;
          case IPE: pmod = CREC_CAIPE; break;
          case CAS: pmod = CREC_CAS; break;
          case CICR: pmod = CREC_CICR; break;
        }
        break;
   case CABUF:
        if (cagparm > -1000 && cagparm <=1000) {
               pval = cagparm;
               pmod = CREC_CABUF;               /* record Ca buf */
	}else
        switch(cagparm) {

        //  case 1:  d5 = popm();
        //           pval = int(d5.val);
        //           pmod = CREC_CABUF;           /* record cai from Ca shell */
        //           break;
        default:  break;
        }
        break;
   case CABUFB:
        if (cagparm > -1000 && cagparm <=1000) {
               pval = cagparm;
               pmod = CREC_CABUFB;               /* record buf bound */
	}else
        switch(cagparm) {

        //  case 1:  d5 = popm();
        //           pval = int(d5.val);
        //           pmod = CREC_CABUFB;         /* record cai from Ca shell */
        //           break;
        default:  break;
        }
        break;
   default: break;
 }
 plotnod[plotnum].pval  = pval;
 plotnod[plotnum].pmod  = pmod;

 if (strcmp(gfrname(),rootframe))                /* if diff than root frame */
   strcpy(plotnod[plotnum].plframe,gfrname());  /* copy current frame name */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int cagparm, node *npnt,  double maxx, double minx)

{
   double t;
 
 plot (param, cagparm, npnt);
 if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
 plotnod[plotnum].pymax  = maxx;
 plotnod[plotnum].pymin  = minx;
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int cagparm, int ename, double maxx, double minx)

{
   double t;
 
 plot (param, cagparm, ename);
 if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
 plotnod[plotnum].pymax  = maxx;
 plotnod[plotnum].pymin  = minx;
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void plot (int param, int ename, double maxx, double minx)

{
   double t;
 
 plot (param, ename);
 if (maxx < minx) { double t; t = maxx; maxx=minx; minx=t; }
 plotnod[plotnum].pymax  = maxx;
 plotnod[plotnum].pymin  = minx;
}


/*---------------------------------------------------------------------*/

#ifdef XSTIM
  #define STIMCODE 1
#else
  #define STIMCODE 0
#endif

void ncinit(int argc, char **argv)

{
   int hlp = 0;

  timepnt = &timebuf;
  startclk = times(timepnt);
  startu = timepnt->tms_utime;
  makestim = STIMCODE;
  prcomps = 0;
 
  setptrs();		/* set pointers to variables that will be set from ncmdlin */
  ncmdlin(argc,argv);	/* set variables from command line */
  infile = progname;
  einfile = infile;
  istat = (int *)signal(SIGINT, SIG_IGN);  /* save original status */
  if (istat != (int*)SIG_IGN)
      (void) signal (SIGINT,onintr);   /* "onintr" at ^C */
  istat = (int *)signal(SIGTERM, SIG_IGN);  /* save original status */
  if (istat != (int*)SIG_IGN)
      (void) signal (SIGTERM,onintr);      /* "onintr" at TERM */
  istat = (int *)signal(SIGPIPE, SIG_IGN);  /* save original status */
  if (istat != (int*)SIG_IGN)
      (void) signal (SIGPIPE,onintr);      /* "onintr" at TERM */
  signal(SIGFPE, fpecatch);

  if (hlp && !silent)
     ncfprintf (stdout,"## %s version %s\n",progname,print_version(ncversion));
  runyet = 0;		/* resets time */
  ncleanup();		/* erase lists  (in "ncmak.c")  */
  initpl(0);		/* init plot stuff, in "ncplot.c" */
  init();		/* init variables, symbol table, in "init.c" */
  if (setdebug) debug  = setdebug;      /* "-y n" overrides "debug=n" */
  if (setdebgz) debugz = setdebgz;      /* "-z n" overrides "debugz=n" */
}

/*---------------------------------------------------------------------*/

void mrun(int runtyp, double step)
{

if (debug & 1) ncfprintf (stderr,"mrun start\n");

  vinit();				/* set symbol table from C variables */

  if (!makestim) {			/* nc running */
    if (!runyet) {
        initchan();                     /* set up sequential-state chans */
    }
    if (!runyet || elempnt) {
        findconnect(0);
    }
    switch (runtyp)
      {

       case RUN:  runsim(0.0,1);
                  runyet = 0;
                  break;

       case STEP: runsim(step,0);
                  break;

      }
   }
   else {				/* makestim running */
    switch (runtyp)
      {
       case RUN:
       case STEP:
                simtime = step + simtime;  /* if XSTIM, set the time, at least */
                break;
      }
   }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void run(void)

{
  mrun (RUN,0);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void step(double step)

{
  mrun (STEP,step);
}


/*---------------------------------------------------------------------*/

void ncexit(void)

{
  if (makestim) {
     ncfprintf(stimout,"## %s version %s\n",progname,print_version(ncversion));
     prtime(stimout);          /* print out time */
     fclose (stimout);
     if (stimout!=stdout) sortstimfile();
  }
  else bexit();

}
/*---------------------------------------------------------------------*/

extern plotfr plotnod[PLOTNODSIZ];  /* holds nodes to be recorded */
extern int plotnum;            /* number of plots to be displayed */
extern int vidmode;             /* =1 -> graphics, =0 -> text */

void mplot (double y, double x, int n, int i, int pflag);
void plotinit(int plotnum);
void plotrst(int plotnum);
void plotpen (int val, int i);
void plotchar (int val,int lines, int i);
void plotcsiz (double val, int plotnum);
void plotfilt (int nfilt, double *timec, int plotnum);
void plotfilt (int nfilt, double *timec);
void plotvpenc  (double (vpen)(int,double,double), int plotnum);
void plotvpenc  (double (vpen)(int,double,double));
void plotname (const char *plname, int plotnum);
void plotname (const char *plname);
void plotn (int pl, int plotnum);
void plotn (int pl);
void plotsize (double plotsiz, int plotnum);
void plotsize (double plotsiz);
void plotmode (int pmod, int plotnum);
void plotmode (int pmod);

/* - - - - - - - - - - - - - - - -*/

void graph(double x, double y1)

/* draw a graph */
{
   mplot (y1,x,1,0,1);
}

void graph(double x, double y1, double y2)
{
   int n;

  n=2;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
}

void graph(double x, double y1, double y2, double y3)
{
   int n;

  n=3;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
  mplot (y3,x,n,2,1);
}

void graph(double x, double y1, double y2, double y3, double y4)
{
   int n;

  n=4;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
  mplot (y3,x,n,2,1);
  mplot (y4,x,n,3,1);
}

void graph(double x, double y1, double y2, double y3, double y4, double y5)
{
   int n;

  n=5;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
  mplot (y3,x,n,2,1);
  mplot (y4,x,n,3,1);
  mplot (y5,x,n,4,1);
}

void graph(double x, double y1, double y2, double y3, double y4, double y5, double y6)
{
   int n;

  n=6;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
  mplot (y3,x,n,2,1);
  mplot (y4,x,n,3,1);
  mplot (y5,x,n,4,1);
  mplot (y6,x,n,5,1);
}

void graph(double x, double y1, double y2, double y3, double y4, double y5, double y6, double y7)
{
   int n;

  n=7;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
  mplot (y3,x,n,2,1);
  mplot (y4,x,n,3,1);
  mplot (y5,x,n,4,1);
  mplot (y6,x,n,5,1);
  mplot (y7,x,n,6,1);
}

void graph(double x, double y1, double y2, double y3, double y4, double y5, double y6, double y7, double y8)
{
   int n;

  n=8;
  mplot (y1,x,n,0,1);
  mplot (y2,x,n,1,1);
  mplot (y3,x,n,2,1);
  mplot (y4,x,n,3,1);
  mplot (y5,x,n,4,1);
  mplot (y6,x,n,5,1);
  mplot (y7,x,n,6,1);
  mplot (y8,x,n,7,1);
}

/* - - - - - - - - - - - - - - - -*/

void graph(double x, double *y, int n)

{
   int i;

  if (y!=NULL) {
    for (i=0; i<n; i++) {
       mplot (y[i],x,n,i,1);
    }
  }
}

/* - - - - - - - - - - - - - - - -*/

void graph_x (double xmax, double xmin)
/* access the next plot */
{
   int np;
   if (plotnum<PLOTNODSIZ-1) np = plotnum+1;
   			                         /* bug: can't use for indiv separate graphs */
        if (!strcmp(gfrname(),rootframe))        /* if root fr */
              plotnum = -1;                /* number of plots for this graph */
        if (xmax < xmin) { double t; t = xmax; xmax = xmin; xmin=t; }
        plotnod[np].pxmax  = xmax;
        plotnod[np].pxmin  = xmin;
        plotnod[np].xrange = xmax - xmin;
}

/* - - - - - - - - - - - - - - - -*/

void graph_y (double ymax, double ymin)

/* make new plot y axis */

{
   if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
   if (strcmp(gfrname(),rootframe))         /* if diff than root frame */
     strncpy(plotnod[plotnum].plframe,gfrname(),PNAMSIZ-1); /* cp frame nm */
   if (ymax < ymin) { double t; t = ymax; ymax = ymin; ymin=t; }
   plotnod[plotnum].pymax  = ymax;
   plotnod[plotnum].pymin  = ymin;
   plotnod[plotnum].yrange = ymax - ymin;
   plotnod[plotnum].pmod  = GREC;
}

/* - - - - - - - - - - - - - - - -*/
void graph_init (void)

{
    plotinit(plotnum+1);
}

/* - - - - - - - - - - - - - - - -*/
void graph_restart(void)

{
    plotrst(plotnum);
}

/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen)
{
   plotpen (pen,0);
}

/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen1, int pen2)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
}
/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen1, int pen2, int pen3)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
   plotpen (pen3,2);
}
/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen1, int pen2, int pen3, int pen4)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
   plotpen (pen3,2);
   plotpen (pen4,3);
}
/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
   plotpen (pen3,2);
   plotpen (pen4,3);
   plotpen (pen5,4);
}

/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5, int pen6)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
   plotpen (pen3,2);
   plotpen (pen4,3);
   plotpen (pen5,4);
   plotpen (pen6,5);
}

/* - - - - - - - - - - - - - - - -*/
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5, int pen6, int pen7)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
   plotpen (pen3,2);
   plotpen (pen4,3);
   plotpen (pen5,4);
   plotpen (pen6,5);
   plotpen (pen7,6);
}
/* - - - - - - - - - - - - - - - -*/

void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5, int pen6, int pen7, int pen8)
{
   plotpen (pen1,0);
   plotpen (pen2,1);
   plotpen (pen3,2);
   plotpen (pen4,3);
   plotpen (pen5,4);
   plotpen (pen6,5);
   plotpen (pen7,6);
   plotpen (pen8,7);
}

/* - - - - - - - - - - - - - - - -*/
void graph_pen (int *val, int n)
{
   int i;

  if (val!=NULL) { 
    for (i=0; i<n; i++ ) {     /* set plot colors before graph init */
      plotpen (val[i],i);
    }
  }
}

/* - - - - - - - - - - - - - - - -*/
void graph_char (char val)
{
  plotchar (val,LINES,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void graph_cchar (char val)
{
  plotchar (val,LINES,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void graph_csiz (double val)
{
   plotcsiz (val,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void graph_set (const char *name, int plnum, double plsize, double plval)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotsize(plsize,plotnum);
    plotval(plval,plotnum);
}
/* - - - - - - - - - - - - - - - -*/

void graph_set (const char *name, int plnum, double plsize)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotsize(plsize,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_pen (int pen)
{
   plotpen (pen,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_param (const char *name, int plnum, double plsize)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotsize(plsize,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_param (const char *name, int pen, int plnum, double plsize)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotsize(plsize,plotnum);
    plotpen(pen,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_param (const char *name, int pen, int plnum, double plsize, int pmod)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotsize(plsize,plotnum);
    plotpen(pen,plotnum);
    plotmode(pmod,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_param (const char *name, int pen, int plnum)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotpen(pen,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_param (const char *name, int pen, int plnum, double plsize, double plval)

{
    plotname(name,plotnum);
    plotn(plnum,plotnum);
    plotsize(plsize,plotnum);
    plotval(plval,plotnum);
    plotpen(pen,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void graph_filt (double *val, int n)

{
#define TTIMEC 10
    if (n >= TTIMEC) n = TTIMEC - 1;
    if (!makestim) plotfilt (n,val,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void graph_vpen (double (*vpen)(int, double, double))
{
    if (!makestim) plotvpenc (vpen,plotnum);
}

/* - - - - - - - - - - - - - - - -*/

void plot_vpen (double (*vpen)(int, double, double))
{
    if (!makestim) plotvpenc (vpen,plotnum);
}

/*------------------------------------------------------------*/

double v(node *npnt)
{
    double rval=0;
    comp *cpnt;

  if (npnt==NULL) return(0.0);
  if (npnt->comptr) rval = npnt->comptr->v;
  return rval;
}

/* - - - - - - - - - - - - - - - -*/

double v(int node1a, int node1b, int node1c, int node1d)

{
    node *npnt;

  if (!(npnt=findnode(node1a, node1b, node1c, node1d, "record V"))) return(0.0);
  return (v(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double v(int node1a, int node1b, int node1c)

{
    node *npnt;

  if (!(npnt=findnode(node1a, node1b, node1c, NULLVAL, "record V"))) return(0.0);
  return (v(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double v(int node1a, int node1b)

{
    node *npnt;

  if (!(npnt=findnode(node1a, node1b, NULLVAL, NULLVAL, "record V"))) return(0.0);
  return (v(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double v(int node1a)

{
    node *npnt;

  if (!(npnt=findnode(node1a, NULLVAL, NULLVAL, NULLVAL, "record V"))) return(0.0);
  return (v(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double v(int elemnum, double fracdist)
{
  if (!makestim) reccable(elemnum, fracdist);
}

/* - - - - - - - - - - - - - - - -*/

double i(node *npnt)

{
    double rval=0;
    comp *cpnt;

  if (npnt==NULL) return(0.0);
  if ((cpnt=npnt->comptr)) {
     if (cpnt->miscfl & VEXT) {		/* comp is v-clamped */
           rval = cpnt->extvi;
	   // fprintf (stderr,"i: %g\n",rval);
     }
     else if (cpnt->miscfl & IEXT) {	/* comp is i-clamped */
           rval = cpnt->exti;
     }
  }
  return rval;
}
/* - - - - - - - - - - - - - - - -*/

double i(int node1a, int node1b, int node1c, int node1d)

{
    node *npnt;

  if (!(npnt=findnode(node1a, node1b, node1c, node1d, 
	"record I"))) return(0.0);
  return (i(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double i(int node1a, int node1b, int node1c)

{
    node *npnt;

  if (!(npnt=findnode(node1a, node1b, node1c, NULLVAL, "record I"))) return(0.0);
  return (i(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double i(int node1a, int node1b)

{
    node *npnt;

  if (!(npnt=findnode(node1a, node1b, NULLVAL, NULLVAL, "record I"))) return(0.0);
  return (i(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double i(int node1a)

{
    node *npnt;

  if (!(npnt=findnode(node1a, NULLVAL, NULLVAL, NULLVAL, "record I"))) return(0.0);
  return (i(npnt));
}

/* - - - - - - - - - - - - - - - -*/

double l(node *npnt)
{
    double rval;
    comp *cpnt;
    photrec *rpnt;


 if (findphotrec (npnt->nodenm1, npnt->nodenm2, npnt->nodenm3, npnt->nodenm4, 
				&rpnt, "record L")) {
    rval = rpnt->iflux / stiminc;
    return rval;
 }
 else return 0;
}

/* - - - - - - - - - - - - - - - -*/

double record_nt (node *npnt, int nt)

/* nt = GLU, AMPA, KAINATE, NMDA, CNQX, GABA, BIC, PTX, GLY
        STRY, CAMP, CGMP
*/

{
    double rval=0;

  if (npnt!=NULL) {
     if (npnt->comptr) {
      rval = getnt (npnt->comptr,nt);
     }
  }
  return rval;
}

/* - - - - - - - - - - - - - - - -*/

double record_ca (node *npnt, int caval, int pval)

/* CACONC, CAV, CAIT, CAIP, CAIE, CAIPE, CABUF, CAS, CICR */

{
    double rval=0;

  if (npnt->comptr) {
         cacomp *capnt;

    if (capnt=npnt->comptr->capnt) {
       switch (caval) {

       case CABUF: 
          if (capnt->cab) {
            if (pval==0) {
	       if (capnt->caos!=NULL) rval = capnt->caos[0];
	       else                   rval = capnt->cao;
	    }
            else if (pval>(capnt->cashell-1))
                 rval = capnt->cab[capnt->cashell-1];
            else rval = capnt->cab[pval-1];
          }
          break;
	case CA:
          if (caval==0) {
	     if (capnt->caos!=NULL) rval = capnt->caos[0];
	     else	            rval = capnt->cao;
          }
          else if (caval>(capnt->cashell-1))
               rval = capnt->cais[capnt->cashell-1];
          else rval = capnt->cais[pval-1];
          break;
       case VVREV:
          rval = capnt->vrev; break;
       case I:
          rval = -capnt->ica;   break;
       case IPE:
          rval = -(capnt->ipump+capnt->iexch); break;
       case IE:
           rval = -capnt->iexch; break;
       case IP:
          rval = -capnt->ipump; break;
       case CAS:
          rval = capnt->cas; break;
       case CICR:
          rval = capnt->cicrflux; break;
      }
   }
   else {
        static int ca_err_printed = 0;

      if (!ca_err_printed) {
          ncfprintf (stderr,"# nc: error, node %s does not have a Ca compartment.\n",
                 prnode(npnt->nodenm1,npnt->nodenm2,
                        npnt->nodenm3,npnt->nodenm4));
          ca_err_printed = 1;
      }
    }
  }
  return rval;
}
/* - - - - - - - - - - - - - - - -*/

double record_ca (node *npnt, int caval)

{
 return record_ca (npnt, caval, 0);
}

/* - - - - - - - - - - - - - - - -*/

double record_chan (int elnum, int nf, int pval)

{
  double rval=0;

 switch (nf) {
   case G:   nf = NRECG;  break;
   case G0:  nf = NRECG0; break;
   case G1:  nf = NRECG1; break;
   case G2:  nf = NRECG2; break;
   case G3:  nf = NRECG3; break;
   case M:   nf = NRECGM; break;
   case H:   nf = NRECGH; break;
   case VVREV:nf = NRECGV; break;
   case V:   nf = NRECGV; break;
   case I:   nf = NRECGI; break;
   case CA:  nf = NRECCA; break;
   default:  nf=0;
 }
#ifndef XSTIM
 rval =  recchan(elnum, nf, pval);
#endif
 return rval;
}

/* - - - - - - - - - - - - - - - -*/

double record_chan (int elnum, int nf)

{
 return record_chan (elnum, nf, 0);
}

/* - - - - - - - - - - - - - - - -*/

double record_chan (elem *epnt, int nf, int pval)

{
  double rval=0;

 switch (nf) {
   case G:   nf = NRECG;  break;
   case G0:  nf = NRECG0; break;
   case G1:  nf = NRECG1; break;
   case G2:  nf = NRECG2; break;
   case G3:  nf = NRECG3; break;
   case M:   nf = NRECGM; break;
   case H:   nf = NRECGH; break;
   case VVREV:nf = NRECGV; break;
   case V:   nf = NRECGV; break;
   case I:   nf = NRECGI; break;
   case CA:  nf = NRECCA; break;
   default:  nf=0;
 }
#ifndef XSTIM
 rval =  recchan(epnt, nf, pval);
#endif
 return rval;
}

/* - - - - - - - - - - - - - - - -*/

double record_chan (chan *chpnt, int nf, int pval)

{
  double rval=0;
  int cnum;

 switch (nf) {
   case G:   nf = NRECG;  break;
   case G0:  nf = NRECG0; break;
   case G1:  nf = NRECG1; break;
   case G2:  nf = NRECG2; break;
   case G3:  nf = NRECG3; break;
   case M:   nf = NRECGM; break;
   case H:   nf = NRECGH; break;
   case VVREV:nf = NRECGV; break;
   case V:   nf = NRECGV; break;
   case I:   nf = NRECGI; break;
   case CA:  nf = NRECCA; break;
   default:  nf=0;
 }
#ifndef XSTIM
 rval = rechan(chpnt, nf, pval, cnum=0);
#endif
 return rval;
}

/* - - - - - - - - - - - - - - - -*/

double record_synapse(int snum, int nf)

{
   double rval=0;
  if (nf==CA) nf=NRECCA;
  else nf = NRECA0 + (nf - FA0);
#ifndef XSTIM
  rval =  recsynapse(snum, nf);
#endif
  return rval;
}

/* - - - - - - - - - - - - - - - -*/

double record_synapse(synapse *spnt, int nf)

{
   synap *sp;
   double rval=0;

  sp = (synap *)spnt->lptr;        /* pointer to synapse from elem */
  if (nf==CA) nf=NRECCA;
  else nf = NRECA0 + (nf - FA0);
#ifndef XSTIM
  rval =  recsynap(sp,nf,spnt->elnum);
#endif
  return (rval);
}

/*---------------------------------------------------------------------*/

/* defined in "modcode.cc" : */

void set_display (int elemtype, int disptype, int narg2,
             node *nd1, node *nd2, node *nexc, int exceptype,
             double zrange1, double zrange2,
             int dcolor, Symbol *vpen, double(*vpenn)(int,int),
             int cmap, double vmax,double vmin,
             double dscale, int hide, int excl, double stime);

void set_disp_rot (double xrot,  double yrot,  double zrot,
             double dxcent, double dycent, double dzcent,
             double rxcent, double rycent, double rzcent, double dsize);

static double zrange1=-LARGENODE;
static double zrange2= LARGENODE;

/*---------------------------------------------------------------------*/

void display_size(double size)
{
  dsize = size;
  set_disp_rot (0, 0, 0, 0, 0, 0, 0, 0, 0, size);
}

/*---------------------------------------------------------------------*/

void display (int disptype, node *nd1, int dcolor, double dscale)

{
    int elemtype, narg2, exceptype, hide, cmap, excl;
    double vmax, vmin, stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype=0, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap=0, vmax=0, vmin=0,
             dscale, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, node *nd1, int dcolor, double dscale)

{
    int narg2, exceptype, hide, cmap, excl;
    double vmax, vmin, stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap=0, vmax=0, vmin=0,
             dscale, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, node *nd1, int dcolor)

{
   int elemtype;
   double dscale;
   display (elemtype=0, disptype, nd1, dcolor, dscale=1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, node *nd1, int dcolor)

{
   double dscale;
   display (elemtype, disptype, nd1, dcolor, dscale=1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, node *nd1, int exceptype, node *nexc, int dcolor, double dscale)

{
    int elemtype,narg2, hide, cmap, excl;
    double vmax, vmin, stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2;

   nd2  = ndt();
   set_display (elemtype=0, disptype, narg2=0,
             nd1, nd2, nexc, exceptype,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap=0, vmax=0, vmin=0,
             dscale, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, 
             node *nd1, int exceptype, node *nexc, int dcolor, double dscale)

{
    int narg2, hide, cmap, excl;
    double vmax, vmin, stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2;

   nd2  = ndt();
   set_display (elemtype, disptype, narg2=0,
             nd1, nd2, nexc, exceptype,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap=0, vmax=0, vmin=0,
             dscale, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, node *nd1, double (*vpenn)(int elnum, int dcolor))

{
    int elemtype,narg2, dcolor, exceptype, hide, cmap, excl;
    double dscale, vmax, vmin, stime;
    Symbol *vpen;
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype=0, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn,
             cmap=0, vmax=0, vmin=0,
             dscale=1, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, node *nd1, double (*vpenn)(int elnum, int dcolor))

{
    int narg2, dcolor, exceptype, hide, cmap, excl;
    double dscale, vmax, vmin, stime;
    Symbol *vpen;
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn,
             cmap=0, vmax=0, vmin=0,
             dscale=1, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, node *nd1, int dcolor, double vmax, double vmin)

{
    int elemtype,narg2, exceptype, hide, cmap, excl;
    double dscale, stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype=0, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap=0, vmax, vmin,
             dscale=1, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, node *nd1, int dcolor, double vmax, double vmin)

{
    int narg2, exceptype, hide, cmap, excl;
    double dscale, stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap=0, vmax, vmin,
             dscale=1, hide=0, excl=0, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, node *nd1, int excl, int dcolor, 
			double vmax, double vmin, int cmap, double dscale)

{
    int elemtype,narg2, exceptype, hide;
    double stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype=0, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap, vmax, vmin,
             dscale, hide=0, excl, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, node *nd1, int excl, int dcolor, 
			double vmax, double vmin, int cmap, double dscale)

{
    int narg2, exceptype, hide;
    double stime;
    Symbol *vpen;
    double (*vpenn)(int,int);
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn=NULL,
             cmap, vmax, vmin,
             dscale, hide=0, excl, stime=0);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, 
             node *nd1,
             int dcolor, double(*vpenn)(int elnum, int color),
             int cmap, double vmax,double vmin,
             double dscale, int hide, int excl, double stime)

{
    int elemtype,narg2, exceptype;
    Symbol *vpen;
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype=0, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn,
             cmap, vmax,vmin,
             dscale, hide, excl, stime);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, 
             node *nd1,
             int dcolor, double(*vpenn)(int elnum, int color),
             int cmap, double vmax,double vmin,
             double dscale, int hide, int excl, double stime)

{
    int narg2, exceptype;
    Symbol *vpen;
    node *nd2, *nexc;

   nd2  = ndt();
   nexc = ndt();
   set_display (elemtype, disptype, narg2=0,
             nd1, nd2, nexc, exceptype=0,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn,
             cmap, vmax,vmin,
             dscale, hide, excl, stime);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int disptype, 
             node *nd1, node *nd2, node *nexc, int exceptype,
             int dcolor, double(*vpenn)(int, int),
             int cmap, double vmax,double vmin,
             double dscale, int hide, int excl, double stime)

{
    int elemtype,narg2;
    Symbol *vpen;

   set_display (elemtype=0, disptype, narg2=1,
             nd1, nd2, nexc, exceptype,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn,
             cmap, vmax,vmin,
             dscale, hide, excl, stime);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int disptype, 
             node *nd1, node *nd2, node *nexc, int exceptype,
             int dcolor, double(*vpenn)(int, int),
             int cmap, double vmax,double vmin,
             double dscale, int hide, int excl, double stime)

{
    int narg2;
    Symbol *vpen;

   set_display (elemtype, disptype, narg2=1,
             nd1, nd2, nexc, exceptype,
             zrange1, zrange2,
             dcolor, vpen=NULL, vpenn,
             cmap, vmax,vmin,
             dscale, hide, excl, stime);
  ndt_free(nd1);
  ndt_free(nd2);
  ndt_free(nexc);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int elnum, 
             int dcolor, double(*vpenn)(int, int),
             int cmap, double vmax,double vmin,
             double dscale, int hide)

{
    int narg2;
    Symbol *vpen=NULL;

  if (!makestim) 
     if (elemtype==ELEMENT) {
        ncdispe (elnum,zrange1,zrange2,dcolor,vpen,vpenn,vmax,vmin,
                 dscale,hide,cmap);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (int elemtype, int elnum, int dcolor, double dscale)

{
    int hide,cmap;
    Symbol *vpen;
    double vmax, vmin;
    double(*vpenn)(int, int);

  if (!makestim) 
    if (elemtype==ELEMENT) {
       ncdispe (elnum,zrange1,zrange2,dcolor,vpen=NULL,vpenn=NULL,vmax=0,vmin=0,
                                  dscale,hide=0,cmap=0);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display (elem *epnt, int dcolor, double dscale)

{
    int hide,cmap;
    Symbol *vpen;
    double vmax, vmin;
    double(*vpenn)(int, int);

  if (!makestim) ncdispe (epnt,dcolor,vpen=NULL,vpenn=NULL,vmax=0,vmin=0,dscale,hide=0,cmap=0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display_z (double zmax, double zmin)

{
  zrange1 = zmin;
  zrange2 = zmax;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display_stim (double sttime, double attime, double dscale, int cmap, double max_flux, double min_flux)

{
//    int narg2, dcolor, exceptype, hide, cmap, excl;
//    double vmax, vmin, stime;
//    Symbol *vpen;
//    double (*vpenn)(int,int);
//    node *nd1, *nd2, *nexc;
//
//   nd1  = ndt();
//   nd2  = ndt();
//   nexc = ndt();
//   set_display (ELEMENT, STIM, narg2=0,
//             nd1, nd2, nexc, exceptype,
//	     zrange1, zrange2,
//             dcolor=0, vpen=NULL, vpenn=NULL,
//             cmap=0, vmax=max_flux, vmin=min_flux,
//             dscale, hide=0, excl=0,stime=attime);
//  ndt_free(nd1);
//  ndt_free(nd2);
//  ndt_free(nexc);
//
//  more direct to call dispstim directly:
//
  dispstim(sttime,attime,dscale,cmap,max_flux,min_flux);

}
/*---------------------------------------------------------------------*/

void display_stim (double sttime, double attime, double dscale, double max_flux, double min_flux)
{
    int cmap;

  display_stim(sttime,attime,dscale,cmap=0,max_flux,min_flux);
}

/*---------------------------------------------------------------------*/

void display_stim (double sttime, double attime, double dscale)
{
   display_stim (sttime, attime, dscale, -LARGENODE, LARGENODE);
}

/*---------------------------------------------------------------------*/

void display_stim (double attime, double dscale, double max_flux, double min_flux)
{
   display_stim (simtime, attime, dscale, max_flux, min_flux);
}

/*---------------------------------------------------------------------*/

void display_stim (double attime, double dscale)
{
   display_stim (simtime, attime, dscale, -LARGENODE, LARGENODE);
}

/*---------------------------------------------------------------------*/

void display_page(void)

{
  gpage();
}
/*---------------------------------------------------------------------*/

void comp_file(const char *filename)

{
#define  CPFILNAM 100
  static char compfile[CPFILNAM]={0};
  FILE *ftemp;
  
  if (*compfile) fclose(compout);
  *compfile = 0;
  strcpy (compfile,filename);
  if (*compfile) {
       if ((ftemp=fopen(compfile,"w")) == NULL) {
         ncfprintf (stderr,"%s: can't open file '%s'\n", progname,compfile);
         return;
       }
       else compout = ftemp;
  }
  else {
        compout = stdout;
  }
  prcomps = 1;
}

/*---------------------------------------------------------------------*/

chattrib *make_chan (elem *epnt, int ctype, int stype)

/* NB: to make Ca channel, need to cast returned value to "cattrib *" */

{
   chattrib *apnt;

 apnt = chanattr(epnt,ctype);
 apnt->stype = stype;
 return (apnt);
}

/* - - - - - - - - - - - - - - - -*/

elem *make_chan (node *npnt, int ctype, int stype)

/* NB: to make Ca channel, need to cast returned value to "cattrib *" */

{
   elem *epnt;
   chattrib *apnt;

 epnt = at(npnt,CHAN); 
 apnt = chanattr(epnt,ctype);
 apnt->stype = stype;
 return (epnt);
}

/* - - - - - - - - - - - - - - - -*/

double get_chan_trconc(chattrib *cpnt)

{
    double trconc;

   if (cpnt->trconc!=NULLVAL) trconc = cpnt->trconc; /* if user set */
   else trconc = getchantype(cpnt->ctype,cpnt->stype)->trconc;
   return (trconc);
}

/* - - - - - - - - - - - - - - - -*/

nattrib *make_chnoise (elem *epnt)

{
   nattrib *npnt;

  npnt = (nattrib *)chanattr(epnt,CCHNOISE);
  return npnt;
}

/* - - - - - - - - - - - - - - - -*/

nattrib *make_vesnoise (elem *epnt)

{
   nattrib *npnt;

  npnt = (nattrib *)chanattr(epnt,VESNOISE);
  return npnt;
}

/* - - - - - - - - - - - - - - - -*/

synapse *make_synapse (node *node1, node *node2)

{
   synapse *spnt; 

  spnt = (synapse *)conn (node1, node2, SYNAPSE);
  return spnt;
}

/* - - - - - - - - - - - - - - - -*/

synapse *make_synapse (node *node1, node *node2, int connum)

{
   synapse *spnt; 

  spnt = (synapse *)conn (node1, node2, SYNAPSE);
  spnt->connum = connum;
  return spnt;
}

/* - - - - - - - - - - - - - - - -*/

gapjunc *make_gj (node *node1, node *node2, double maxcond)

{
   gapjunc *g; 

  g = (gapjunc *)conn (node1, node2, GJ);
  g->gmax = maxcond;
  return (gapjunc *)g;
}

/* - - - - - - - - - - - - - - - -*/

gapjunc *make_pnx (node *node1, node *node2, double maxcond)

{
   gapjunc *g; 

  g = (gapjunc *)conn (node1, node2, PNX);
  g->gmax = maxcond;
  return (gapjunc *)g;
}

/* - - - - - - - - - - - - - - - -*/

vbuf *make_vbuf (node *node1, node *node2, double offset, double gain, double tau, double delay, int stype)

{
   vbuf *b; 

  b = (vbuf *)conn (node1, node2, BUF);
  b->offset = offset;
  b->gain = gain;
  b->delay = delay;
  b->tau = tau;
  b->lphp = stype;
  return (vbuf *)b;
}

/* - - - - - - - - - - - - - - - -*/

vbuf *make_vbuf (node *node1, node *node2, double offset, double gain, double tau, double delay)
{
   return make_vbuf (node1, node2, offset, gain, tau, delay, LP);
}

/* - - - - - - - - - - - - - - - -*/

nbuf *make_nbuf (node *node1, node *node2, double offset, double gain, double ntoffset, int ntrans)
{
   nbuf *n; 

  n = (nbuf *)conn (node1, node2, NBUF);
  n->offset = offset;
  n->ntoffset = ntoffset;
  n->gain = gain;
  n->ntrans = ntrans;
  return (nbuf *)n;
}

/* - - - - - - - - - - - - - - - -*/

double *make_filt(double val)

{
    double *filt;

  filt = (double *)emalloc(sizeof(double));
  filt[0] = val;
  return filt;
}
  
/* - - - - - - - - - - - - - - - -*/

double *make_filt(double val1,double val2)

{
    double *filt;

  filt = (double *)emalloc(2*sizeof(double));
  filt[0] = val1;
  filt[1] = val2;
  return filt;
}
  
/* - - - - - - - - - - - - - - - -*/

double *make_filt(double val1,double val2,double val3)

{
    double *filt;

  filt = (double *)emalloc(3*sizeof(double));
  filt[0] = val1;
  filt[1] = val2;
  filt[2] = val3;
  return filt;
}
  
/* - - - - - - - - - - - - - - - -*/

double *make_filt(double val1,double val2,double val3,double val4)

{
    double *filt;

  filt = (double *)emalloc(4*sizeof(double));
  filt[0] = val1;
  filt[1] = val2;
  filt[2] = val3;
  filt[3] = val4;
  return filt;
}

/* - - - - - - - - - - - - - - - -*/

double *make_filt(double val1,double val2,double val3,double val4,double val5)

{
    double *filt;

  filt = (double *)emalloc(5*sizeof(double));
  filt[0] = val1;
  filt[1] = val2;
  filt[2] = val3;
  filt[3] = val4;
  filt[4] = val5;
  return filt;
}
 
/* - - - - - - - - - - - - - - - -*/

elem *foreach (elem *ept, int etype)

{
    int efilt;
    register elem *epnt=ept;

  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SPHERE:
     case SYNAPSE:
     case CHAN:
     case ROD:
     case CONE:
     case CHR:
     case VTRANSDUCER:
     case ITRANSDUCER:
     case GJ:
     case PNX:
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:	 efilt = 2;
                        break;

     default:           efilt = 0;
                        break;
  }
  if (efilt==2) 
    for (; epnt; epnt=epnt->next) { 			/* look for elem type */
      if (epnt->ctype==etype) break;
    }
  return epnt;
}

/* - - - - - - - - - - - - - - - -*/

elem *foreach (elem *ept, int etype, double radius, 
			double (distfunc)(elem *e1, node *n2), node *npntw)
{
    int efilt, match;
    register elem *epnt=ept;

  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SPHERE:
     case SYNAPSE:
     case CHAN:
     case ROD:
     case CONE:
     case CHR:
     case VTRANSDUCER:
     case ITRANSDUCER:
     case GJ:
     case PNX:
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:	 efilt = 2;
                        break;

     default:           efilt = 0;
                        break;
  }
  if (efilt==2) 
    for (; epnt; epnt=epnt->next) { 			/* look for elem type */
      if (epnt->ctype!=etype) continue;
      if (epnt->modif) continue;
      match = 1;
      if (distfunc(epnt,npntw) > radius) continue;
    }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na1, int na2, int na3, int na4,
				      int nb1, int nb2, int nb3, int nb4, int xxx)

/* find elements that match 2 nodes */

{
    int efilt;
    node *npnt1, *npnt2;
 
  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SYNAPSE:
     case GJ:
     case PNX:
     case RESISTOR:
     case CAP:
     case BATT:
     case BUF:		efilt = 2;

     default:           efilt = 0;
                        break;
  }
  for (; epnt; epnt=epnt->next) { 			/* look for elem type */
    switch (epnt->ctype) {

      case ROD:
      case CONE:
      case CHR:
      case VTRANSDUCER:
      case ITRANSDUCER:
      case CHAN:
      case SPHERE:
      case LOAD:
      case GNDCAP:
      case GNDBATT:
      case CACOMP: continue;
                   break;
      default:     break;
    }
    if (efilt==2) if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    if ((npnt1=epnt->nodp1)!=NULL && (npnt2=epnt->nodp2)!=NULL) { 

      if (na2>=0 && na2!=npnt1->nodenm2) continue;
      if (na1>=0 && na1!=npnt1->nodenm1) continue;
      if (na3>=0 && na3!=npnt1->nodenm3) continue;
      if (na4>=0 && na4!=npnt1->nodenm4) continue;

      if (nb2>=0 && nb2!=npnt2->nodenm2) continue;
      if (nb1>=0 && nb1!=npnt2->nodenm1) continue;
      if (nb3>=0 && nb3!=npnt2->nodenm3) continue;
      if (nb4>=0 && nb4!=npnt2->nodenm4) continue;
      break;
    }
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na1, int na2, int na3, int nb1, int nb2, int nb3, int xxx) 

{
   int na4, nb4;
   return (foreach (epnt, etype, na1, na2, na3, na4=-1, nb1, nb2, nb3, nb4=-1, xxx));
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na1, int na2, int nb1, int nb2, int xxx) 

{
   int na3, na4, nb3, nb4;
   return (foreach (epnt, etype, na1, na2, na3=-1, na4=-1, nb1, nb2, nb3=-1, nb4=-1, xxx));
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *ept, int etype, int a, int b, int c, int d) 
{
    int efilt, match;
    register node *npnt;
    register elem *epnt=ept;
    register int na=a, nb=b, nc=c, nd=d;

  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SPHERE:
     case SYNAPSE:
     case CHAN:
     case ROD:
     case CONE:
     case CHR:
     case VTRANSDUCER:
     case ITRANSDUCER:
     case GJ:
     case PNX:
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:	 efilt = 2;
                        break;

     default:           efilt = 0;
                        break;
  }
  for (; epnt; epnt=epnt->next) { 				/* look for elem type */
    if (efilt==2) if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    if (npnt=epnt->nodp1) {					/* test first node */
      match = 1;
      if      (nb>=0 && nb!=npnt->nodenm2) match=0;
      else if (na>=0 && na!=npnt->nodenm1) match=0;
      else if (nc>=0 && nc!=npnt->nodenm3) match=0;
      else if (nd>=0 && nd!=npnt->nodenm4) match=0;
    }
    if (match) break; 
    else {
      if (npnt=epnt->nodp2) {					/* test second node */
        if      (nb>=0 && nb!=npnt->nodenm2) continue;
        else if (na>=0 && na!=npnt->nodenm1) continue;
        else if (nc>=0 && nc!=npnt->nodenm3) continue;
        else if (nd>=0 && nd!=npnt->nodenm4) continue;
        break;
      }
    } 
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na, int nb, int nc) 

{
   int nd;
   return (foreach (epnt, etype, na, nb, nc, nd=-1));
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na, int nb) 

{
   int nc,nd;
   return (foreach (epnt, etype, na, nb, nc=-1, nd=-1));
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *ept, int etype, int a, int b, int c, int d, 
				      int *pa, int *pb, int *pc, int *pd)
{
    int efilt, match;
    register node *npnt;
    register elem *epnt=ept;
    register int na=a, nb=b, nc=c, nd=d;

  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SPHERE:
     case SYNAPSE:
     case CHAN:
     case ROD:
     case CONE:
     case CHR:
     case VTRANSDUCER:
     case ITRANSDUCER:
     case GJ:
     case PNX:
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:	 efilt = 2;
                        break;

     default:           efilt = 0;
                        break;
  }
  for (; epnt; epnt=epnt->next) { 				/* look for elem type */
    if (efilt==2) if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    match = 0;
    if (npnt=epnt->nodp1) {					/* test first node */
      match = 1;
      if      (nb>=0 && nb!=npnt->nodenm2) match=0;
      else if (na>=0 && na!=npnt->nodenm1) match=0;
      else if (nc>=0 && nc!=npnt->nodenm3) match=0;
      else if (nd>=0 && nd!=npnt->nodenm4) match=0;
    }
    if (match) {
        if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
        if (pb!=NULL) *pb=npnt->nodenm2;
        if (pc!=NULL) *pc=npnt->nodenm3;
        if (pd!=NULL) *pd=npnt->nodenm4;
        break;
    }
    else {
      if (npnt=epnt->nodp2) {					/* test second node */
        if      (nb>=0 && nb!=npnt->nodenm2) continue;
        else if (na>=0 && na!=npnt->nodenm1) continue;
        else if (nc>=0 && nc!=npnt->nodenm3) continue;
        else if (nd>=0 && nd!=npnt->nodenm4) continue;
                    
        if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
        if (pb!=NULL) *pb=npnt->nodenm2;
        if (pc!=NULL) *pc=npnt->nodenm3;
        if (pd!=NULL) *pd=npnt->nodenm4;
	break;
      }
    } 
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na, int nb, int nc, int *pa, int *pb, int *pc)

{
   int nd;
   int *pd;

  return (foreach (epnt, etype, na, nb, nc, nd=-1, pa, pb, pc, pd=NULL));
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *epnt, int etype, int na, int nb, int *pa, int *pb)

{
   int  nc,  nd;
   int *pc, *pd;

  return (foreach (epnt, etype, na, nb, nc=-1, nd=-1, pa, pb, pc=NULL, pd=NULL));
}

/* - - - - - - - - - - - - - - - - */

elem *foreachn (elem *ept, int etype, int a, int b, node **pn)

/* Use hash index of (ct,cn) to find next matching element for connect_cell(). */
/*  Epnt must be set to the entry in enhashtab[] corresponding to (ct,cn) */
/* Have removed testing nc,nd, pa,pb because these aren't used by connect_cell(). */
/* Assumes na,nb are non-negative. */

{
   register elem *epnt=ept;
   register int na=a, nb=b;

//    int match;
 

//  if (nc >= 0 || nd >= 0) 			/* if nc or nd are specified, revert to std foreach() */
//        return foreach (epnt, etype, na, nb, nc, pc);

  for (; epnt; epnt=epnt->hnnext) { 		/* look for elem type */
    if (epnt->ctype!=etype) continue;
    // match = 1;
    if      (nb!=epnt->node1b) continue;	/* test first node */
    else if (na!=epnt->node1a) continue;
    if (epnt->modif) continue;
 
    *pn = epnt->nodp1;				/* return matching node */
    break;
//    else {						/* test second node */
//        if      (nb!=epnt->node2b) continue;
//        else if (na!=epnt->node2a) continue;
//                    
//	*pn = epnt->nodp2;				/* return matching node */
//	break;
//    } 
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

int foreachn_etype;
int foreachn_a;
int foreachn_b;
node **foreachn_pn;

elem *setforeachn2 (int etype, int a, int b, node **pn)

/* save args for foreachn(elem *) below. */
/*  ct = a, cn = b */

{
  foreachn_etype = etype;
  foreachn_a = a;
  foreachn_b = b;
  foreachn_pn = pn;
  return getepnt(a,b);
}

/* - - - - - - - - - - - - - - - - */

elem *foreachn2 (elem *ept)

/* Use args from foreachn(elem *) above to find next matching elem. */

{
   register elem *epnt=ept;
   register int etype = foreachn_etype;
   register int na = foreachn_a;
   register int nb = foreachn_b;

  for (; epnt; epnt=epnt->hnnext) { 		/* look for elem type */
    if (epnt->ctype!=etype) continue;
    if      (nb!=epnt->node1b) continue;	/* test first node */
    else if (na!=epnt->node1a) continue;
    if (epnt->modif) continue;
 
    *foreachn_pn = epnt->nodp1;			/* return matching node */
    break;
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *ept, int etype, int a, int b, int c, int d, 
				      int *pa, int *pb, int *pc, int *pd,
					double radius, 
					double (distfunc)(node *n1, node *n2), node *npntw)
{
    int efilt, match;
    register node *npnt;
    register elem *epnt=ept;
    register int na=a, nb=b, nc=c, nd=d;

 
  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SPHERE:
     case SYNAPSE:
     case CHAN:
     case ROD:
     case CONE:
     case CHR:
     case VTRANSDUCER:
     case ITRANSDUCER:
     case GJ:
     case PNX:
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:	 efilt = 2;
                        break;

     default:           efilt = 0;
                        break;
  }
  for (; epnt; epnt=epnt->next) { 				/* look for elem type */
    if (efilt==2) if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    match = 0;
    if (npnt=epnt->nodp1) {					/* test first node */
      match = 1;
      if      (nb>=0 && nb!=npnt->nodenm2) match=0;
      else if (na>=0 && na!=npnt->nodenm1) match=0;
      else if (nc>=0 && nc!=npnt->nodenm3) match=0;
      else if (nd>=0 && nd!=npnt->nodenm4) match=0;
    }
    if (distfunc(npnt,npntw) > radius) match = 0;
    if (match) {
        if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
        if (pb!=NULL) *pb=npnt->nodenm2;
        if (pc!=NULL) *pc=npnt->nodenm3;
        if (pd!=NULL) *pd=npnt->nodenm4;
	break;
    }
    else {
      if (npnt=epnt->nodp2) {					/* test second node */
        if      (nb>=0 && nb!=npnt->nodenm2) continue;
        else if (na>=0 && na!=npnt->nodenm1) continue;
        else if (nc>=0 && nc!=npnt->nodenm3) continue;
        else if (nd>=0 && nd!=npnt->nodenm4) continue;
      
        if (distfunc(npnt,npntw) <= radius)
        {
            if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
            if (pb!=NULL) *pb=npnt->nodenm2;
            if (pc!=NULL) *pc=npnt->nodenm3;
            if (pd!=NULL) *pd=npnt->nodenm4;
	    break;
        }
      }
    } 
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

elem *foreachn (elem *ept, int etype, int a, int b,
					double radius, 
					double (distfunc)(node *n1, node *n2), node *npntw, node **pn)
{
   register elem *epnt=ept;
   register int na=a, nb=b;

  // if (nc >= 0 || nd >= 0) 			/* if nc or nd are specified, revert to std foreach() */
  //       return foreach (epnt, etype, na, nb, pc, radius, distfunc, npntw);

  for (; epnt; epnt=epnt->hnnext) { 				/* look for elem type */
    if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    if      (nb!=epnt->node1b) continue;				/* test first node */
    else if (na!=epnt->node1a) continue;
//    else {
//       if      (nb!=epnt->node2b) continue;			/* test second node */
//       else if (na!=epnt->node2a) continue;
//       {
//	    *pn = epnt->nodp2;
//       }
//    } 
    if (distfunc(epnt->nodp1,npntw) > radius) continue;
    *pn = epnt->nodp1;
    break;
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

int foreachn3_etype;
int foreachn3_a;
int foreachn3_b;
double (*foreachn3_distfunc)(node *n1, node*n2);
double foreachn3_radius;
node *foreachn3_npntw;
node **foreachn3_pn;

elem *setforeachn3 (int etype, int a, int b, double radius, 
			double (distfunc)(node *n1, node *n2), node *npntw, node **pn)

/* save args for foreachn(elem *) below. */
/*  ct = a, cn = b */

{
  foreachn3_etype = etype;
  foreachn3_a = a;
  foreachn3_b = b;
  foreachn3_radius = radius;
  foreachn3_distfunc = distfunc;
  foreachn3_npntw = npntw;
  foreachn3_pn = pn;
  return getepnt(a,b);
}

/* - - - - - - - - - - - - - - - - */

elem *foreachn3 (elem *ept)

/* use args from setforeachn3() above. */

{
    register elem *epnt=ept;
    register int etype = foreachn3_etype;
    register int na=foreachn3_a;
    register int nb=foreachn3_b;
    register node *npntw = foreachn3_npntw;
    register double radius = foreachn3_radius;

  for (; epnt; epnt=epnt->hnnext) { 				/* look for elem type */
    if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    if      (nb!=epnt->node1b) continue;				/* test first node */
    else if (na!=epnt->node1a) continue;
    if (foreachn3_distfunc(epnt->nodp1,npntw) > radius) continue;
    *foreachn3_pn = epnt->nodp1;
    break;
  }
  return epnt;
}


/* - - - - - - - - - - - - - - - - */

elem *foreach (elem *ept, int etype, int a, int b, int c, int d, 
				      int *pa, int *pb, int *pc, int *pd,
					double radius, 
					double (distfunc)(elem *e1, node *n2), node *npntw)
{
    int efilt, match;
    register node *npnt;
    register elem *epnt = ept;
    register int na=a, nb=b, nc=c, nd=d;

  switch (etype) {

     case ELEMENT:
                        efilt = 1;
                        break;
     case CABLE:
     case SPHERE:
     case SYNAPSE:
     case CHAN:
     case ROD:
     case CONE:
     case CHR:
     case VTRANSDUCER:
     case ITRANSDUCER:
     case GJ:
     case PNX:
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:	 efilt = 2;
                        break;

     default:           efilt = 0;
                        break;
  }
  for (; epnt; epnt=epnt->next) { 				/* look for elem type */
    if (efilt==2) if (epnt->ctype!=etype) continue;
    if (epnt->modif) continue;
    match = 0;
    if (npnt=epnt->nodp1) {					/* test first node */
      match = 1;
      if (na>=0 && na!=npnt->nodenm1) match=0;
      if (nb>=0 && nb!=npnt->nodenm2) match=0;
      if (nc>=0 && nc!=npnt->nodenm3) match=0;
      if (nd>=0 && nd!=npnt->nodenm4) match=0;
    }
    if (distfunc(epnt,npntw) > radius) match = 0;
    if (match) {
      if (npnt) {
        if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
        if (pb!=NULL) *pb=npnt->nodenm2;
        if (pc!=NULL) *pc=npnt->nodenm3;
        if (pd!=NULL) *pd=npnt->nodenm4;
	break;
      }
    }
    else {
      match = 1;
      if (npnt=epnt->nodp2) {					/* test second node */
        if (na>=0 && na!=npnt->nodenm1) match=0;
        if (nb>=0 && nb!=npnt->nodenm2) match=0;
        if (nc>=0 && nc!=npnt->nodenm3) match=0;
        if (nd>=0 && nd!=npnt->nodenm4) match=0;
      }
      if (distfunc(epnt,npntw) > radius) match = 0;
      if (match) {
        if (npnt) {
          if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
          if (pb!=NULL) *pb=npnt->nodenm2;
          if (pc!=NULL) *pc=npnt->nodenm3;
          if (pd!=NULL) *pd=npnt->nodenm4;
	  break;
        }
      }
    } 
  }
  return epnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb, int nc, int nd, 
				      int *pa, int *pb, int *pc, int *pd)
{
    register node *npnt = npt;

  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (nc>=0 && nc!=npnt->nodenm3) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    if (nd>=0 && nd!=npnt->nodenm4) continue;

    if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
    if (pb!=NULL) *pb=npnt->nodenm2;
    if (pc!=NULL) *pc=npnt->nodenm3;
    if (pd!=NULL) *pd=npnt->nodenm4;
    break;
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */


node *foreach (node *npt, int na, int nb, int nc, int *pa, int *pb, int *pc)
{
    register node *npnt = npt;

  for (; npnt; npnt=npnt->next) {
    if (nc>=0 && nc!=npnt->nodenm3) continue;
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;

    if (pa!=NULL) *pa=npnt->nodenm1;
    if (pb!=NULL) *pb=npnt->nodenm2;
    if (pc!=NULL) *pc=npnt->nodenm3;
    break;
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreachn (node *npt, int na, int nb, int *pc)
{
    register node *npnt=npt;

  for (; npnt; npnt=npnt->hcnext) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    if (pc!=NULL) *pc=npnt->nodenm3;
    break; 
  }
  return npnt;
}


/* - - - - - - - - - - - - - - - - */
/*
node *foreach (node *npt, int na, int nb, int nc, int *pa, int *pb, int *pc)
{
    int match;
    register node *npnt = npt;

  for (; npnt; npnt=npnt->next) {
    match = 1;
    if      (nc>=0 && nc!=npnt->nodenm3) match=0;
    else if (nb>=0 && nb!=npnt->nodenm2) match=0;
    else if (na>=0 && na!=npnt->nodenm1) match=0;
    if (match) {
       if (pa!=NULL) *pa=npnt->nodenm1;
       if (pb!=NULL) *pb=npnt->nodenm2;
       if (pc!=NULL) *pc=npnt->nodenm3;
       break;
    }
  }
  return npnt;
}
*/
/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb, int *pa, int *pb)
{
    register int match;
    register node *npnt = npt;

  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;

    if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
    if (pb!=NULL) *pb=npnt->nodenm2;
    break;
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb, int nc, int nd, 
		      int *pa, int *pb, int *pc, int *pd,
			double radius, 
			double (distfunc)(node *n1, node *n2), node *npntw)
{
    register node *npnt=npt;

  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    if (nc>=0 && nc!=npnt->nodenm3) continue;
    if (nd>=0 && nd!=npnt->nodenm4) continue;
    if (distfunc(npnt,npntw) > radius) continue;

    if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
    if (pb!=NULL) *pb=npnt->nodenm2;
    if (pc!=NULL) *pc=npnt->nodenm3;
    if (pd!=NULL) *pd=npnt->nodenm4;
    break;
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb, int nc, 
		      int *pa, int *pb, int *pc, 
			double radius, 
			double (distfunc)(node *n1, node *n2), node *npntw)
{
    register node *npnt=npt;
 
  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    if (nc>=0 && nc!=npnt->nodenm3) continue;
    if (distfunc(npnt,npntw) > radius) continue;

    if (pa!=NULL) *pa=npnt->nodenm1;			/* store node number if requested */
    if (pb!=NULL) *pb=npnt->nodenm2;
    if (pc!=NULL) *pc=npnt->nodenm3;
    break;
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb, int nc, int nd)
{
    register node *npnt=npt;

  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    if (nc>=0 && nc!=npnt->nodenm3) continue;
    if (nd>=0 && nd!=npnt->nodenm4) continue;
    break; 
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb, int nc)
{
    register node *npnt=npt;

  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    if (nc>=0 && nc!=npnt->nodenm3) continue;
    break; 
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreach (node *npt, int na, int nb)
{
    register node *npnt=npt;

  for (; npnt; npnt=npnt->next) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    break; 
  }
  return npnt;
}

/* - - - - - - - - - - - - - - - - */

node *foreachn(node *npt, int na, int nb)
{
    register node *npnt=npt;

  for (; npnt; npnt=npnt->hcnext) { 				/* look for node */
    if (nb>=0 && nb!=npnt->nodenm2) continue;
    if (na>=0 && na!=npnt->nodenm1) continue;
    break; 
  }
  return npnt;
}

/*-----------------------------------------------------------*/

chan *findchan (conlst *cpnt, int ctype, int stype)

/* find the first channel of a given type (ctype, stype) located at a compartment */

{
  register conlst *c;
  register chan *p;

  for (p=NULL,c=cpnt; c!=NULL; c=c->next) {
    if ((p=(chan *)c->conpnt)!=NULL) {
       // fprintf (stderr,"chan %d\n",p->ctype);
       if (p->ctype==ctype && (stype<0 || p->stype==stype)) {
         break;
       }
    }
  }
  if (c==NULL) return (chan *)c;
  else return p;
}


/*-----------------------------------------------------------*/

#include "notinit.cc"

/*-----------------------------------------------------------*/

int streq (const char *str1, const char *str2)
{
    return (!strcmp(str1,str2));
}

/*-----------------------------------------------------------*/

// #define CBUFSIZ 16384
// 
// char *system(char *str)
// 
// {
//   int i,c;
//   FILE *cstdout;
//   char *cbuf;
// 
//  if (!(cbuf=emalloc(CBUFSIZ))) {
//     ncfprintf (stderr,"can't allocate char buf\n");
//     return NULL;
//  }
// 
//  if (cstdout=popen(str,"r")) {
//    for (i=0; i<10; i++) cbuf[i] = 0;
//    for (i=0;(c = getc(cstdout)) != EOF && i<CBUFSIZ-1; i++) {
//      cbuf[i] = c;               /* copy output of pipe into buffer */
//    }
//    cbuf[i++] = (char)NULL;
// 
//    pclose(cstdout);
//    return (cbuf);
//  }
// }

/*-----------------------------------------------------------*/

double set_int_val (double val)

{
    double logval, tval;

  logval = floor(log10(val));
  tval = int (val / exp(logval*M_LN10) + 0.5);
  if      (tval > 0.8 && tval <= 1.2) tval = 1.0;
  else if (tval > 1.2 && tval <= 3.1) tval = 2.0;
  else if (tval > 3.1 && tval <= 8.0) tval = 5.0;
  else tval = 10.0;
  val = tval * exp(logval*M_LN10);
  return val;
}


