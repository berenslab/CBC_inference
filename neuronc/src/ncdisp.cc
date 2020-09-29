/* Segment ncdisp in Program nc */

/* draws anatomy of neuronal circuits */

/*      Nov 90                  R.G. Smith */

extern "C" {

#include <stdio.h>
#include <stdlib.h>
#ifdef CPML
#include <cpml.h>
#endif
#include <math.h>
#include <string.h>
#include "stdplt.h"

}

#include "nc.h"
#include "ncsub.h"
#include "ndef.h"
#include "ncelem.h"
#include "ncomp.h"
#include "y.tab.h"
#include "control.h"
#include "colors.h"
#include "gprim.h"
#include "ncio.h"

#ifndef PI
#define PI	3.14159265358979323846264
#endif

#define DEBUG

#ifdef DEBUG
#include "ncdebug.h"
#endif

double *rcolors=NULL;

extern double ncmat[4][4];		/* scale, rotation matrix for display */
extern double ncmatc[4][4];		/* scale, rotation matrix for camera */
double xyscal=1.0;			/* scale for drawing pictures 0 -> 1 */
double zscal=0;				/* z scale for drawing pictures */
extern double linewidth;		/* set by ncmain */
extern int disp_ray;			/* =1 -> render 3D ray-tracing */

extern elem *elempnt;			/* pointer to element list */
extern node *nodepnt;			/* pointer to node list */
extern recstim *recspnt;		/* pointer to recstim list */
extern recstim *recsend;		/* pointer to end of recstim list */
extern photrec *recpnt;			/* pointer to receptor list */
extern char *progname;			/* name of program */

extern int ncerror;                     /* error flag set by "warning()" */
extern int runyet;			/* step or run statement has run */

char *findsym(int num);
double ncabs(double x);
void ncdraw (elem *epnt, int dcolor, Symbol *vpen, double(*vpenn)(int elnum,int color),
		double vmax, double vmin, double dscale, int hide, int cmap);
void purge (void);
void drnode(node *npnt, double dscale, int color, double (*mat)[4]);
photrec *makr(photorec *epnt, int ctype);
int readstim(double stoptime);
void execerror(const char *s, const char*t);
void drsphere (double x, double y, double z, double dia, double dscale,
	int color, int hide, double (*mat)[4], int fill);
void drcable (double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	double dia, double dia2, double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4]);
void drsynap (double x1, double y1, double z1, 
	double x2, double y2, double z2, synapse *epnt,
	double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4]);
void drconn (double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	int ctype, double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4]);
void drload (double x, double y, double z, int ctype, 
	double n1dia, double n2dia, int color, int hide, double (*mat)[4]);
void drphotrec (double x, double y, double z, 
	     photorec *epnt, double dscale,
	     double n1dia,double n2dia,
	     int color, int hide, double (*mat)[4]);
void drelec (double x, double y, double z, 
             electrode *epnt, double dscale,
             int color, int hide, double (*mat)[4]);
void raysphere (double x, double y, double z, double dia, double dscale, 
	int color, int hide, double (*mat)[4], int fill);
void raycable (double x1, double y1, double z1, 
		double x2, double y2, double z2, 
		double dia, double dia2, double dscale, 
		double n1dia, double n2dia, 
		int color, int hide, double (*mat)[4]);
void raysynap (double x1, double y1, double z1, 
	double x2, double y2, double z2, synapse *epnt,
	double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4]);
void rayconn (double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	int ctype, double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4]);
void rayload (double x, double y, double z, int ctype, 
	double n1dia, double n2dia, int color, int hide, double (*mat)[4]);
void rayphotrec (double x, double y, double z, 
		     photorec *epnt, double dscale,
		     double n1dia,double n2dia,
		     int color, int hide, double (*mat)[4]);
void rayelec (double x, double y, double z, 
		     electrode *epnt, double dscale,
		     int color, int hide, double (*mat)[4]);
void raynode(node *npnt, double dscale, int color, double (*mat)[4]);

void transf (double x,double y,double z,
	      double *tx,double *ty,double *tz, double (*mat)[4]); 
void ppurge();
void elec_dr (int color, double dscale, double dia, double length, 
		double foreshorten, int hide, int elnum);
void phot_dr (int type, int pigm, int color, double dscale,
		double dia, double foreshorten, int hide);
void synapse_dr (synapse *, int color, double vrev, double dscale, double dia, 
		double length, double foreshorten, int hide);
void gapjunc_dr (int color, double length, double dia, 
		double foreshorten, int hide);
double callfunc(Symbol *funcp, int npar, double par1, double par2, 
					 double par3, double par4);
void callfunc8(Symbol *funcp, int npar, double par1, double par2, 
			double par3, double par4, double par5, 
			double par6, double par7, double par8);
int loccheck (double x, double y, double z);	/* check if valid location */
void raycalib(double x, double y, double len, double size, int color);
char *prnode(int n1, int n2, int n3, int n4);
void delphotrec (photrec *rpnt);
void efree(void *ptr);
char *emalloc(unsigned int size);
double rechan(conn *elpnt, int nf, int pval, int cnum);
double recsynap(synap *spnt, int nf, int snum);
elem *findelem(int num);
void delrstim(recstim *rspnt);
elem *getepnt(nodeint nodenm1, nodeint nodenm2);
elem *foreachn (elem *epnt, int etype, int na, int nb, node **ppn);
node *foreachn (node *npnt, nodeint nodenm1, nodeint nodenm2);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, const char *str);
node *findnode(nodeint node1a, nodeint node1b, const char *str);
node *getnpnt(nodeint nodenm1, nodeint nodenm2);

/*------------------------------------*/

/* Definitions for user-defined "draw icon" functions". */

static Symbol *elec_drpnt = NULL;
static void (*rhd)(int,double,double,double,double,int,int) = NULL;
static Symbol *phot_drpnt = NULL;
static void (*phd)(int,int,int,double,double,double,int) = NULL;
static Symbol *syn_drpnt = NULL;
static void (*syd)(synapse*,int,double,double,double,double,double,int) = NULL;
static Symbol *gapj_drpnt = NULL;
static void (*gjd)(int,double,double,double,int) = NULL;

/*------------------------------------*/

int checknodez(node *npnt,double zmin,double zmax)
{
    int nomatch;
    double z;

  nomatch = 0;
  if (npnt) {
    z = npnt->zloc;
    if (z != LARGENODE) {
      if (zmin < zmax) {
        if ((zmin > z) || (z > zmax)) nomatch = 1;
      }
      else {
        if ((zmin > z) && (z > zmax)) nomatch = 1;
      }
    }
  }
  return nomatch;
}

/*------------------------------------*/

int checkelemz(elem *epnt,double z1,double z2)
{
     int nomatch;
     double z;	

  nomatch = 0;
  if (checknodez(epnt->nodp1,z1,z2)) nomatch = 1; 
  if (checknodez(epnt->nodp2,z1,z2)) nomatch = 1; 
  return nomatch;
}

/*------------------------------------*/

void ncdisp(int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype, 
	int na, int nb, int nc, int nd, double zrange1, double zrange2, 
	int color, Symbol *vpen, double(*vpenn)(int elnum,int color),
	double vmax, double vmin, double dscale, int hide, int excl, int cmap)

/* Draw neural circuit elements that connect to 
   node1. If any of the node numbers == NULLVAL,
   then ignore this dimension.  If "excl"=1, display
   only elements that connect exclusively to given node.
   If exceptype is set or if the "except" node number is
   set then include appropriate exceptions.
*/

{
   int match,nomatch,exceptf;
   elem *epnt;
   node *npnt;

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		 ncfprintf (stderr,"ncdisp start n1a %d n1b %d n1c %d n1d %d excl %d\n",
				n1a,n1b,n1c,n1d,excl);
#endif

/* ncfprintf (stderr,"ncdisp n1a %d n1b %d n1c %d n1d %d excl %d\n",
					n1a,n1b,n1c,n1d,excl);/* */

   exceptf = (na!=NULLVAL || nb!=NULLVAL || nc!=NULLVAL || nd!=NULLVAL);

/* -  -  -  -  -  -  -  -  -  -  -  -  - */

   if (excl==0 && exceptf==0 && 
		 n1a >= 0 && 
		 n1b >= 0 && 
		 n1c <= 0 && 
		 n1d==NULLVAL) {	/* run faster loop */

  for (epnt=getepnt(n1a,n1b); epnt=foreachn(epnt,elemtype,n1a,n1b,&npnt); epnt = epnt->hnnext) {

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2))
	ncfprintf (stderr,"ncdisp fast loop %d ct %d cn %d n %d ct2 %d cn2 %d\n",
			epnt->elnum, epnt->node1a, epnt->node1b, epnt->node1c, epnt->node2a, epnt->node2b);
#endif

     if (elemtype  && elemtype!=ELEMENT && epnt->ctype!=elemtype) continue;
     match = 1;
     if (n1b != NULLVAL)			/* check node 1, cn first */
        if (epnt->node1b!=n1b) match=0;
     if (n1a != NULLVAL) 
        if (epnt->node1a!=n1a) match=0;
     if (n1c != NULLVAL)
        if (epnt->node1c!=n1c) match=0;

     if (!match) {
	 match = 1;				 /* else check node 2 */
         if (n1b != NULLVAL)
            if (epnt->node2b!=n1b) continue;
         if (n1a != NULLVAL)
            if (epnt->node2a!=n1a) continue;
         if (n1c != NULLVAL)
            if (epnt->node2c!=n1c) continue;
     } 

     if (!match) continue;
     if (checkelemz(epnt,zrange1,zrange2)) match = 0;
     if (!match) continue;
     ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);
   }


 } else {	// else do full display with excl, n1d

/* -  -  -  -  -  -  -  -  -  -  -  -  - */

   for (epnt=elempnt; epnt; epnt=epnt->next) {	/* search all elems */

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2))
	ncfprintf (stderr,"ncdisp long loop %d ct %d cn %d n %d ct2 %d cn2 %d\n",
			epnt->elnum, epnt->node1a, epnt->node1b, epnt->node1c, epnt->node2a, epnt->node2b);
#endif

     if (elemtype  && elemtype!=ELEMENT && epnt->ctype!=elemtype) continue;
     match = 1;
     if (n1a != NULLVAL) 			/* check node 1 first */	
        if (epnt->node1a!=n1a) match=0;
     if (n1b != NULLVAL)
        if (epnt->node1b!=n1b) match=0;
     if (n1c != NULLVAL)
        if (epnt->node1c!=n1c) match=0;
     if (n1d != NULLVAL)
        if (epnt->node1d!=n1d) match=0;

     if ((!excl && !match) || (excl && (epnt->node2a!=NULLVAL
				    ||  epnt->node2b!=NULLVAL
				    ||  epnt->node2c!=NULLVAL
				    ||  epnt->node2d!=NULLVAL))) {
	 if (!excl) match = 1;				 /* else check node 2 */
         if (n1a != NULLVAL)
            if (epnt->node2a!=n1a) match=0;
         if (n1b != NULLVAL)
            if (epnt->node2b!=n1b) match=0;
         if (n1c != NULLVAL)
            if (epnt->node2c!=n1c) match=0;
         if (n1d != NULLVAL)
            if (epnt->node2d!=n1d) match=0;
     } 

     if (!match) continue;
     nomatch = 0;
     if (exceptf) {
       nomatch = 1;
       if (na!=NULLVAL && epnt->node1a!=na) nomatch=0;     
       if (nb!=NULLVAL && epnt->node1b!=nb) nomatch=0;     
       if (nc!=NULLVAL && epnt->node1c!=nc) nomatch=0;     
       if (nd!=NULLVAL && epnt->node1d!=nd) nomatch=0;     

       if (!nomatch) {
       nomatch = 1;
         if (na!=NULLVAL && epnt->node2a!=na) nomatch=0;     
         if (nb!=NULLVAL && epnt->node2b!=nb) nomatch=0;     
         if (nc!=NULLVAL && epnt->node2c!=nc) nomatch=0;     
         if (nd!=NULLVAL && epnt->node2d!=nd) nomatch=0;     
	}
       if (nomatch) 
           if (exceptype && epnt->ctype!=exceptype) nomatch=0;
      }
      else if (exceptype && epnt->ctype==exceptype) nomatch=1;

     if (checkelemz(epnt,zrange1,zrange2)) nomatch = 1;
     if (nomatch) continue;
     ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);
   }
 }

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdisp end\n");
#endif
}

/*------------------------------------*/

void ncdispc(int n1a, int n1b, int n1c, int n1d,
	int n2a, int n2b, int n2c, int n2d, int elemtype, int exceptype, 
	int na, int nb, int nc, int nd, double zrange1, double zrange2, 
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
	double vmax, double vmin, double dscale, int hide, int cmap)

/* Draw anatomy of neural circuit that connects to 
   both node1 and node2. If any of the node
   numbers == NULLVAL, then ignore this dimension.
 */

{
   elem *epnt;
   int n1pri, n2pri, compat,exceptf;
   int n1_used, match, nomatch;
   int qn1a, qn1b, qn1c, qn1d, qn2a, qn2b, qn2c, qn2d;	/* nodes from query */

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		 ncfprintf (stderr,"ncdispc n1a %d n1b %d n1c %d n1d %d\n",
			n1a,n1b,n1c,n1d);
#endif

/* ncfprintf (stderr,
	"ncdisp n1a %d n1b %d n1c %d n1d %d n2a %d n2b %d n2c %d n2d %d\n",
				n1a,n1b,n1c,n1d,n2a,n2b,n2c,n2d); /* */
  n1pri = 0;
  if (n1a!=NULLVAL) n1pri++;
  if (n1b!=NULLVAL) n1pri++;
  if (n1c!=NULLVAL) n1pri++;
  if (n1d!=NULLVAL) n1pri++;

  n2pri = 0;
  if (n2a!=NULLVAL) n2pri++;
  if (n2b!=NULLVAL) n2pri++;
  if (n2c!=NULLVAL) n2pri++;
  if (n2d!=NULLVAL) n2pri++;

  compat = 1;				/* check if nodes are same */
  if      (n1a != NULLVAL && n2a != NULLVAL && n1a != n2a) compat = 0;
  else if (n1b != NULLVAL && n2b != NULLVAL && n1b != n2b) compat = 0;
  else if (n1c != NULLVAL && n2c != NULLVAL && n1c != n2c) compat = 0;
  else if (n1d != NULLVAL && n2d != NULLVAL && n1d != n2d) compat = 0;

  qn1a = n1a;
  qn1b = n1b;
  qn1c = n1c;
  qn1d = n1d;
  qn2a = n2a;
  qn2b = n2b;
  qn2c = n2c;
  qn2d = n2d;
  if (compat) 				/* if nodes are partially same */
    if (n1pri < n2pri) {		/* reverse order */
      qn1a = n2a;
      qn1b = n2b;
      qn1c = n2c;
      qn1d = n2d;
      qn2a = n1a;
      qn2b = n1b;
      qn2c = n1c;
      qn2d = n1d;
    }

   exceptf = (na!=NULLVAL || nb!=NULLVAL || nc!=NULLVAL || nc!=NULLVAL);

   for (epnt=elempnt; epnt; epnt=epnt->next) {	/* search all elems */

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2))
	ncfprintf (stderr,"ncdispc long loop %d ct %d cn %d n %d ct2 %d cn2 %d\n",
			epnt->elnum, epnt->node1a, epnt->node1b, epnt->node1c, epnt->node2a, epnt->node2b);
#endif

     if (elemtype && elemtype!=ELEMENT && epnt->ctype!=elemtype) continue;
     n1_used = 0;
     match = 1;
     if (qn1a != NULLVAL) 			/* check node 1 first */	
        if (epnt->node1a!=qn1a) match=0;
     if (qn1b != NULLVAL)
        if (epnt->node1b!=qn1b) match=0;
     if (qn1c != NULLVAL)
        if (epnt->node1c!=qn1c) match=0;
     if (qn1d != NULLVAL)
        if (epnt->node1d!=qn1d) match=0;

     if (match) n1_used = 1;
     else {					/* else check node 2 */
	 match = 1;
         if (qn1a != NULLVAL)
            if (epnt->node2a!=qn1a) match=0;
         if (qn1b != NULLVAL)
            if (epnt->node2b!=qn1b) match=0;
         if (qn1c != NULLVAL)
            if (epnt->node2c!=qn1c) match=0;
         if (qn1d != NULLVAL)
            if (epnt->node2d!=qn1d) match=0;
     }

     if (match) {				/* if query node 1 matches */
      if (!n1_used) {				/*   if node 1 not used yet */
        if (qn2a != NULLVAL) 			/* does qn2 match node 1 ? */	
           if (epnt->node1a!=qn2a) match=0;
        if (qn2b != NULLVAL)
           if (epnt->node1b!=qn2b) match=0;
        if (qn2c != NULLVAL)
           if (epnt->node1c!=qn2c) match=0;
        if (qn2d != NULLVAL)
           if (epnt->node1d!=qn2d) match=0;
      }
      else {				     /* else if node 2 not used yet */
        if (qn2a != NULLVAL)	     		/* does qn2 match node 2 ? */	
           if (epnt->node2a!=qn2a) match=0;
        if (qn2b != NULLVAL)
           if (epnt->node2b!=qn2b) match=0;
        if (qn2c != NULLVAL)
           if (epnt->node2c!=qn2c) match=0;
        if (qn2d != NULLVAL)
           if (epnt->node2d!=qn2d) match=0;
      }   /* else (n1_used) */
     }   /* if (match) */

     if (!match) continue;
     nomatch = 0;
     if (exceptf) {
       nomatch = 1;
       if (na!=NULLVAL && epnt->node1a!=na) nomatch=0;     
       if (nb!=NULLVAL && epnt->node1b!=nb) nomatch=0;     
       if (nc!=NULLVAL && epnt->node1c!=nc) nomatch=0;     
       if (nd!=NULLVAL && epnt->node1d!=nd) nomatch=0;     

       if (!nomatch) {
       nomatch = 1;
         if (na!=NULLVAL && epnt->node2a!=na) nomatch=0;     
         if (nb!=NULLVAL && epnt->node2b!=nb) nomatch=0;     
         if (nc!=NULLVAL && epnt->node2c!=nc) nomatch=0;     
         if (nd!=NULLVAL && epnt->node2d!=nd) nomatch=0;     
	}
       if (nomatch) 
           if (exceptype && epnt->ctype!=exceptype) nomatch=0;
     }
     else if (exceptype && epnt->ctype==exceptype) nomatch=1;

     if (checkelemz(epnt,zrange1,zrange2)) nomatch = 1;
     if (nomatch) continue;
     ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);
   }

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdispc end\n");
#endif
}

/*------------------------------------*/

void ncdispn(int n1a, int n1b, int n1c, int n1d, 
	int n2a, int n2b, int n2c, int n2d, 
	int elemtype, int exceptype, 
	int na, int nb, int nc, int nd, 
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
	double vmax, double vmin, double dscale, int hide, int excl, int cmap)

/* Draw neural circuit elements that connect to
   the range of nodes from node1 to node2.
   If any of the node numbers == NULLVAL,
   then ignore this dimension. 
   If "excl"=1, display only elements that connect 
   exclusively to given node.
*/

{
   int swap;
   int qn1a, qn1b, qn1c, qn1d, qn2a, qn2b, qn2c, qn2d;	/* nodes from query */
   int exceptf;
   elem *epnt;
   int match,nomatch;

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdispn start...\n");
#endif

  if      (n1a > n2a) swap = 1;
  else if (n1b > n2b) swap = 1;
  else if (n1c > n2c) swap = 1;
  else if (n1d > n2d) swap = 1;
  else                swap = 0;

  qn1a = n1a;
  qn1b = n1b;
  qn1c = n1c;
  qn1d = n1d;
  qn2a = n2a;
  qn2b = n2b;
  qn2c = n2c;
  qn2d = n2d;
  if (swap) {			/* make qn1 = min, qn2 = max */
      qn1a = n2a;
      qn1b = n2b;
      qn1c = n2c;
      qn1d = n2d;
      qn2a = n1a;
      qn2b = n1b;
      qn2c = n1c;
      qn2d = n1d;
  }

   exceptf = (na!=NULLVAL || nb!=NULLVAL || nc!=NULLVAL || nc!=NULLVAL);

   for (epnt=elempnt; epnt; epnt=epnt->next) {	/* search all elems */

     if (elemtype && elemtype!=ELEMENT && epnt->ctype!=elemtype) continue;
     match = 1;
     if (qn1a != NULLVAL) {
        if (epnt->node1a<qn1a) match=0;
	if (qn2a != NULLVAL)
           if (epnt->node1a>qn2a) match=0;
     }
     if (qn1b != NULLVAL) {
        if (epnt->node1b<qn1b) match=0;
	if (qn2b != NULLVAL)
           if (epnt->node1b>qn2b) match=0;
     }
     if (qn1c != NULLVAL) {
        if (epnt->node1c<qn1c) match=0;
	if (qn2c != NULLVAL)
           if (epnt->node1c>qn2c) match=0;
     }
     if (qn1d != NULLVAL) {
        if (epnt->node1d<qn1d) match=0;
	if (qn2d != NULLVAL)
           if (epnt->node1d>qn2d) match=0;
     }

     if ((!excl && !match) || (excl && (epnt->node2a!=NULLVAL
				    ||  epnt->node2b!=NULLVAL
				    ||  epnt->node2c!=NULLVAL
				    ||  epnt->node2d!=NULLVAL))) {
       if (!excl) match = 1;
       if (qn1a != NULLVAL) {
          if (epnt->node2a<qn1a) match=0;
          if (qn2a != NULLVAL) 
             if (epnt->node2a>qn2a) match=0;
       }
       if (qn1b != NULLVAL) {
          if (epnt->node2b<qn1b) match=0;
          if (qn2b != NULLVAL)
              if (epnt->node2b>qn2b) match=0;
       }
       if (qn1c != NULLVAL) {
          if (epnt->node2c<qn1c) match=0;
          if (qn2c != NULLVAL)
             if (epnt->node2c>qn2c) match=0;
       }
       if (qn1d != NULLVAL) {
          if (epnt->node2d<qn1d) match=0;
          if (qn2d != NULLVAL)
             if (epnt->node2d>qn2d) match=0;
       }
     }		/* if (!match) */

     if (!match) continue;
     nomatch = 0;
     if (exceptf) {
       nomatch = 1;
       if (na!=NULLVAL && epnt->node1a!=na) nomatch=0;     
       if (nb!=NULLVAL && epnt->node1b!=nb) nomatch=0;     
       if (nc!=NULLVAL && epnt->node1c!=nc) nomatch=0;     
       if (nd!=NULLVAL && epnt->node1d!=nd) nomatch=0;     

       if (!nomatch) {
       nomatch = 1;
         if (na!=NULLVAL && epnt->node2a!=na) nomatch=0;     
         if (nb!=NULLVAL && epnt->node2b!=nb) nomatch=0;     
         if (nc!=NULLVAL && epnt->node2c!=nc) nomatch=0;     
         if (nd!=NULLVAL && epnt->node2d!=nd) nomatch=0;     
	}
       if (nomatch) 
           if (exceptype && epnt->ctype!=exceptype) nomatch=0;
     }
     else if (exceptype && epnt->ctype==exceptype) nomatch=1;

     if (checkelemz(epnt,zrange1,zrange2)) nomatch = 1;
     if (nomatch) continue;
     ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);
   }

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr," end\n");
#endif
}

/*------------------------------------*/

void ncdispe (int elemnum, double zrange1, double zrange2,
			int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
			double vmax, double vmin, double dscale, int hide, int cmap)

/* Draw a neural element, given its element number.
*/

{
   elem *epnt;

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdispe start...\n");
#endif

//   for (epnt=elempnt; epnt; epnt=epnt->next) {	/* search all elems */
//      if (epnt==NULL) {
//         ncfprintf (stderr,"ncdispe: can't find element %d\n",elemnum);
//         return;  
//      }
//      if (epnt->elnum == elemnum) {
//         if (checkelemz(epnt,zrange1,zrange2)) continue;
//	 ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);
//	 break;
//      }
//   }

    if ((epnt=findelem(elemnum))==NULL) {
        ncfprintf (stderr,"ncdispe: can't find element %d\n",elemnum);
	return;
    }
    if (checkelemz(epnt,zrange1,zrange2)) return;
    ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr," end\n");
#endif
}

/*------------------------------------*/

void ncdispe (elem *epnt, int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
			double vmax, double vmin, double dscale, int hide, int cmap)
{
#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdispe start...\n");
#endif

    ncdraw(epnt,color,vpen,vpenn,vmax,vmin,dscale,hide,cmap);

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr," end\n");
#endif
}

/*------------------------------------*/

void ncdrnod (int n1a, int n1b, int n1c, int n1d, 
	int elemtype, int exceptype, int na, int nb, int nc, int nd,
	double zrange1, double zrange2,
	int color, Symbol *vpen, double(*vpenn)(int elnum, int color),
	double dscale, int drlabel)

/* Draw node number of node1. If any of the node
   numbers == NULLVAL, then ignore this dimension.
   If "excl"=1, display only elements that connect 
   exclusively to given node.
 */

{
   node *npnt;
   int match,nomatch;

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
	 ncfprintf (stderr,"ncdrnod n1a %d n1b %d n1c %d n1d %d exceptype %d\n",
		n1a,n1b,n1c,n1d,exceptype);
#endif

/* ncfprintf (stderr,"ncdrnod n1a %d n1b %d n1c %d n1d %d\n", 
	n1a,n1b,n1c,n1d);/* */

/* -  -  -  -  -  -  -  -  -  -  -  -  - */

 if (n1a >= 0 && 
     n1b >= 0 && 
     n1c <= 0 && 
     n1d < 0 && 
     exceptype==0) {	/* Run fast loop without n1c or n1d: */
			/*  we're usually displaying all a cell's nodes */

  for (npnt=getnpnt(n1a,n1b); npnt=foreachn(npnt,n1a,n1b); npnt=npnt->hcnext) {

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2))
		ncfprintf (stderr,"ncdrnod fast loop %d %d %d\n",npnt->nodenm1,npnt->nodenm2,npnt->nodenm3);
#endif
    if ((n1c == 0) && (npnt->nodenm3 != n1c)) continue;			// if displaying e.g. soma only
    if (checknodez(npnt,zrange1,zrange2)) return;
    if (drlabel) {
	    if (npnt->label==NULLVAL) return;		// skip if not label 
            else color = npnt->label;
    }
    if (disp_ray) raynode(npnt,dscale,color,ncmatc);
    else          drnode (npnt,dscale,color,ncmat);
  }
 } else {

/* -  -  -  -  -  -  -  -  -  -  -  -  - */

 for (npnt=nodepnt; npnt; npnt=npnt->next) {    /* search all nodes */

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2))
		ncfprintf (stderr,"ncdrnod long loop %d %d\n",npnt->nodenm1,npnt->nodenm2);
#endif
     match = 1;
     if (n1a != NULLVAL) 			/* check node 1 first */	
        if (npnt->nodenm1!=n1a) match=0;
     if (n1b != NULLVAL)
        if (npnt->nodenm2!=n1b) match=0;
     if (n1c != NULLVAL)
        if (npnt->nodenm3!=n1c) match=0;
     if (n1d != NULLVAL)
        if (npnt->nodenm4!=n1d) match=0;

     if (!match) continue;
     nomatch = 0;
     if (na!=NULLVAL || nb!=NULLVAL || nc!=NULLVAL || nc!=NULLVAL) {
       nomatch = 1;
       if (na!=NULLVAL && npnt->nodenm1!=na) nomatch=0;     
       if (nb!=NULLVAL && npnt->nodenm2!=nb) nomatch=0;     
       if (nc!=NULLVAL && npnt->nodenm3!=nc) nomatch=0;     
       if (nd!=NULLVAL && npnt->nodenm4!=nd) nomatch=0;     

       if (nomatch) 
           if (exceptype && npnt->ctype!=exceptype) nomatch=0;
      }
      else if (exceptype && npnt->ctype==exceptype) nomatch=1;

     if (checknodez(npnt,zrange1,zrange2)) nomatch = 1;
     if (nomatch) continue;
     if (drlabel) {
	    if (npnt->label==NULLVAL) continue;		// skip if not label 
            else color = npnt->label;
     }
     if (disp_ray) raynode(npnt,dscale,color,ncmatc);
     else          drnode (npnt,dscale,color,ncmat);
   }
 }
#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdrnod end\n");
#endif
}


/*------------------------------------*/

int findrcp(int num1, int num2, int num3, int num4, photrec **rpnt, const char *str)

/* find PHOTREC (temporary photoreceptor) for "dispstim()" below.  */

{
   photrec *pnt;

  for (pnt=recpnt; pnt; pnt=(photrec *)pnt->next) {
    if ((pnt->ctype==PHOTREC) &&
	(pnt->recnm1 == num1) && (pnt->recnm2 == num2) &&
        (pnt->recnm3==num3) && (pnt->recnm4==num4)) {
        *rpnt = pnt;
        return 1;
    }
  }
  if (str) {
     if (num4 != NULLVAL)
         ncfprintf(stderr,"\n%s: can't find photrec %d %d %d %d\n",
                        str,num1,num2,num3,num4);
     else if (num3 != NULLVAL)
         ncfprintf(stderr,"\n%s: can't find photrec %d %d %d\n",str,num1,num2,num3);
     else if (num2 != NULLVAL)
         ncfprintf(stderr,"\n%s: can't find photrec %d %d\n",str,num1,num2);
     else
         ncfprintf(stderr,"\n%s: can't find photrec %d\n",str,num1);
   }
  return 0;
}

/*------------------------------------*/

/* colors defined in "colors.h", implemented in "pl/mprint?.c. */

#define NUMCOLS1 32

static int color_lut1[NUMCOLS1]=
		{BLACK,BLUE,GREEN,CYAN,RED,MAGENTA,BROWN,WHITE,
		 GRAY,LTBLUE,LTGREEN,LTCYAN,LTRED,LTMAG,YELLOW,BRTWHT,
		 VCOL1,VCOL2,VCOL3,VCOL4,VCOL5,VCOL6,VCOL7,VCOL8,
		 VCOL9,VCOL10,VCOL11,VCOL12,VCOL13,VCOL14,VCOL15,VCOL16};

#define NUMCOLS2 16
#define OFFSET_COLMAP2 16
static int *color_lut2 = &color_lut1[OFFSET_COLMAP2];

#define NUMCOLS3 32
static int color_lut3[NUMCOLS3]=
		{NCOL1,NCOL2,NCOL3,NCOL4,NCOL5,NCOL6,NCOL7,NCOL8,
                 NCOL9,NCOL10,NCOL11,NCOL12,NCOL13,NCOL14,NCOL15,NCOL16,
		 NCOL17,NCOL18,NCOL19,NCOL20,NCOL21,NCOL22,NCOL23,NCOL24,
                 NCOL25,NCOL26,NCOL27,NCOL28,NCOL29,NCOL30,NCOL31,NCOL32};

#define NUMCOLS4 16
static int color_lut4[NUMCOLS4]=
		{PCOL1,PCOL2,PCOL3,PCOL4,PCOL5,PCOL6,PCOL7,PCOL8,
                 PCOL9,PCOL10,PCOL11,PCOL12,PCOL13,PCOL14,PCOL15,PCOL16};

#define NUMCOLS5 16
static int color_lut5[NUMCOLS5]=
		{RCOL1,RCOL2,RCOL3,RCOL4,RCOL5,RCOL6,RCOL7,RCOL8,
                 RCOL9,RCOL10,RCOL11,RCOL12,RCOL13,RCOL14,RCOL15,PCOL16};

#define NUMCOLS6 16
static int color_lut6[NUMCOLS6]=
		{DCOL1,DCOL2,DCOL3,DCOL4,DCOL5,DCOL6,DCOL7,DCOL8,
                 DCOL9,DCOL10,DCOL11,DCOL12,DCOL13,DCOL14,DCOL15,DCOL16};

#define NUMCOLS7 16
static int color_lut7[NUMCOLS7]=
		{BCOL1,BCOL2,BCOL3,BCOL4,BCOL5,BCOL6,BCOL7,BCOL8,
                 BCOL9,BCOL10,BCOL11,BCOL12,BCOL13,BCOL14,BCOL15,BCOL16};


/* See ccols.txt for the variety of colors: (from /usr/lib/X11/rgb.txt) */

#define NUMCOLS8 100
static int color_lut8[NUMCOLS8]=

{GCOL0, GCOL1, GCOL2, GCOL3, GCOL4, GCOL5, GCOL6, GCOL7, GCOL8, GCOL9, GCOL10,
GCOL11, GCOL12, GCOL13, GCOL14, GCOL15, GCOL16, GCOL17, GCOL18, GCOL19, 
GCOL20, GCOL21, GCOL22, GCOL23, GCOL24, GCOL25, GCOL26, GCOL27, GCOL28,
GCOL29, GCOL30, GCOL31, GCOL32, GCOL33, GCOL34, GCOL35, GCOL36, GCOL37,
GCOL38, GCOL39, GCOL40, GCOL41, GCOL42, GCOL43, GCOL44, GCOL45, GCOL46,
GCOL47, GCOL48, GCOL49, GCOL50, GCOL51, GCOL52, GCOL53, GCOL54, GCOL55,
GCOL56, GCOL57, GCOL58, GCOL59, GCOL60, GCOL61, GCOL62, GCOL63, GCOL64,
GCOL65, GCOL66, GCOL67, GCOL68, GCOL69, GCOL70, GCOL71, GCOL72, GCOL73,
GCOL74, GCOL75, GCOL76, GCOL77, GCOL78, GCOL79, GCOL80, GCOL81, GCOL82,
GCOL83, GCOL84, GCOL85, GCOL86, GCOL87, GCOL88, GCOL89, GCOL90, GCOL91,
GCOL92, GCOL93, GCOL94, GCOL95, GCOL96, GCOL97, GCOL98, GCOL99};

#define NUMCOLS9 100
static int color_lut9[NUMCOLS9]=

{CCOL0, CCOL1, CCOL2, CCOL3, CCOL4, CCOL5, CCOL6, CCOL7, CCOL8, CCOL9, CCOL10,
CCOL11, CCOL12, CCOL13, CCOL14, CCOL15, CCOL16, CCOL17, CCOL18, CCOL19, 
CCOL20, CCOL21, CCOL22, CCOL23, CCOL24, CCOL25, CCOL26, CCOL27, CCOL28,
CCOL29, CCOL30, CCOL31, CCOL32, CCOL33, CCOL34, CCOL35, CCOL36, CCOL37,
CCOL38, CCOL39, CCOL40, CCOL41, CCOL42, CCOL43, CCOL44, CCOL45, CCOL46,
CCOL47, CCOL48, CCOL49, CCOL50, CCOL51, CCOL52, CCOL53, CCOL54, CCOL55,
CCOL56, CCOL57, CCOL58, CCOL59, CCOL60, CCOL61, CCOL62, CCOL63, CCOL64,
CCOL65, CCOL66, CCOL67, CCOL68, CCOL69, CCOL70, CCOL71, CCOL72, CCOL73,
CCOL74, CCOL75, CCOL76, CCOL77, CCOL78, CCOL79, CCOL80, CCOL81, CCOL82,
CCOL83, CCOL84, CCOL85, CCOL86, CCOL87, CCOL88, CCOL89, CCOL90, CCOL91,
CCOL92, CCOL93, CCOL94, CCOL95, CCOL96, CCOL97, CCOL98, CCOL99};

struct color_map {
	int ncolors;
	int *lut;
	};

#define MAXNCMAP 20
#define    NCMAP 10

int ncmap = NCMAP;
int *new_cmap=NULL;

static color_map
colormaps[MAXNCMAP] =	   {NUMCOLS1, color_lut1,	/* for plots        */
			    NUMCOLS2, color_lut2,	/* for vcolor, etc. */
			    NUMCOLS3, color_lut3,	/* for vcolor, psps, spikes */
			    NUMCOLS4, color_lut4,	/* vcolor blue-red */
			    NUMCOLS5, color_lut5,	/* vcolor black-red */
			    NUMCOLS6, color_lut6,	/* vcolor black-green  */
			    NUMCOLS7, color_lut7,	/* vcolor black-blue  */
			    NUMCOLS8, color_lut8,	/* gray dis pstim   */
			    NUMCOLS9, color_lut9,	/* variety of colors */
			    	   0, 0			/* for user-defined map */
		          };

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setcmap(int *arrp, int cmapsize) 

/* set a user-defined color map */

{
     int i;

   if (ncmap >= MAXNCMAP-1) {
       ncfprintf (stderr,"%s: too many colormaps %d\n",progname,ncmap+1);
       return 0;
   }
   new_cmap = (int*)emalloc(cmapsize*sizeof(int));
   for (i=0; i<cmapsize; i++) {
     new_cmap[i] = arrp[i];
   }
   colormaps[ncmap].ncolors=cmapsize;
   colormaps[ncmap].lut=new_cmap;
   ncmap++;
   return ncmap;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setcmap(datum d1)

/* set a user-defined color map */

{
     int i,cmapsize;

   if (ncmap >= MAXNCMAP-1) {
       ncfprintf (stderr,"%s: too many colormaps %d\n",progname,ncmap+1);
       return 0;
   }
   cmapsize = d1.arrp->dim[0];
   new_cmap = (int*)emalloc(cmapsize*sizeof(int));
   for (i=0; i<cmapsize; i++) {
     new_cmap[i] = int(d1.arrp->arr[i]);
   }
   colormaps[ncmap].ncolors=cmapsize;
   colormaps[ncmap].lut=new_cmap;
   ncmap++;
   return ncmap;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void dispstim(double start_time, double stop_time, double dia, int cmap, double max_flux, double min_flux)

/* display stimulus from start_time through stop_time */

/* Set max_flux to -LARGENODE for automatic scaling */

{
//#define RTIMERND   0.00001		/* down-round value */
#define RTIMERND   1e-10		/* down-round value */
#define DSTIME     0.0001		/* time incr for reading in stims */

   recstim *rspnt, *next;
   photrec *rpnt,*rnext;
   elem *epnt;
   static double start,stop, mask;
   double flux, flux_range, flux_val, norm_val;
   double flux_range_max, flux_range_min, flux_dev, flux_avg;
   double x,y,tdia,width,dist;
   int i,n,color_val,col_indx,*color_lut,ncolors,fill=1;
   char nbuf[30];
   static double oldmax=0, oldmin=0, oldtime=0;
   static int made_recs = 0;
   static int old_made_recs = 0;

#ifdef DEBUG
  if (debug & NCDISP && debugz & 1)  ncfprintf (stderr,"dispstim interval %g %g\n",start_time,stop_time);
#endif

 if (makestim) return;

 if (!runyet) { 

   /* If the "step" or "run" statements haven't been run yet, */
   /*  no photoreceptors have been created so before we can   */
   /*  read out their light fluxes, we need to make dummy     */
   /*  photoreceptors and stimulate with all the stimulus     */
   /*  events.  */

  if (!made_recs) {
    for (epnt=elempnt; epnt; epnt=epnt->next) {  	  /* check all elements */
      if (epnt->ctype == ROD || epnt->ctype == CONE || epnt->ctype == CHR ||
  	epnt->ctype == VTRANSDUCER || epnt->ctype == ITRANSDUCER) {
       makr((photorec*)epnt, PHOTREC);	  /* make a dummy photoreceptor */
      }
    }
    made_recs = 1;
  }
	/* now read in the stimulus events in increments of DSTIME: */

  for (start=start_time-RTIMERND; start<stop_time; start+=DSTIME) {
    stop = start + DSTIME;
    if (stop > stop_time) continue;

#ifdef DEBUG
  if (debug & NCDISP && debugz & 1)  ncfprintf (stderr,"dispstim checking time %g %g %g\n",start, stop, stop_time);
#endif

    if (!readstim(stop)) {  /* read new photrec inputs */
         ncfprintf (stderr,"dispstim: invalid stimulus at time %g\n",stop_time);
	 execerror ("Missing stimulus: "," stopping... ");
    }   

	/* check new stimulus events and possibly use them */

   for (rspnt=recspnt; rspnt; ) {

/*  ncfprintf (stderr,"stim num %d type '%c' time %g rtime %g start %g stop %g\n",
              rspnt->recnm1, rspnt->ctype,rspnt->time, stop_time,start,stop);  /* */

   if ((rspnt->time >= start) && (rspnt->time <stop)) {
      switch (rspnt->ctype) {

       case 'a':                /* absolute additive */
          if (findrcp (rspnt->recnm1, rspnt->recnm2,
			 rspnt->recnm3,rspnt->recnm4,
					&rpnt, "dispstim")) {
             rpnt->aflux[rspnt->stimchan] = rspnt->val;
	  }
#ifdef DEBUG
          if (debug & NCDISP && debugz & 1)  ncfprintf (stderr,"abs stim %d %g %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel);
#endif
         break;

     case 'b':                /* delta masking */
          if (findrcp (rspnt->recnm1, rspnt->recnm2,
			 rspnt->recnm3, rspnt->recnm4,
                                        &rpnt, "dispstim")) {
            rpnt->mflux[rspnt->stimchan] += rspnt->val;
            rpnt->mask[rspnt->stimchan] = rspnt->mask;
          }
#ifdef DEBUG
          if (debug & NCDISP && debugz & 1)  ncfprintf (stderr,"del mask stim %d %g %g %g %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel,rpnt->mflux,rpnt->mask);
#endif
         break;
       case 'c':                /* absolute masking */
          if (findrcp (rspnt->recnm1, rspnt->recnm2,
			 rspnt->recnm3,rspnt->recnm4,
					&rpnt, "dispstim")) {
             rpnt->mflux[rspnt->stimchan] = rspnt->val;
             rpnt->mask[rspnt->stimchan] = rspnt->mask;
	  }
#ifdef DEBUG
          if (debug & NCDISP && debugz & 1)  ncfprintf (stderr,"abs mask stim %d %g %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel);
#endif
         break;

     case 'd':                /* delta additive */
          if (findrcp (rspnt->recnm1, rspnt->recnm2,
			 rspnt->recnm3, rspnt->recnm4,
                                        &rpnt, "dispstim")) {
            rpnt->aflux[rspnt->stimchan] += rspnt->val;
          }
#ifdef DEBUG
          if (debug & NCDISP && debugz & 1)  ncfprintf (stderr,"del stim %d %d %g %g %g %12.12g\n",
                                        rspnt->recnm1,rspnt->recnm2,rspnt->val,rspnt->wavel,rpnt->aflux,rspnt->time);
#endif
         break;
      } /* switch */	

		/* now delete the stimulus event if we've used it: */

      next = rspnt->next;            /* delete and patch pointers */
      delrstim(rspnt);
      rspnt = next;
   }        /* if */

   else {   /* no stimulus here */
     rspnt = rspnt->next;
   }
  }  /* for (rspnt=recspnt;;) */
  
 }  /* for (start=;;) */

} /* if (!runyet) */


 if (max_flux == -LARGENODE) { /* If max_flux is not specified, */
			  /* find max and min stimulus vals */
   for (n=0,flux_avg=0,max_flux=-LARGENODE,min_flux=LARGENODE,rpnt=recpnt; 
			rpnt; rpnt=(photrec*)rpnt->next){
	flux = 0;
	for (i=0; i<numstimchan; i++) {		/* add up all disp channels */
           mask = rpnt->mask[i];
	   if (mask>1) mask = 1;
	   if (mask<0) mask = 0;
   	   if (rpnt->mflux[i] <= unmaskthr) mask = 0;       /* intensity masking */
   	   flux += rpnt->aflux[i]*(1-mask) + rpnt->mflux[i]*mask; /* masking */
	}
   	flux_val = flux;
	if (flux_val > max_flux)
		max_flux = flux_val;
	if (flux_val < min_flux)
		min_flux = flux_val;
	flux_avg += flux_val;
	n++;
    }
    if (n==0) n=1;
    flux_avg /= n; 
  } else {	/* max_flux is specified */
	flux_avg = (max_flux + min_flux) * 0.5;
  }
  flux_range_max = max_flux - flux_avg;
  flux_range_min = flux_avg - min_flux;
  if (flux_range_max > flux_range_min) flux_dev = flux_range_max;
  else                                 flux_dev = flux_range_min;
  if (flux_range == 0)  flux_range = 1;

  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  if (cmap==0) cmap=7;
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next){
	flux = 0;
	for (i=0; i<numstimchan; i++) {		/* add up all disp channels */
           mask = rpnt->mask[i];
	   if (mask>1) mask = 1;
	   if (mask<0) mask = 0;
   	   if (rpnt->mflux[i] <= unmaskthr) mask = 0;       /* intensity masking */
   	   flux += rpnt->aflux[i]*(1-mask) + rpnt->mflux[i]*mask; /* masking */
	   // if (mask>0) fprintf (stderr,"mask %g mflux %g\n", mask,rpnt->mflux[i]);
	}
	if (flux_dev <= 0) flux_dev = 1;
	norm_val = (flux - flux_avg) / (2*flux_dev) + 0.5;
	if (norm_val < 0) norm_val = 0;
	else if (norm_val >= 1) norm_val = 0.999;
	col_indx = (int)(norm_val * ncolors);
	color_val = color_lut[col_indx];
        if (! loccheck (rpnt->xloc,rpnt->yloc,0.0)) {/* check if valid loc */
	  execerror ("nc: dispstim: invalid photoreceptor location","stopping.");
        }
 	drsphere (rpnt->xloc,rpnt->yloc,0.0,dia,1.0,color_val,0,ncmat,fill); 
  }
		/* now delete the dummy photoreceptors we made above */

//   Old version deleted dummy photorecs every call,
//    now we leave them, with memories of stimuli to save time

  if (start_time==simtime && old_made_recs) {
    for (rpnt=recpnt; rpnt; ) {
	rnext = (photrec*)rpnt->next;
	if (rpnt->ctype==PHOTREC) {	/* delete dummy photorec */
		delphotrec(rpnt);
	}	
	rpnt = rnext;
    }
  }
  old_made_recs = made_recs;

#define COLBARLEN 0.40
#define XBAR  0.30
#define YBAR  0.94

					/* generate color bar */
  tdia = 0.02;
  dist = COLBARLEN/ncolors;
  x = XBAR; y = YBAR;

   gframe ("Col_bar");
   gorigin (x,y);
   for (width= -tdia/2.0; width<tdia/2.0; width+= linewidth) {
     gmove (0, width);
     for (i=0; i<ncolors; i++){
	gpen(color_lut[i]);
	gmove (i*dist+.001, width);
        gdraw ((i+1)*dist, width);
     }
  }
  gpen(0);
  gmove(-0.1,-0.03);
  sprintf (nbuf,"%.3g",oldmin);
  gtext(nbuf);

  gmove(COLBARLEN+0.05,-0.03);
  sprintf (nbuf,"%.3g",oldmax);
  gtext(nbuf);

  gmove(COLBARLEN+0.20,-0.0);
  sprintf (nbuf,"%.4f",oldtime);
  gtext(nbuf);

  gmove(-0.1,-0.002);
  //gpen(color_lut[0]);
  gpen(15);
  gtext("MIN");
  gmove(-0.1,-0.03);
  sprintf (nbuf,"%.3g",min_flux);
  gtext(nbuf);
  gmove(COLBARLEN+0.05,-0.002);
  //gpen(color_lut[ncolors-1]);
  gpen(15);
  gtext("MAX");
  gmove(COLBARLEN+0.05,-0.03);
  sprintf (nbuf,"%.3g",max_flux);
  gtext(nbuf);
  gmove(COLBARLEN+0.20,-0.0);
  sprintf (nbuf,"%.4f",stop_time);
  gtext(nbuf);

  ppurge();
  oldmax = max_flux;
  oldmin = min_flux;
  oldtime = stop_time;

  gframe ("..");
  grmframe ("Col_bar");
}

/*------------------------------------*/

void dispstim(double stop_time, double dia, int cmap, double max_flux, double min_flux)

{
  dispstim(simtime, stop_time, dia, cmap, max_flux, min_flux);
}

/*------------------------------------*/

void dispstim(double start_time, double stop_time, double dia, int cmap)

{
  dispstim(start_time, stop_time, dia, cmap, -LARGENODE, LARGENODE);
}

/*------------------------------------*/

void dispstim(double stop_time, double dia, double max_flux, double min_flux )

{ 
     int cmap;
  dispstim(simtime, stop_time, dia, cmap=1, max_flux, min_flux);
}

/*------------------------------------*/

void dispstim(double stop_time, double dia, int cmap)

{
  dispstim(simtime, stop_time, dia, cmap, -LARGENODE, LARGENODE);
}

/*------------------------------------*/

int dispcolor (elem *epnt, int dcolor, int cmap)

/* find the color of an element from the color map. */

{
    int color;
    int *color_lut, ncolors,ncol;

  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  //if (cmap==0) cmap = 1;
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  ncol = ncolors - 1;
  if (dcolor==NULLVAL) {
    switch (epnt->ctype) {
      case SPHERE: dcolor = BLUE; break; 
      case CABLE:  dcolor = GREEN; break; 
      case SYNAPSE: dcolor = MAGENTA; break; 
      case GJ: dcolor = CYAN; break; 
      default: dcolor = YELLOW; break;
    }
  }
  color = limit (dcolor,ncol,0);
  return color_lut[color];
}

/*------------------------------------*/

int dispvcolor (elem *epnt, int dcolor, double vmax, double vmin, int cmap)

/* find the color of an element from its voltage. */

{
    double range,v,t;
    node *npnt;
    comp *cpnt;
    int color,colorn,vfound;
    int *color_lut, ncolors, ncol;

  vfound = 0;
  if (vmax < vmin) { t = vmax; vmax = vmin; vmin = t; }
  if (vmax==NULLVAL) vmax = plmax;
  if (vmin==NULLVAL) vmin = plmin;
  range = vmax - vmin;
  range = abs(range);
  if (npnt=epnt->nodp1)
    if (cpnt=npnt->comptr) {
       v = cpnt->v;
       vfound = 1;
     }
  if (!vfound) v = vmin;

  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  color = (int) ((v-vmin) / range * ncolors);
  ncol = ncolors - 1;
  colorn = limit (color,ncol,0);
  return color_lut[colorn];
}

/*------------------------------------*/

int disprcolor (elem *epnt, int dcolor, double rmax, double rmin, int cmap)

/* find the color of an element from its cell region. */

{
    int color,region=0;
    int *color_lut, ncolors, ncol;

  if (epnt!=NULL) region = epnt->region;
  color = 0;
  if (rcolors!=NULL) {
    color = int(rcolors[region]);
  }
  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  if (cmap==0) cmap = 0;			/* default color = base colors */
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  if (color<=0) {
    color = region+1;
    ncol = ncolors - 1;
    color = limit (color,ncol,0);
  }
  return color_lut[color];
}

/*------------------------------------*/

int displcolor (elem *epnt, int dcolor, double lmax, double lmin, int cmap)

/* find the color of an element from its light flux. */

{
    double range,flux,t;
    photrec *rpnt;
    comp *cpnt;
    int color,lfound;
    int *color_lut, ncolors,ncol;

  lfound = 0;
  if (lmax < lmin) { t = lmax; lmax = lmin; lmin = t;}
  if (lmax==NULLVAL) lmax = 1e6;
  if (lmin==NULLVAL) lmin = 0;
  range = lmax - lmin;
  range = abs(range);
  if (rpnt=(photrec*)epnt->lptr) {
     if (epnt->ctype == ROD 
	|| epnt->ctype == CONE
	|| epnt->ctype == CHR
	|| epnt->ctype==VTRANSDUCER
	|| epnt->ctype==ITRANSDUCER) {

       flux = rpnt->iflux;
   /* ncfprintf (stderr,"flux %g\n",flux); /* */
       lfound = 1;
     }
  }
  if (!lfound) flux = lmin;
  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  if (cmap==0) cmap = 7;			/* default color = gray */
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  if (cmap == 1) {
    color = (int) (flux-lmin);
  } else {
    color = (int) ((flux-lmin) / range * ncolors);
  }
  ncol = ncolors - 1;
  color = limit (color,ncol,0);
  return color_lut[color];
}

/*------------------------------------*/

int dispcacolor (elem *epnt, int dcolor, double camax, double camin, int cmap)

/* find the color of an element from its Ca conc. */

{
    double range,ca,t;
    node *npnt;
    comp *cpnt;
    cacomp *capnt;
    int color,cafound;
    int *color_lut, ncolors,ncol;

  cafound = 0;
  if (camax < camin) { t = camax; camax = camin; camin = t;}
  if (camax==NULLVAL) camax = 5e-6;
  if (camin==NULLVAL) camin = 0;
  range = camax - camin;
  range = abs(range);
  camin = min(camax,camin);
  if (npnt=epnt->nodp1)
    if (cpnt=npnt->comptr) {
       if (capnt=cpnt->capnt) {
         ca = capnt->cais[0];
         cafound = 1;
       }
     }
  if (!cafound) ca = camin;
  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  if (cmap==0) cmap = 3;			/* default color = blue - red */
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  color = (int) ((ca-camin) / range * ncolors);
  ncol = ncolors - 1;
  color = limit (color,ncol,0);
  return color_lut[color];
}

/*------------------------------------*/

int dispccolor (elem *epnt, int dcolor, double fmax, double fmin, int cmap)

/* find the color of a state from a channel connected to an element. */

{
    int n;
    double range,val,t;
    node *npnt;
    comp *pnt;
    conn *cpnt;
    chan *chpnt;
    int color,found;
    int *color_lut, ncolors,ncol;
    conlst *lpnt;

#define NAASTATE  6		/* state 6 is active state */
#define NAISTATE1 7		/* state 7 is inactive state */
#define NAISTATE2 8		/* state 8 is inactive state */
#define NAISTATE3 9		/* state 9 is inactive state */
#define KASTATE   5		/* state 4 is active state */

  found = 0;
  if (fmax < fmin) { t = fmax; fmax = fmin; fmin = t;}
  if (fmax==NULLVAL) fmax = 1;
  if (fmin==NULLVAL) fmin = 0;
  range = fmax - fmin;
  range = abs(range);
  fmin = min(fmax,fmin);
  if (npnt=epnt->nodp1)
    if (pnt=npnt->comptr) {
      for (n=0,lpnt=pnt->clst; lpnt; lpnt=lpnt->next, n++) {  /* find connections */
        if (!(cpnt=lpnt->conpnt)) continue;
	switch (cpnt->ctype) {
          case NA: 
		if (cpnt->stype == 2) {		/* if Na type 2 */
		  switch (dcolor) {
		   case NAACOLOR: 
		     val = rechan (cpnt, NRECG, NAASTATE,0); 
	             found = 1;
		     break;
		   case NAICOLOR: 
		     val = rechan (cpnt,  NRECG, NAISTATE1,0); 
		     val += rechan (cpnt, NRECG, NAISTATE2,0); 
		     val += rechan (cpnt, NRECG, NAISTATE3,0); 
	             found = 1;
		     break;
		  }
		}
	        break;
          case K: 
		if (cpnt->stype == 2 || cpnt->stype==3) {	/* if K type 1,6 */
		  switch (dcolor) {
		   case KACOLOR: 
		     val = rechan (cpnt, NRECG, NAASTATE,0); 
	             found = 1;
		     break;
		  }
		}
	        break;
          case SYNAPSE: 
		  switch (dcolor) {
		   case SRCOLOR: 
		     val = recsynap (((synap*)cpnt), NRECA9,0); 
	             found = 1;
		     break;
		   case SGCOLOR: 
		     val = rechan (cpnt, NRECG0,0,0); 
	             found = 1;
		     break;
		  }
	        break;
	}
      }  /* for lpnt= */
    }
  if (!found) val = fmin;
  cmap--;
  if (cmap<0) cmap = 0;
  if (cmap>=ncmap) cmap = ncmap-1;
  if (cmap==0) cmap = 3;			/* default color = blue - red */
  color_lut = colormaps[cmap].lut;
  ncolors   = colormaps[cmap].ncolors;
  if (color_lut==NULL) {
       color_lut = colormaps[0].lut;
       ncolors   = colormaps[0].ncolors;
  }
  color = (int) ((val-fmin) / range * ncolors);
  ncol = ncolors - 1;
  color = limit (color,ncol,0);
  return color_lut[color];
}

/*------------------------------------*/

void ncdraw (elem *epnt, int ecolor, Symbol *vpen, double(*vpenn)(int elnum,int color),
		double vmax, double vmin,double dscale, int hide, int cmap)

/* draw a neural element on the screen */

{
   node *np1,*np2;
   conlst *lpnt;
   elem *nepnt;
   double n1dia,n2dia;
   double x1,y1,z1,x2,y2,z2,xd,yd;
   int dcolor;

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr,"ncdraw start...\n");
#endif

  if (!epnt) return;

  np1 = epnt->nodp1;
  np2 = epnt->nodp2;

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2)) {
	ncfprintf (stderr,"ncdraw element type: '%s' ",findsym(epnt->ctype));
	ncfprintf (stderr,"dia %g dscale %g color %d\n",((cable *)epnt)->dia,
			dscale,ecolor);
  }
#endif

  if (!np1) {
#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 2))
		ncfprintf (stderr,"ncdraw: returning, nodp1 not set\n");
#endif
     return;
     }

	/* find diameter of any spheres connected to node: */

  if (np1) for (n1dia=0.0,lpnt=np1->elemlst; lpnt; lpnt=lpnt->next) {
      nepnt = ((elem *)lpnt->conpnt);
      if (nepnt->ctype==SPHERE) {
        if (n1dia<((sphere*)nepnt)->dia) n1dia = ((sphere*)nepnt)->dia;
      } 
    }
  if (np2) for (n2dia=0.0,lpnt=np2->elemlst; lpnt; lpnt=lpnt->next) {
      nepnt = ((elem *)lpnt->conpnt);
      if (nepnt->ctype==SPHERE) {
        if (n2dia<((sphere*)nepnt)->dia) n2dia = ((sphere*)nepnt)->dia;
      } 
    }

  if (np1) { 
	x1 = np1->xloc; y1 = np1->yloc; z1 = np1->zloc;
	if (!loccheck(x1,y1,0.0)) { 
	  ncfprintf (stderr,
	  "# Invalid display location (%g %g %g) for node 1 '%s' of %s # %d\n",
		x1,y1,z1,
		prnode(np1->nodenm1,np1->nodenm2,np1->nodenm3,np1->nodenm4),
		findsym(epnt->ctype), epnt->elnum);
        }
	if (!loccheck(0.0,0.0,z1)) { 
		z1 = 0.0;
	}
  }
  if (np2) { 
	x2 = np2->xloc; y2 = np2->yloc; z2 = np2->zloc; 
	if (!loccheck(x2,y2,0.0)) {
	  ncfprintf (stderr,
	  "# Invalid display location (%g %g %g) for node 2 '%s' of %s # %d\n",
		x2,y2,z2,
		prnode(np2->nodenm1,np2->nodenm2,np2->nodenm3,np2->nodenm4),
		findsym(epnt->ctype), epnt->elnum);
	}
	if (!loccheck(0.0,0.0,z2)) { 
		z2 = 0.0;
	}
  }
				/* set color from voltage */

  switch (ecolor) {
	 case   VCOLOR: dcolor = dispvcolor (epnt,ecolor,vmax,vmin,cmap); break;
         case   LCOLOR: dcolor = displcolor (epnt,ecolor,vmax,vmin,cmap); break;
         case  CACOLOR: dcolor = dispcacolor(epnt,ecolor,vmax,vmin,cmap); break;
         case   RCOLOR: dcolor = disprcolor (epnt,ecolor,vmax,vmin,cmap); break;
	 case  SRCOLOR:
	 case  SGCOLOR:
	 case NAICOLOR:
	 case NAACOLOR: dcolor = dispccolor (epnt,ecolor,vmax,vmin,cmap); break;
	       default: dcolor = dispcolor  (epnt,ecolor,          cmap); break;
  }

  if (vpen) {
     dcolor = int(callfunc (vpen,2,(double)epnt->elnum,(double)ecolor,0.0,0.0));
  } 
  if (vpenn) {
     dcolor = int(vpenn(epnt->elnum,ecolor));
  } 

 if (disp_ray) {

  switch (epnt->ctype) {

   case CABLE:	raycable (x1,y1,z1,x2,y2,z2, 
		((cable*)epnt)->dia,((cable*)epnt)->dia2,
		dscale,n1dia,n2dia,dcolor,hide,ncmatc);
		break;

   case SPHERE:	raysphere (x1,y1,z1, ((sphere*)epnt)->dia, dscale,
			dcolor,hide,ncmatc,1);
		break;

   case SYNAPSE: raysynap (x1,y1,z1,x2,y2,z2,(synapse *)epnt, dscale,n1dia,n2dia,
			dcolor,hide,ncmatc);
		 break;

   case RESISTOR:
   case AXIALRES:
   case GJ:	rayconn (x1,y1,z1,x2,y2,z2,epnt->ctype,dscale,
			n1dia,n2dia,dcolor,hide,ncmatc);
		 break;

   case VTRANSDUCER:	 
   case ITRANSDUCER:	 
   case CHR:	 
   case ROD:	 
   case CONE:   rayphotrec (x1,y1,z1,
		   (photorec *)epnt,dscale,n1dia,n2dia,dcolor,hide,ncmatc);
		 break;

   case GNDCAP:
   case LOAD:	 rayload (x1,y1,z1,
			  epnt->ctype,n1dia,n2dia,dcolor,hide,ncmatc);
		break;
   case ELECTRODE:  rayelec (x1,y1,z1,
		   (electrode *)epnt,dscale,dcolor,hide,ncmatc);
		 break;
   }
  }  
 else /* display 2D */

  switch (epnt->ctype) {

   case CABLE:  drcable (x1,y1,z1,x2,y2,z2, 
			((cable*)epnt)->dia, ((cable*)epnt)->dia2,
			dscale,n1dia,n2dia,dcolor,hide,ncmat);
		break;

   case SPHERE:	drsphere (x1,y1,z1, ((sphere*)epnt)->dia, dscale,
			dcolor,hide,ncmat,1);
		break;

   case SYNAPSE: drsynap (x1,y1,z1,x2,y2,z2, (synapse *)epnt, dscale,n1dia,n2dia,
			dcolor,hide,ncmat);
		 break;

   case RESISTOR:
   case AXIALRES:
   case GJ:	drconn (x1,y1,z1,x2,y2,z2,epnt->ctype,dscale,
			n1dia,n2dia,dcolor,hide,ncmat);
		 break;

   case VTRANSDUCER:	 
   case ITRANSDUCER:	 
   case CHR:	 
   case ROD:	 
   case CONE:   drphotrec (x1,y1,z1,
		   (photorec *)epnt,dscale,n1dia,n2dia,dcolor,hide,ncmat);
		 break;

   case GNDCAP:
   case LOAD:	 drload (x1,y1,z1,epnt->ctype,n1dia,n2dia,dcolor,hide,ncmat);
		break;

   case ELECTRODE:   drelec (x1, y1, z1,
		   (electrode *)epnt,dscale,dcolor,hide,ncmat);
		break;
  }

  ppurge();		/* flush the output stream */

#ifdef DEBUG
  if ((debug & NCDISP) && (debugz & 1))
		ncfprintf (stderr," ncdraw end\n");
#endif
}

/*--------------------------------------*/

void ppurge (void)

{
  fflush (stdout);
}

/*--------------------------------------*/

void drcomp(double xloc, double yloc, double zloc, double dia, 
		int color, double vmax, double vmin)

/* draw a compartment */

{
   int dcolor;

 if (color) dcolor = color;
 else       dcolor = BLUE;
 if (disp_ray)
  raysphere (xloc,yloc,zloc,dia,1.0,dcolor,0,ncmatc,0);
 else
  drsphere (xloc,yloc,zloc,dia,1.0,dcolor,0,ncmat,0);
 ppurge();
}

/*--------------------------------------*/

void dr_node (node *npnt, double dscale)

{
    double xoffs,charsiz;
    static char numbuf[20];

  if (dscale<0) {
     dscale = -dscale;
     int d = int(dscale);
     dscale = (dscale - d) * 10;
     if (dscale==0) dscale = 1.0;
     if (npnt->labeltext != NULL) {
         sprintf (numbuf,"%s",npnt->labeltext);
     } else {
       switch (d) {
        case 1: sprintf (numbuf,"%d",npnt->nodenm1); break;
        case 2: sprintf (numbuf,"%d",npnt->nodenm2); break;
        case 3: sprintf (numbuf,"%d",npnt->nodenm3); break;
        case 4: sprintf (numbuf,"%d",npnt->nodenm4); break;
        case 5: sprintf (numbuf,"%d %d",npnt->nodenm1,npnt->nodenm2); break;
        case 6: sprintf (numbuf,"%d %d",npnt->nodenm2,npnt->nodenm3); break;
        case 7: sprintf (numbuf,"%d %d %d",npnt->nodenm1,npnt->nodenm2,npnt->nodenm3);
  	      break;
      }
    }
  }
  else if (npnt->nodenm3 != NULLVAL)
      sprintf (numbuf,"%d %d %d",npnt->nodenm1,npnt->nodenm2,npnt->nodenm3);
  else if (npnt->nodenm2 != NULLVAL)
      sprintf (numbuf,"%d %d",npnt->nodenm1,npnt->nodenm2);
  else
      sprintf (numbuf,"%d",npnt->nodenm1);
  charsiz = .010 * dscale;
  // xoffs = strlen(numbuf) / 2.0;		/* find center of number */
  xoffs = charsiz * .6;
  grmove (-xoffs, -charsiz*.3);	/* label node numbers */
  gcwidth (charsiz);
  gtext (numbuf);
  ppurge();
}

/*--------------------------------------*/

void drnode(node *npnt, double dscale, int color, double (*mat)[4])

/* display node numbers in 3D location */

{
    double x,y,z,tx,ty,tz;

  x = npnt->xloc;
  y = npnt->yloc;
  z = npnt->zloc;
  if (z==LARGENODE) z = 0.0;
  transf (x,y,z,&tx,&ty,&tz,mat); 
  if (color) gpen (color);
  else       gpen (WHITE);
  gmove (tx, ty);
  dr_node (npnt,dscale);
}

/*--------------------------------------*/

void drcable (double x1, double y1, double z1, 
		double x2, double y2, double z2, 
		double dia, double dia2, 
		double dscale, double n1dia, double n2dia, 
		int color, int hide, double (*mat)[4])
                
/*
   Draw a cable element in 2D projection of 3D rotation.  First, 
transform the element's 2 end points by a 3D rotation matrix, and
then draw the 2D projection of the element in a subframe, in a 
standard position (0,0) and orientation (horizontal).  The subframe is 
translated ("origin()") to one endpoint of the element,
and rotated ("rotate()") by an amount equal to the original
element's 3D orientation angle projected into 2D. 

*/

{
    double tx1,ty1,tz1, tx2,ty2,tz2;
    double width,width2,tdia,tdia2;
    double dx,dy,dz,dist,disto,dratio,theta,r1,r2;
    int fill;

  if (dia2<0) dia2 = dia;
  if (z1==LARGENODE) z1 = 0.0;
  if (z2==LARGENODE) z2 = 0.0;
  if (dscale<0) {
     dscale = -dscale;
     // if (dia < dscale) dia = dscale;	/* make diameter at least = dscale */
     tdia = dia * xyscal;
     // if (dia2 < dscale) dia2 = dscale;	/* make diameter at least = dscale */
     tdia2 = dia2 * xyscal;
  }
  else {
     tdia  = dia *  dscale * xyscal;
     tdia2 = dia2 * dscale * xyscal;
  }
  r1 = n1dia * dscale * xyscal / 2.0;
  r2 = n2dia * dscale * xyscal / 2.0;
  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  disto = sqrt(dx*dx + dy*dy + dz*dz) * xyscal;
  if (disto==0.0) disto = 1e-6;
  transf (x1,y1,z1,&tx1,&ty1,&tz1,mat);
  transf (x2,y2,z2,&tx2,&ty2,&tz2,mat);
  dy = ty2-ty1;
  dx = tx2-tx1;
  dist = sqrt(dx*dx + dy*dy);
  dratio = dist / disto;		/* amount of foreshortening */
  r1 *= dratio;
  r2 *= dratio;
  /* ncfprintf (stderr,"dratio %g\n",dratio); /* */
  if (ncabs(dx)<1e-8) dx = .0001;
  theta = atan(dy/dx);
  if (dx<0) theta += PI;
  gframe ("cable");
  gorigin (tx1,ty1);
  grotate (theta/DEG);  
  if (color!=NULLVAL) gpen (color);
  else       gpen (GREEN);
  if (hide) {
    if (r1+r2 < dist) {
      gmove (r1,     -tdia/2.0);
      gdraw (r1,      tdia/2.0);
      gdraw (dist-r2, tdia2/2.0);
      gdraw (dist-r2,-tdia2/2.0);
      gdraw (r1,     -tdia/2.0);
    } 
  }
  else {
     double x1,y1,x2,y2;

   if (r1+r2 < dist) {          /* if cable not inside sphere */
     width  = tdia  / 2.0;
     width2 = tdia2 / 2.0;
     fill = 1;
     if (width > width2) {	/* if tapered */
        x1 = r1;
        y1 = -width2;
        x2 = dist-r2;
        y2 = -width;
        gtri(x1,y1,x2,y1,x1,y2,fill);
        y1 = width2;
        y2 = width;
        gtri(x1,y1,x2,y1,x1,y2,fill);
     }
     else if (width < width2) {
        x1 = dist-r2;
        y1 = -width;
        x2 = r1;
        y2 = -width2;
        gtri(x1,y1,x2,y1,x1,y2,fill);
        y1 = width;
        y2 = width2;
        gtri(x1,y1,x2,y1,x1,y2,fill);
     }
     if (width > width2) width = tdia2 / 2.0;
     x1 = r1;
     y1 = -width;
     x2 = dist-r2;
     y2 = width;
     grect(x1,y1,x2,y1,x2,y2,x1,y2,fill);
   }
  }
  gframe ("..");
  grmframe ("cable");
/*
  ncfprintf (stderr,"x1 %g y1 %g x2 %g y2 %g theta %g\n",tx1,ty1,tx2,ty2,theta);
*/
}

/*--------------------------------------*/

void drsphere (double x, double y, double z, double dia, double dscale, 
	int color, int hide, double (*mat)[4],int fill)
{
    double tx,ty,tz;

  if (dscale<0) dscale = -dscale;
  if (z==LARGENODE) z = 0.0;
  transf (x,y,z,&tx,&ty,&tz,mat); 
  if (color!=NULLVAL) gpen (color);
  else       gpen (BLUE);
  gmove (tx, ty);
  gcirc (dia*dscale*xyscal/2.0,fill);
}

/*--------------------------------------*/

/*  Don't use this subroutine, assume graphics driver
    has its own circle and rectangle routines. */

/*
void ncdrcirc (double x, double y, double rad, int fill)
{ 
   double theta, chord;
   double plXpos, plYpos, sin(double t), cos(double t);
   
#define rmove(dx,dy)	gmove(plXpos + (dx), plYpos + (dy))
#define rpmove(dr,t)	rmove((dr)*cos((double)t), (dr)*sin((double)t))
#define rdraw(dx,dy)	gdraw(plXpos + (dx), plYpos + (dy))
#define rpdraw(dr,t)	rdraw((dr)*cos((double)t), (dr)*sin((double)t))

 gmove (x,y);
 plXpos = x;
 plYpos = y;
 rpmove (rad,-PI/2);
 if      (rad > .2) chord = .05;
 else if (rad > .1) chord = .1;
 else               chord = .2;
 rpmove (rad*chord/2,-PI);
 for (theta=0.0; theta < PI*2; theta += chord)
   rpdraw (rad*chord,theta);

}
*/

/*--------------------------------------*/

void drsynap (double x1, double y1, double z1, 
	double x2, double y2, double z2, synapse* epnt,
	double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4])
{
    double dx,dy,tx1,ty1,tz1, tx2,ty2,tz2,theta,length,dist;

  if (color==NULLVAL) color=RED;
  if (z1==LARGENODE) z1 = 0.0;
  if (z2==LARGENODE) z2 = 0.0;
  transf (x1,y1,z1,&tx1,&ty1,&tz1,mat);
  transf (x2,y2,z2,&tx2,&ty2,&tz2,mat);
  dy = ty2-ty1;
  dx = tx2-tx1;
  dist = sqrt(dx*dx + dy*dy) / xyscal;
  length = dist;
  if (ncabs(dx)<1e-8) dx = .0001;
  theta = atan(dy/dx);
  if (dx<0) theta += PI;
  gframe ("syn");
  gorigin (tx1,ty1);
  gsize (xyscal);
  if (ncabs(theta)>1e-4) grotate (theta/DEG);  
  gmove (0,0);
  if ((xyscal*dist)>.0002) 
	(*syd)(epnt, color, epnt->vrev, dscale, 1.0, length, dist, hide);
  gframe ("..");
  grmframe ("syn");
}

/*--------------------------------------*/

void syn_draw (synapse *spnt, int color, double vrev, double dscale, double dia, 
		double length, double foreshorten, int hide)
{
    int fill=0;

  dia *= dscale;			/* draw circle with line */
  if (dia < 0) dia = -dia;
  if (color < 0) {
    if (vrev < -.04) color = RED;
    else             color = CYAN;
  }
  gpen (color);
  if (length > 1e-3) {
     gmove (length/2.0,0.0);
     if (dia > 0.001) gcirc (dia/2.0,fill);
     else             gcirc (0.001,fill);
     gmove (0,0);
     gdraw (length,0);
  }
  else                gcirc (0.001,fill);
}

/*--------------------------------------*/

void drcconn (double x1, double y1, double z1, double x2, double y2, double z2)
{
 if (disp_ray)
   rayconn (x1,y1,z1,x2,y2,z2, AXIALRES, 1.0, 0.0, 0.0, YELLOW, 0, ncmatc);
 else
   drconn (x1,y1,z1,x2,y2,z2, AXIALRES, 1.0, 0.0, 0.0, YELLOW, 0, ncmat);
}

/*--------------------------------------*/

void drconn (double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	int ctype, double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4])
{
    double dx,dy, tx1,ty1,tz1, tx2,ty2,tz2, dist, length, theta;

  if (z1==NULLVAL) z1 = 0.0;
  if (z2==NULLVAL) z2 = 0.0;
  if (color==NULLVAL) color = YELLOW;
  transf (x1,y1,z1,&tx1,&ty1,&tz1,mat);
  transf (x2,y2,z2,&tx2,&ty2,&tz2,mat);
  dy = ty2-ty1;
  dx = tx2-tx1;
  dist = sqrt(dx*dx + dy*dy) / xyscal;
  if (ncabs(dx)<1e-8) dx = .0001;
  theta = atan(dy/dx);
  if (dx<0) theta += PI;
  gframe ("gj");
  gorigin (tx1,ty1);
  gsize (xyscal);
  grotate (theta/DEG);  
  gmove (0,0);
  length = 1.0;

  if ((xyscal*dist)>.0002) 
    switch (ctype){

    case GJ: (*gjd) (color, dscale, length, dist, hide);
      break;

    case RESISTOR:
    case AXIALRES:
      gpen (color);
      gmove (0,0);
      gdraw (dist,0);
      break;
   } 
  gframe ("..");
  grmframe ("gj");
}

/*--------------------------------------*/

void  gj_draw (int color, double dscale, double length, 
			double foreshorten, int hide)
{
   double dia;
   int fill;

 fill = 0;
 gpen (color);
 if (dscale>=0) {
   dia = 1.0;
   if (length > 1e-4) gmove (foreshorten*length/2.0, 0);
   gcirc (dia/2.0,fill);
   gmove (0,0);
   if (length > 1e-4) gdraw (foreshorten*length,0);
 }
 else {
   dia = 1.0 * -dscale;
   if (length > 1e-4) gmove (length/2.0, 0);
   gcirc (dia/2.0,fill);
 }
}

/*--------------------------------------*/

void drphotrec (double x, double y, double z, 
		     photorec *epnt, double dscale,
		     double n1dia,double n2dia,
		     int color, int hide, double (*mat)[4])

/* see "drcable()" above for explanation of
   frame, origin() and rotation.  The screen
   rotation in this case is derived from translating
   two endpoints of a line drawn along the z axis.
   The subframe "photrec" is rotated in relation to
   the root frame, so that photoreceptors are always
   drawn (created) along the horizontal (X) axis.
 */
 
{
    double tx1,ty1,tz1,tx2,ty2,tz2;
    double dia,len,dist,dx,dy,theta;
    int ctype,pigm;
    static int color_trans[TOTREC] = {MAGENTA,RED,GREEN,BLUE};

  ctype = epnt->ctype;
  if ((dia=epnt->dia)==NULLVAL) dia = 2 * dscale;
  if (dia < 0) dia = -dia;
  if (color==NULLVAL) {
      if ((pigm=(int)(epnt->pigm)) == NULLVAL)
         pigm = (ctype==CONE ? 1 : 0);
      pigm = limit(pigm,0,TOTREC-1);
      color= color_trans[pigm];
  }
  len = 1.0;				/* length of photorec in um */
  if (dia<=0) dia = 2.0;
  if (z==NULLVAL) z = 0.0;
  transf (x,y,z,&tx1,&ty1,&tz1,mat);
  transf (x,y,z+len,&tx2,&ty2,&tz2,mat);
  dy = ty2-ty1;
  dx = tx2-tx1;
  if (ncabs(dy)<1e-8) dy = .0005;
  dist = sqrt(dx*dx + dy*dy) / xyscal;
  if (ncabs(dx)<1e-8) dx = .0005;
  theta = atan(dy/dx);
  if (dx<0) theta += PI;
  gframe ("photorec");
  gorigin (tx1,ty1);
  gsize (xyscal);
  grotate (theta/DEG);  
  gmove (0,0);
  gdraw (dist*0.01,0);
  gmove (0,0);
  if ((xyscal*dist)>.0002) 
  (*phd)(ctype,pigm,color,dscale,dia,dist,hide);
  gframe ("..");
  grmframe ("photorec");
}

/*--------------------------------------*/

void phot_draw(int ctype, int pigm, int color, 
		double dscale, double dia, double foreshorten, int hide)

/* Default subroutine to draw a photoreceptor */
/* Given that the pen has already been set at the origin (0,0),
   draw the photoreceptor's icon.  Icon is oriented horizontally
   to the right (zero degrees).
*/

{
  double width,dist;

  dist = foreshorten * 10.0;		/* photoreceptor extends 10 um */
  dia *= dscale;
  gpen (color);
  switch (ctype) {

  case ROD:
     width = dia*.75;
     if (hide) {
        gdraw (dist*.25,0);
        gmove (dist*.25,-width*.5);
        gdraw (dist*.25,width*.5);
        gdraw (dist,width*.5);
        gdraw (dist,-width*.5);
        gdraw (dist*.25,-width*.5);
     }
     else {
          int fill=1;
          double x1,y1,x2,y2;

       gdraw (dist*.25,0);
       x1 = dist*.25;
       y1 = -width*.5;
       x2 = dist;
       y2 = width*.5;
       grect(x1,y1,x2,y1,x2,y2,x1,y2,fill);
     }
     break;

  case VTRANSDUCER:
  case ITRANSDUCER:
  case CHR:
  case CONE:
	width = dia;
	gdraw (dist*.25,0);
	gmove (dist*.25,width/2.);
	gdraw (dist*.25,-width/2.);
	gdraw (dist, 0.0);
	gdraw (dist*.25,width/2.);
     break;

  }	/* switch */
}

/*--------------------------------------*/

void drload (double x, double y, double z, int ctype, 
	double n1dia, double n2dia, int color, int hide, double (*mat)[4])
                              
/* see "drphotrec()" above for explanation of
   frame, origin() and rotation. */
 
{
    double tx1,ty1,tz1,tx2,ty2,tz2;
    double len,dist,dx,dy,theta;

  len = 10.;				/* length of load in um */
  if (z==NULLVAL) z = 0.0;
  transf (x,y,z,&tx1,&ty1,&tz1,mat);
  transf (x,y,z+len,&tx2,&ty2,&tz2,mat);
  dy = ty2-ty1;
  dx = tx2-tx1;
  dist = sqrt(dx*dx + dy*dy);
  if (ncabs(dx)<1e-8) dx = .0001;
  theta = atan(dy/dx);
  if (dx<0) theta += PI;
  gframe ("load");
  gorigin (tx1,ty1);
  grotate (theta/DEG);  
  if (color!=NULLVAL) gpen (color);
  else       gpen (BLUE);
  switch (ctype) {

  case LOAD:
        gmove (0, 0);
        gdraw (dist*.25,0);
        gmove (dist*.25,-dist*.15);
        gdraw (dist*.25,dist*.15);
        gdraw (dist*.75,dist*.15);
        gdraw (dist*.75,-dist*.15);
        gdraw (dist*.25,-dist*.15);
        gmove (dist*.75,0);
        gdraw (dist, 0);
     break;

  case GNDCAP:
     break;

  }	/* switch */
	
  gframe ("..");
  grmframe ("load");

}

/*--------------------------------------*/

void elec_draw(int color, double dscale, double dia, double length, 
		double foreshorten, int hide, int elnum)

/* Default subroutine to draw a recording electrode */
/* Given that the pen has already been set at the origin (0,0),
   draw the recording electrode's icon.  Icon is oriented 
   horizontally to the right (zero degrees).
   If dscale has been set, then draw the node number according to
   "dr_node()".
*/

{
   elem *epnt;
   node *np;

  dia *= dscale;
  gpen (color);
  if (color!=NULLVAL) gpen (color);
  else       gpen (MAGENTA);
  gcirc (dia/2.0,hide=1);
  gmove (0.02,0);
  epnt = findelem(elnum);
  if (dscale!=1.0 && ((np=epnt->nodp1)!=NULL)) 
     dr_node(np,dscale);
}

/*--------------------------------------*/

void drelec (double x, double y, double z, 
                     electrode *epnt, double dscale,
                     int color, int hide, double (*mat)[4])

/* Draw recording electrode, taken from drphotrec. 
   See "drcable()" above for explanation of frame, 
   origin() and rotation.  The screen rotation in this 
   case is derived from translating two endpoints of 
   a line drawn along the z axis. The subframe 
   "photrec" is rotated in relation to the root frame, 
   so that recording electrodes are always drawn (created) 
   along the horizontal (X) axis.
 */
 
{
    double tx1,ty1,tz1,tx2,ty2,tz2;
    double dia,len,dist,dx,dy,theta;
    int ctype,pigm;
    static int color_trans[TOTREC] = {MAGENTA,RED,GREEN,BLUE};

  ctype = epnt->ctype;
  if ((dia=epnt->dia)==NULLVAL) dia = 2;
  if (dia < 0) dia = -dia;
  if (dia<=0) dia = 2.0;
  len = epnt->length;             /* length of recording elec in um */
  if (len<=0) len = 10.0;
  if (z==NULLVAL) z = 0.0;
  transf (x,y,z,&tx1,&ty1,&tz1,mat);
  transf (x,y,z+len,&tx2,&ty2,&tz2,mat);
  dy = ty2-ty1;
  dx = tx2-tx1;
  dist = sqrt(dx*dx + dy*dy) / xyscal;
  if (ncabs(dx)<1e-8) dx = .0001;
  theta = atan(dy/dx);
  if (dx<0) theta += PI;
  gframe ("electrode");
  gorigin (tx1,ty1);
  gsize (xyscal);
  grotate (theta/DEG);  
  gmove (0,0);
  gdraw (dist*0.01,0);
  gmove (0,0);
  if (xyscal>.0002) 
   (*rhd)(color,dscale,dia,len,dist,hide,epnt->elnum);
  gframe ("..");
  grmframe ("electrode");
}

/*--------------------------------------*/

void drcalib(double x, double y, double len, double size, int color)
{
    double dist;
    static char numbuf[20];

 if (disp_ray) raycalib (x, y-.03, len, size, color);
  else {
   dist = len * xyscal;
   if (color!=NULLVAL) gpen(color);
   else       	       gpen(CYAN);
   gmove (x,y);
   gdraw (x-dist,y);
   gmove (x-dist/2-.03,y-.03);
   sprintf (numbuf,"%g um",len);  
   gtext (numbuf);  
   gpurge();
  }
}

/*--------------------------------------*/

int loccheck (double x, double y, double z)

/* check whether point is at valid location */

{

   int locok;

  locok = 1;
  if (x==NULLVAL) locok = 0;
  if (y==NULLVAL) locok = 0;
  if (z==NULLVAL) locok = 0;
  if (x>=LARGNOD) locok = 0;
  if (y>=LARGNOD) locok = 0;
  if (z>=LARGNOD) locok = 0;
  return (locok);
}

/*--------------------------------------*/

void set_icons(void)

{
  Symbol *f,*lookup(const char *name);

  if (rhd==NULL) {
    if ((f=lookup("elec_dr"))==NULL) rhd = elec_draw;
    else {
      elec_drpnt = f;	/* set pointer to user function */
      rhd = elec_dr;      /* set pointer to call user function */
    }
  }
  if (phd==NULL) {
    if ((f=lookup("phot_dr"))==NULL) {
      phd = phot_draw;
    }
    else {
      phot_drpnt = f;	/* set pointer to user function */
      phd = phot_dr;      /* set pointer to call user function */
    }
  }
  if (syd==NULL) {
    if ((f=lookup("synapse_dr"))==NULL) syd = syn_draw;
    else {
      syn_drpnt = f;	/* set pointer to user function */
      syd = synapse_dr;      /* set pointer to call user function */
    }
  }
  if (gjd==NULL) {
    if ((f=lookup("gapjunc_dr"))==NULL) gjd = gj_draw;
    else {
      gapj_drpnt = f;	/* set pointer to user function */
      gjd = gapjunc_dr;      /* set pointer to call user function */
    }
  }
}

/*--------------------------------------*/

void set_elec_dr (void (*elec_drpnt)(int,double,double,double,double,int,int))
{
  if (elec_drpnt) rhd = elec_drpnt;
}

/*--------------------------------------*/

void set_phot_dr (void (*phot_drpnt)(int,int,int,double,double,double,int))
{
  if (phot_drpnt) phd = phot_drpnt;
}

/*--------------------------------------*/

void set_synapse_dr (void (*synapse_drpnt)(synapse *,int,double,double,double,double,double,int))
{
  if (synapse_drpnt) syd = synapse_drpnt;
}

/*--------------------------------------*/

void set_gapjunc_dr (void (*gapjunc_drpnt)(int,double,double,double,int))
{
  if (gapjunc_drpnt) gjd = gapjunc_drpnt;
}

/*--------------------------------------*/

void elec_dr (int color, double dscale, double dia, double length, 
		double foreshorten, int hide, int elnum)

/* Call a user-defined subroutine for drawing a recording point. */

{
   double x=0;

 callfunc8
(elec_drpnt,6,(double)color,dscale,dia,length,foreshorten,(double)hide,(double)elnum,x);
}

/*--------------------------------------*/

void phot_dr (int type, int pigm, int color, double dscale, double dia, 
		double foreshorten, int hide)

/* Call a user-defined subroutine for drawing a photoreceptor icon. */

{
   double x=0;

   callfunc8 (phot_drpnt,7,(double)type, (double)pigm, (double)color, 
		dscale, dia, foreshorten, (double)hide, x);
}

/*--------------------------------------*/

void synapse_dr (synapse *, int color, double vrev, double scale, double dia, 
	double length, double foreshorten, int hide)

/* Call a user-defined subroutine for drawing a synapse. */
/* Arg1, a pointer to the synapse, is used in an alternate syn_dr() function */

{
   double x=0;

   callfunc8 (syn_drpnt,7, (double)color, vrev, scale, dia, length, 
		foreshorten,(double)hide,x);
}

/*--------------------------------------*/

void gapjunc_dr (int color, double dscale, double length, 
		double foreshorten, int hide)

/* Call a user-defined subroutine for drawing a photoreceptor icon. */

{
   double x=0;

   callfunc8 (gapj_drpnt,5,(double)color, dscale, length, foreshorten,   
		(double)hide,x,x,x);
}

/*--------------------------------------*/
