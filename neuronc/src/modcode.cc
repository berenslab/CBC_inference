/* modcode.cc */

/* This module contains all subroutines called by "nc.y"
   that create neural element data structures.  This
   module is compiled for use by both "nc" and "stim".
   For the "stim" version, conditional compilation statements
   remove or modify some of the subroutines so that they
   do only things associated with making and controlling
   stimuli.  Thus, both "stim" and "nc" interpret the same
   language but "stim" ignores all the "compartment" stuff and
   deals only with neural elements and stimuli.  After
   compilation, the "stim" version is "mv"ed to "modcodes.o.
*/

extern "C" {

//#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "stdplt.h"
#undef scale
}

#include "nc.h"
#include "y.tab.h"
#include "ndef.h"
#include "control.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "ncplot.h"
#include "stim.h"
#include "gprim.h"
#include "vec.h"


#include "ncio.h"

/* #define DEBUG 	/* */

#ifdef DEBUG
#include "ncdebug.h"
#endif

/* To compile "stim" to totally ignore neural elements, define XSTIMM.  

   In a complex circuit, it is possible that the arrangement
   of photoreceptors and their stimuli need to be determined from
   information only available from some already-defined neural
   elements.  

   To compile stim to understand neural elements define XSTIM.

*/

#ifdef XSTIMM		/* define XSTIMM when stim should ignore neur elems */
  #define XSTIM		/* define XSTIM for stim  */
#endif

#ifndef XSTIM
#define XMOD			/* define XMOD when not XSTIM */
#endif
		
/* variables for stim */

int recepnode = 0;                      /* receptor node from conn() */
int recepn2 = 0;
int recepn3 = 0;
int stimfflg = 0;		/* =1 -> stimulus comes from file */
int stimflg = 0;		/* =1 -> stim functions, create stim file */

/* interpreter innards: */

extern datum *stack;
extern int interp;			/* running the interpreter */

#define pushm(d) stack[stackp++] = (d)		/*  */
#define popm()  stack[--stackp]		/* function pop still needed */

/* #define popm() pop()			/*  */

#define NSTACK 256

extern char *progname;
extern char *runfile;
extern int vidmode;
extern int debug;
extern int debugz;
extern int setdebug;
extern int setdebgz;
extern int setdisp;
extern int disp;
extern int disp_ray;
extern int stimhlp;		/* in stimsub.c */
extern int stimelem;		/* = 1 => stim understands neural elements */
extern int nolabels;    	/* no labels in morphology display */

extern int stopping;
extern int returning;
extern int breaking;
extern int continuing;

extern  int   stackp;
extern	Inst	*pc;		/* program counter during execution */

elem *elpnt;			/* current element pointer */

int curelnum=0;			/* current element being read in */
extern node *nodepnt;		/* pointer to node list */
extern photrec *recpnt;		/* pointer to photrec list */
extern recnod *reclist;		/* pointer to photrec stimulus list */
extern elem *elempnt;		/* pointer to current element */
extern int cumelem;
extern int cumnode;
extern int cumrec;

extern double simtime;		/* time variable set in "ncsub" (control.h) */
extern Symbol *timeptr;		/* pointer to "time" for setting time */


#ifdef __cplusplus
extern "C" {
#endif

extern char *strncat(char *dest, const char *src, size_t n);
int mkstemp(char *pnt);

#ifdef __cplusplus
}
#endif

#ifdef XMOD		/* variables for modcode */

#define STIMCODE 0

extern plotfr plotnod[PLOTNODSIZ];  /* holds nodes to be recorded */
extern int plotnum;		/* number of plots to be displayed */
extern int runyet;
			/* functions for modcode */

recstim *makvstim(double start, int nodnm1, int nodnm2, int nodnm3, int nodnm4,
			double value, int bgend, const char *action);
recstim *makrstim(double start, int recnm1, int recnm2, int recnm3, int recnm4,
		double inten, double wavel, double mask, int stimchan, int bgend, const char *action);

int inrect(double x, double y, double maxx, double minx, 
			double maxy, double miny, double orient);
int incirc(double x, double y, double circx, 
			double circy, double rad);

extern double ncmat[4][4];
void transf (double x,double y,double z, double *tx,double *ty,double *tz, 
							double (*mat)[4]);
/* functions for nc */

void makcheckerboardx(double width, double height, int xn, int yn, double orient,
		double xloc, double yloc, double tfreq, double inten, double contrast,
		const char *action, double start, double dur, double **rndarr, int *nfr, int ckrseed);

void makgaborx(double speriod, double sphase, double orient, double xloc, double yloc,
		double tfreq, double drift, double inten, double contrast, double wavel,
		double xenv, double yenv, int makenv, int sq, double mask, int stimchan,
		const char *action, double start, double dur);

void maksineannx(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten, double contrast,
		double wavel, double renv, int makenv, int sq, double mask, int stimchan,
		const char *action, double start, double dur);

void makwindmillx(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten, double contrast,
		double wavel, double renv, int makenv, int sq, double mask, int stimchan,
		const char *action, double start, double dur);

void sector_mask(double orient, double width, double val, double stimtime, double dur);

#else			/* variables for makestim */

#define STIMCODE 1

plotfr plotnod[PLOTNODSIZ];  	/* holds nodes to be recorded */
int plotnum;			/* number of plots to be displayed */
int runyet;

			/* functions for makestim */
void stcomment (void);
int makblur(double rad, double blur_ampl, double scale, double scatter_rad, 
			double scatter_pow, double scatter_ampl);
void find_photarrsiz(double *xmax,double *xmin,double *ymax,double *ymin);
int initc(int num, int size);
void recback(int array, double inten, double wavel);
void abslist(double ratio, double time, double mask, int stimchan);
void makrect(int arr, int arrsize, double width, double length, 
		double xoff, double yoff, double sscale, double orient, 
		double inten, double wavel);
void recpat(int array, int arrsize, double xoff, double yoff, double sscale);
void stimlist(double delratio, double time, double mask, int stimchan);

/* functions for stim */
void makcheckerboard(int arr, int arrsiz, double xsize, double ysize, int xn, int yn, 
		double xloc, double yloc, double xcent, double ycent,
		double scale, double orient, double tfreq, 
		double inten, double contrast, double start, double dur,
		double **stim_rndarr, int *stim_nfr, int ckrseed);

void makgabor(int arr, int arrsize, double contrast, double inten,
                double xloc, double yloc, double xcent, double ycent, 
		double xenv, double yenv, int makenv, int sq, 
		double tfreq, double drift, double speriod, double sphase,
                double scale, double orient, double wavel,
                double mask, double start, double dur);

void maksineann(int arr, int arrsize, double contrast, double inten,
                double xloc, double yloc, double xcent, double ycent, 
		double renv, int makenv, int sq,
                double tfreq, double drift, double speriod, double sphase,
                double scale, double wavel, double mask,  
		double start, double dur);

void makwindmill(int arr, int arrsize, double contrast, double inten,
                double xloc, double yloc, double xcent, double ycent, 
		double renv, int makenv, int sq,
                double tfreq, double drift, double speriod, double sphase,
                double scale, double wavel, double mask, 
		double start, double dur);

void makimage(int arr, int arrsize, double width, double length,
                double xoff, double yoff, double orient, double scale,
                double inten, int type, const char *imgfil);

void makspot(int arr, int arrsize, double dia, double xloc, double yloc,
		double xcent, double ycent, double sscale, double inten, 
		double start, double dur, double wavel, double mask, int stimchan, int invert);
void makbar(int arr, int arrsize, double width, double length, double xloc, double yloc, 
		double xcent, double ycent, double scale, double orient, double inten,
                double start, double dur, double wavel, double mask, int stimchan);

void delc(int num);

#ifdef __cplusplus
extern "C" {
#endif

/* void mktemp(char *f);	/* */
char *mktemp(char *f);	/* */

#ifdef __cplusplus
}
#endif

#endif				/* makestim */

int setstimchan(int stimchan);

extern char *rootframe;

double xrot,yrot,zrot;
double dxcent=LARGENODE,dycent=LARGENODE,dzcent=LARGENODE;
double rxcent=0,rycent=0,rzcent=0;
double dsize = 200;
int cmap=0;
int calibcolor = NULLVAL; 
static int synnum=0;		/* temporary store for synaptic "ename" */
extern int nocond;		/* no condense, in ncsub */

char *stfile=0;			/* stimulus file for plotinit */

// extern FILE *stdout;
// extern FILE *stderr;

FILE *stimin=0;			/* stimulus file */
FILE *stimout=stdout;           /* stimulus output file, used only by stim */
static char stimfile[80]={0};

#ifdef __cplusplus
extern "C" {
#endif

double sqrt(double);
double sin(double);

#include "gr.h"

FILE *fopen(const char *path, const char *mode);

#ifdef __cplusplus
}
#endif

node *maknod(nodeint nodea, nodeint nodeb, nodeint nodec,
	nodeint noded);
elem *makelem(int etype, elem *epnt);
elem *tmpelem();
elem *eninstall (elem *epnt);   
double getval(const char *str);
double mindist(double dist1, double dist2, double dist3);
char *prnode(int n1, int n2, int n3, int n4);
char *prnode(node *npnt);
void execerror(const char *s, const char *t);
void ncleanup();
void savemodel(const char *filnam);
void restoremodel(const char *filnam);
Symbol *getvar(pdatum *p);
void erasrecnod(elem *epnt);	/* erase possible rec stim */
void checknod(elem *elpnt);
void checkelemnode(elem *elpnt);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, nodeint node1d, const char *s);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, const char *s);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c);
void setnodlst(node *npnt, elem *cpnt);
void unsetnodlst(node *npnt, elem *cpnt);
char *prnum(char *adj, double csiz, char *fmt, double v1, double v2, double v3);
int gausnn (double mean, double stdev, double density, double ms, int grseed,
	double framex, double framey, double xcent, double ycent, 
	int numcells, double **xarr, double **yarr, int first_center,
	int filout, int textfl, int printfl);
double *darr2(Symbol *sp, int narg);
void varcopy();
void varcopyu();
int checknum(datum d1);
int checkstr(datum d1);
double *makfiltarr (int newsiz, int oldsiz, double *oarr, double val);
void plotinit(int plotnum);
void findconnect(int donefl);
void initcomp(void);
void condense(void);
void execute(Inst *pc);
int forexec(Inst *pc);
void maklst(conlst **head, conn *cpnt);
void initchan();
void runsim(double time, int run);
int smaxmin(double val);
char *findsym(int num);
int numdyad(dyadlst *head);
double norm_line (vec a, vec b, vec p, vec q, double *t, double *s,
			vec *lp1, vec *lp2);
void efree(void *ptr);
void delnode(node *npnt);
void delelem(elem *epnt);
elem *findelem(int num);
void copyattr(attrib *src, attrib *dest);
attrib *getattrib(int elnum);
void plotfunc(Symbol *func, int nplot);
void drcalib (double x, double y, double length, double size, int color);
void setrot (double xrot, double yrot, double zrot, 
	double xcent, double ycent, double zcent, 
	double rxcent, double rycent, double rzcent, double scal);
void initray (double xrot, double yrot, double zrot,
        double xcent, double ycent, double zcent,
        double rxcent, double rycent, double rzcent, double scal);
void set_icons(void);
void dispstim(double stime, double dscale, int cmap, double flux_max, double flux_min);
void dispstim(double sttime, double stime, double dscale, int cmap, double flux_max, double flux_min);
void dcomp(int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype, 
   	int na, int nb, int nc, int nd, 
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color), 
	double vmax, double vmin,
    	double dscale, int hide, int excl);
void ncdrnod (int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype, 
	int na, int nb, int nc, int nd, 
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color), 
	double dscale, int drlabel);
void ncdisp(int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype, 
	int na, int nb, int nc, int nd, 
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color), 
	double vmax, double vmin, 
	double dscale, int hide, int excl, int cmap);
void ncdispc(int n1a, int n1b, int n1c, int n1d, int n2a, int n2b, int n2c, 
	int n2d, int elemtype, int exceptype, int na, int nb, int nc, int nd,
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color), 
	double vmax, double vmin, double dscale, int hide, int cmap);
void ncdispn(int n1a, int n1b, int n1c, int n1d, int n2a, 
	int n2b, int n2c, int n2d, int elemtype, int exceptype, 
	int na, int nb, int nc, int nd,
	double zrange1, double zrange2, int color, 
	Symbol *vpen, double (*vpenn)(int elnum, int color), 
	double dscale, double vmax, double vmin, 
	int hide, int excl, int cmap);
void ncdispe (int elemnum, double zrange1, double zrange2, 
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color), 
	double vmax, double vmin, 
	double dscale, int hide, int cmap);
void set_display (int elemtype, int disptype, int narg2, 
	     int n1a,int n1b,int n1c,int n1d,
	     int n2a,int n2b,int n2c,int n2d,
	     int na, int nb, int nc, int nd, int exceptype,
	     double zrange1, double zrange2,
	     int dcolor, Symbol *vpen, double(*vpenn)(int, int), 
	     int cmap, double vmax,double vmin, 
	     double dscale, int hide, int excl, double stime);
void set_display (int elemtype, int disptype, int narg2, 
	     int n1a,int n1b,int n1c,int n1d,
	     int n2a,int n2b,int n2c,int n2d,
	     int na, int nb, int nc, int nd, int exceptype,
	     double zrange1, double zrange2,
	     int dcolor, Symbol *vpen, double(*vpenn)(int, int), 
	     int cmap, double vmax,double vmin, 
	     double dscale, int hide, int excl, double stime);
chantype *getchantype(int ctype, int cnum);
int setcmap(datum d1);
void break_fixup(int before, int after);	/* in code.cc */

/*------------------------------------------*/

/* Always remember to check a new subroutine added to "modcode.c" for
   its effect on the corresponding program "stim". 
   The stim subroutines should be called in exactly 
   the same way as their modcode counterparts (exception: stimulus subroutines),
   and they should always leave the stack (popm) and program counter
   (pc) correct.  However if XSTIMM is defined, they will
   ignore all neural elements to save time and memory.
*/ 

/*------------------------------------------*/

void eramod(void)

/* erase all of the model stuff, to allow
   more space to run something else. 
*/

{
  ncleanup();				/* erase all lists, reset pointers */
}
/*------------------------------------------*/

void savemod(void)

/* Save all state variables in the simulation

*/

{
     datum d1;
     const char *filnam;

  if (!checkstr(d1=popm())) {
      return;
  }
  filnam = d1.str;		/* file to be read */
  savemodel(filnam);		/* save compartment v, ca, nt, etc. */
}

/*------------------------------------------*/

void restoremod(void)

/* Restore all state variables in the simulation
   from a previous "save".
*/

{
     datum d1;
     const char *filnam;

  if (!checkstr(d1=popm())) {
      return;
  }
  filnam = d1.str;		/* file to be read */
  restoremodel(filnam);		/* erase all lists, reset pointers */
}

/*------------------------------------------*/

chattrib *chanattr(elem *elpnt, short int ctype) 

/* set up attributes for a channel,
   and link them with the parent element.
*/

{
#ifndef XSTIMM
   register attrib *pnt,*apnt;
   chattrib *makchattr(void);
   cattrib *makcattr(void);
   nattrib *maknattr(void);

  switch (ctype) {

  case GLU:
  case KAINATE:
  case GABA:
  case GLY:
  case CHRC:
  case NA:
  case  K:
  case  KCa:
  case  ClCa:
     apnt = (attrib *)makchattr();	/* make new NA/K channel attribute */
     break;
  case AMPA:
  case NMDA:
  case CGMP:
  case SYN2:
  case CACOMP:
  case CA:
     apnt = (attrib *)makcattr();	/* make new Ca chan attribute */
     break;
  case VESNOISE:
  case CCHNOISE:
  case NUNIT:
     apnt = (attrib *)maknattr();	/* make new noise attribute */
     break;
  } /* switch */

  apnt->ctype = ctype;			/* set type of new attrib */
  if (!elpnt->attpnt) {			/* if element has no attributes */
     elpnt->attpnt = apnt;		/* add directly to element */
  }
  else {				/* else add to element's attr list */
     for (pnt=elpnt->attpnt; pnt->attpnt; pnt=pnt->attpnt) ;
     pnt->attpnt = apnt;
  }
  return ((chattrib *)apnt);
#else
   chattrib *makchattr(void);
   chattrib *nullattrib = makchattr();  /* make null channel attribute */

  return (nullattrib);
#endif
}

/*------------------------------------------*/

void getnod (datum *d)

/* Get up to 4 node dimensions off the stack. */

{
    int i,narg;

  narg = (long int)*pc++;
  if (narg>MAXNODIM) {
     ncfprintf (stderr,"getnod: more than 4 dimensions\n");
     execerror ("Bad stack", "stopping");
  }
  if (narg > 3) checknum(d[3]=popm());	/* get node number */
  else {d[3].val = NULLVAL; d[3].vtype=NUMBER;}
  if (narg > 2) checknum(d[2]=popm());
  else {d[2].val = NULLVAL; d[2].vtype=NUMBER;}
  if (narg > 1) checknum(d[1]=popm());
  else {d[1].val = NULLVAL; d[1].vtype=NUMBER;}
  if (narg > 0) checknum(d[0]=popm());
  else {d[0].val = NULLVAL; d[0].vtype=NUMBER;}
  for (i=0; i<MAXNODIM; i++) 
    d[i].val = smaxmin(d[i].val);
}

/*------------------------------------------*/

elem *modify (int elnum)

/* set "modify" flag in an element */

{

#ifndef XSTIMM
  elpnt = tmpelem();
  elpnt->modif = elnum;		/* save original elemnum as ref. */
#endif
  return elpnt;
}

/*------------------------------------------*/

void xmod(void)

/* modify an element */

{
  datum d1;
  Symbol *param;
  pdatum p={0};

  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case MODIFY: if (checknum(d1 = popm())) {
#ifndef XSTIMM
    		    elpnt = tmpelem();			/* Make new element; */
  		    elpnt->modif = (int)(d1.val);	/* mark it modified. */
#endif
		 }
		 break;

    case ENAME:  getvar(&p);	
		 if (synnum) {		/* check for synaptic ename */
		   *p.val = synnum;	/* save element number in var */
		    synnum = 0;
		 }
		 else *p.val = cumelem;	/* save element number in var */
		 if (*p.type==UNDEF) *p.type = VAR;
		 *p.vtype = NUMBER;
		 elpnt->nocondens = 1;
		 break;

    case ELABL:  if (checkstr(d1=popm())) {
#ifndef XSTIMM
		   if (d1.str)
  		     elpnt->elabl = (char *)d1.str;
		   else 
  		     elpnt->elabl = (char *)"";
#endif
		 }
		 break;
   }
}

/*------------------------------------------*/

node *erasenode(node *npt)

/* Delete a node and any elements that connect to it. */

/* For each element on the node's element list, remove the
   element pointer in all the nodes connected to that element,
   then delete the element.  Finally delete the original node.
*/

{
   int found, erased;
   int nod1, nod2, nod3, nod4;
   conlst *lpnt,*lnext;
   elem *epnt;
   node *nnext;
   register node *npnt=npt;

#ifndef XSTIMM
 if (npnt==NULL) return NULL;
 nnext = npnt->next;
 erased = 0;
 if ((lpnt=npnt->elemlst)) {	/* if any elems point to node */
   for (lpnt=npnt->elemlst; lpnt; lpnt=lnext) { /* delete element */
      lnext = lpnt->next;
      epnt = (elem *)lpnt->conpnt;
      if (epnt->nodp1) {
	unsetnodlst(epnt->nodp1,epnt);     /* free the nodpnt to this elem */
      }
      if (epnt->nodp2) {
	unsetnodlst(epnt->nodp2,epnt);     /* free the nodpnt to this elem */
      }
      if (makestim) erasrecnod(epnt);	   /* erase possible rec stim */
      delelem(epnt);                       /* erase the element */
   }
 }
 if (!erased) delnode(npnt);
#endif
 return (nnext);
}

/*------------------------------------------*/

void eranode(void)

/* Delete a node and any elements that connect to it. */

/* For each element on the node's element list, remove the
   element pointer in all the nodes connected to that element,
   then delete the element.  Finally delete the original node.
*/

{
   int nod1, nod2, nod3, nod4;
   datum d[MAXNODIM];
   node *npnt;

 getnod(d);					/* get node, check dims */

#ifndef XSTIMM
 nod1 = (int)d[0].val;
 nod2 = (int)d[1].val;
 nod3 = (int)d[2].val;
 nod4 = (int)d[3].val;
 if (!(npnt=findnode(nod1,nod2,nod3,nod4,"Eranode"))) {
   execerror ("Missing node; ","stopping...");
   return; 
 }
 npnt = erasenode(npnt);
#endif
}

/*------------------------------------------*/

void eraseelem(elem *)

/* Delete an element and any nodes that connect only to it. */

{
  datum d1;
  elem *epnt;
  int elnum;

#ifndef XSTIMM
      checkelemnode(epnt);		/* set pointers if not done */
      if (epnt->nodp1) {
	unsetnodlst(epnt->nodp1,epnt);      /* free the nodpnt to this elem */
	if (!epnt->nodp1->elemlst) delnode(epnt->nodp1);
      }
      if (epnt->nodp2) {
	unsetnodlst(epnt->nodp2,epnt);      /* free the nodpnt to this elem */
	if (!epnt->nodp2->elemlst) delnode(epnt->nodp2);
      }
      if (makestim) erasrecnod(epnt);	  /* erase the rec stim */
      delelem(epnt);                      /* erase the element */
#endif
}

/*------------------------------------------*/

void erelem(void)

/* Delete an element and any nodes that connect only to it. */

{
  datum d1;
  elem *epnt;
  int elnum;

if (checknum(d1 = popm())) {

#ifndef XSTIMM
   elnum = int(d1.val);
   if (!(epnt=findelem(elnum))) {
        ncfprintf (stderr,"erelem: can't find element %d\n",elnum);
        execerror ("Missing element: ","stopping...");
        return;
    }
    else {				/* first find nodes that connect */
      eraseelem(epnt);
    }
#endif
  }
}

/*------------------------------------------*/

void xcable(void)
{
  datum d1;
  Symbol *param;
  cable *epnt;

  param  = (Symbol *)*pc++;
  checknum(d1 = popm());
  epnt = (cable *)makelem(CABLE,elpnt);		/* make new cable */
  elpnt = (elem *)epnt;				/* for membtyp() below */

  switch (param->type)
   {
    case LENGTH: epnt->length = d1.val; 
		break;
    case DIA:    epnt->dia    = d1.val;
		 if (epnt->dia2==NULLVAL) epnt->dia2 = epnt->dia;
		break;
    case DIA2:   epnt->dia2   = d1.val;
		break;
    case CPLAM:  epnt->cplam = d1.val;
		break;
    
   }
/* ncfprintf (stdout,"cable #%d %s %6.3f\n",
			elpnt->elnum,param->name,d1.val);  /* */
}

/*------------------------------------------*/

void xsphere(void)
{
  datum d1;
  Symbol *param;
  sphere *epnt;

  param = (Symbol *)*pc++;
  checknum(d1 = popm());
  epnt = (sphere *)makelem(SPHERE,elpnt);	/* make new sphere */
  elpnt = (elem *)epnt;				/* for membtyp() below */

   if (param==NULL) {
	epnt->dia = d1.val;
   }
   else switch (param->type)
    {
	case DIA: epnt->dia = d1.val;
		break;
	case RADIUS: epnt->dia = d1.val * 2;
		break;
    }
/*  ncfprintf (stdout,"sphere %s %g\n",param->name,d1.val); */
}

/*------------------------------------------*/

void xelec(void)
{
  datum d1;
  Symbol *param;
  electrode *epnt;

  param = (Symbol *)*pc++;
  epnt = (electrode *)makelem(ELECTRODE,elpnt);	/* make new electrode */
  elpnt = (elem *)epnt;

   if (param!=NULL) {
    checknum(d1 = popm());
    switch (param->type)
    {
	case DIA:     epnt->dia = d1.val;
		break;
	case LENGTH:  epnt->length = d1.val;
		break;
	case RS:      epnt->r = d1.val;
		break;
	case CAP:     epnt->c = d1.val;
		break;
	case VREST:   epnt->vrest = d1.val;
		break;
    }
  }
/*  ncfprintf (stdout,"electrode %s %g\n",param->name,d1.val); */
}

/*------------------------------------------*/

void xcachan(void)
{
  datum d1;
  Symbol *param;
  static cattrib *apnt; 
  int npops;
  static int pump=0;

  param  = (Symbol *)*pc++;
  elpnt = makelem(CHAN,elpnt);			/* make dummy elem "chan"  */
  switch (param->type)
   {
    case CA:
    case SYN2:
    case CGMP:
    case AMPA:
    case NMDA:
    case CAPUMP:
    case CAEXCH:
    case CABUF:
    case CICR:
    case IP3:
    	 	 npops=0; 
		 break;
    case CACOMP: 
    	 	 npops=0; 
		 break;
    case TYPE:  
    case VVREV:  
    case TAU:
    case TAUM:
    case TAUH:
    case TAUN:
    case TAUA:
    case TAUB:
    case TAUC:
    case TAUD:
    case TAUE:
    case TAUF:
    case OFFSET:
    case OFFSETM:
    case OFFSETH:
    case MAXCOND:
    case CAPERM:
    case CAKD:
    case CAHC:
    case DENSITY:
    case NDENSITY:
    case UNIT:
    case N:

    case CAO:
    case CAI:
    case TCAI:
    case CBOUND:
    case CASHELL:
    case CAOSHELL:
    case KEX:
    case VMAX:
    case KM:
    case MG:
    case KD:
    case BTOT:
    case BTOTI:

    case CAS2:		/* IP3 params */
    case IP3I:
    case BIP3:
    case VIP3:

    case CAS:		/* CICR params */
    case VM2:
    case VM3:
    case M:
    case P:
    case C1:
    case KA:
    case KF:
    case KR:
    case K2:

		 npops=1;
		break;
   }
 
  if (npops==1) checknum(d1 = popm());

  switch (param->type)
   {
    case CACOMP: 
		 if (elpnt->ctype==SYNAPSE) { 
		   apnt->caflg = 1;         /* don't make sep attr for cacomp */
		   break;
		 }
		 else /* if element has no attribute, it must be a cacomp */
		   if (!elpnt->attpnt) elpnt->ctype=CACOMP;

    case NMDA:
    case AMPA:
    case SYN2:
    case CGMP:					/* link attribs to elem */
    case CA:    apnt = (cattrib *)chanattr(elpnt,param->type); 
		pump = 0;
		break;
    case CAPUMP:apnt->pump = 1; pump = param->type;
		break;
    case CAEXCH:apnt->exch = 1; pump = param->type;
		break;
    case CABUF: apnt->cabuf = 1; pump = param->type;
		break;
    case CICR:  apnt->cicr = 1; pump = param->type;
		break;
    case IP3:  apnt->ip3 = 1; pump = param->type;
		break;
    case TYPE:  if (apnt) apnt->stype = (int)(d1.val);
		break;
    case VVREV:  if (apnt) apnt->vrev = d1.val;
		break;
    case TAU:  if (apnt) {
		  if (apnt->taua!= NULLVAL) apnt->taua *= d1.val;
		  else 			    apnt->taua  = d1.val;
		  if (apnt->taub!= NULLVAL) apnt->taub *= d1.val;
		  else 			    apnt->taub  = d1.val;
		  if (apnt->tauc!= NULLVAL) apnt->tauc *= d1.val;
		  else 			    apnt->tauc  = d1.val;
		  if (apnt->taud!= NULLVAL) apnt->taud *= d1.val;
		  else 			    apnt->taud  = d1.val;
		}
		break;
    case TAUM:  if (apnt) {
		  if (apnt->taua!= NULLVAL) apnt->taua *= d1.val;
		  else 			    apnt->taua  = d1.val;
		  if (apnt->taub!= NULLVAL) apnt->taub *= d1.val;
		  else 			    apnt->taub  = d1.val;
		}
		break;
    case TAUH:  if (apnt) {
		  if (apnt->tauc!= NULLVAL) apnt->tauc *= d1.val;
		  else 			    apnt->tauc  = d1.val;
		  if (apnt->taud!= NULLVAL) apnt->taud *= d1.val;
		  else 			    apnt->taud  = d1.val;
		}
		break;
    case TAUA:  if (apnt) {
		  if (apnt->taua!= NULLVAL) apnt->taua *= d1.val;
		  else 			    apnt->taua  = d1.val;
		}
		break;
    case TAUB:  if (apnt) {
		  if (apnt->taub!= NULLVAL) apnt->taub *= d1.val;
		  else 			    apnt->taub  = d1.val;
		}
		break;
    case TAUC:  if (apnt) {
		  if (apnt->tauc!= NULLVAL) apnt->tauc *= d1.val;
		  else 			    apnt->tauc  = d1.val;
		}
		break;
    case TAUD:  if (apnt) {
		  if (apnt->taud!= NULLVAL) apnt->taud *= d1.val;
		  else 			    apnt->taud  = d1.val;
		}
		break;
    case TAUE:  if (apnt) {
		  if (apnt->taue!= NULLVAL) apnt->taue *= d1.val;
		  else 			    apnt->taue  = d1.val;
		}
		break;
    case TAUF:  if (apnt) {
		  if (apnt->tauf!= NULLVAL) apnt->tauf *= d1.val;
		  else 			    apnt->tauf  = d1.val;
		}
		break;
    case OFFSET: if (apnt) {
		   apnt->voffsm  = d1.val;
		   apnt->voffsh  = d1.val;
		}
		break;
    case OFFSETM: if (apnt) apnt->voffsm  = d1.val;
		break;
    case OFFSETH: if (apnt) apnt->voffsh  = d1.val;
		break;
    case MAXCOND: if (apnt) apnt->maxcond = d1.val;
		break;
    case CAPERM: if (apnt) apnt->caperm = d1.val;
		break;
    case CAKD:  if (apnt) apnt->cakd = d1.val;
		break;
    case CAHC:  if (apnt) apnt->cahc = d1.val;
		break;
    case DENSITY: if (apnt) apnt->density = d1.val;
		break;
    case NDENSITY: if (apnt) apnt->ndensity = d1.val;
		break;
//    case N:    if (apnt) apnt->n = int(d1.val);
//		break;
    case N: 	  if (apnt) {
		    switch (pump) {
			case 0: apnt->n  = d1.val; break; // n unitary
			case CICR: apnt->ncicr  = d1.val; break;
		    }
		}
		break;
    case M: 	  if (apnt) {
		    switch (pump) {
			case CICR: apnt->mcicr  = d1.val; break;
		    }
		}
		break;
    case P: 	  if (apnt) {
		    switch (pump) {
			case CICR: apnt->pcicr  = d1.val; break;
		    }
		}
		break;
    case C1:   if (apnt) {
		  switch (pump) {
                      case CICR: apnt->c1cicr = d1.val; break;	  
       		  }
	       }
	       break;
    case UNIT: if (apnt) apnt->unitary = d1.val;
		break;
    case CAO: 	  if (apnt) apnt->cao = d1.val;
		break;
    case CAI: 	  if (apnt) apnt->cai = d1.val;
		break;
    case TCAI: 	  if (apnt) apnt->tcai = d1.val;
		break;
    case KEX: 	  if (apnt) apnt->kex   = d1.val;
		break;
    case BTOT:	  if (apnt) apnt->btot = d1.val;
		break;
    case BTOTI:	  if (apnt) apnt->btoti = d1.val;
		break;
    case KD:	  if (apnt) apnt->bkd = d1.val;
		break;
    case VMAX: 	  if (apnt) {
		    switch (pump) {
			case CAPUMP: apnt->vmax  = d1.val; break;
			case CABUF:  apnt->bmax  = d1.val; break;
		    }
		}
		break;
    case KM: 	  if (apnt) {
		    switch (pump) {
			case CAPUMP: apnt->pkm = d1.val; break; 
		    }
		  }
		break;
    case CBOUND:  if (apnt) apnt->cabnd = d1.val;
		break;
    case CASHELL:  if (apnt) apnt->cashell = (short int)d1.val;
		break;
    case CAOSHELL:  if (apnt) apnt->caoshell = (short int)d1.val;
		break;
    case MG: if (apnt) apnt->mg = d1.val;
		break;

    case CAS: 	  if (apnt) { 					/* CICR params */
		    switch (pump) {
			case CICR: apnt->cas  = d1.val; break;
		    }
		}
		break;
    case VM2:   if (apnt) {
		    switch (pump) {
			case CICR: apnt->vm2 = d1.val; break; 
		    }
		  }
		break;
    case VM3:   if (apnt) {
		    switch (pump) {
			case CICR: apnt->vm3 = d1.val; break; 
		    }
		  }
		break;
    case KA:   if (apnt) {
		    switch (pump) {
			case CICR: apnt->kacicr = d1.val; break; 
		    }
		  }
		break;
    case KF:   if (apnt) {
		    switch (pump) {
			case CICR: apnt->kfcicr = d1.val; break; 
		    }
		  }
		break;
    case KR:   if (apnt) {
		    switch (pump) {
			case CICR: apnt->krcicr = d1.val; break; 
		    }
		  }
		break;
    case K1:    if (apnt) {
	            switch (pump) {
			case CICR: apnt->k1cicr = d1.val; break;
	            }
	          }
                break;    
    case K2:   if (apnt) {
		    switch (pump) {
			case CICR: apnt->k2cicr = d1.val; break; 
		    }
		  }
		break;

    case CAS2: 	  if (apnt) {
		    switch (pump) {
			case IP3: apnt->cas2  = d1.val; break;
		    }
		}
		break;
    case IP3I: 	  if (apnt) {
		    switch (pump) {
			case IP3: apnt->ip3i  = d1.val; break;
		    }
		}
		break;
    case BIP3:   if (apnt) {
		    switch (pump) {
			case IP3: apnt->bip3 = d1.val; break; 
		    }
		  }
		break;
    case VIP3:   if (apnt) {
		    switch (pump) {
			case IP3: apnt->vip3 = d1.val; break; 
		    }
		  }
		break;

    case OI:    if (apnt) {
                  switch (pump) {
                      case IP3: apnt->oip3 = d1.val; break;
                  }
               }
               break;

    case A3:    if (apnt)  {
                   switch (pump)  {
                       case IP3: apnt->a3ip3 = d1.val; break;
                   }
                }
                break;

    case D1:    if (apnt)  {
                   switch (pump)  {
                       case IP3: apnt->d1ip3 = d1.val; break;
                   }
                }
                break;

    case D2:    if (apnt)  {
                   switch (pump)  {
                       case IP3: apnt->d2ip3 = d1.val; break;
                   }
                }
                break;

    case D3:    if (apnt)  {
                   switch (pump)  {
                       case IP3: apnt->d3ip3 = d1.val; break;
                   }
                }
                break;

    case D4:    if (apnt)  {
                   switch (pump)  {
                       case IP3: apnt->d4ip3 = d1.val; break;
                   }
                }
                break;

    case A2:    if (apnt)  {
                   switch (pump)  {
                       case IP3: apnt->a2ip3 = d1.val; break;
                   }
                }
                break;

   case B2:    if (apnt) {
                  switch (pump) {
                      case IP3: apnt->b2ip3  = d1.val; break;
                  }
               }
               break;

   case K3:    if (apnt) {
                   switch (pump) {
                       case IP3: apnt->k3ip3  = d1.val; break;
                   }
               }
               break;

   case V2:    if (apnt)  {
                  switch (pump)  {
                      case IP3: apnt->v2ip3 = d1.val; break;
                  }
               }
               break;

   case V3:    if (apnt) {
                  switch (pump) {
                      case IP3: apnt->v3ip3  = d1.val; break;
                  }
               }
               break;

   case V4:    if (apnt) {
                  switch (pump) {
                      case IP3: apnt->v4ip3  = d1.val; break;
                  }
               }
               break;

   case MTYPE: if (apnt) {
                  switch (pump) {
                      case IP3: apnt->mtypeip3  = d1.val; break;
                  }
               }
               break;

   case MIP3:  if (apnt) {
                  switch (pump) {
                      case IP3: apnt->mip3  = d1.val; break;
                  }
               }
               break;

   case HIP3:  if (apnt) {
                  switch (pump) {
                      case IP3: apnt->hip3  = d1.val; break;
                  }
               }
               break;

   }

/* ncfprintf (stderr,"ca chan %s %g\n",param->name,d1.val); /* */ 
}

/*------------------------------------------*/

void xchan(void)
{
  datum d1;
  Symbol *param;
  static chattrib *apnt; 
  int npops;

  param  = (Symbol *)*pc++;
  elpnt = makelem(CHAN,elpnt);		/* make dummy elem "chan"  */
  switch (param->type)
   {
    case GLU:
    case GABA:
    case GLY:
    case CHRC:
    case AMPA:
    case KAINATE:
    case NMDA:
    case CGMP:
    case SYN2:
    case K:
    case KCa:
    case ClCa:
    case NA:     npops=0; 
		 break;
    case TYPE:  
    case VVREV:  
    case TAU:
    case TAUM:
    case TAUH:
    case TAUN:
    case TAUA:
    case TAUB:
    case TAUC:
    case TAUD:
    case TAUE:
    case TAUF:
    case D1:
    case D2:
    case K1:
    case K2:
    case MG:
    case N:
    case UNIT:
    case THRESH:
    case OFFSET:
    case OFFSETM:
    case OFFSETH:
    case TRCONC:
    case MAXCOND:
    case CAPERM:
    case CAKD:	
    case CAHC:	
    case NDENSITY:
    case DENSITY:npops=1;
		 break;
   }
 
  if (npops==1) checknum(d1 = popm());

  switch (param->type)
   {
    case K:
    case KCa:
    case ClCa:
    case NA:				/* link attribs to elem */
    case SYN2:  
    case GLU:  
    case AMPA:  
    case KAINATE:  
    case NMDA:  
    case GABA:  
    case GLY:  
    case CHRC:  
    case CGMP:  apnt = (chattrib *)chanattr(elpnt,param->type); 
		break;
    case TYPE:  if (apnt) apnt->stype = (int)(d1.val);
		break;
    case VVREV:  if (apnt) apnt->vrev = d1.val;
		break;
    case TAU:  if (apnt) {
		  if (apnt->taua!= NULLVAL) apnt->taua *= d1.val;
		  else 			    apnt->taua  = d1.val;
		  if (apnt->taub!= NULLVAL) apnt->taub *= d1.val;
		  else 			    apnt->taub  = d1.val;
		  if (apnt->tauc!= NULLVAL) apnt->tauc *= d1.val;
		  else 			    apnt->tauc  = d1.val;
		  if (apnt->taud!= NULLVAL) apnt->taud *= d1.val;
		  else 			    apnt->taud  = d1.val;
		}
		break;
    case TAUN:
    case TAUM:  if (apnt) {
		  if (apnt->taua!= NULLVAL) apnt->taua *= d1.val;
		  else 			    apnt->taua  = d1.val;
		  if (apnt->taub!= NULLVAL) apnt->taub *= d1.val;
		  else 			    apnt->taub  = d1.val;
		}
		break;
    case TAUH:  if (apnt) {
		  if (apnt->tauc!= NULLVAL) apnt->tauc *= d1.val;
		  else 			    apnt->tauc  = d1.val;
		  if (apnt->taud!= NULLVAL) apnt->taud *= d1.val;
		  else 			    apnt->taud  = d1.val;
		}
		break;
    case TAUA:  if (apnt) {
		  if (apnt->taua!= NULLVAL) apnt->taua *= d1.val;
		  else 			    apnt->taua  = d1.val;
		}
		break;
    case TAUB:  if (apnt) {
		  if (apnt->taub!= NULLVAL) apnt->taub *= d1.val;
		  else 			    apnt->taub  = d1.val;
		}
		break;
    case TAUC:  if (apnt) {
		  if (apnt->tauc!= NULLVAL) apnt->tauc *= d1.val;
		  else 			    apnt->tauc  = d1.val;
		}
		break;
    case TAUD:  if (apnt) {
		  if (apnt->taud!= NULLVAL) apnt->taud *= d1.val;
		  else 			    apnt->taud  = d1.val;
		}
		break;
    case TAUE:  if (apnt) {
		  if (apnt->taue!= NULLVAL) apnt->taue *= d1.val;
		  else 			    apnt->taue  = d1.val;
		}
		break;
    case TAUF:  if (apnt) {
		  if (apnt->tauf!= NULLVAL) apnt->tauf *= d1.val;
		  else 			    apnt->tauf  = d1.val;
		}
		break;
    case TRCONC: if (apnt) apnt->trconc  = d1.val;
		break;
    case THRESH: 
    case OFFSET: if (apnt) {
		  apnt->voffsm  = d1.val;
		  apnt->voffsh  = d1.val;
		}
		break;
    case OFFSETM: if (apnt) apnt->voffsm  = d1.val;
		break;
    case OFFSETH: if (apnt) apnt->voffsh  = d1.val;
		break;
    case MAXCOND: if (apnt) apnt->maxcond = d1.val;
		break;
    case CAPERM: if (apnt) apnt->caperm = d1.val;
		break;
    case CAKD:  if (apnt) apnt->cakd = d1.val;
		break;
    case CAHC:  if (apnt) apnt->cahc = d1.val;
		break;
    case DENSITY: if (apnt) apnt->density = d1.val;
		break;
    case NDENSITY: if (apnt) apnt->ndensity = d1.val;
		break;
    case N: if (apnt) apnt->n = int(d1.val);
		break;
    case UNIT: if (apnt)apnt->unitary = d1.val;
		break;
    case D1: if (apnt) apnt->d1 = d1.val;
		break;
    case D2: if (apnt) apnt->d2 = d1.val;
		break;
    case K1: if (apnt) apnt->k1 = d1.val;
		break;
    case K2: if (apnt) apnt->k2 = d1.val;
		break;
   }

/* ncfprintf (stderr,"chan %s %g\n",param->name,d1.val); /* */ 
}

/*------------------------------------------*/

void initrgen(void)

{
  datum d1;
  Symbol *param;
  static int ngen=0;
  static int rsd=0;
  void initrand(int ngen, int rsd);

  param  = (Symbol *)*pc++;
  if (param==0) {
     checknum(d1 = popm()); 
     ngen = (int)d1.val;
  }
  else switch (param->type)
    {
     case RSD:  checknum(d1 = popm()); 
		rsd = (int)d1.val;
		break;
     case INITRAND: 
		initrand(ngen,rsd);
		ngen = rsd = 0;
		break; 
    }
}  

/*------------------------------------------*/

void noise(void)
{
  datum d1;
  Symbol *param;
  static nattrib *apnt=0; 
  int npops;
  double chnoise, vesnoise;

  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case VESNOISE:
    case CCHNOISE:
    case N:	 
    case RSD:	
    case COV:	
    case REFR:	
    case UNIT:	
    case TAUF:	
    case VSIZE:
    case VCOV:	  npops = 1; 
		 break;
   }

  if (npops==1) checknum(d1 = popm());

#ifndef XSTIMM
  switch (param->type)
   {	 					/* link attributes to elem */
    case VESNOISE:				 /* quantal vesicle noise */
		 vesnoise = d1.val;
		 if (vesnoise)
		   apnt = (nattrib *)chanattr(elpnt,param->type); 
		 else apnt = (nattrib *)NULL;
		 break;

    case CCHNOISE:				 /* quantal channel noise */
		 chnoise = d1.val;
		 if (chnoise)
		   apnt = (nattrib *)chanattr(elpnt,param->type);
		 else apnt = (nattrib *)NULL;
		 break;

    case N:	 if (apnt) apnt->n = (int)d1.val;
		 break;

    case COV:	 if (apnt) apnt->cov = d1.val;
		 break;

    case REFR:	 if (apnt) apnt->refr = d1.val;
		 break;

    case RSD:	 if (apnt) apnt->rseed = (int)d1.val;
		 break;

    case UNIT:	 if (apnt) apnt->unitary = d1.val;
		 break;

    case TAUF:	 if (apnt) apnt->tauf = d1.val;
		 break;

    case VSIZE:	if (apnt) {
		  if (apnt->ctype==VESNOISE)
			 apnt->vsize = d1.val;
		  else 
	          ncfprintf (stderr,"Snoise: vsize is illegal with chnoise.\n");
		}
		 break;

    case VCOV:	if (apnt) {
		  if (apnt->ctype==VESNOISE)
			 apnt->vcov = d1.val;
		  else 
	          ncfprintf (stderr,"Snoise: vcov is illegal with chnoise.\n");
		}
		 break;

   }
#endif
}

/*------------------------------------------*/

void xsynapse(void)
{
  datum d1;
  Symbol *param;
  int i,npops,nfilts;
  static int xfilt=0;
  synapse *epnt;
  static nattrib *napnt=(nattrib*)NULL;
#define TTIMEC 10
  double ttimec[TTIMEC];

  param  = (Symbol *)*pc++;
  epnt = (synapse *)makelem(SYNAPSE,elpnt);	/* make new synapse */
  if (elpnt->ctype==ELEMENT) napnt = (nattrib *)NULL; /* reset chan attr */
  elpnt = (elem *)epnt;

  switch (param->type)
   {
    case OPEN:  
    case CLOSE: 
    case RESP: 
    case MESGOUT: 
    case GLU: 
    case AMPA: 
    case KAINATE: 
    case NMDA: 
    case CNQX: 
    case GABA: 
    case BIC: 
    case PTX: 
    case GLY: 
    case STRY: 
    case CAMP: 
    case CGMP: 
    case SYN2: 
    case V: 
    case CA: 
		npops=0;
		break;
    case VVREV:  
    case THRESH:
    case NFILT1:		/* old syntax included for compatibility */
    case NFILT1H:
    case NFILT2:
    case NFILT3: 
    case HGAIN: 
    case OFFSET: 
    case TFALL2:
    case TFALL3:
    case CGAIN: 
    case CKD: 
    case CHC: 
    case EGAIN:
    case VGAIN:
    case MAXCOND:
    case MAXSRATE:
    case MRRPOOL:
    case RRPOOL:
    case RRPOOLG:
    case N:
    case UNIT:
    case TAUF:
    case TRCONC:
    case MESGCONC:
    case KD:	
    case HCOF:	
    case EXPON:
    case LINEAR:
    case DYAD:	
    case SPOST:	
    case TYPE:  npops = 1;
		break;
    case TIMEC1:
    case TIMEC1H:
    case TIMEC2:
    case TIMEC3: npops=0;
		nfilts = (long int)*pc++; /* # of tau's actually entered */
		if (nfilts >= TTIMEC) nfilts = TTIMEC - 1;
		for (i=0; i<TTIMEC; i++) {
		   ttimec[i] = 0;
		}
		for (i=nfilts-1; i>=0; i--) {
		   checknum(d1 = popm());
		   ttimec[i] = d1.val;
                }
		break;
   }

  if (npops==1) checknum(d1 = popm());

  switch (param->type) {
    case RESP:   synnum = cumelem;
		break;
    case N: 				/* make temporary attrib */
    case UNIT: 
    case TAUF:  if (!napnt) napnt=(nattrib*)chanattr(elpnt,NUNIT);
		break;
     default:   break;
  }
  switch (param->type)
   {
    case OPEN:  epnt->ntact = param->type; 
		epnt->secmsg = 0;
		break;
    case CLOSE: epnt->ntact = param->type; 
		epnt->secmsg = 1;
		break;
    case CA:  	
    case V:  	epnt->sens  = param->type;
		break;
    case VVREV:  epnt->vrev  = d1.val;
		break;
    case THRESH: epnt->thresh= d1.val;
		break;
    case HGAIN: epnt->filt1hg= d1.val;
		break;
    case OFFSET: epnt->filt1ho= d1.val;
		break;
    case NFILT1:  epnt->nfilt1 = (short int) d1.val; xfilt = 1; break;
    case NFILT1H: epnt->nfilt1h= (short int) d1.val; xfilt = 1; break;
    case NFILT2:  epnt->nfilt2 = (short int) d1.val; xfilt = 1; break;
    case NFILT3:  epnt->nfilt3 = (short int) d1.val; xfilt = 1; break;

    case TIMEC1:  if (!xfilt) epnt->nfilt1 = nfilts;
		  epnt->timec1=makfiltarr(epnt->nfilt1,0,epnt->timec1,0);
		  if (xfilt) {		/* if nfilts was given separately */
		     for (i=0; i<epnt->nfilt1; i++) 
                       epnt->timec1[i]= ttimec[0];
		  }
		  else {		/* if nfilts from num of tau's */
		     for (i=0; i<nfilts; i++) 
                       epnt->timec1[i]= ttimec[i];
		  }
		  xfilt = 0;
		break;
    case TIMEC1H:  if (!xfilt) epnt->nfilt1h = nfilts;
		  epnt->timec1h=makfiltarr(epnt->nfilt1h,0,epnt->timec1h,0);
		   if (xfilt) {		/* if nfilts was given separately */
		     for (i=0; i<epnt->nfilt1h; i++) 
                       epnt->timec1h[i]= ttimec[0];
		  }
		  else {		/* if nfilts from num of tau's */
		     for (i=0; i<nfilts; i++) 
                       epnt->timec1h[i]= ttimec[i];
		  }
		  xfilt = 0;
		break;
    case TIMEC2:  if (!xfilt) epnt->nfilt2 = nfilts;
		  epnt->timec2=makfiltarr(epnt->nfilt2,0,epnt->timec2,0);
		  if (xfilt) {		/* if nfilts was given separately */
		     for (i=0; i<epnt->nfilt2; i++) 
                       epnt->timec2[i]= ttimec[0];
		  }
		  else {		/* if nfilts from num of tau's */
		     for (i=0; i<nfilts; i++) 
                       epnt->timec2[i]= ttimec[i];
		  }
		  xfilt = 0;
		break;
    case TFALL2: epnt->tfall2 = d1.val;
		break;
    case TIMEC3:  if (!xfilt) epnt->nfilt3 = nfilts;
		  epnt->timec3=makfiltarr(epnt->nfilt3,0,epnt->timec3,0);
	          if (xfilt) {	/* if nfilts was given separately */
		     for (i=0; i<epnt->nfilt3; i++) 
                       epnt->timec3[i]= ttimec[0];
		  }
		  else {	/* if nfilts was given by num of tau's */
		     for (i=0; i<nfilts; i++) 
                       epnt->timec3[i]= ttimec[i];
		  }
		  xfilt = 0;
		break;
    case TFALL3: epnt->tfall3 = d1.val;
		break;
    case CGAIN: epnt->cgain   = d1.val;
		break;
    case COFF: epnt->coff   = d1.val;
		break;
    case CKD:   epnt->ckd     = d1.val;
		break;
    case CHC:   epnt->chc     = d1.val;
		break;
    case MAXCOND: epnt->maxcond = d1.val;
		break;
    case MAXSRATE: epnt->maxsrate = d1.val;
		break;
    case MRRPOOL: epnt->mrrpool = d1.val;
		break;
    case RRPOOL: epnt->rrpool = d1.val;
		break;
    case RRPOOLG: epnt->rrpoolg = d1.val;
		break;
    case TRCONC: epnt->trconc = d1.val;
		break;
    case MESGCONC: epnt->mesgconc = d1.val;
		break;
    case KD:	epnt->nkd   = d1.val;
		break;
    case EGAIN:	epnt->caegain = d1.val;
		break;
    case VGAIN:	epnt->vgain = d1.val;
		break;
    case HCOF:	epnt->npow  = d1.val;
		break;
    case LINEAR:if (epnt->curve == NULLVAL) {
		 epnt->curve = LINEAR;
    		 epnt->ngain = d1.val;
		}
		break;
    case EXPON:	if (epnt->curve == NULLVAL) {
		 epnt->curve = EXPON;
    		 epnt->ngain = d1.val;
		}
		break;
    case DYAD:	if (epnt->dyadelem <= 0) {
    		 epnt->dyadelem = int(d1.val);
		}
		break;
    
    case SPOST:	if (epnt->spost == NULLVAL) {
		 epnt->spost = int(d1.val);
		}
		break;
    
    case N:  	napnt->n = int(d1.val);	/* channel attrib made above */
		break;
    case UNIT:  napnt->unitary = d1.val;
		break;
    case TAUF:  napnt->tauf = d1.val;
		break;
    case GLU: 
    case AMPA: 
    case KAINATE: 
    case NMDA: 
    case CNQX: 
    case GABA: 
    case BIC: 
    case PTX: 
    case GLY: 
    case STRY: 
    		epnt->mesg1 = param->type; break;
    case CAMP: 
    case CGMP:  epnt->mesg2 = param->type; break;

    default:	break;
   }
/*ncfprintf(stderr,"synapse %s %g curv %g\n",
		param->name,d1.val,elpnt->curve); /* */ 
}

/*------------------------------------------*/

void xgj(void)
{
  datum d1,d[MAXNODIM];
  Symbol *param;
  gapjunc *epnt;
  static int revflag=0;
 
  param  = (Symbol *)*pc++;
  epnt = (gapjunc *)makelem(GJ,elpnt);	/* make new gap junction */
  elpnt = (elem *)epnt;			/* for membtyp() below */
  if (param==NULL) {			/* conductance */
	checknum(d1=popm());
	epnt->gmax = d1.val;
  }
  else {
   switch (param->type) {
    case GJ:
    case OPEN:
    case CLOSE:
    case MESGIN:
    case REV:
		break;
    case CAMP:
    case CGMP:	getnod(d);		/* get node, check dims */
	       break;
    case AREA:
    case DIA:
    case VGAIN:
    case OFFSET:
    case TAUN:
    case GMAX:
    case GNV: checknum(d1 = popm());
	      break;
   }

   switch (param->type) {
    case AREA: epnt->area  = d1.val;
		break;
    case DIA:  epnt->area  = MPI * d1.val * d1.val * 0.25;
		break;
    case OFFSET:if (revflag)  epnt->rvoff = d1.val;
    		else epnt->voff = d1.val;
		break;
    case VGAIN: if (revflag) epnt->rvgain = d1.val;
		else epnt->vgain  = d1.val;
		break;
    case TAUN: if (revflag) epnt->rtaun= d1.val;
    		else epnt->taun= d1.val;
		break;
    case REV:  revflag = 1;
		break;
    case GMAX: epnt->gmax= d1.val;
		break;
    case GNV: epnt->gnv= d1.val;
		break;
    case CAMP: 
    case CGMP:	epnt->nodeca = (int)d[0].val;       /* set node number */ 
		epnt->nodecd = (int)d[1].val;
		epnt->nodecd = (int)d[2].val;
		epnt->nodecd = (int)d[3].val;
		epnt->modtyp = param->type; 
		break;
    case OPEN: 
		epnt->sign = 1;
		break;
    case CLOSE:
		epnt->sign = 0;
		break;
    case GJ:
		revflag = 0;
		break;
   }
 }
}

/*------------------------------------------*/

void rload(void)
{
  datum d1;
  loadelem *epnt;

  checknum(d1 = popm());
  epnt = (loadelem *)makelem(LOAD,elpnt);	/* make new load */
  elpnt = (elem *)epnt;				/* for membtyp() below */
  epnt->r = d1.val;
}

void rcap(void)
{
  datum d1;
  capac *epnt;

  checknum(d1 = popm());
  epnt = (capac *)makelem(CAP,elpnt);		/* make new cap */
  epnt->c = d1.val;
}

void xgcap(void)
{
  datum d1;
  capac *epnt;
 
  checknum(d1 = popm());
  epnt = (capac *)makelem(GNDCAP,elpnt);	/* make new gnd cap */
  elpnt = (elem *)epnt;				/* for membtyp() below */
  epnt->c = d1.val;
}

void xresistor(void)
{
  datum d1;
  resistor *epnt;

  checknum(d1 = popm());
  epnt = (resistor *)makelem(RESISTOR,elpnt);	/* make new resistor */
  epnt->r = d1.val;
}

void xdiode(void)
{
  datum d1;
  diode *epnt;

  checknum(d1 = popm());
  epnt = (diode *)makelem(DIODE,elpnt);		/* make new diode */
  epnt->r = d1.val;
}

void rbatt(void)
{
  datum d1;
  batt *epnt;

  checknum(d1 = popm());
  epnt = (batt *)makelem(BATT,elpnt);		/* make new batt */
  epnt->v = d1.val;
}

void xgbatt(void)
{
  datum d1;
  batt *epnt;

  checknum(d1 = popm());
  epnt = (batt *)makelem(GNDBATT,elpnt);	/* make new gndbatt */
  epnt->v = d1.val;
}

void membtyp(void)
{
  datum d1;
  Symbol *param;
  int ctype;

  param  = (Symbol *)*pc++;
  checknum(d1 = popm());
  ctype = elpnt->ctype;
  switch (param->type)
   {
    case RM: switch (ctype) {
		case CABLE:  ((cable *)elpnt)->Rm = d1.val; break;
		case SPHERE: ((sphere *)elpnt)->Rm = d1.val; break;
	     } break;
    case RI: switch (ctype) {
		case CABLE:  ((cable *)elpnt)->Ri = d1.val; break;
	     } break;
    case RG: switch (ctype) {		/* for gj only */
		case GJ:  ((gapjunc *)elpnt)->specres = d1.val; break;
	     } break;
    case CM: switch (ctype) {
		case CABLE:  ((cable *)elpnt)->Cm = d1.val; break;
		case SPHERE: ((sphere *)elpnt)->Cm = d1.val; break;
	     } break;
    case JNOISE: elpnt->jnoise = d1.val; break;
    case RSD:    elpnt->rsd = int(d1.val); break;

    case VREST: switch (ctype) {
		case CABLE:  ((cable *)elpnt)->vrest = d1.val; break;
		case SPHERE: ((sphere *)elpnt)->vrest = d1.val; break;
		case GNDCAP:   ((capac *)elpnt)->vrest = d1.val; break;
		case LOAD:   ((loadelem *)elpnt)->vrest = d1.val; break;
	     } break;
    case VVREV: switch (ctype) {
		case CABLE:  ((cable *)elpnt)->vrev = d1.val; break;
		case SPHERE: ((sphere *)elpnt)->vrev = d1.val; break;
		case LOAD: ((loadelem *)elpnt)->vrev = d1.val; break;
	     } break;
    default: break;
   }
}

void xvbuf(void)
{
  vbuf *epnt;

  epnt = (vbuf *)makelem(BUF,elpnt);	/* make new vbuf */
  elpnt = (elem *)epnt;			/* for params below */
  epnt->delay = 0;			/* no delay */
  epnt->offset = 0;			/* no offset */
  epnt->gain = 1.0;			/* gain = vbuf */
}

void xvbufd(void)
{
   datum d1;

  checknum(d1 = popm());
  ((vbuf *)elpnt)->delay = d1.val;	/* delay in sec */
}

void xvbufo(void)
{
   datum d1;

  checknum(d1 = popm());
  ((vbuf *)elpnt)->offset = d1.val;	/* offset in volts */
}

void xvbufg(void)
{
   datum d1;

  checknum(d1 = popm());
  ((vbuf *)elpnt)->gain = d1.val;	/* voltage gain */
}

/*------------------------------------------*/

void xnbuf(void)
{
  nbuf *epnt;

  epnt = (nbuf *)makelem(NBUF,elpnt);	/* make new nbuf */
  elpnt = (elem *)epnt;			/* for params below */
  epnt->ntrans = 0;			/* null ntrans */
}

void xnbufo(void)
{
   datum d1;

  checknum(d1 = popm());
  ((nbuf *)elpnt)->offset = d1.val;	/* offset in volts */
}

void xnbuft(void)
{
   datum d1;

  checknum(d1 = popm());
  ((nbuf *)elpnt)->ntoffset = d1.val;	/* offset in nt */
}

void xnbufg(void)
{
   datum d1;

  checknum(d1 = popm());
  ((nbuf *)elpnt)->gain = d1.val;	/* voltage gain */
}

void xnbufn(void)
{
   datum d1;

  checknum(d1 = popm());
  ((nbuf *)elpnt)->ntrans = d1.val;	/* ntrans */
}

/*------------------------------------------*/

int smaxmin (double val)

{ 
   nodeint ival;
   int vali;

   static nodeint minnode=   (nodeint)(1L << (BITS(nodeint)-1));
   static nodeint maxnode= (nodeint)(~(1L << (BITS(nodeint)-1)));
   //static nodeint minnode= -32768;
   //static nodeint maxnode= 32767;


  vali = int(val);  
  if (vali > maxnode) vali = maxnode;
  if (vali < minnode) vali = minnode;
  ival = (nodeint) vali;

  return (ival);
}


/*------------------------------------------*/

void conn1(void)
{
  datum d[MAXNODIM];

  getnod(d);			/* get node, check dims */
  elpnt = tmpelem();		/* make new element */
  elpnt->node1a = (int)d[0].val;	/* set node number */ 
  elpnt->node1b = (int)d[1].val;
  elpnt->node1c = (int)d[2].val;
  elpnt->node1d = (int)d[3].val;
  elpnt->node2a = NULLVAL;
  elpnt->node2b = NULLVAL;
  elpnt->node2c = NULLVAL;
  elpnt->node2d = NULLVAL;
  eninstall (elpnt);            /* install new element in hash node table (req node1a,1b) */
  checknod(elpnt);		/* make node associated with element */
}

/*------------------------------------------*/

void conn1l(void)
{
  datum d1,d2,d3;
  int narg;
  node *npnt;

  narg = (long int)*pc++;		/* get loc off stack first */
  if (narg > 2) d3 = popm();
  else d3.val = LARGENODE;
  if (narg > 1) d2 = popm();
  else d2.val = LARGENODE;
  checknum(d1 = popm());
  conn1();

#ifndef XSTIMM
  checknod(elpnt);		/* make node associated with element */
  if (! (npnt=elpnt->nodp1)) {
	ncfprintf (stderr,"conn1L: can't find node %s\n",
	    prnode(elpnt->node1a,elpnt->node1b,elpnt->node1c,elpnt->node1d));
	return;
  }
  npnt->xloc = d1.val;
  npnt->yloc = d2.val;
  npnt->zloc = d3.val;
#endif
}

/*------------------------------------------*/

double dist3d(node *n1, node *n2)
{
   double dist,xdist,ydist,zdist;

  if (!n1 || !n2) return 0.0;
  if (n1->xloc >= LARGNOD) {
   ncfprintf(stderr,"\ndist3d: location of node 1: %s has not been defined.\n",
		 prnode(n1->nodenm1,n1->nodenm2,n1->nodenm3,n1->nodenm4));
   execerror ("Missing location: ","stopping...");
  }  
  if (n2->xloc >= LARGNOD) {
   ncfprintf(stderr,"\ndist3d: location of node 2: %s has not been defined.\n",
		 prnode(n2->nodenm1,n2->nodenm2,n2->nodenm3,n1->nodenm4));
   execerror ("Missing location: ","stopping...");
  }  
/*
  if (n1->xloc >= LARGNOD || n2->xloc >= LARGNOD)
    xdist = 0; 
  else
    xdist = n1->xloc - n2->xloc;
  if (n1->yloc >= LARGNOD || n2->yloc >= LARGNOD)
    ydist = 0; 
  else
    ydist = n1->yloc - n2->yloc;
*/
  if (n1->zloc >= LARGNOD || n2->zloc >= LARGNOD)
    zdist = 0; 

  xdist = n1->xloc - n2->xloc;
  ydist = n1->yloc - n2->yloc;
  zdist = n1->zloc - n2->zloc;

  dist = sqrt (xdist*xdist + ydist*ydist + zdist*zdist);
/* ncfprintf (stderr,"%g %g %g %g %g %.20g\n",
	n1->xloc,n1->yloc,n1->zloc,n2->xloc,n2->yloc,n2->zloc); /* */
  return (dist);
}


/*------------------------------------------*/

double dist2d(node *n1, node *n2)
{
   double dist,xdist,ydist;

  if (!n1 || !n2) return 0.0;

/*  if (n1->xloc >= LARGNOD) {
   ncfprintf(stderr,"\ndist2d: location of node %s has not been defined.\n", prnode(n1));
   execerror ("Missing location: ","stopping...");
  }  
  if (n2->xloc >= LARGNOD) {
   ncfprintf(stderr,"\ndist2d: location of node %s has not been defined.\n", prnode(n2));
   execerror ("Missing location: ","stopping...");
  }  
*/

/*
  if (n1->xloc >= LARGNOD || n2->xloc >= LARGNOD)
    xdist = 0; 
  else
    xdist = n1->xloc - n2->xloc;
  if (n1->yloc >= LARGNOD || n2->yloc >= LARGNOD)
    ydist = 0; 
 */
  xdist = n1->xloc - n2->xloc;
  ydist = n1->yloc - n2->yloc;
  dist = sqrt (xdist*xdist + ydist*ydist);
/* ncfprintf (stderr,"%g %g %g %g %g %.20g\n",
	n1->xloc,n1->yloc,n1->zloc,n2->xloc,n2->yloc,n2->zloc); /* */
  return (dist);
}

/*------------------------------------------*/

double distzd(node *n1, node *n2)
{
   double zdist;

  if (!n1 || !n2) return 0.0;
  if (n1->zloc >= LARGNOD) {
   ncfprintf(stderr,"\ndistzd: Z location of node %s has not been defined.\n",
		 prnode(n1->nodenm1,n1->nodenm2,n1->nodenm3,n1->nodenm4));
   execerror ("Missing location: ","stopping...");
  }  
  if (n2->zloc >= LARGNOD) {
   ncfprintf(stderr,"\ndist2d: Z location of node %s has not been defined.\n",
		 prnode(n2->nodenm1,n2->nodenm2,n2->nodenm3,n1->nodenm4));
   execerror ("Missing location: ","stopping...");
  }  
/*
  if (n1->zloc >= LARGNOD || n2->zloc >= LARGNOD)
    zdist = 0; 
  else
    zdist = n1->zloc - n2->zloc;
 */
  zdist = n1->zloc - n2->zloc;
/* ncfprintf (stderr,"%g %g %g %g %g %.20g\n",
	n1->xloc,n1->yloc,n1->zloc,n2->xloc,n2->yloc,n2->zloc); /* */
  return (zdist);
}

/*------------------------------------------*/

void conn1m(void)

/* "at node: elem offset=frac put node" */

/* make a new node along a cable element */

{
  datum d[MAXNODIM],d5,d6;
  int elnum,found;
  nodeint nod1a,nod1b,nod1c,nod1d;
  node *npnt,*npnt2,*npnt3;
  elem *epnt;
  attrib *apnt,*napnt;
  conlst *lpnt;
  double xdist,ydist,zdist,distfrac,len,ndia;

  getnod(d);			/* get node, check dims */
  nod1a=(int)d[0].val;
  nod1b=(int)d[1].val;
  nod1c=(int)d[2].val;
  nod1d=(int)d[3].val;

  checknum(d5 = popm());	/* get fractional distance */
  checknum(d6 = popm());	/* get cable element # rel to parent */

  conn1();		/* this gets parent node number off stack, */
			/*  makes new temp. element, and connects parent node */
			/*  at one end but does not set type of element, */
			/*  nor its loc */

#ifndef XSTIMM
  if (!(npnt=findnode(elpnt->node1a,elpnt->node1b,
			elpnt->node1c,elpnt->node1d,"Atmake"))) {
     ncfprintf (stderr,"conn1m: can't find node %s\n",
			 prnode(elpnt->node1a,elpnt->node1b,
				elpnt->node1c,elpnt->node1d));
    execerror ("Missing parent node: ","stopping...");
  }
				
  elnum = (int)(d6.val); 		/* find previously made element */
  if (!(epnt=findelem(elnum))) {
      ncfprintf (stderr,"Edist: can't find element %d\n",elnum);
      return;
  }
  else checkelemnode(epnt);

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

  distfrac = d5.val;
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

#endif     /* XSTIMM */
}

/*------------------------------------------*/

void conn2s(void)
{
  datum d[MAXNODIM];

  getnod(d);				/* get node, check dims */
  elpnt->node1a = (int)d[0].val;
  elpnt->node1b = (int)d[1].val;
  elpnt->node1c = (int)d[2].val;
  elpnt->node1d = (int)d[3].val;
  checknod(elpnt);		/* make nodes associated with element */
}

/*------------------------------------------*/

void conn2sl(void)
{
  datum d1,d2,d3;
  int narg;
  node *npnt;

  narg = (long int)*pc++;	/* get loc off stack first */
  if (narg > 2) d3 = popm();
  else d3.val = LARGENODE;
  if (narg > 1) d2 = popm();
  else d2.val = LARGENODE;
  checknum(d1 = popm());
  conn2s();

#ifndef XSTIMM
  checknod(elpnt);		/* make nodes associated with element */
  if (! (npnt=elpnt->nodp1)) {
	ncfprintf (stderr,"conn2s: can't find node %s\n",
	   prnode(elpnt->node1a,elpnt->node1b,elpnt->node1c,elpnt->node1d));
	return;
  }
  npnt->xloc = d1.val;
  npnt->yloc = d2.val;
  npnt->zloc = d3.val;
#endif
}

/*------------------------------------------*/

void conn2d(void)
{
  datum d[MAXNODIM];

  getnod(d);				/* get node, check dims */
  elpnt = tmpelem();			/* make new element */
  elpnt->node2a = (int)d[0].val;
  elpnt->node2b = (int)d[1].val;
  elpnt->node2c = (int)d[2].val;
  elpnt->node2d = (int)d[3].val;
  eninstall (elpnt);                   /* install new element in hash node table (req node1a,1b) */
}

/*------------------------------------------*/

void conn2dl(void)
{
  datum d1,d2,d3;
  int narg;
  node *npnt;

  narg = (long int)*pc++;
  if (narg > 2) d3 = popm();
  else d3.val = LARGENODE;
  if (narg > 1) d2 = popm();
  else d2.val = LARGENODE;
  checknum(d1 = popm());
  conn2d();
#ifndef XSTIMM
  checknod(elpnt);		/* make nodes associated with element */
  if (! (npnt=elpnt->nodp2)) {
	ncfprintf (stderr,"conn2dl: can't find node %s\n",
	    prnode(elpnt->node2a,elpnt->node2b,elpnt->node2c,elpnt->node2d));
	return;
  }
  npnt->xloc = d1.val;
  npnt->yloc = d2.val;
  npnt->zloc = d3.val;
#endif
}

/*------------------------------------------*/

void gettype(void)

/* Return numerical type for element. */

{
   datum d={0};
   Symbol *field;
 
  field = (Symbol *)*pc++;
  if (field) d.val= field->type;
  d.vtype = NUMBER;
  d.type = VAR;
  pushm (d);
}

/*------------------------------------------*/

double get_length (elem *epnt)
{
   conlst *lpnt;
   node *npnt;
   sphere *sepnt;
   cable *cepnt;
   double len, mindia1, mindia2;
   int found;

  checknod(epnt);

 found = 0;
 if ((npnt=epnt->nodp1)) {           /* check for spheres at node 1 */
   mindia1 = 1e10;
   for (lpnt=npnt->elemlst; lpnt; lpnt=lpnt->next) {
      if ((sepnt=(sphere *)lpnt->conpnt)) {
        if (sepnt->ctype==SPHERE) {
          found = 1;
          if (mindia1 > sepnt->dia) mindia1 = sepnt->dia;
        }
      }
   }  /* for (lpnt;;) */
   if (!found) mindia1 = 0.0;
 }  /* nodp1 */

 found = 0;
 if ((npnt=epnt->nodp2)) {           /* check for spheres at node 2 */
   mindia2 = 1e10;
   for (lpnt=npnt->elemlst; lpnt; lpnt=lpnt->next) {
      if ((sepnt=(sphere *)lpnt->conpnt)) {
        if (sepnt->ctype==SPHERE) {
          found = 1;
          if (mindia2 > sepnt->dia) mindia2 = sepnt->dia;
        }
      }
   }  /* for (lpnt;;) */
   if (!found) mindia2 = 0.0;
 }  /* nodp2 */

   cepnt=(cable*)epnt;
   if ((len=cepnt->length)==NULLVAL) {           /* if length not set */
     if (!epnt->nodp1 && !epnt->nodp2) len = 10.0;/* default length */
     else len = dist3d(epnt->nodp1,epnt->nodp2);   /* calc length */
     len -= (mindia1 + mindia2) * 0.5;
     if (len < 0.0) len = 1e-3;
   }
   else len = ((cable*)epnt)->length;
   return len;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

const char *get_elabel (elem *epnt, int field)

{
  const char *rval;

  switch (field) {
   case TYPE:
     rval = findsym(epnt->ctype);
     break;
   case ELABL:
     rval = epnt->elabl;
     break;
   default:
     rval = NULL;
     break;
  }
  return rval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

double get_efield (elem *epnt, int field)

{
  double rval,rm,ri,cm;
  static int err;

  if (epnt==NULL) return 0;
  rval = 0;
#ifndef XSTIMM
  switch (field) {
   case ELNUM:
	rval = epnt->elnum;
	break;
   case NTYPE:
	rval = epnt->ctype;
	break;
   case REGION:
	rval = epnt->region;
	break;
   case LENGTH:
	if (epnt->ctype != CABLE) {
	  rval = 0.0;
	  break;
	}
	if (((cable*)epnt)->length != NULLVAL) {
	  rval = ((cable*)epnt)->length;
	}
	else {			/* find length from location */
  	  rval = get_length(epnt);
	}
	break;
   case N3DIST:
	checknod(epnt);
	rval = dist3d(epnt->nodp1,epnt->nodp2);
	break;
   case N2DIST:
	checknod(epnt);
	rval = dist2d(epnt->nodp1,epnt->nodp2);
	break;
   case DIA:
	rval = ((cable *)epnt)->dia;
	break;
   case DIA2:
	rval = ((cable *)epnt)->dia2;
	break;
   case RM:
	rm = ((cable *)epnt)->Rm;
	if (rm==NULLVAL) rm = drm / exp(log(dqrm) * (tempcel-dbasetc)/10.0);
	rval = rm;
	break;
   case RI:
	ri = ((cable *)epnt)->Ri;
	if (ri==NULLVAL) ri = dri / exp(log(dqrm) * (tempcel-dbasetc)/10.0);
	rval = ri;
	break;
   case CM:
	cm = ((cable *)epnt)->Cm;
	if (cm==NULLVAL) cm = dcm;
	rval = cm;
	break;
   case CPLAM:
	rval = ((cable *)epnt)->cplam;	/* only for cable */
	break;
   case NODE1A:
	rval = epnt->node1a;
	break;
   case NODE1B:
	rval = epnt->node1b;
	break;
   case NODE1C:
	rval = epnt->node1c;
	break;
   case NODE1D:
	rval = epnt->node1d;
	break;
   case NODE2A:
	rval = epnt->node2a;
	break;
   case NODE2B:
	rval = epnt->node2b;
	break;
   case NODE2C:
	rval = epnt->node2c;
	break;
   case NODE2D:
	rval = epnt->node2d;
	break;
   case MODIFY:
	rval = epnt->modif;
	break;
   case MAXCOND:
	if (epnt->ctype == SYNAPSE) {
	    synap *spnt;

	  if (spnt=(synap *)epnt->lptr) {
	    rval = spnt->maxcond;
	  } else 
	    rval = ((synapse*)epnt)->maxcond;
	}
	else if (epnt->ctype == CHAN) {
	     attrib *apnt;
	     chan *chpnt;
	  if (chpnt=(chan *)epnt->lptr) {
	    rval = chpnt->maxcond;
	  } else {
	    if (!(apnt=epnt->attpnt)) {
	      if (!err) {
	        ncfprintf (stderr,"Chanfield: can't find elem %d\n",epnt->elnum);
	        execerror ("Missing element; ","stopping...");
		err=1;
	      }
	      return (0);
	    }
	    rval = ((mattrib *)apnt)->maxcond;
	  }
	}
	break;
   case SCURVE:
	if (epnt->ctype == SYNAPSE) {
	   synap *spnt;
          if (spnt=(synap *)epnt->lptr)    /* get ptr to synapse from elem */
	    rval = spnt->curve;
	  else
	    rval = ((synapse *)epnt)->curve;
	 }
	break;
   case SDYAD:				/* return elem number of original */
	if (epnt->ctype == SYNAPSE)
	  rval = ((synapse *)epnt)->dyadelem; /* elem # of dyad (orig. synap) */
	  if (rval < 0) rval = 0;
	break;
   case DYAD:				/* return elem # of first dyad */
	if (epnt->ctype == SYNAPSE) {
	   synap *spnt;
	 if ((spnt=(synap *)epnt->lptr)) {   /* get ptr to synapse from elem */
	  /* if (spnt->dyad) 
	     rval = spnt->dyad->sdyad; */
		/* need to scan all synapse elements, find one that */
		/*  has lptr same as sdyad above, then find elem number */
	  }
	 }
	break;
   case NUMDYAD:			/* return number of dyads */
			/* Note: this only works well after "step" or "run" */
			/*  so that low-lev synap structures are set up. */
			/*  It is very inefficient searching through elems. */

	if (epnt->ctype == SYNAPSE) {
	   synap *spnt;
	 if ((spnt=(synap *)epnt->lptr))    /* get ptr to synapse from elem */
	   rval = numdyad(spnt->spost);
	 else {				  /* else check all elements */
	    int n;
	    elem *espnt;

	   for (n=0,espnt=elempnt; espnt; espnt=espnt->next) { 
	     if (espnt->ctype == SYNAPSE && ((synapse *)espnt)->dyadelem>0)
		if (((synapse*)espnt)->ngain==epnt->elnum) n++;
	   }
	   rval = n;
	  }
	 }
	break;
  
   case NUMSPRE:	/* return number of presynaptic synaptic elements */
			/* Note: this only works well after "step" or "run" */
			/*  so that low-lev synap structures are set up. */
			/*  It is very inefficient searching through elems. */

	if (epnt->ctype == SYNAPSE) {
	   synap *spnt;
	 if ((spnt=(synap *)epnt->lptr))    /* get ptr to synapse from elem */
	   rval = numdyad(spnt->spre);
	 else {				  /* else check all elements */
	    int n;
	    elem *espnt;

	   for (n=0,espnt=elempnt; espnt; espnt=espnt->next) { 
	     if (espnt->ctype == SYNAPSE && 
		((synapse *)espnt)->spost!=NULLVAL)

		if (((synapse*)espnt)->spost==epnt->elnum) n++;
	   }
	   rval = n;
	  }
	 }
	break;
  }
#endif
  return rval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

double get_efield (int elnum, int field)

{
  elem *epnt;
  static int err=0;

 if (!(epnt=findelem(elnum))) {
      if (!err) {
        ncfprintf (stderr,"Elemfield: can't find elem %d\n",elnum);
        execerror ("Missing element; ","stopping...");
        err=1;
      }
      return 0;
 }
 return (get_efield(epnt,field));
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

attrib *get_chanattr (elem *epnt)

/* return a channel attrib if found */

{
  int found;
  attrib *apnt;

  if (epnt==NULL) return NULL;
  if (epnt->ctype!=CHAN && epnt->ctype!=SYNAPSE) return NULL;
  for (found=0,apnt=epnt->attpnt; apnt; apnt=apnt->attpnt) {
    switch (apnt->ctype) {

	case GLU:
	case GABA:
	case GLY:
	case CHRC:
	case NA:
	case  K:
	case  KCa:
	case  ClCa:
	case AMPA:
	case KAINATE:
	case NMDA:
	case CGMP:
	case SYN2:
	case CACOMP:
	case CA:   
		found = 1;
     		break;
	case VESNOISE:
	case CCHNOISE:
	case NUNIT:
		break;
    }
    if (found) break;
  }
  if (found) return (apnt);
  else	    return (NULL);
} 

/* - - - - - - - - - - - - - - - - - - - - - - - - */

node *get_elemnode1 (elem *epnt)

{
 if (!(epnt)) return NULL;
 if (epnt->nodp1==NULL) checkelemnode(epnt);
 return (epnt->nodp1);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

node *get_elemnode2 (elem *epnt)

{
 if (!(epnt)) return NULL;
 if (epnt->nodp2==NULL) checkelemnode(epnt);
 return (epnt->nodp2);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

void efield(void)

/* push elem field on stack */

{
    datum d1={0},d2={0};
    Symbol *field;
    elem *epnt;
    int elnum;
    static int err=0;
    double ri, rm, cm;

 field = (Symbol *)*pc++;
 checknum(d1 = popm());
 d2.val = 0;
 d2.type = VAR;
 d2.vtype = NUMBER;
#ifndef XSTIMM
 elnum = (int)(d1.val);
 if (!(epnt=findelem(elnum))) {
      if (!err) {
        ncfprintf (stderr,"Elemfield: can't find elem %d\n",elnum);
        //execerror ("Missing element; ","stopping...");
	err = 1;
      }
      d2.val = 0;
 }
 else {
  switch (field->type) {
   case TYPE:
	d2.str = get_elabel(epnt,field->type);
 	d2.vtype = STRING;
	break;
   case ELABL:
	d2.str = get_elabel(epnt,field->type);
 	d2.vtype = STRING;
	break;
   default:
	d2.val = get_efield(epnt,field->type); 
	break;
  }
 }
#endif
 d2.type = VAR;
 pushm(d2);
}

/*------------------------------------------*/

int get_chan_nstate(int elnum)

{
   attrib *apnt;
   chantype *chtyp;
   int retval;

 if (!(apnt=getattrib(elnum))) {
      fprintf (stderr,"Get_chan_numstate: can't find elem %d\n",elnum);
      return(0);
  }
  chtyp = getchantype(apnt->ctype,apnt->stype);
  if (chtyp) retval = chtyp->numstate;
  else       retval = 0;
 return retval;
}


/*------------------------------------------*/

void cfield(void)

/* push chan field on stack */

{
    datum d1={0},d2={0};
    Symbol *field;
    attrib *apnt;
    int elnum;
    static int err=0;
    chantype *chtyp;
    chan *chpnt;
    elem *epnt;
   
 field = (Symbol *)*pc++;
 checknum(d1 = popm());
 d2.val = 0;
 d2.type = VAR;
 d2.vtype = NUMBER;
#ifndef XSTIMM
 elnum = (int)(d1.val);
 if (!(apnt=getattrib(elnum))) {
      ncfprintf (stderr,"Chanfield: can't find elem %d\n",elnum);
      execerror ("Missing element; ","stopping...");
      return;
 }
 else {
  switch (field->type) {
   case TYPE:
	d2.str = findsym(apnt->ctype);
 	d2.vtype = STRING;
	break;
   case NTYPE:
	d2.val = apnt->ctype;
	break;
   case STYPE:
	d2.val = apnt->stype;
	break;
   case NSTATE:
        chtyp = getchantype(apnt->ctype,apnt->stype);
	if (chtyp) d2.val = chtyp->numstate;
        else       d2.val = 0;
	break;
   case MAXCOND:
        if (!(epnt=findelem(elnum))) {
	    if (!err) {
              ncfprintf (stderr,"Elemfield: can't find elem %d\n",elnum);
              //execerror ("Missing element; ","stopping...");
	      err=1;
	    }
            return;
        }
	if (chpnt= (chan *)epnt->lptr) {
	    d2.val = chpnt->maxcond;
	  } else {
	    d2.val = ((mattrib *)apnt)->maxcond;
	}
	break;
  }
 }
#endif
 d2.type = VAR;
 pushm(d2);
}

/*------------------------------------------*/

double get_nfield(node *npnt, int field, int elnum)

{
    int i,found;
    double rval;
    conlst *lpnt=(conlst*)NULL;
    elem *epnt=(elem*)NULL;
    attrib *apnt;

  switch (field) {
   case ELNUM:
      for (lpnt=npnt->elemlst,i=0; lpnt && i<elnum; i++)	/* find the element */
           lpnt=lpnt->next;
      if (lpnt) epnt=(elem *)lpnt->conpnt;
      if (!lpnt || !epnt) {
        ncfprintf (stderr,"Nodefield: can't find elem %d relative to node %s\n",
				elnum,prnode(npnt));
        execerror ("Missing element; ","stopping..."); /* */
        rval = -1;			
      }
      else rval = epnt->elnum;				/* element's number */
	break;
   case EXIST:
	rval = 1;
	break;
   case NUMCONN:
        for (lpnt=npnt->elemlst,i=0; lpnt; i++)		/* how many connections */
          lpnt=lpnt->next; 
	rval = i;
	break;
   case NUMSYN:
        for (lpnt=npnt->elemlst,i=0; lpnt; i++) {	/* how many synapses */
	  if (lpnt) epnt=(elem *)lpnt->conpnt;
          lpnt=lpnt->next; 
	  if (epnt) { if (epnt->ctype!=SYNAPSE) continue; }
	}
	rval = i;
	break;
   case SYNAPSE:			/* find synapse elem rel to node */
   	for (lpnt=npnt->elemlst,i=0; lpnt && i<elnum; ) { /* find the element */
	  if (lpnt) epnt=(elem *)lpnt->conpnt;
	  if (epnt) { if (epnt->ctype==SYNAPSE) i++; }
          lpnt=lpnt->next;
   	}
   	if (i != elnum) {
     ncfprintf (stderr,"Nodefield: can't find synapse %d relative to node %s\n",
				elnum,prnode(npnt));
	execerror ("Missing element; ","stopping..."); /* */
	rval = -1;			
	}
   	else rval = epnt->elnum;			/* element's number */
   break;
 
   case XLOC:
	rval = npnt->xloc;
	break;
   case YLOC:
	rval = npnt->yloc;
	break;
   case ZLOC:
	rval = npnt->zloc;
	break;
   case CACOMP:
        rval = 0;
	if (runyet) {
	  if (npnt->comptr) {
	     if (npnt->comptr->capnt) { rval = 1; }
	  }
	}
	else {		/* look through node's elements for a Ca channel */
	  for (lpnt=npnt->elemlst; lpnt; lpnt=lpnt->next) {
	    if (epnt=(elem*)lpnt->conpnt) {
	      switch (epnt->ctype) {
		 case CACOMP: found = 1; break;
		 case CABLE:
		 case SPHERE:
		 case CHAN:
		   for (found=0,apnt=epnt->attpnt; apnt; apnt=apnt->attpnt) {
		      switch (apnt->ctype) {
			case NMDA:
			case AMPA:
			case SYN2:
			case CGMP:
			case CA:
			     if (((mattrib *)apnt)->maxcond > 0 ||
			         ((mattrib *)apnt)->density > 0 ||
			         ((mattrib *)apnt)->ndensity > 0) found=1;
			     break;
			 default:
			     break;
		      }
		      if (found) break;
		   }
		  break;
		default:
		  break;
	      }
	    }
	    if (found) break;
	  }
	  if (found) rval = 1;
	  else       rval = 0;
	}
	break;
  }
  return rval;
}

/*------------------------------------------*/

double get_nfield(node *npnt, int field)

{
  return get_nfield(npnt, field, 0);
}

/*------------------------------------------*/

void nfield(void)

/* push node field on stack */

{
    datum d[MAXNODIM]={0},d5={0},d6={0};
    Symbol *field;
    int i, elnum, nod1, nod2, nod3, nod4;
    int found;
    node *npnt;

 field = (Symbol *)*pc++;			/* get rel elem number */
 if (field == (Symbol *)0 || field == (Symbol*)SYNAPSE) {
    d5 = popm(); 
    elnum = (int)d5.val;
 }
 else elnum = 0;
 getnod(d);					/* get node, check dims */

#ifndef XSTIMM
 nod1 = (int)d[0].val;
 nod2 = (int)d[1].val;
 nod3 = (int)d[2].val;
 nod4 = (int)d[3].val;
 if (!(npnt=findnode(nod1,nod2,nod3,nod4,NULL))) {
   if (!field || (field->type != EXIST)) {
     execerror ("Nfield: nissing node; ","stopping...");
     return; 
   }
   else d6.val = 0;		/* node -> exist = 0 */
 }
 else if (field == (Symbol *)0) {			/* "node -> elemnum" */
   d6.val = get_nfield(npnt,ELNUM,elnum);
 }
 else {
   d6.val = get_nfield(npnt,field->type,elnum);
 }
#endif
 d6.vtype = NUMBER;
 d6.type = VAR;
 pushm(d6);
}

/*------------------------------------------*/

double mindist (double dist1, double dist2, double dist3)

/* find minimum distance between a point and
   a line segment. If the point is off the end
   of the line segment, return the distance to
   that end.
*/

{
    double mdist,dist12,dist22,dist32,tmp1,tmp;

				/* check to see if node is off the end: */
   if (dist1 > dist2 && (dist1*dist1 - dist2*dist2) > dist3*dist3)
		 mdist = dist2;
   else
   if (dist1 < dist2 && (dist2*dist2 - dist1*dist1) > dist3*dist3)
		 mdist = dist1;

				/* formula for distance between a point
				/* and a line segment in three dimensions: */
   else {
      if (dist3==0.0) dist3 = 1e-20;
      dist12 = dist1 * dist1;
      dist22 = dist2 * dist2;
      dist32 = dist3 * dist3;
      tmp1 = (dist32 + dist12 - dist22) / (2 * dist3);
      if ((tmp=dist12-tmp1*tmp1) <= 0) mdist = 0;
      else mdist = sqrt (tmp);
   }
  return mdist;
}

/*------------------------------------------*/

void eedist(void)

/* return 3D distance between 2 elements */

{
    datum d1,d2;
    int elnum1,elnum2;
    int inside1, inside2;
    double dist,dist1,dist2,dist3,distc;
    node *np11,*np12,*np21,*np22;
    elem *epnt1,*epnt2;
    vec p11, p12, p21, p22, p1, p2, pt1, pt2;
    double t,s;

 checknum(d1 = popm());
 checknum(d2 = popm());
 elnum1 = (int)(d1.val);
 elnum2 = (int)(d2.val);
 dist = 0;

#ifndef XSTIMM

 if (!(epnt1=findelem(elnum1))) {
      ncfprintf (stderr,"E3dist: can't find element %d\n",elnum1);
      execerror ("Missing element: "," stopping... ");
      return;  
 }
 if (!(epnt2=findelem(elnum2))) {
      ncfprintf (stderr,"E3dist: can't find element %d\n",elnum2);
      execerror ("Missing element: "," stopping... ");
      return;  
 }
 np11 = epnt1->nodp1;
 np12 = epnt1->nodp2;
 np21 = epnt2->nodp1;
 np22 = epnt2->nodp2;
 if (!np11) {
       ncfprintf (stderr,"E3dist: can't find first node for element %d\n",elnum1);
       return;  
 }
 if (!np21) {
       ncfprintf (stderr,"E3dist: can't find first node for element %d\n",elnum2);
       return;  
 }

 p11.x = np11->xloc;
 p11.y = np11->yloc;
 p11.z = np11->zloc;
 if (np11->xloc >= LARGNOD) p11.x = 0;
 if (np11->yloc >= LARGNOD) p11.y = 0;
 if (np11->zloc >= LARGNOD) p11.z = 0;

 p12.x = np12->xloc;
 p12.y = np12->yloc;
 p12.z = np12->zloc;
 if (np12->xloc >= LARGNOD) p12.x = 0;
 if (np12->yloc >= LARGNOD) p12.y = 0;
 if (np12->zloc >= LARGNOD) p12.z = 0;

 p21.x = np21->xloc;
 p21.y = np21->yloc;
 p21.z = np21->zloc;
 if (np21->xloc >= LARGNOD) p21.x = 0;
 if (np21->yloc >= LARGNOD) p21.y = 0;
 if (np21->zloc >= LARGNOD) p21.z = 0;

 p22.x = np22->xloc;
 p22.y = np22->yloc;
 p22.z = np22->zloc;
 if (np22->xloc >= LARGNOD) p22.x = 0;
 if (np22->yloc >= LARGNOD) p22.y = 0;
 if (np22->zloc >= LARGNOD) p22.z = 0;

 switch (epnt1->ctype) {

  case CABLE: 
 
    if (!np12) {
      ncfprintf (stderr,"E2dist: can't find second node for element %d\n",elnum1);
      return;  
    }
   break;

  case SPHERE: 
  default:				/* element is not cable */
    
    p12.x = p11.x;
    p12.y = p11.y;
    p12.z = p11.z;
    break;

  }

 switch (epnt2->ctype) {

  case CABLE: 
 
    if (!np22) {
      ncfprintf (stderr,"E2dist: can't find second node for element %d\n",elnum2);
      return;  
    }
   break;

  case SPHERE: 
  default:				/* element is not cable */
    
    p22.x = p21.x;
    p22.y = p21.y;
    p22.z = p21.z;
    break;

  }
		/* find closest points */

 distc = norm_line (p11,vsub(p12,p11),p21,vsub(p22,p21),&t,&s,&p1,&p2);

 if (distc < 0) {		/* if line segs are parallel */
    inside1 = inside2 = 0;
    if (vdist (p11,p12) > vdist(p21,p22)) {     /* use line 1 for dist3 */
	dist3 = vlen(vsub(p12,p11));
	if ((vdist (p21,p11)+vdist(p21,p12))<(vdist (p22,p11)+vdist(p22,p12))) {
          dist1 = vlen(vsub(p21,p11));
          dist2 = vlen(vsub(p21,p12));
	}
	else {
          dist1 = vlen(vsub(p22,p11));
          dist2 = vlen(vsub(p22,p12));
	}
    }
    else {					/* use line 2 for dist3 */
        dist3 = vlen(vsub(p22,p21));
	if ((vdist (p11,p21)+vdist(p11,p22))<(vdist (p12,p21)+vdist(p12,p22))) {
          dist1 = vlen(vsub(p11,p21));
          dist2 = vlen(vsub(p11,p22));
	} else {
          dist1 = vlen(vsub(p12,p21));
          dist2 = vlen(vsub(p12,p22));
	}
    }
  /* ncfprintf (stderr,"dist1 %g dist2 %g dist3 %g\n",dist1,dist2,dist3); */

    dist = mindist (dist1,dist2,dist3) - ((cable*)epnt1)->dia*0.5 - 
					 ((cable*)epnt2)->dia*0.5;
 }

 else {			/* line segs not parallel */
  inside1 = (0 <= t && t <= 1);	     /* point of intersection is inside seg 1 */
  inside2 = (0 <= s && s <= 1);      /* point of intersection is inside seg 2 */

  if (inside1) {
    if (inside2) {				/* if inside both */
       dist = distc - ((cable*)epnt1)->dia*0.5 - ((cable*)epnt2)->dia*0.5;
    }
    else {					/* if inside 1 but not 2 */
     if (vdist(p2,p21) < vdist(p2,p22)) {	/* 21 is closer */
        
        dist1 = vlen(vsub(p21,p11));
        dist2 = vlen(vsub(p21,p12));
     }
     else {					/* 22 is closer */
        dist1 = vlen(vsub(p22,p11));
        dist2 = vlen(vsub(p22,p12));
     }
     dist3 = vlen(vsub(p11,p12));
     dist = mindist (dist1,dist2,dist3) - ((cable*)epnt1)->dia*0.5;

    }
  }
  else if (inside2){				/* if inside2 but not 1 */
     if (vdist(p1,p11) < vdist(p1,p12)) {	/* 11 is closer */
        
        dist1 = vlen(vsub(p11,p21));
        dist2 = vlen(vsub(p11,p22));
     }
     else {					/* 12 is closer */
        dist1 = vlen(vsub(p12,p21));
        dist2 = vlen(vsub(p12,p22));
     }
     dist3 = vlen(vsub(p21,p22));
     dist = mindist (dist1,dist2,dist3) - ((cable*)epnt2)->dia*0.5;
  }

  else {					/* outside both */
     if (vdist(p1,p11) < vdist(p1,p12)) pt1 = p11;
     else                               pt1 = p12;
     if (vdist(p2,p21) < vdist(p2,p22)) pt2 = p21;
     else                               pt2 = p22;
     if (vdist(p1,pt1) < vdist(p2,pt2)) {
        dist1 = vlen(vsub(pt2,p11));
        dist2 = vlen(vsub(pt2,p12));
        dist3 = vlen(vsub(p11,p12));
        dist = mindist (dist1,dist2,dist3) - ((cable*)epnt1)->dia*0.5;
     }
     else {
        dist1 = vlen(vsub(pt1,p21));
        dist2 = vlen(vsub(pt1,p22));
        dist3 = vlen(vsub(p21,p22));
        dist = mindist (dist1,dist2,dist3) - ((cable*)epnt2)->dia*0.5;
     }
  }
 }
#endif
  dist = max(dist,0);			/* no negative distances */
  d1.val = dist;
  d1.vtype = NUMBER;
  d1.type = VAR;
  pushm(d1);
}

/*------------------------------------------*/

double endist3d(elem *epnt, node *npnt)

/* return the distance between an element and a node */

{
  node *np1,*np2;
  double dist,dist1,dist2,dist3;

 dist = 0;
 np1 = epnt->nodp1;
 np2 = epnt->nodp2;
 if (!np1) {
     //ncfprintf (stderr,"Endist3d: can't find node for element %d\n",epnt->elnum);
     return (LARGENUM);  
 }

 dist1 = dist3d (npnt,np1);

 switch (epnt->ctype) {

#define max(x, y)	(((x) < (y)) ? (y) : (x))

  case CABLE: 
 
    if (!np2) {
       //ncfprintf (stderr,"E2dist: can't find node for element %d\n",epnt->elnum);
       return(LARGENUM);  
    }

   dist2 = dist3d (npnt,np2);
   dist3 = dist3d (np1,np2);
   dist = mindist (dist1,dist2,dist3);
   dist = max(dist,0);			/* no negative distances */
   break;

  case SPHERE: 
    dist = dist1 - ((sphere*)epnt)->dia*0.5;
    dist = max(dist,0);			/* no negative distances */
    break;

  default:				/* element is not cable */
    dist = dist1;
    break;
  }
  return dist;
}

/*------------------------------------------*/

double endist2d(elem *epnt, node *npnt)

/* return the distance between an element and a node */
{
  node *np1,*np2;
  double dist,dist1,dist2,dist3;

 np1 = epnt->nodp1;
 np2 = epnt->nodp2;
 if (!np1) {
     ncfprintf (stderr,"Endist2d: can't find node for element %d\n",epnt->elnum);
     return (0);  
 }

 dist1 = dist2d (npnt,np1);

 switch (epnt->ctype) {

  case CABLE: 
 
   if (!np2) {
      ncfprintf (stderr,"E2dist: can't find node for element %d\n",epnt->elnum);
      return(0);  
   }

   dist2 = dist2d (npnt,np2);
   dist3 = dist2d (np1,np2);
   dist = mindist (dist1,dist2,dist3);
   dist = max(dist,0);			/* no negative distances */
   break;

  case SPHERE: 
    dist = dist1 - ((sphere*)epnt)->dia*0.5;
    dist = max(dist,0);			/* no negative distances */
    break;

  default:				/* element is not cable */
    dist = dist1;
    break;
  }
  return dist;
}

/*------------------------------------------*/

double endistzd(elem *epnt, node *npnt)

/* Return the Z distance between an element and a node. */
/*  Calculate 2D and 3D min distances, then solve for Z */

{
  node *np1,*np2;
  double dist,dist1,dist2,dist3;
  double distx,dist1x,dist2x,dist3x;
  double zdist, zdistsq;

 np1 = epnt->nodp1;
 np2 = epnt->nodp2;
 if (!np1) {
     ncfprintf (stderr,"Endistzd: can't find node for element %d\n",epnt->elnum);
     return (0);  
 }

 dist1  = dist2d (npnt,np1);
 dist1x = dist3d (npnt,np1);

 switch (epnt->ctype) {

  case CABLE: 
 
   if (!np2) {
      ncfprintf (stderr,"E2dist: can't find node for element %d\n",epnt->elnum);
      return(0);  
   }

   dist2  = dist2d (npnt,np2);
   dist3 = dist2d (np1,np2);
   dist = mindist (dist1,dist2,dist3);
   dist = max(dist,0);			/* no negative distances */

   dist2x = dist3d (npnt,np2);
   dist3x = dist3d (np1,np2);
   distx = mindist (dist1x,dist2x,dist3x);
   distx = max(distx,0);			/* no negative distances */
   break;

  case SPHERE: 
    dist = dist1 - ((sphere*)epnt)->dia*0.5;
    dist = max(dist,0);			/* no negative distances */
    distx = dist1x - ((sphere*)epnt)->dia*0.5;
    distx = max(distx,0);			/* no negative distances */
    break;

  default:				/* element is not cable */
    dist = dist1;
    distx = dist1x;
    break;
  }
  zdistsq = distx*distx - dist*dist;
  if (zdistsq < 0) zdistsq = 0;
  zdist = sqrt(zdistsq);
  return zdist;
}

/*------------------------------------------*/

elem *get_elempnt(int elnum);

double endist2d(int elnum, node *npnt)

{
  elem *epnt;

  if ((epnt=get_elempnt(elnum))==NULL) return 0;
  else return endist2d(epnt,npnt);
}

/*------------------------------------------*/

double endist3d(int elnum, node *npnt)

{
  elem *epnt;

  if ((epnt=get_elempnt(elnum))==NULL) return 0;
  else return endist3d(epnt,npnt);
}

/*------------------------------------------*/

double endistzd(int elnum, node *npnt)

{
  elem *epnt;

  if ((epnt=get_elempnt(elnum))==NULL) return 0;
  else return endistzd(epnt,npnt);
}

/*------------------------------------------*/

void e3dist(void)

/* return 3D distance between a node and an element */

{
    datum d1={0},d2={0},d[MAXNODIM]={0};
    int nod1a, nod1b, nod1c, nod1d, elnum;
    node *npnt;
    elem *epnt;

 checknum(d2 = popm());
 elnum = (int)(d2.val);
 getnod(d);					/* get node, check dims */

#ifndef XSTIMM
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;

 if (!(npnt=findnode(nod1a,nod1b,nod1c,nod1d,"E3dist"))) return; 

 if (!(epnt=findelem(elnum))) {
      ncfprintf (stderr,"E3dist: can't find element %d\n",elnum);
      execerror ("Missing element: "," stopping... ");
      return;  
 }
#endif
  d1.val = endist3d(epnt,npnt);
  d1.vtype = NUMBER;
  d1.type = VAR;
  pushm(d1);
}

/*------------------------------------------*/

void e2dist(void)

/* return 2D distance between a node and an element */

{
    datum d1={0},d2={0},d[MAXNODIM]={0};
    int nod1a, nod1b, nod1c, nod1d, elnum;
    node *npnt;
    elem *epnt;

 checknum(d2 = popm());
 elnum = (int)(d2.val);
 getnod(d);					/* get node, check dims */

#ifndef XSTIMM
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;

 if (!(npnt=findnode(nod1a,nod1b,nod1c,nod1d,"E2dist"))) return; 

 if (!(epnt=findelem(elnum))) {
      ncfprintf (stderr,"E2dist: can't find element %d\n",elnum);
      execerror ("Missing element: "," stopping... ");
      return;  
 }

#endif
  d1.val = endist2d(epnt,npnt);
  d1.vtype = NUMBER;
  d1.type = VAR;
  pushm(d1);
}

/*------------------------------------------*/

void ezdist(void)

/* return 1D Z distance between a node and an element */

{
    datum d1={0},d2={0},d[MAXNODIM]={0};
    int nod1a, nod1b, nod1c, nod1d, elnum;
    node *npnt;
    elem *epnt;

 checknum(d2 = popm());
 elnum = (int)(d2.val);
 getnod(d);					/* get node, check dims */

#ifndef XSTIMM
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;

 if (!(npnt=findnode(nod1a,nod1b,nod1c,nod1d,"E2dist"))) return; 

 if (!(epnt=findelem(elnum))) {
      ncfprintf (stderr,"Ezdist: can't find element %d\n",elnum);
      execerror ("Missing element: "," stopping... ");
      return;  
 }

#endif
  d1.val = endistzd(epnt,npnt);
  d1.vtype = NUMBER;
  d1.type = VAR;
  pushm(d1);
}


/*------------------------------------------*/

double enfrac(elem *epnt, node *npnt)

/* return point in an element that is the closest to a node */

{
    double dist, dist1, dist2, dist3, frac;
    node *np1,*np2;

 if (epnt==NULL) return 0;

 np1 = epnt->nodp1;
 np2 = epnt->nodp2;
 
 dist1 = dist3d (npnt,np1);

 if (!np1) {
         ncfprintf (stderr,"Fdist: can't find node for element %d\n",epnt->elnum);
         return 0;  
 }
 if (epnt->ctype==CABLE) {
    if (!np2) {
         ncfprintf (stderr,"Edist: can't find node for element %d\n",epnt->elnum);
         return 0;  
    }

    dist2 = dist3d (npnt,np2);
    dist3 = dist3d (np1,np2);
    if (dist3 == 0.0) dist3 = 1e-20;
    dist = mindist (dist1,dist2,dist3);
    if      (dist == dist1) frac = 0.0;
    else if (dist == dist2) frac = 1.0;
    else frac = sqrt (dist1*dist1 - dist*dist) / dist3;
 }

 else {				/* element is not cable */
    frac = 0.0;
 }
 return frac;
}

/*------------------------------------------*/

void efrac(void)

/* return point in an element that is the closest to a node */

{
    datum d1={0},d2={0},d[MAXNODIM]={0};
    int nod1a, nod1b, nod1c, nod1d, elnum;
    double frac,dist,dist1,dist2,dist3;
    node *npnt;
    elem *epnt;

 checknum(d2 = popm());
 elnum = (int)(d2.val);
 frac = 0;
 getnod(d);					/* get node, check dims */

#ifndef XSTIMM
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;

 if (!(npnt=findnode(nod1a,nod1b,nod1c,nod1d,"Fdist"))) return; 

 if (!(epnt=findelem(elnum))) {
      ncfprintf (stderr,"Fdist: can't find element %d\n",elnum);
      execerror ("Missing element: "," stopping... ");
      return;
 }
 frac = enfrac(epnt,npnt);
#endif
  d1.val = frac;
  d1.vtype = NUMBER;
  d1.type = VAR;
  pushm(d1);
}

/*------------------------------------------*/

void n3dist(void)

/* return distance between two nodes */

{
    datum d1={0}, d[MAXNODIM];
    int nod1a, nod1b, nod1c, nod1d;
    int nod2a, nod2b, nod2c, nod2d;
    node *npnt1,*npnt2;

 getnod(d);					/* get node, check dims */
 nod2a = (int)d[0].val;
 nod2b = (int)d[1].val;
 nod2c = (int)d[2].val;
 nod2d = (int)d[3].val;

 getnod(d);					/* get node, check dims */
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;

#ifndef XSTIMM
 if (!(npnt1=findnode(nod1a,nod1b,nod1c,nod1d,"N3dist"))) return; 
 if (!(npnt2=findnode(nod2a,nod2b,nod2c,nod2d,"N3dist"))) return; 

 d1.val = dist3d (npnt1,npnt2);
#else
 d1.val = 0;
#endif
 d1.vtype = NUMBER;
 d1.type = VAR;
 pushm(d1);
}

/*------------------------------------------*/

void n2dist(void)

/* Return 2D distance between two nodes: */
/*  ignore z distance. */

{
    datum d1={0}, d[MAXNODIM];
    int nod1a, nod1b, nod1c, nod1d;
    int nod2a, nod2b, nod2c, nod2d;
    node *npnt1,*npnt2;

 getnod(d);					/* get node, check dims */
 nod2a = (int)d[0].val;
 nod2b = (int)d[1].val;
 nod2c = (int)d[2].val;
 nod2d = (int)d[3].val;

 getnod(d);					/* get node, check dims */
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;
 
#ifndef XSTIMM
 if (!(npnt1=findnode(nod1a,nod1b,nod1c,nod1d,"N2dist"))) return; 
 if (!(npnt2=findnode(nod2a,nod2b,nod2c,nod2d,"N2dist"))) return; 

 d1.val = dist2d (npnt1,npnt2);
#else
 d1.val = 0;
#endif
 d1.vtype = NUMBER;
 d1.type = VAR;
 pushm(d1);
}

/*------------------------------------------*/

void nzdist(void)

/* Return Z distance between two nodes */

{
    datum d1={0}, d[MAXNODIM];
    int nod1a, nod1b, nod1c, nod1d;
    int nod2a, nod2b, nod2c, nod2d;
    node *npnt1,*npnt2;

 getnod(d);					/* get node, check dims */
 nod2a = (int)d[0].val;
 nod2b = (int)d[1].val;
 nod2c = (int)d[2].val;
 nod2d = (int)d[3].val;

 getnod(d);					/* get node, check dims */
 nod1a = (int)d[0].val;
 nod1b = (int)d[1].val;
 nod1c = (int)d[2].val;
 nod1d = (int)d[3].val;
 
#ifndef XSTIMM
 if (!(npnt1=findnode(nod1a,nod1b,nod1c,nod1d,"Nzdist"))) return; 
 if (!(npnt2=findnode(nod2a,nod2b,nod2c,nod2d,"Nzdist"))) return; 

 d1.val = distzd (npnt1,npnt2);
#else
 d1.val = 0;
#endif
 d1.vtype = NUMBER;
 d1.type = VAR;
 pushm(d1);
}

/*------------------------------------------*/

void xgausnn(void)

/* Make random array of cells with defined regularity. */

{
  int i, first_center;
  Symbol *param,*var;
  datum d1={0},d2={0},d3={0};
  static int narg,numcells=0,numpts=0,grseed=0,ginfo=0;
  static double mean=0,stdev=0,density=0,reg=0;
  static double xframe=0,yframe=0,xcenter=0,ycenter=0;
  static double *x=0,*y=0;
  double *dpnt;


  param = (Symbol *)*pc++;
  switch (param->type) {

    case N: 	checknum(d1 = popm());
		numcells = int(d1.val);
		break;

    case NND:	checknum(d1 = popm());
		mean = d1.val;
		break;

    case NNSTD:	checknum(d1 = popm());
		stdev = d1.val;
		break;

    case RSD:  checknum(d1 = popm());
		grseed = int(d1.val);
		break;

    case DENSITY: checknum(d1 = popm());
		density = d1.val;
		break;

    case REG:   checknum(d1 = popm());
		reg = d1.val;
		break;

    case GINFO: checknum(d1 = popm());
		ginfo = int(d1.val);
		break;

    case SIZE:  narg = (long int) *pc++;
                if (narg > 2) d3 = popm();
                else d3.val = 0;
                if (narg > 1) d2 = popm();
                else d2.val = 0;
 		checknum(d1 = popm());
		xframe = d1.val;
		yframe = d2.val;
		break;

    case CENTER: narg = (long int) *pc++;
                if (narg > 2) d3 = popm();
                else d3.val = 0;
                if (narg > 1) d2 = popm();
                else d2.val = 0;
 		checknum(d1 = popm());
		xcenter = d1.val;
		ycenter = d2.val;
		break;

    case GAUSNN:var = (Symbol *)(*pc++);          /* array to be created */
#ifndef XSTIMM

  		if (setdisp>=0) disp = setdisp;	/* "-d n" overrides "disp=n" */
 		vidmode = (int)getval("vidmode");
    		numpts=gausnn(mean,stdev,density,reg,grseed,xframe,yframe,
			xcenter,ycenter,numcells,&x,&y,first_center=0,0,!vidmode,ginfo);

     		d1.val = numpts;	/* set up nc array dimensions */
     		d1.type = VAR;
     		d1.vtype = NUMBER;
     		d2.val = 2;
     		d2.type = VAR;
     		d2.vtype = NUMBER;
		pushm(d1);
		pushm(d2);
		dpnt = darr2(var,2);	     /* allocate space for array */
		for (i=0; i<numpts; i++) {
		  *dpnt++ = x[i];	     /* make x and y values be */
		  *dpnt++ = y[i];	     /*  successive values in array */
		}
		var->vtype = NUMBER;
		mean = stdev = density = reg = 0.0;
		xframe = yframe = xcenter = ycenter = 0;
		numcells = grseed = ginfo = 0;
		if (x) efree (x);	/* erase the subroutine's arrays */
		if (y) efree (y);
#else
     		d1.val = 0;	
     		d1.type = VAR;
     		d1.vtype = NUMBER;
#endif
		pushm (d1); 	/* finally, return n */
		break; 
  } /* switch */
}

/*------------------------------------------*/

void xphotrec(void)

/* read in data from rod or cone statement:

		1	x position
		2	y position  (optional)
*/

{
  datum d1,d2,d3;
  int narg;
  Symbol *param;
  node *npnt;
  double xpos,ypos,zpos;
  photorec *epnt;

#ifdef XSTIM
  recnod *maksnod(),*rpnt;
#endif

  param = (Symbol *)*pc++;			/* receptor type */
  narg  = (long int)*pc++;
  if (narg > 2) d3 = popm();			/* get z loc */
  else d3.val = LARGENODE;
  if (narg > 1) d2 = popm();			/* get y loc */
  else d2.val = LARGENODE;
  if (narg > 0) checknum(d1 = popm());		/* get x loc */
  else d1.val = LARGENODE;
  
  zpos  = d3.val;				/* z val */
  ypos  = d2.val;				/* y val */
  xpos  = d1.val;				/* x val */

  epnt = (photorec *)makelem(param->type,elpnt); /* make new photoreceptor */
  elpnt = (elem *)epnt;				/* save for recparm() below */

  if (narg > 2) 
     epnt->zpos  = zpos;			/* z val */
  if (narg > 1) 
     epnt->ypos  = ypos;			/* y val */
  epnt->xpos     = xpos;			/* x val */

						/* no z val in * struct photrec */
  if (!epnt->modif) {				/* never modify node */
#ifndef XSTIMM
    if (! (npnt=epnt->nodp1)) {
	ncfprintf (stderr,"photrec: can't find node %s\n",
		prnode(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d));
	return;
    }
  if (npnt->xloc>=LARGNOD) npnt->xloc = xpos;
  if (npnt->yloc>=LARGNOD) npnt->yloc = ypos;
#endif

#ifdef XSTIM
      rpnt = maksnod();
      rpnt->recnm1 = elpnt->node1a;     /* node numbers are defined by conn() */
      rpnt->recnm2 = elpnt->node1b;
      rpnt->recnm3 = elpnt->node1c;
      rpnt->recnm4 = elpnt->node1d;
      if (xpos!=LARGENODE) rpnt->xpos = xpos;
      else                 rpnt->xpos = 0;
      if (ypos!=LARGENODE) rpnt->ypos = ypos;
      else                 rpnt->ypos = 0;
#endif    
  }   /* if (!elpnt->modif) */
}

/*------------------------------------------*/

void xrecparm(void)

/* get parameters for rods and cones */

{
  datum d1;
  Symbol *param;
  int npops;
  photorec * ppnt;

  param  = (Symbol *)*pc++;
  
  switch (param->type) {
    case DIA:   
    case MAXCOND: 
    case PIGM:  
    case PATHL: 
    case ATTF:  
    case FILT:  
    case TIMEC1:
    case LOOPG:
    case DARKNOISE:
    case CHANNOISE:
    case PHOTNOISE: 
    case UNIT:
    case LINIT:
    case STIMCHAN:
    case RSD:   npops = 1;
		break;
    case SAVE:  
    case RESTORE:
		npops = 0;
		break;
   }
  
  if (npops==1) checknum(d1 = popm());

#ifndef XSTIMM

		/* Here, elpnt is a photrec, set by photrec() above */

  ppnt = (photorec *)elpnt; 
  switch (param->type) {
    case DIA:   ppnt->dia = d1.val;
		break;
    case MAXCOND: ppnt->maxcond = d1.val; 
		break;
    case PIGM:  ppnt->pigm  = d1.val; 
		break;
    case PATHL: ppnt->pathl = d1.val; 
		break;
    case ATTF:  ppnt->attf = d1.val; 
		break;
    case FILT:  ppnt->filt = (int)d1.val; 
		break;
    case TIMEC1:ppnt->timec1 = d1.val; 
		break;
    case LOOPG: ppnt->loopg = d1.val; 
		break;
    case RSD:	ppnt->phrseed = (int)d1.val; 
		break;
    case SAVE:  ppnt->save = 1; 
		break;
    case RESTORE:
		ppnt->restore = 1; 
		break;
    case PHOTNOISE:
		ppnt->photnoise = (int)d1.val; 
		break;
    case CHANNOISE:
		ppnt->channoise = d1.val; 
		break;
    case UNIT:
		ppnt->unitary = d1.val; 
		break;
    case DARKNOISE:
		ppnt->darknoise = d1.val; 
		break;
    case LINIT:
		ppnt->linit = d1.val; 
		break;
    case STIMCHAN:
		ppnt->stimchan = d1.val; 
		break;
   }
#endif 		/* XSTIMM */
}

/*------------------------------------------------------*/

double axy (double x, double y)

/* return angle in radians of stimulus at (x,y) point */


{
    double t;
#define SMALLVAL 1e-6

  if (x > SMALLVAL) {
    if (y > 0) t = atan(y/x);
    else       t = 2*MPI + atan(y/x);
  }
  else if (x < -SMALLVAL) {
    if (y>0)   t = MPI + atan(y/x);
    else       t = MPI + atan(y/x);
  }

  else {  /* x == 0 */

    if (y > SMALLVAL)       t =  MPI/2; 
    else if (y < -SMALLVAL) t = - MPI/2;
    else t = 0;
  };
  return t;
};

/********************************************************/

extern double dsintinc;		/* default sine time incr in control.h */
extern double dsintres;		/* default sine time time res in control.h */

int blurfl=0;			/* blur array done */
int cumstim=0;			/* number of stimulus reports */
int arraymade=0;		/* stimulus array already allocated */
int stimarrsize;		/* size of stimulus array */


#ifdef XSTIM
					/* this part for "stim", not "nc" */
void xstim()

/* read in data from stimulus statement:
*/

{
  datum d1,d2,d3,d4,d[MAXNODIM];
  Symbol *param;
  int stimarr, backarr;
  int narg, spec;
  FILE *ftemp;

static double start=0.0;
static double dur=0.1;
static double stop=0.1;
static double wavel=1.0;
static double inten=100.0;
static double backgr=0.0;
static double blurrad=0.0,tblur=0.0;
static double blur_ampl=1.0;
static double scatter_rad=0.0;
static double scatter_pow=0.0;
static double scatter_ampl=0.0;
static double tscatter_rad=0.0;
static double tscatter_pow=0.0;
static double tscatter_ampl=0.0;
static double speriod = 20;
static double contrast = 0.5;
static double tfreq = 0.0;
static double sscale = 1.0;
static double xenv = 10.0;
static double yenv = 10.0;
static double orient = 0.0;
static double sphase = 0.0;
static double puffconc = 0.0;
static double drift=0;
static double mask=0;
static int stimchan=0;
static int tmod=0;
static int xn,yn;
static int makenv=0;
static int setxenv=0;
static int setyenv=0;
static int sq=0;		/* make sine into square wave */
static int invert=0;		/* invert spot for mask */

static int stype;
/* static int rnum1,rnum2,rnum3,rnum4;  /* not used for "stim" */
static double xloc=0, yloc=0;
static double xcent=0, ycent=0;		/* center coords for stim array */
static int centerfl=0;                  /* center specified */
static double xmax,xmin,ymax,ymin;	/* borders of of photrec array */
static double xrange, yrange;		/* size of photrec array */
static int stimflag=0;                  /* start stimulus */
extern int blursize;			/* radius of blur array from makblur()*/
double toffset,tperiod,tincr,stime,stimend;
const char *imgfil = NULL;		/* image file, raw format */
double *stim_rndarr;			/* stimulus array for makcheckerboard()*/
int stim_nfr;				/* stimulus number of frames */
extern double stimdia;                     /* defined in "stimsub.c" */
extern double stimydist;                   /* defined in "stimsub.c" */

  param  = (Symbol *)*pc++;
  
  switch (param->type)
   {
    case SIMAGE:
    		checkstr(d1 = popm());
		imgfil = d1.str;	
		break;
    case SECTOR:
    case BAR: 
    case SPOT:  checknum(d1 = popm());
                stimdia   = d1.val;
                stype = param->type; 
                break;

    case WINDMILL:  
    case SINEANN:  
    case GABOR:  
    case SINE:  checknum(d1 = popm());
                speriod   = d1.val;
                stype = param->type; 
                break;
    case CHECKERBOARD:
		narg = (long int)*pc++;
  		if (narg > 3) d4 = popm();
		else d4.val = 25;
  		if (narg > 2) d3 = popm();
		else d3.val = 25;
		xn = d3.val; 
		yn = d4.val;
  		if (narg > 1) d2 = popm();
		else d2.val = 1.0;
 		checknum(d1 = popm());
		stimdia = d1.val;
		stimydist = d2.val;
		stype = param->type; 
		break;
    case RECT:  narg = (long int)*pc++;
  		if (narg > 1) d2 = popm();
		else d2.val = 1.0;
 		checknum(d1 = popm());
		stimdia = d1.val;
		stimydist = d2.val;
		stype = param->type; 
		break;
    case NODE:
		getnod(d);			/* get node num, check dims */
    /*          rnum1 = d[0].val;		/* get node number */
    /*          rnum2 = d[1].val;
                rnum3 = d[2].val;
                rnum4 = d[3].val;
    */
                break;
    case VTRANSDUCER: 
    case ITRANSDUCER: 
    case ROD: 
    case CONE:  
    case CHR:  
		getnod(d);			/* get node num, check dims */
    /*          rnum1 = d[0].val;		/* get node number */
    /*          rnum2 = d[1].val;
                rnum3 = d[2].val;
                rnum4 = d[3].val;
    */
                stype = param->type; 
                break;
    case LOC:   narg  = (long int)*pc++;
                if (narg > 0) {
                  if (narg > 1) {
                    checknum(d2 = popm());	/* read in y loc */
                    yloc = d2.val;
                  }
		  else yloc = 0;
 		checknum(d1 = popm());
                xloc = d1.val;
               }
                break;
    case CENTER:narg  = (long int)*pc++;
                if (narg > 0) {
                  if (narg > 1) {
                    if (narg > 2) 
                       d3 = popm();
                    d2 = popm();                /* read in y center */
                    ycent = d2.val;
                  }
		  else ycent = 0;
 		checknum(d1 = popm());
                xcent = d1.val;
		centerfl = 1;
                }
                break;
    case SSCALE: checknum(d1 = popm());
                sscale   = d1.val;
                break;
    case START: checknum(d1 = popm());
                start = d1.val;
                stimflag = 1;
                break;
    case STOP:  checknum(d1 = popm());
                stop = d1.val;
		dur = start - stop;
                stimflag = 1;
                break;
    case DUR:   checknum(d1 = popm());
                dur   = d1.val;
                break;
    case CONTRAST: checknum(d1 = popm());
                contrast   = d1.val;
                break;
    case TFREQ: checknum(d1 = popm());
                tfreq = d1.val;
                break;
    case DRIFT: checknum(d1 = popm());
                drift = int(d1.val);
                break;
    case ORIENT: checknum(d1 = popm());
                orient = d1.val;
                break;
    case SPHASE: checknum(d1 = popm());
                sphase = d1.val;
                break;
    case XENV:	checknum(d1 = popm());
                xenv = d1.val;
		setxenv = 1;
                break;
    case YENV:	checknum(d1 = popm());
                yenv = d1.val;
		setyenv = 1;
                break;
    case SQ:    checknum(d1 = popm());
                sq = d1.val;
                break;
    case WAVEL: spec = (long int)*pc++;
		switch (spec) {
		  case 0: 
 		    checknum(d1 = popm());
		    wavel = d1.val;
		    break;
		  case XENON:
                    wavel = 0;
                    break;
		  case SUN:
                    wavel = 1;
                    break;
		  case TUNGSTEN:
                    wavel = 2;
                    break;
		}
		break;

    case INTEN: checknum(d1 = popm());
                inten = d1.val;
                break;

    case MASK: checknum(d1 = popm());
                mask = d1.val;
                break;

    case STIMCHAN: checknum(d1 = popm());
                stimchan = setstimchan(int(d1.val));
                break;

    case BACKGR: checknum(d1 = popm());
                backgr = d1.val;
                backarr = 1;
                stype = param->type; 
                break;

    case VCLAMP:
    case CCLAMP: checknum(d1 = popm());
                inten = d1.val;
                stype = param->type; 
                break;

    case PUFF:  checknum(d1 = popm());
                puffconc = d1.val;
		pc++;
                break;

    case BLUR:  checknum(d1 = popm());		 /* get blur diameter */
                tblur = d1.val/2;		 /* divide dia by 2 = radius */
		if (blurfl) {
		   if (blurrad!=tblur) blurfl = 0;  /* remake blur function */
                }
		blurrad = tblur;
                break;

    case SCATTER:
		narg = (long int)*pc++;
                if (narg > 0) {
                  if (narg > 1) {
                    if (narg > 2) 
                       checknum(d3=popm());		/* falloff power */
                    checknum(d2=popm());                /* get dia */
                  }
 		 checknum(d1 = popm());			/* get amplitude */
                 tscatter_ampl = d1.val;
                 tscatter_rad  = d2.val/2;
                 tscatter_pow  = d3.val;
		}
		else {					/* use cat optics */
                 tscatter_ampl = .15;
                 tscatter_rad  = 2.05;
                 tscatter_pow  = 2.5;
		}
		if (blurfl) {
		 if (scatter_ampl!=tscatter_ampl &&
                    scatter_rad !=tscatter_rad  &&
		    scatter_pow !=tscatter_pow) blurfl=0; /* remake blur func */
		}
                scatter_ampl = tscatter_ampl;
                scatter_rad  = tscatter_rad;
                scatter_pow  = tscatter_pow;
                break;
    case SFILE: if (*stimfile) fclose (stimout);
                *stimfile = 0;
		if (!checkstr(d1=popm())) break; 	 /* get filename */
                strcpy (stimfile,d1.str);  		/* copy filename */
                if (stimhlp) ncfprintf (stderr,"stimfil '%s'\n",stimfile);
                if (*stimfile) {
                        stimfflg = 1; 
                        if ((ftemp=fopen(stimfile,"w")) == NULL) {
                          ncfprintf (stderr,"%s: can't open file '%s'\n",
						progname,stimfile);
                          break;
                        }
                        else stimout = ftemp;
                }
                else {
                        stimfflg = 0;
                        stimout = stdout;
                }
		stcomment();			/* print comment in stim file */
                break;

    case STIM:  
        backarr = BACKARR;
        stimarr = STIMARR;
					/* find photoreceptor array size */

    if (interp) varcopy();
    if (setdebug) debug = setdebug;
    if (setdebgz) debugz = setdebgz;

     switch (stype) {
       case BAR:
       case RECT:
       case CHECKERBOARD:
       case SPOT:
       case SINE:
       case SINEANN:
       case GABOR:
       case BACKGR:
        if (!blurfl) {
	   blursize=makblur(blurrad,blur_ampl,sscale,  /* make a blur function */
		scatter_rad, scatter_pow, scatter_ampl); 
           blurfl = 1;
	   ncfprintf (stderr,"## %d %-6s blur     array size %d\n",
				++cumstim,findsym(stype),blursize*2); /* */
	   if (arraymade) {			/* if backgr array made */
		 delc(2);			/* delete existing arrays */
		 arraymade = 0;
	   }
        }
	if (!arraymade) {
	   find_photarrsiz(&xmax,&xmin,&ymax,&ymin);  /* find maxmin of photrecs */
	   // if (!centerfl) {
	     xcent = (xmax+xmin) * 0.5;
	     ycent = (ymax+ymin) * 0.5;
	   // }
	   xrange = xmax-xmin;
	   yrange = ymax-ymin;
	   if (sscale <=0) sscale = 1.0;
	   stimarrsize = (long int)(max(xrange,yrange) / sscale);
	   if (stimarrsize==0) stimarrsize = 1;
	   if (stimarrsize&1) stimarrsize++;
	   if (blurrad>0) stimarrsize += blursize*2;

	   ncfprintf (stimout,"## %d %-6s stimulus array size %d\n",
				++cumstim,findsym(stype),stimarrsize);	/* */

	   if (initc(2,stimarrsize)<=0) {		/* make stimulus to fit */
	       ncfprintf (stderr,"Stimulus array is too big = %d\n",stimarrsize);
	       ncfprintf (stderr,"Blur array dia %d\n",blursize*2);
	       ncfprintf (stderr,"Stimulus xmax %g xmin %g ymax %g ymin %g\n",
			xmax,xmin,ymax,ymin);
    	       execerror ("can't continue"," stopping...");
	   }
	   arraymade = 1;
	}
	 break;
       default:
	  break;
      } /* switch stype */

        switch (stype) {

        case 0:
        default: break;

/* Vclamp, cclamp, rod and cone are done in xstim() for "nc" below.
   They are not performed by (the program) "stim" because they don't
   require blur, and sometimes need to be done in addition to a
   predefined, blurred stimulus.
*/

        case VCLAMP:
/*              vclist(start,    rnum1,rnum2,rnum3,rnum4,inten,wavel,"vcl");
                vclist(start+dur,rnum1,rnum2,rnum3,rnum4,inten,wavel,"eff");
*/              break;

        case CCLAMP:
/*              vclist(start,    rnum1,rnum2,rnum3,rnum4,inten,wavel,"ccl");
                vclist(start+dur,rnum1,rnum2,rnum3,rnum4,inten,wavel,"iof");
*/              break;

        case VTRANSDUCER:
        case ITRANSDUCER:
        case CHR:
        case ROD:
        case CONE: 				/* flash one photrec */
/*              vclist(start,    rnum1,rnum2,rnum3,rnum4, inten,wavel,"del");
                vclist(start+dur,rnum1,rnum2,rnum3,rnum4,-inten,wavel,"del"); 
*/              break;

        case PUFF: 				/* puff neurotrans/2nd msng */
                  break;

        case BACKGR:
                  recback(stimarr,0.0,wavel);   /* zero photrec stim inten */
                  recback(backarr,backgr,wavel); /* set photrec backgr inten */
                  abslist(0.0, start, mask, stimchan);   /* make action list */
                  break;
        case BAR:
		makrect(stimarr, stimarrsize, stimdia, (double)CONVSIZE,
			xloc-xcent,yloc-ycent,sscale,orient,inten,wavel);
                if (stimflag) {
			double xr;
                  stimflag = 0;
		  if (int(stimdia/sscale+0.5) & 1) xr = 0; else xr = 0.5;
                  recpat(stimarr, stimarrsize,
			xcent+xr*sscale,ycent,sscale);   /* set photrec inten */
                  stimlist(1.0, start, mask, stimchan);  /* make an action list */
                  stimlist(-1.0, start+dur, mask, stimchan);
                  recback(stimarr,0.0,1.0);   /* zero photrec stim inten */
                }
                break;

        case CHECKERBOARD:
	    	makcheckerboard(stimarr,stimarrsize, stimdia, stimydist, xn, yn, 
				xloc, yloc, xcent, ycent, 
				sscale, orient, tfreq, inten, contrast,start,dur, 
				&stim_rndarr, &stim_nfr, rseed);  
                break;
        case RECT:
		makrect(stimarr, stimarrsize, stimdia, stimydist,
			xloc-xcent,yloc-ycent,sscale,orient,inten,wavel);
                if (stimflag) {
                  stimflag = 0;
                  recpat(stimarr,stimarrsize,xcent,ycent,sscale); /* set photrec inten */
                  stimlist(1.0, start,mask,stimchan);             /* make an action list */
                  stimlist(-1.0, start+dur,mask,stimchan);
                  recback(stimarr,0.0,1.0);   /* zero photrec stim inten */
                }
                break;

        case SPOT:
                makspot(stimarr,stimarrsize,stimdia,xloc,yloc,xcent,ycent,
						sscale,inten,start,dur,wavel,mask,stimchan,invert=0);
                break;

        case SECTOR:				/* sector mask */
		double w,h,xl,yl,th1,th2;
		double mask;
#define BARLENGTH 1000
		if (stimdia < 0) stimdia = -stimdia;
		th1 = orient - stimdia/2;
		if (th1 <  0)   th1 += 360;
		if (th1 >= 360) th1 -= 360;
		th2 = th1 + stimdia + 180;
		if (th2 > 360) th2 -= 360;
		w = h = BARLENGTH;

		makbar(stimarr, stimarrsize, w, h, xl=-w*0.5*cos(th1/180.0*MPI)+xloc, 
			yl=w*0.5*sin(th1/180.0*MPI)+yloc, xcent=0, ycent=0, 
			sscale=1, th1, inten, start, dur, wavel, mask=1, stimchan);
		makbar(stimarr, stimarrsize, w, h, xl=-w*0.5*cos(th2/180.0*MPI)+xloc, 
			yl=w*0.5*sin(th2/180.0*MPI)+yloc, xcent=0, ycent=0, 
			sscale=1, th2, inten, start, dur, wavel, mask=1, stimchan);
#undef BARLENGTH

                break;

        case SINE:
		if (!setxenv) xenv = 1e30;
		if (!setyenv) yenv = 1e30;
		makgabor(stimarr, stimarrsize, contrast, inten, 
			xloc, yloc, xcent, ycent,
                	xenv, yenv, makenv=0, sq, tfreq, drift, speriod, sphase,
                	sscale, orient, wavel, mask, start, dur);
		setxenv = setyenv = 0;
	  	break;

        case GABOR:
		makgabor(stimarr, stimarrsize, contrast, inten, 
			xloc, yloc, xcent, ycent,
                	xenv, yenv, makenv=1, sq, tfreq, drift, speriod, sphase,
                	sscale, orient, wavel, mask, start, dur);
	  	break;

        case SINEANN:
		maksineann(stimarr, stimarrsize, contrast, inten, 
			xloc, yloc, xcent, ycent,
                	xenv, makenv=1, sq, tfreq, drift, speriod, sphase, 
			sscale, wavel, mask, start, dur);
	  	break;  

        case WINDMILL:
		makwindmill(stimarr, stimarrsize, contrast, inten, 
			xloc, yloc, xcent, ycent,
                	xenv, makenv=1, sq, tfreq, drift, speriod, sphase, 
			sscale, wavel, mask, start, dur);
	  	break;  

        case SIMAGE:
		makimage(stimarr, stimarrsize, xenv, yenv, 
			xloc-xcent, yloc-ycent, orient, sscale,
			inten, 0, imgfil);
                if (stimflag) {
                  stimflag = 0;
                  recpat(stimarr,stimarrsize,xcent,ycent,sscale); /* set photrec inten */
                  stimlist(1.0, start,mask,0);         /* make an action list */
                  stimlist(-1.0, start+dur,mask,0);
                  recback(stimarr,0.0,1.0);   /* zero photrec stim inten */
                }

        }       /* switch (stype) */

	stype = 0;		/* reset stype for next stim statement */
	orient = 0;
	scatter_ampl = 0;		 /* get scatter amplitude */
	scatter_pow = 0;		 /* get scatter amplitude */
        break;  /* case STIM */

  }     /* switch */
}

#endif 			/* end of xstim for generating stimulus file */

/*------------------------------------------------------*/

void find_photarrsiz(double *xmax,double *xmin,double *ymax,double *ymin)

/* Find max, min for the array of photoreceptors
    described by the user's program.
*/

{
  recnod *npt;
  double xma,xmi,yma,ymi;

  xma = yma = 0.0;
  xmi = ymi = 1e10;
  for (npt=reclist; npt; npt=npt->next) {
    xma  = max(xma,npt->xpos);
    xmi  = min(xmi,npt->xpos);
    yma  = max(yma,npt->ypos);
    ymi  = min(ymi,npt->ypos);
  }
  *xmax = xma;
  *xmin = xmi;
  *ymax = yma;
  *ymin = ymi;
}

/********************************************************/

#ifdef XMOD

void xstim(void)		/* this part for nc (modcode.o): */

/* read in data from stimulus statement: */

{
#define  STFILNAM 40
#define  NCHECKS 25
  datum d1,d2,d3,d4,d[MAXNODIM];
  Symbol *param;
  int narg, spec;
  double xpos, ypos, val;
  static double dia, ydist, halfx, halfy, radius;
  static double orientrad, sphasenorm, sino, coso;
  double xr, yr, xcent, ycent;
  static double start,dur;
  static double inten,sinten,backgr,xsine;
  static double xedgemin, xedgemax, yedgemin, yedgemax;
  static int stype,rnum1,rnum2,rnum3,rnum4,filuse=0;
  static const char *filename;
  static double xloc = 0;
  static double yloc = 0;
  static double ymins = 0;
  static double ymaxs = 0;
  static double xenv = 50;
  static double yenv = 50;
  static double wavel = 1.0;
  static double speriod = 20;
  static double contrast = 0.5;
  static double tfreq = 0.0;
  static double orient = 0.0;
  static double sphase = 0.0;
  static double sscale = 1.0;
  static char stimfile[STFILNAM]={0};
  static char stimfilz[STFILNAM]={0};
  static double puffconc = 0.0;
  static int puffmsg=0;			/* CAMP or CGMP */
  static double drift = 0;
  static int makenv = 0;
  static int setxenv = 0;
  static int setyenv = 0;
  static int sq=0;		/* make sine into square wave */
  static int xn,yn;
  static int stimchan=0;
  static double mask = 0;
  static const char *action;
  static double *stim_rndarr;	/* stimulus array */
  static int stim_nfr;		/* stimulus number of frames */
  photorec *ept;
  photrec *rpnt;
  char *imgfil = NULL;		/* image file, raw format */
  FILE *ftemp, *fimg;

  param  = (Symbol *)*pc++;

  switch (param->type)
   {
    case CHECKERBOARD:  
		narg = (long int)*pc++;
  		if (narg > 3) d4 = popm();
		else d4.val = NCHECKS;
  		if (narg > 2) d3 = popm();
		else d3.val = NCHECKS;
		xn = d3.val; 
		yn = d4.val;
  		if (narg > 1) d2 = popm();
		else d2.val = 1.0;
 		checknum(d1 = popm());
		dia   = d1.val;
		ydist = d2.val;
		stype = param->type; 
		break;
    case RECT:  
		narg = (long int)*pc++;
  		if (narg > 1) d2 = popm();
		else d2.val = 1.0;
 		checknum(d1 = popm());
		dia   = d1.val;
		ydist = d2.val;
		stype = param->type; 
		break;
    case SIMAGE: 
		break;
    case SECTOR:
    case BAR: 
    case SPOT:  checknum(d1 = popm());
		dia   = d1.val;
		stype = param->type; 
		break;
    case WINDMILL:  
    case SINEANN:  
    case GABOR:  
    case SINE:  checknum(d1 = popm());
		speriod  = d1.val;
		stype = param->type; 
		break;
    case NODE:
		getnod(d);			/* get node num, check dims */
		rnum1 = (int)d[0].val;		/* get node number */
		rnum2 = (int)d[1].val;
		rnum3 = (int)d[2].val;
		rnum4 = (int)d[3].val;
		break;

    case VTRANSDUCER: 
    case ITRANSDUCER: 
    case CHR: 
    case ROD: 
    case CONE:
		getnod(d);			/* get node num, check dims */
		rnum1 = (int)d[0].val;		/* get node number */
		rnum2 = (int)d[1].val;
		rnum3 = (int)d[2].val;
		rnum4 = (int)d[3].val;
		stype = param->type; 
		break;
    case LOC:	narg  = (long int)*pc++;
  		if (narg > 0) {
  		  if (narg > 1) {
     	 	    d2 = popm();		/* read in y loc */
		    yloc = d2.val;
  		  }
                else yloc = 0;
 		checknum(d1 = popm());   	/* x loc */
		xloc = d1.val;
		}
		break;
    case CENTER: narg  = (long int)*pc++;
  		if (narg > 0) {
  		  if (narg > 1) {
  		    if (narg > 2) 
     	 	       d3 = popm();
     	 	    d2 = popm();		/* read in y loc */
  		  }
 		checknum(d1 = popm());	 	/* x loc */
		}
		break;
    case START:	checknum(d1 = popm());
		start = d1.val;
		break;
    case DUR:	checknum(d1 = popm());
		dur   = d1.val;
		break;
    case CONTRAST: checknum(d1 = popm());
                contrast   = d1.val;
                break;
    case TFREQ:	checknum(d1 = popm());
		tfreq  = d1.val;
		break;
    case DRIFT:	checknum(d1 = popm());
		drift  = int(d1.val);
		break;
    case ORIENT:checknum(d1 = popm());
		orient  = d1.val;
		break;
    case SPHASE:checknum(d1 = popm());
		sphase  = d1.val;
		break;
    case XENV:	checknum(d1 = popm());
		xenv  = d1.val;
		setxenv = 1;
		break;
    case YENV:	checknum(d1 = popm());
		yenv  = d1.val;
		setyenv = 1;
		break;
    case SQ:	checknum(d1 = popm());
		sq  = d1.val;
		break;
    case WAVEL: spec = (long int)*pc++;	
		switch (spec) {
		  case 0:
 		    checknum(d1 = popm());
		    wavel = d1.val;
		    break;
		  case SUN:
		    wavel = 0;
		    break;
		  case XENON:
		    wavel = 1;
		  case TUNGSTEN:
                    wavel = 2;
                    break;
		  break;
		}
		break;
    case SSCALE:checknum(d1 = popm());
		sscale = d1.val;
		break;
    case INTEN:	checknum(d1 = popm());
		inten = d1.val;
		break;
    case BACKGR: checknum(d1 = popm());
		backgr = d1.val;
		stype = param->type; 
		break;
    case MASK:	checknum(d1 = popm());
		mask = d1.val;
		break;
    case STIMCHAN: checknum(d1 = popm());
		stimchan = setstimchan(int(d1.val));
		break;
    case VCLAMP:
    case CCLAMP: checknum(d1 = popm());
		inten = d1.val;			/* get node number */
		stype = param->type; 
		break;
    case PUFF:  checknum(d1 = popm());
                puffconc = d1.val;
		stype = param->type;
  		param = (Symbol *)*pc++;
  		puffmsg = param->type;
                break;

    case BLUR:	checknum(d1 = popm());
		/* blurrad = d1.val; */
		break;
    case SCATTER:
		narg = (long int)*pc++;
                if (narg > 0) {
                  if (narg > 1) {
                    if (narg > 2) 
                       checknum(d3=popm());		/* falloff power */
                    checknum(d2=popm());                /* get dia */
                  }
 		 checknum(d1 = popm());			/* get amplitude */
                /* scatter_ampl = d1.val;
                 scatter_rad  = d2.val/2;
                 scatter_pow  = d3.val;
		*/
		}
                break;
    case SFILE:	if (!checkstr(d1=popm())) break; 	 /* get filename */
		filename = d1.str;
		if (stimfflg && *stimfile) fclose(stimin); /* close previous */
                strcpy (stimfile,filename);  		/* copy filename */
		stfile = stimfile;			/* for plotinit */
		if (! *filename) { filuse = 0; break; }
		if ((ftemp=fopen(stimfile,"r"))==NULL) { /* open file */
                  strcpy(stimfilz,stimfile);
                  strncat(stimfilz,".gz",STFILNAM-5);
		  if ((ftemp=fopen(stimfilz,"r"))==NULL){ /* look for .gz file */
		     filuse = 0;		/* if not, don't use file */
		     break;
		  }
		  else {				/* .gz file exists */
		     fclose(ftemp);
                     strcpy(stimfilz,"zcat ");		/* check for .gz file */
                     strncat(stimfilz,stimfile,STFILNAM-5);
                     strncat(stimfilz,".gz",STFILNAM-5);
                     if ((ftemp=popen(stimfilz,"r"))==NULL) {
			ncfprintf (stderr,"%s: can't open file '%s'\n",
					progname,stimfile);
			filuse = 0;
			break;
		     }
		  }
		}
		stimin = ftemp;
		filuse = 1;
		break;

    case STIM:	 			/* real time stimulus without blur */

					/* set masking light stimulus */
	if (mask>0) action = "b";	/* b -> mask, d -> additive */
	else        action = "d";

	switch (stype) {

	case VCLAMP:
		makvstim(start, rnum1,rnum2,rnum3,rnum4,inten,0,"evclamp");
		makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-inten,1,"evclamp");
		break;

	case CCLAMP:
		makvstim(start, rnum1,rnum2,rnum3,rnum4,inten,0,"gclamp");
		makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-inten,1,"gclamp");
		break;

        case PUFF: 				/* puff neurotrans/2nd msng */
	   switch (puffmsg) {
		case GLU:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"i");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"i");
		break;
		case AMPA:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"j");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"j");
		break;
		case KAINATE:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"k");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"k");
		break;
		case NMDA:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"l");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"l");
		break;
		case CNQX:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"m");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"m");
		break;
		case GABA:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"n");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"n");
		break;
		case BIC:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"o");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"o");
		break;
		case PTX:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"p");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"p");
		break;
		case GLY:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"q");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"q");
		break;
		case STRY:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"r");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"r");
		break;
		case CAMP:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"s");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"s");
		break;
		case CGMP:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"t");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"t");
		break;
		case CA:
		 makvstim(start,    rnum1,rnum2,rnum3,rnum4, puffconc,0,"u");
		 makvstim(start+dur,rnum1,rnum2,rnum3,rnum4,-puffconc,1,"u");
		break;
	    }
           break;
	case ROD:			/* add to inten on specific receptor */
	case CONE: 
	case CHR: 
	case VTRANSDUCER: 
	case ITRANSDUCER: 
	 makrstim(start+dur-RNDUP,rnum1,rnum2,rnum3,rnum4,-inten,wavel,-mask,stimchan,1,action);
	 makrstim(start,    rnum1,rnum2,rnum3,rnum4, inten,wavel,mask,stimchan,0,action);
	      break;	

       }

      if (filuse) stimfflg = 1;
      else {

        stimfflg = 0;

	switch (stype) {

	case VCLAMP:			/* don't make these stim twice */ 
	case CCLAMP:
	case PUFF:
	case ROD:
	case CONE: 
	case CHR: 
	case VTRANSDUCER: 
	case ITRANSDUCER: 
		  break;

	case BACKGR: 					/* abs inten */

	       if (mask) action = "c";	/* c -> mask, a -> additive */
	       else      action = "a";
	       			/* look for existing photoreceptors: */
		if (runyet) {
	  	 for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
		  if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
		      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
		    makrstim (start,rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				rpnt->recnm4,backgr,wavel,mask,stimchan,0,action);
		  }
                 }
				/* look for new rod, code elements, */
			        /*  even if run yet: */

	  	 for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
		  if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
		      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
		    makrstim (start,ept->node1a,ept->node1b,ept->node1c,
				ept->node1d,backgr,wavel,mask,stimchan,0,action);
		 }
		break;


	case BAR:			/* stim for possibly all receptors */
	  halfx = dia / 2;
	  halfy = 1000;
	  if (runyet) {
	   for (rpnt=recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
	    if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
		      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
	    if (inrect(rpnt->xloc,rpnt->yloc, xloc+halfx, xloc-halfx,
			yloc+halfy, yloc-halfy, orient)) {
		makrstim(start,    rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,inten,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,-inten,wavel,-mask,stimchan,1,action);
	     }
	   }
	  }
	  else {
	   for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
	    if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
		      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
	    if (inrect(ept->xpos,ept->ypos,
			xloc+halfx, xloc-halfx,
			yloc+halfy, yloc-halfy, orient)) {
		makrstim(start,    ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,inten,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,-inten,wavel,-mask,stimchan,1,action);
	    }
	   }
	  }
	  break;

	case SINE:			/* stim for possibly all receptors */
	  if (!setxenv) xenv = 1e30;
	  if (!setyenv) yenv = 1e30;
	  makgaborx(speriod, sphase, orient, xloc, yloc, tfreq, drift,
                inten, contrast, wavel, xenv, yenv, makenv=0, sq,
		mask, stimchan, action, start, dur);
          setxenv = setyenv = 0;
	  break;

	case GABOR:			/* stim for possibly all receptors */
	    makgaborx(speriod, sphase, orient, xloc, yloc, tfreq, drift,
                inten, contrast, wavel, xenv, yenv, makenv=1, sq,
		mask, stimchan, action, start, dur);
	  break;

	case SINEANN:			/* stim for possibly all receptors */

	  maksineannx(speriod, sphase, xloc, yloc, tfreq, drift, inten, 
		contrast, wavel, xenv, makenv=1, sq, mask, stimchan, action, start, dur);

	  break;

	case WINDMILL:			/* stim for possibly all receptors */

	  makwindmillx(speriod, sphase, xloc, yloc, tfreq, drift, inten, 
		contrast, wavel, xenv, makenv=1, sq, mask, stimchan, action, start, dur);

	  break;

	case CHECKERBOARD:		/* stim for possibly all receptors */
	    makcheckerboardx(dia, ydist, xn, yn, orient, xloc, yloc, 
			    tfreq, inten, contrast,action,start,dur,
			    &stim_rndarr,&stim_nfr, rseed);  
	  break;

	case SECTOR:		/* sector_mask */
	     sector_mask(orient, dia, inten, start, dur);
	  break;

	case RECT:			/* stim for possibly all receptors */
	  halfx = dia / 2;
	  halfy = ydist / 2;
          if (runyet) {
	  for (rpnt=(photrec*)recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
	    if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR || 
		      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
	    if (inrect(rpnt->xloc,rpnt->yloc, xloc+halfx, xloc-halfx,
			yloc+halfy, yloc-halfy, orient)) {
		makrstim(start,    rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,inten,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,-inten,wavel,-mask,stimchan,1,action);
	    }
	   }
          }
	  else {
	   for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
	    if (ept->ctype == ROD || ept->ctype == CONE ||  ept->ctype == CHR ||
		      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
	    if (inrect(ept->xpos,ept->ypos, xloc+halfx, xloc-halfx,
			yloc+halfy, yloc-halfy, orient)) {
		makrstim(start,    ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,inten,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,-inten,wavel,-mask,stimchan,1,action);
	    }
	  }
	 }
	  break;

	case SIMAGE:			/* stim for possibly all receptors */

	  if ((fimg=fopen(imgfil,"r"))==NULL) { /* open file */
	    ncfprintf (stderr,"%s: can't open image file '%s'\n", progname,imgfil);
	    break;
	  }

	  halfx = dia / 2;
	  halfy = ydist / 2;
          if (runyet) {
	  for (rpnt=(photrec*)recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
	    if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
		      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER)
	    if (inrect(rpnt->xloc,rpnt->yloc, xloc+halfx, xloc-halfx,
			yloc+halfy, yloc-halfy, orient)) {
	        fread(&val,1,1,fimg);
		makrstim(start,    rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,inten+val,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,-(inten+val),wavel,-mask,stimchan,1,action);
	    }
	   }
          }
	  else {
	   for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
	    if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
		      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER)
	    if (inrect(ept->xpos,ept->ypos, xloc+halfx, xloc-halfx,
			yloc+halfy, yloc-halfy, orient)) {
	        fread(&val,1,1,fimg);
		makrstim(start,    ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,inten+val,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,-(inten+val),wavel,-mask,stimchan,1,action);
	    }
	  }
	 }
	  break;

	case SPOT: 
	  radius = dia / 2;
	  if (runyet) {			/* look for existing photoreceptors */
	   for (rpnt=(photrec*)recpnt; rpnt; rpnt=(photrec*)rpnt->next) {
	    if (rpnt->ctype == ROD || rpnt->ctype == CONE || rpnt->ctype == CHR ||
		      rpnt->ctype == VTRANSDUCER || rpnt->ctype == ITRANSDUCER) {
	    if (incirc(rpnt->xloc,rpnt->yloc, xloc, yloc,radius)) {
		makrstim(start,    rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,inten,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,rpnt->recnm1,rpnt->recnm2,rpnt->recnm3,
				   rpnt->recnm4,-inten,wavel,-mask,stimchan,1,action);
	     }
	    }
	   }
	  }
	  else {
	   for (ept=(photorec*)elempnt; ept; ept=(photorec*)ept->next) {
	    if (ept->ctype == ROD || ept->ctype == CONE || ept->ctype == CHR ||
		      ept->ctype == VTRANSDUCER || ept->ctype == ITRANSDUCER) {
	    if (incirc(ept->xpos,ept->ypos, xloc, yloc,radius)) {
		makrstim(start,    ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,inten,wavel,mask,stimchan,0,action);
		makrstim(start+dur-RNDUP,ept->node1a,ept->node1b,ept->node1c,
				   ept->node1d,-inten,wavel,-mask,stimchan,1,action);
	    }
	   }
 	  }
	 }
	  break;
	 } 	/* switch (stype) */

	 start = 0.0;
	}	/* else  (!filuse) */

	stype = 0;		/* reset stype for next stim statement */
        orient = 0; 
        tfreq = 0; 
        drift = 0; 
        mask = 0;
        stimchan = 0;

	break;	/* case STIM */


  }	/* switch  (param->type) */
}

#endif

/*----------------------------------------------------*/

void sortstimfile(void)

/* Sort the stimulus file after all stimulus events 
   have been put into it.

   This routine was allows stimuli to be defined out of sequence.
This is essential because a stimulus consists of "on" and "off"
times.  When 2 stimuli are defined that overlap in time, the
"off" time of one stimulus may be defined before the "on" time of the
next -- but the "on" times may in fact coincide.

*/

{
   static char ftemplate[] = "stimXXXXXX";
   static char stemp[30],tbuf[120];

  if (makestim) {  
    strcpy (stemp, ftemplate);
    mkstemp (stemp);
    sprintf (tbuf,"sort -g -k 11 -k 1 -k 2 -k 3 -k 4 -k 5 -k 7 -k 8 -k 9 %s > %s && mv %s %s",
					stimfile,stemp,stemp,stimfile);
    system (tbuf);
  }
}

/*----------------------------------------------------*/

void vplot(void)

/* read in data from "plot V[]"  with default max min
*/

{
  datum d1,d2,d3,d4,d5;
  int narg,type,pmod,pval,cagparm;
  Symbol *param;
  

 param  = (Symbol *)*pc++;
 narg = (long int)*pc++;

 if (narg > 3) checknum (d4 = popm());
 else d4.val = NULLVAL;
 if (narg > 2) checknum (d3 = popm());
 else d3.val = NULLVAL;
 if (narg > 1) checknum(d2 = popm());
 else d2.val = NULLVAL;
 if (narg > 0) d1 = popm();			/* can't use getnod() here */
 else d1.val = NULLVAL;				/* because of S: below */
  
#ifndef XSTIMM
 if (++plotnum >= PLOTNODSIZ) plotnum = PLOTNODSIZ-1;
 if (narg>0) {
   if ((((long int)param)<2000) ||
      ((((long int)param)>2000) && (param->type!=S))) {
       plotnod[plotnum].cnod1 = (int)d1.val;
       plotnod[plotnum].cnod2 = (int)d2.val;
       plotnod[plotnum].cnod3 = (int)d3.val;
       plotnod[plotnum].cnod4 = (int)d4.val;
     }
 }
 pval = pmod = 0;

 if ((type=(long int)param) < 2000) { /* if Ca (G) parameter */

   cagparm = (long int)*pc++;
   switch (type) {
   case G:
	switch(cagparm) {

	  case 0:  pmod = NRECG0; break;	/* record conductance */
	  case 1: d5 = popm();
		   pval = int(d5.val);
		   pmod = NRECG;		/* record state fraction */
		   break;
	  case VVREV:pmod = NRECGV; break;
	  case I:   pmod = NRECGI; break;
	  case M:   pmod = NRECGM; break;
	  case H:   pmod = NRECGH; break;
	  case CA:  pmod = NRECGC; break;
	}
	break;

   case CA:
	switch(cagparm) {

	  case 0:  pmod = CREC_CACONC; 
		   pval = 1; break;	/* record cai ([Ca] inside) */
	  case 1:  d5 = popm();
		   pval = int(d5.val);
		   pmod = CREC_CACONC;		/* record cai from Ca shell */
		   break;
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
	switch(cagparm) {

	  case 1:  d5 = popm();
		   pval = int(d5.val);
		   pmod = CREC_CABUF;		/* record buf from Ca shell */
		   break;
	default:  break;
	}
	break;
   case CABUFB:
	switch(cagparm) {

	  case 1:  d5 = popm();
		   pval = int(d5.val);
		   pmod = CREC_CABUFB;		/* record buf bound from Ca shell */
		   break;
	default:  break;
	}
	break;
   default: break;
   }
 }
 else 
  switch (param->type) {
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
	pmod  = NRECA0 + (param->type - FA0);
	break;
   case G0: case G1: case G2: case G3:
	pmod  = NRECG0 + (param->type - G0); 
	break;
   case S:
 	plotnod[plotnum].spnt  = (char *)d1.sym;
	pmod  = SREC;
	break;
   case CABLE:
	pmod  = CABLE;
 	plotnod[plotnum].cnod2 = (int)(d2.val * CABLFRAC);
	break;
   case GLU:
   case AMPA:
   case KAINATE:
   case NMDA:
   case CNQX:
   case GABA:
   case PTX:
   case BIC:
   case GLY:
   case CHRC:
   case STRY:
   case CGMP:
   case CAMP: pmod  = param->type - GLU + CREC_GLU; break;
   case VAR:
   case FUNCTION: if (interp) varcopyu();

#ifndef XSTIM
	plotfunc(param,plotnum);
#endif
	pmod = FREC;
	break;
 }
 plotnod[plotnum].pval  = pval;
 plotnod[plotnum].pmod  = pmod;

 if (strcmp(gfrname(),rootframe))		 /* if diff than root frame */
   strcpy(plotnod[plotnum].plframe,gfrname());  /* copy current frame name */
#endif
}

/*------------------------------------------------------*/

void vplotm(void)

/* read in max, min data from "plot V() max min" statement
*/

{
  datum d1,d2;

 checknum(d2 = popm());
 checknum(d1 = popm());
 
#ifndef XSTIMM
 if (d1.val < d2.val) { double t; t = d1.val; d1.val=d2.val; d2.val=t; }
 plotnod[plotnum].pymax  = d1.val;
 plotnod[plotnum].pymin  = d2.val;
#endif
}

/*------------------------------------------------------*/

void xrecord(void)

/* allow user to record from nodes directly */

{
    datum d1,d2,d3,d4,d5;
    Symbol *param;
    int nod1, nod2, nod3, nod4;
    int narg, caparm, pmod, pval, type;
    double record(int cnod1, int cnod2, int cnod3, int cnod4, int pmod,int pval);
    double reccable(int elemnum, double fracdist);

#ifdef DEBUG
 ncfprintf (stderr,"xrecord\n");
#endif

 param  = (Symbol *)*pc++;
 narg = (long int) *pc++;

 if (narg > 3) d4 = popm();			/* can't use getnod() */
 else d4.val = NULLVAL;			/* because of d2.val below */
 if (narg > 2) d3 = popm();
 else d3.val = NULLVAL;
 if (narg > 1) d2 = popm();
 else d2.val = NULLVAL;
 checknum(d1 = popm());

 pmod = pval = 0;
 if ((type=(long int)param) < 2000) {

   caparm = (long int)*pc++;
   switch (type) {
   case G:
	switch(caparm) {

	  case 0:  pmod = NRECG; break;		/* record conductance */
	  case 1: d5 = popm();
		   pval = int(d5.val);
		   pmod = NRECG;		/* record state fraction */
		   break;
	  case VVREV:pmod = NRECGV; break;
	  case I:   pmod = NRECGI; break;
	  case M:   pmod = NRECGM; break;
	  case H:   pmod = NRECGH; break;
	  case CA:  pmod = NRECGC; break;
	}
	break;

   case CA:
	switch(caparm) {

	  case 0:  pmod = CREC_CACONC; break;	/* record cao ([Ca] outside) */
	  case 1: d5 = popm();
                   pval = int(d5.val);
		   pmod = CREC_CACONC;		/* record cai from Ca shell */
		   break;
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
	switch(caparm) {
	  case 1: d5 = popm();
                   pval = int(d5.val);
		   pmod = CREC_CABUF;		/* record buf from Ca shell */
		   break;
	  default: break;	
	}
	break;
   case CABUFB:
	switch(caparm) {
	  case 1: d5 = popm();
                   pval = int(d5.val);
		   pmod = CREC_CABUFB;		/* record bound buf from Ca shell */
		   break;
	  default: break;	
	}
	break;
   default: break;
   }
 }
 else 
 switch (param->type) {
   case V:
	pmod = VREC;
	break;
   case I:
	pmod = IREC;
	break;
   case IM:
	pmod = MREC;
	break;
   case L:
	pmod = LREC;
	break;
   case VM:
	pmod = WREC;
	break;
   case FA0: case FA1: case FA2: case FA3: case FA4: case FA8: case FA9:
   case FH0: case FH1: case FH2: case FH3: case FH4:
   case FB0: case FB1: case FB2: case FB3: case FB4:
   case FC0: case FC1: case FC2: case FC3: case FC4: case FC9:
   case G0: case G1: case G2: 
	pmod = NRECA0 + (param->type - FA0);
	break;
   case CABLE:
	pmod = CABLE;
	break;

   case GLU:
   case AMPA:
   case KAINATE:
   case NMDA:
   case CNQX:
   case GABA:
   case PTX:
   case BIC:
   case GLY:
   case CHRC:
   case STRY:
   case CGMP:
   case CAMP:
	pmod = param->type - GLU + CREC_GLU;
	break;
  }
  if (pmod==CABLE) {
	double fracdist;
	int elemnum;

     fracdist = d2.val;
     elemnum = (int)d1.val;
     d1.val = reccable(elemnum,fracdist);
   }
   else {
#ifndef XSTIMM
     nod1 = (int)d1.val;
     nod2 = (int)d2.val;
     nod3 = (int)d3.val;
     nod4 = (int)d4.val;
 #ifndef XSTIM
     d1.val = record (nod1,nod2,nod3,nod4,pmod,pval);
 #else
     d1.val = 0.0;
 #endif
#else
     d1.val = 0.0;
#endif
   }
   d1.vtype = NUMBER;
   d1.type = VAR;
   pushm(d1);
}

/*------------------------------------------------------*/

void setcmap (int ncmap)

{
  cmap = ncmap;
}

/*------------------------------------------------------*/

static double xcalib=.95,ycalib=.05;

void disp_calib (double xcalib,double ycalib,double cline,double dsize,int dcolor)

{
if (!makestim) {
   if ((disp & DISP) || (disp & DMOVIE)) {	/* if "disp" variable has been set */
       if (cline>0) {
 	 if (calibcolor != NULLVAL) dcolor = calibcolor;
	 drcalib (xcalib,ycalib,cline,dsize,dcolor);
       }
    }
  }
}
/*------------------------------------------------------*/

void display_calib (double cline,int dcolor)
{
   double size;
   disp_calib (xcalib,ycalib,cline,size=1.0,dcolor);
}

/*------------------------------------------------------*/

void set_disp_rot ( double xrot,  double yrot,  double zrot, 
	     	    double dxcent,double dycent,double dzcent,
		    double rxcent,double rycent,double rzcent,double dsize) 

{
   if (!makestim) {
   setrot(xrot,yrot,zrot,-dxcent,-dycent,-dzcent,rxcent,rycent,rzcent,dsize); 

   if (disp_ray) {				     /* move camera, not objects */
     initray(xrot,yrot,zrot,-dxcent,-dycent,-dzcent,
			rxcent,rycent,rzcent,dsize); /* initialize ray-tracer */
   }
 }
}

/*------------------------------------------------------*/

void display_rot ( double xrot,  double yrot,  double zrot)

{
  if (dxcent==LARGENODE && dycent==LARGENODE) {
      dxcent = dycent = dzcent = 0;
  }

  set_disp_rot (-xrot, yrot, zrot, dxcent, dycent, dzcent,
		    rxcent, rycent, rzcent, dsize);
}

/*------------------------------------------------------*/

void display_center ( double xcent, double ycent, double zcent)

{
  dxcent = xcent;
  dycent = ycent;
  dzcent = zcent;
}

/*------------------------------------------------------*/

static int dcmap=0;

void set_display (int elemtype, int disptype, int narg2, 
	     int n1a,int n1b,int n1c,int n1d,
	     int n2a,int n2b,int n2c,int n2d,
	     int na, int nb, int nc, int nd, int exceptype,
	     double zrange1, double zrange2,
	     int dcolor, Symbol *vpen, double(*vpenn)(int, int), 
	     int cmap, double vmax,double vmin, 
	     double dscale, int hide, int excl, double stime)
{
 if (!makestim) {
    if (n1a<0) n1a = NULLVAL;
    if (n1b<0) n1b = NULLVAL;
    if (n1c<0) n1c = NULLVAL;
    if (n1d<0) n1d = NULLVAL;

    if (n2a<0) n2a = NULLVAL;
    if (n2b<0) n2b = NULLVAL;
    if (n2c<0) n2c = NULLVAL;
    if (n2d<0) n2d = NULLVAL;

    if (na<0) na = NULLVAL;
    if (nb<0) nb = NULLVAL;
    if (nc<0) nc = NULLVAL;
    if (nd<0) nd = NULLVAL;

   set_icons();				/* set user-defined neural elem icons */

   if (cmap==0) cmap = dcmap;

   if (disp & DSTIM) 
    if (disptype == STIM) {
        //if (stime<=0) stime = simtime;
	dispstim(stime,dscale,cmap,vmax,vmin);
        stime=0;
    } 
  if (elemtype==COMPS && (disp & (DCOMP|DCONN))) {
  //if (disp & (DCOMP|DCONN)) {
      if (hide) hidstart();
      initcomp();         /* make compartments, don't erase them */
#ifdef DEBUG
      if (debug & NOCOND) nocond = 1;
#endif
      if (!nocond) condense();
      switch(disptype) {
       case MATCHING:
       case CONNECT: dcomp (n1a,n1b,n1c,n1d,elemtype,exceptype,
					na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,vmax,vmin,
					dscale,hide,excl);

/*	 else  dcomp2 (n1a,n1b,n1c,n1d,n2a,n2b,n2c,n2d,elemtype,exceptype,
					na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,dscale,hide); */
		 break;

      }  /* switch */
    }  /* elemtype==COMPS) */

   if ((elemtype==NODE) && (disp & DNODE)) {
        ncdrnod (n1a,n1b,n1c,n1d,elemtype,exceptype,
					na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,dscale,0);
   }

   if (elemtype==LABEL && !nolabels) {
        ncdrnod (n1a,n1b,n1c,n1d,elemtype,exceptype,
					na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,dscale,1);
   }

   if ((disp & DISP) || (disp & DMOVIE)) {	/* if "disp" variable has been set */

    if (hide) hidstart();

    if (elemtype!=NODE)		/* elemtype is CABLE, SPHERE, SYNAPSE, etc. */
    switch (disptype) {

     case MATCHING:
     case CONNECT: if (narg2==0) ncdisp (n1a,n1b,n1c,n1d,elemtype,exceptype,
					na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,vmax,vmin,
					dscale,hide,excl,cmap);
		 else ncdispc (n1a,n1b,n1c,n1d,n2a,n2b,n2c,n2d,elemtype,
					exceptype,na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,vmax,vmin,
					dscale,hide,cmap);
		 break;
     case RANGE:    ncdispn (n1a,n1b,n1c,n1d,n2a,n2b,n2c,n2d,elemtype,
					exceptype,na,nb,nc,nd,
					zrange1,zrange2,
					dcolor,vpen,vpenn,vmax,vmin,
					dscale,hide,excl,cmap);
		 break;
     case ELEMENT: ncdispe (n1a,	zrange1,zrange2,
					dcolor,vpen,vpenn,vmax,vmin,
					dscale,hide,cmap);
		 break; 
    }
   }		/* if (disp & DISP) */
 }
}

/*------------------------------------------------------*/

void set_display (int elemtype, int disptype, int narg2,
             node *nd1, node *nd2, node *nexc, int exceptype,
             double zrange1, double zrange2,
             int dcolor, Symbol *vpen, double(*vpenn)(int, int),
             int cmap, double vmax,double vmin,
             double dscale, int hide, int excl, double stime)

{
  if (nd1!=NULL && nd2 !=NULL && nexc != NULL)
     set_display (elemtype, disptype, narg2,
             nd1->nodenm1, nd1->nodenm2, nd1->nodenm3, nd1->nodenm4,
             nd2->nodenm1, nd2->nodenm2, nd2->nodenm3, nd2->nodenm4,
             nexc->nodenm1, nexc->nodenm2, nexc->nodenm3, nexc->nodenm4, exceptype,
             zrange1, zrange2,
             dcolor, vpen, vpenn,
             cmap, vmax,vmin,
             dscale, hide, excl, stime);
}


/*------------------------------------------------------*/

void dispnod(void)
{
  datum d1,d2,d3,d4,d[MAXNODIM];
  static int narg1,narg2;
  static int n1a,n1b,n1c,n1d,n2a,n2b,n2c,n2d;
  static int na,nb,nc,nd;
  Symbol *param,nulparm;
  static int disptype = 0, dcolor = NULLVAL;
  static double vmax=NULLVAL,vmin=NULLVAL;
  static int excl=0, hide=0;
  static double cline=0,dscale=1.0,stime=0;
  static int elemtype=0;
  static int exceptype=0,except=0;
  static double zrange1=-LARGENODE,zrange2=LARGENODE;
  static Symbol *vpen=NULL;
  static double(*vpenn)(int,int)=NULL;

  param = (Symbol *)*pc++;
  if (param==0) {
  	nulparm.type = CONNECT;
	param = &nulparm;
  }
#ifdef XMOD
  if (interp) varcopy();
  if (setdebug) debug = setdebug;
  if (setdebgz) debugz = setdebgz;
  if (!runyet || elempnt) {
	findconnect(1);
  }
#endif

  switch (param->type) {

  case ONLY:  excl = 1;
	      break;

  case MATCHING:
  case ELEMENT:
  case RANGE:
  case CONNECT:
 
  disptype = param->type; 
  narg2 = (long int)*pc++;
  if (narg2) {
    pc--;
    getnod(d);				/* args for node2 are first */
    					/*  because they're on top of stack */
    n2a = (int)d[0].val;
    n2b = (int)d[1].val;
    n2c = (int)d[2].val;
    n2d = (int)d[3].val;
  }
  else n2a=n2b=n2c=n2d=NULLVAL;

  narg1 = (long int)*pc++;
  if (narg1) {
    pc--;
    getnod(d);
    if (except) {
       na = (int)d[0].val;
       nb = (int)d[1].val;
       nc = (int)d[2].val;
       nd = (int)d[3].val;
    }
    else {
       n1a = (int)d[0].val;
       n1b = (int)d[1].val;
       n1c = (int)d[2].val;
       n1d = (int)d[3].val;
       if (n1a<0) n1a = NULLVAL;
       if (n1b<0) n1b = NULLVAL;
       if (n1c<0) n1c = NULLVAL;
       if (n1d<0) n1d = NULLVAL;
       na= nb= nc=nd= NULLVAL;		/*  clear except node */
    }
  }
  else {				/* if node not specified: */
    if (except) na= nb= nc= nd= NULLVAL; /*  clear except node */
    else n1a=n1b=n1c=n1d=NULLVAL;   	/*  clear n1 when not except */
  }

  break;

  case COMPS:
  case CABLE:
  case SPHERE:
  case SYNAPSE:
  case CHAN:
  case CACOMP:
  case ROD:
  case CONE:
  case CHR:
  case VTRANSDUCER:
  case ITRANSDUCER:
  case ELECTRODE:
  case GJ:
  case LOAD:
  case RESISTOR:
  case CAP:
  case GNDCAP:
  case BATT:
  case GNDBATT:
  case BUF:
  case NODE:
    if (except) exceptype = param->type;
    else         elemtype = param->type;
    break;

  case EXCEPT:  except = 1;
	      break;

  case XROT:				/* get rotation for display */
    checknum(d1 = popm());
    xrot = -d1.val;
    break;

  case YROT:
    checknum(d1 = popm());
    yrot = d1.val;
    break;

  case ZROT:
    checknum(d1 = popm());
    zrot = d1.val;
    break;

  case SIZE:				/* get size of display window */
    checknum(d1 = popm());
    dsize = d1.val;
    break;

  case DSCALE:				/* get display scale */
    checknum(d1 = popm());
    dscale = d1.val;
    break;

  case CMAP:				/* get color map (0-3) */
    checknum(d1 = popm());
    if (d1.type==ARRAY) {
      cmap=setcmap(d1);
    }
    else cmap = int(d1.val);
    break;

  case NEWPAGE:				/* make new page for movie */
     gpage();
    break;

  case COLOR:				/* get color */
    narg1 = (long int)*pc++;
    if (narg1) { 
       d3 = popm();			/* get possible "max, min" */
       d2 = popm();
       vmax = d2.val;
       vmin = d3.val;
    }
    checknum(d1 = popm());
    dcolor = (int)d1.val;
    break;

  case VPEN:
     vpen = (Symbol *)*(pc++);
     if (interp) varcopyu();
     break;

  case CALIBLIN:			/* get calib line */
    narg1 = (long int)*pc++;
    if (narg1) {
      if (narg1 > 2) d3 = popm();		/* get possible "loc" */
      if (narg1 > 1) d2 = popm();
      checknum(d1 = popm());
      xcalib = d1.val; 
      ycalib = d2.val; 
    }
    checknum(d1 = popm());
    cline = d1.val;
    break;

  case STIM:				/* get stimulus time  */
    narg1 = (long int)*pc++;
    if (narg1==1) {
       checknum(d1 = popm());
       stime = d1.val;
    }
    else stime = simtime;
    disptype = param->type;
    break;

  case HIDE:				/* make picture with hidden lines */
    hide = 1;
    break;

  case CENTER:
    narg1 = (long int)*pc++;
    if (narg1) {
      if (narg1 > 2) d3 = popm();
      else d3.val = 0;
      if (narg1 > 1) d2 = popm();
      else d2.val = 0;
      checknum(d1 = popm());
    }
    dxcent = d1.val;
    dycent = d2.val;
    dzcent = d3.val;
    break;

  case Z:
    narg1 = (long int)*pc++;
    if (narg1) {
       d2 = popm();			/* get possible "max, min" */
       d1 = popm();
    }
    zrange1 = d2.val;
    zrange2 = d1.val;
    break;

  case RMOVE:
    narg1 = (long int)*pc++;
    if (narg1) {
      if (narg1 > 2) d3 = popm();
      else d3.val = 0;
      if (narg1 > 1) d2 = popm();
      else d2.val = 0;
      checknum(d1 = popm());
    }
    rxcent = d1.val;
    rycent = d2.val;
    rzcent = d3.val;
    break;

  case WINDOW:
      checknum(d4 = popm());
      checknum(d3 = popm());
      checknum(d2 = popm());
      checknum(d1 = popm());
    break;
 
  case DISPLAY:
	if (!disptype) {
	  dcmap = cmap;		/* if cmap on line by itself, save */
	  calibcolor = dcolor;
	  return;
	}

#ifdef XMOD 
	set_disp_rot(xrot,  yrot,  zrot, dxcent,dycent,dzcent,
	         rxcent,rycent,rzcent,dsize); 

    set_display (elemtype, disptype, narg2, 
	     n1a,n1b,n1c,n1d,
	     n2a,n2b,n2c,n2d,
	     na, nb, nc, nd, exceptype,
	     zrange1, zrange2,
	     dcolor, vpen, vpenn,
	     cmap, vmax, vmin,
	     dscale, hide, excl, stime);

    if (cline > 0) {
      disp_calib (xcalib,ycalib,cline,dsize,dcolor);
      cline = 0;
    }

    disptype = 0;
    elemtype = 0;
    except = exceptype = 0;
    rxcent=rycent=rzcent=0;
    na=nb=nc=nd=NULLVAL;
    dscale = 1.0;
    dcolor = NULLVAL;
    calibcolor = NULLVAL;
    vmax=vmin=NULLVAL;
    vpen=NULL;
    excl = 0;
    cmap = 0;
    if (hide) hidstop();

#endif		/* XMOD */
  break;     /* case DISPLAY */
 } 	   /* switch (param) */
}


/*------------------------------------------------------*/

void xtransf (datum &x, datum &y, datum &z)

{
   datum d,d2;
   double *p,tx=0,ty=0,tz=0;
   Symbol sp;
   char tnam[40];
   char xbuf[150];


   sprintf (tnam,"transf x,y array\n");
   sp.name = tnam;
   sp.type = ARRAY;

   d2.val = 3;
   d2.type = CONST;
   d2.vtype = NUMBER;
   pushm(d2);
   p=darr2(&sp,1);
   if (p==NULL) {
       sprintf (xbuf,"# Can't allocate 'transf()' array\n");
       execerror ("stopping,",xbuf);
   }
#ifdef XMOD
   transf (x.val, y.val, z.val, &tx, &ty, &tz, ncmat);
#endif
   p[0] = tx;
   p[1] = ty;
   p[2] = tz;
   d.arrp = sp.arrp;
   d.type = CONST;
   d.vtype=NUMBER;
   pushm(d);
}

/*------------------------------------------------------*/


#define NUMDIM MAXNODIM			/* number of node dimensions */

        Symbol *eparam;

void foreacode(void)			/* foreach statement */

/* execute a set of statements for one of three conditions (options):

    1) For every element, or each element of a certain type (e.g. "cable")

    2) For every element, or each element of a certain type (e.g. "cable")
         that matches a node specification.  

    3) for each node that matches a node specification.
*/
         
{
	datum d;
	Inst *savepc = pc;
	int i,narg,narg2,match,efilt,fl3d;
	int arg[NUMDIM];
	int val[NUMDIM];
	nodeint node1,node2,node3,node4;
	double radius,x,y,z;
	pdatum p1={0},p2[NUMDIM]={0};
        node *npnt,*tnpnt,*npnt2=NULL;
	elem *epnt,*tepnt,*lepnt;
   	datum dn[NUMDIM];

    eparam = (Symbol *)*pc++;		/* get element type */ 
    efilt = 0;
 /* if (eparam) ncfprintf (stderr,"foreacode efilt name '%s'\n",eparam->name); */
    if (eparam) switch (eparam->type) {

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
     case LOAD:
     case RESISTOR:
     case CAP:
     case GNDCAP:
     case BATT:
     case GNDBATT:
     case BUF:
     case CACOMP:
			efilt = 2;
			break;	

     default:		efilt = 0;
			break;	

    }

    if (efilt) {
       getvar(&p1);  /* var pointer for element */
       if (*p1.type==UNDEF) *p1.type = VAR;
       *p1.vtype = NUMBER;
    }
    
    narg = (long int)*pc++;		/* number of snode args to expect */

   if (narg) {
    if (narg>NUMDIM) narg = NUMDIM;
    for (i=0; i<NUMDIM; i++) {
       p2[i].val = (double *)NULL;
       val[i]  = NULLVAL;
       arg[i]  = (long int)*pc++;
    }
#ifndef XSTIMM
    execute (savepc+9);			/* evaluate node descriptors */
#endif

   narg2 = (long int)*(savepc+6);	/* get number of "within" node args */
   fl3d = (narg2<0);			/* if neg, 3d calculation of dist */
   if (fl3d) narg2 = -narg2;

   if (narg2 > 0) {			/* if there is a "within" clause */

     if (narg2 > 3) {
	checknum(d=popm()); node4 = (nodeint)d.val;
     }
     else node4 = NULLVAL;
     if (narg2 > 2) {
	checknum(d=popm()); node3 = (nodeint)d.val;
     }
     else node3 = NULLVAL;
     if (narg2 > 1) {
	checknum(d=popm()); node2 = (nodeint)d.val;
     }
     else node2 = NULLVAL;
     checknum(d=popm()); node1 = (nodeint)d.val;
     if (!(npnt2=findnode(node1,node2,node3,node4,"Foreach"))) {
        execerror ("Missing node; ","stopping...");
     }
     checknum(d = popm());		/* get radial distance */	
     radius = d.val;
   }

   for (i=narg-1; i>=0; i--) {
         if (arg[i]) {
	    getvar(&p2[i]);  		/* var pointer for node dim */
            if (*p2[i].type==UNDEF) *p2[i].type = VAR;
            *p2[i].vtype = NUMBER;
	 }
	 else {
           checknum(d = popm()); 	/*else get value of node dim */
	   val[i] = (int)d.val;
         }
   }
  }   /* if (narg) */


#ifndef XSTIMM
 if (efilt) {					/* search all elems */
   for (epnt=elempnt; epnt; epnt=tepnt) {
     tepnt = epnt->next;		/* pre-increment in case deleted */
     lepnt = epnt->last;
     if (efilt==2) if (epnt->ctype!=eparam->type) continue;

     npnt = (node *)NULL;
     if (narg) {
       npnt = epnt->nodp1;		/* test first node */
       if (npnt) 
        for (match=1,i=0; i<narg; i++) {
 	 if (!arg[i]) 			/* if dimension was given a value */
	  switch (i) {
	   case 0: if ((val[i]>=0) &&(val[i]!=npnt->nodenm1)) match = 0; break; 
	   case 1: if ((val[i]>=0) &&(val[i]!=npnt->nodenm2)) match = 0; break; 
	   case 2: if ((val[i]>=0) &&(val[i]!=npnt->nodenm3)) match = 0; break; 
	   case 3: if ((val[i]>=0) &&(val[i]!=npnt->nodenm4)) match = 0; break; 
	  }
        }
       if (narg2) {
	  if (fl3d) { 
	    if    (dist3d(npnt,npnt2) > radius) match = 0;
	  }
	  else if (dist2d(npnt,npnt2) > radius) match = 0;
       }
       if (match) {			/* set first node */
       *p1.val = epnt->elnum;		/* set the element number */
       if (npnt)
        for (i=0; i<narg; i++) {
	 if (arg[i]) {			/* If dimension was given blank var, */
	   switch (i) {			/*   set variable for this iteration */
	   case 0: *p2[i].val = npnt->nodenm1; break; 
	   case 1: *p2[i].val = npnt->nodenm2; break; 
	   case 2: *p2[i].val = npnt->nodenm3; break; 
	   case 3: *p2[i].val = npnt->nodenm4; break; 
	   }	
         }
       }		/* run loop code if elem and node match requested */

       if (!forexec(savepc+(long int)*(savepc+7))) break;	

      }    /* if (match) */

      else {					/* test second node */
      if (lepnt && (lepnt->next==epnt)) {	/* if not deleted yet */
       npnt = epnt->nodp2;		/* test second node */
       if (npnt) 
        for (match=1,i=0; i<narg; i++) {
	 if (!arg[i]) 			/* if dimension was given a value */
	  switch (i) {
	   case 0: if ((val[i]>=0) &&(val[i]!=npnt->nodenm1)) match = 0; break; 
	   case 1: if ((val[i]>=0) &&(val[i]!=npnt->nodenm2)) match = 0; break; 
	   case 2: if ((val[i]>=0) &&(val[i]!=npnt->nodenm3)) match = 0; break; 
	   case 3: if ((val[i]>=0) &&(val[i]!=npnt->nodenm4)) match = 0; break; 
	   }	
        }
        if (narg2) {
	  if (fl3d) {
	    if      (dist3d(npnt,npnt2) > radius) match = 0;
	  }
	  else if   (dist2d(npnt,npnt2) > radius) match = 0;
        }
        if (!match) continue;

       if (match) {			/* set second node */
       *p1.val = epnt->elnum;		/* set the element number */
       if (npnt)
        for (i=0; i<narg; i++) {
	 if (arg[i]) {			/* If dimension was given blank var, */
	   switch (i) {			/*   set variable for this iteration */
	   case 0: *p2[i].val = npnt->nodenm1; break; 
	   case 1: *p2[i].val = npnt->nodenm2; break; 
	   case 2: *p2[i].val = npnt->nodenm3; break; 
	   case 3: *p2[i].val = npnt->nodenm4; break; 
	   }	
         }
       }		/* run loop code if elem and node match requested */

       if (!forexec(savepc+(long int)*(savepc+7))) break;	

       }    /* if (match) */
      }    /* if (epnt) (second node match) */
     }    /* else second node */
    }    /* if (narg) */

    else {		/* run loop code if only elem match requested */

        if (narg2) {
	  if (fl3d) {
	    if    (endist3d(epnt,npnt2) > radius) continue;
	  }
	  else if (endist2d(epnt,npnt2) > radius) continue;
        }

        *p1.val = epnt->elnum;		/* set the element number */
        if (!forexec(savepc+(long int)*(savepc+7))) break;	
    }

  }  /* for (epnt;;) */
 }  /* if (efilt) */



 else 		/* no element, just node to match */
   for (npnt=nodepnt; npnt; npnt=tnpnt) {	/* search all nodes */
      tnpnt = npnt->next;
      for (match=1,i=0; i<narg; i++) {
	if (!arg[i]) 			/* if dimension was given a value */
	  switch (i) {
	   case 0: if ((val[i]>=0) &&(val[i]!=npnt->nodenm1)) match = 0; break; 
	   case 1: if ((val[i]>=0) &&(val[i]!=npnt->nodenm2)) match = 0; break; 
	   case 2: if ((val[i]>=0) &&(val[i]!=npnt->nodenm3)) match = 0; break; 
	   case 3: if ((val[i]>=0) &&(val[i]!=npnt->nodenm4)) match = 0; break; 
	  }	
     }
     if (narg2) {
	  if (fl3d) {
	    if    (dist3d(npnt,npnt2) > radius) match = 0;
	  }
	  else if (dist2d(npnt,npnt2) > radius) match = 0;
     }
     if (!match) continue;

     for (i=0; i<narg; i++) {
	if (arg[i])			/* If dimension was given blank var, */
	  switch (i) {			/*   set variable for this iteration */
	  case 0: *p2[i].val = npnt->nodenm1; break; 
	  case 1: *p2[i].val = npnt->nodenm2; break; 
	  case 2: *p2[i].val = npnt->nodenm3; break; 
	  case 3: *p2[i].val = npnt->nodenm4; break; 
	  }	
     }

    if (!forexec(savepc+(long int)*(savepc+7))) break;	

    }  /* for (npnt;;) */

#endif
  if (!returning) {
	pc = savepc+(long int)*(savepc+8);	/* next stmt */
  }
}

/*------------------------------------------------------*/

int forexec(Inst *body)

{
     int savestackp = stackp;

  execute(body);		/* body of element, node loop */
  if (stopping) {
	if (returning)
		return(0);
	stopping = 0;
	if (breaking) {
		break_fixup(savestackp,stackp);
		breaking = 0;
		return (0);
	}
		continuing = 0;
   }
   return (1);
}

/*------------------------------------------------------*/

static double elim_xmin = -LARGENODE;	/* Window limit for elements */
static double elim_xmax =  LARGENODE;	/*  set by elimit() */
static double elim_ymin = -LARGENODE;
static double elim_ymax =  LARGENODE;
static double elim_zmin = -LARGENODE;
static double elim_zmax =  LARGENODE;

int eclip (node *npt1, node *npt2)

/* Establish a window beyond which elements are clipped.
   If either node (one in the case of single-node elements)
   are outside the window, return 1, otherwise return 0.
*/

{
   int c1,c2;

   if (!npt1 && !npt2) return 1;  		/* pointers not found */
   c1 = c2 = 0;
   if (npt1) {
      if      (npt1->xloc < elim_xmin) { c1 = 1; }
      else if (npt1->xloc > elim_xmax) { c1 = 1; }
      if      (npt1->yloc < elim_ymin) { c1 = 1; }
      else if (npt1->yloc > elim_ymax) { c1 = 1; }
      if      (npt1->zloc < elim_zmin) { c1 = 1; }
      else if (npt1->zloc > elim_zmax) { c1 = 1; }
   }
   if (npt2) {
      if      (npt2->xloc < elim_xmin) { c2 = 1; }
      else if (npt2->xloc > elim_xmax) { c2 = 1; }
      if      (npt2->yloc < elim_ymin) { c2 = 1; }
      else if (npt2->yloc > elim_ymax) { c2 = 1; }
      if      (npt2->zloc < elim_zmin) { c2 = 1; }
      else if (npt2->zloc > elim_zmax) { c2 = 1; }
   }
  return (c1||c2);
}

/*------------------------------------------------------*/

void elimit (int param, double emax, double emin)

{
  switch (param) {
 	case X: 
	  	if (emax < emin) 
			{ double t; t = emax; emax=emin; emin=t; }
		elim_xmax = emax;
		elim_xmin = emin;
		break;
  	case Y:
	  	if (emax < emin) 
			{ double t; t = emax; emax=emin; emin=t; }
		elim_ymax = emax;
		elim_ymin = emin;
		break;
  	case Z:
	  	if (emax < emin) 
			{ double t; t = emax; emax=emin; emin=t; }
		elim_zmax = emax;
		elim_zmin = emin;
		break;
  }
}

/*------------------------------------------------------*/

void elimit (elem *epnt)
{
  findconnect(1);		  /* set pointers if not done */
  if (eclip (epnt->nodp1,epnt->nodp2)) {

	/* first free the nodes' pointers to this elem */
    if (epnt->nodp1) {
	unsetnodlst(epnt->nodp1,epnt);
	if (!epnt->nodp1->elemlst) delnode(epnt->nodp1);
    }
    if (epnt->nodp2) {
	unsetnodlst(epnt->nodp2,epnt);
	if (!epnt->nodp2->elemlst) delnode(epnt->nodp2);
    }
    delelem(epnt);                     /* erase the element */
  }
}

/*------------------------------------------------------*/

void elimit (int elnum)

{
   elem *epnt;

#ifndef XSTIMM
   	if (!(epnt=findelem(elnum))) {
          ncfprintf (stderr,"elimit: can't find element %d\n",elnum);
          execerror ("Missing element: ","stopping...");
    	}
    	else {				  /* first find nodes that connect */
	  elimit (epnt);
	}
#endif
}

/*------------------------------------------------------*/

void elimit (void)

/* Set "soft" window limits on element location.
   Delete elements entirely or partially outside window. 
*/

{
	datum d1,d2;
  	elem *epnt;
  	int elnum;
	Symbol *param;

  param = (Symbol *)*pc++;		/* get element type */ 
  switch (param->type) {
  	case X: 
  	case Y:
  	case Z:
		checknum(d2 = popm());
		checknum(d1 = popm());
  		elimit (param->type,d1.val,d2.val);
		break;
  	case ELEMENT: 
		checknum(d1 = popm());
  		elimit (param->type,d1.val,d2.val);
		break;
  } 
}

/*------------------------------------------------------*/

void findconnect (int donefl)

/* Find and/or make all nodes that connect to all elements.
   Do only once, except do always when donefl = 0.
*/

{
   elem *epnt;
   static int done = 0;

#ifndef XSTIMM
 if (elempnt==NULL) done = 0;			/* reset at beginning */
 if (!done || !donefl)
   for (epnt=elempnt; epnt; epnt=epnt->next) { /* check nodes for element*/
     checkelemnode(epnt);
   }
 done = 1;
#endif

/*
if (debug & NCBASE)
 for (npnt=nodepnt; npnt; npnt=npnt->next)
  ncfprintf (stderr,"node %d %d conn %d\n",npnt->nodenm1,npnt->nodenm2, 
		(npnt->elemlst ? npnt->elemlst->num : NULLVAL));

 for (epnt=elempnt; epnt; epnt=epnt->next)
  ncfprintf (stderr,"element %d node1 %d %d %d %d node2 %d %d %d\n",
		epnt->elnum,epnt->node1a,epnt->node1b,epnt->node1c,
		epnt->node1d,epnt->node2a,epnt->node2b);

	/* */
}

/*------------------------------------------------------*/

void checknod (elem *epnt)

/* Find the nodes that an element connects to;
   set up the element's node pointers.
   If node doesn't exist, create a new node.
   Don't set node's element pointer because element
   has not been fully defined yet.
*/

{
   node  *npnt,*newnod;
   int found;

#ifndef XSTIMM
    if (!epnt) return;
    found = 0;
    if (epnt->node1a != NULLVAL) {
      if ((npnt=
	findnode(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d,
			(char *)NULL))) {
	   epnt->nodp1 = npnt;			/* set node pointer */
/*	   setnodlst(npnt,epnt); 		/* set elem pointer */
	   found = 1;
      }
      if (! found) {				/* if not, then make it */
 	    newnod=maknod(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d);
/* 	    setnodlst(newnod=maknod(
	      epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d),epnt); /* */
	    epnt->nodp1 = newnod;		/* set node pointer */
      }
    }	/* if (epnt->node1a != NULLVAL) */

    found = 0;
    if (epnt->node2a == NULLVAL) return;	/* ignore if no second node */
      if ((npnt=
	findnode(epnt->node2a,epnt->node2b,epnt->node2c,epnt->node2d,
			(char *)NULL))) {
	   epnt->nodp2 = npnt;			/* set node pointer */
/*	   setnodlst(npnt,epnt);	 	/* set elem pointer */
	   found = 1;
      }
    if (! found) {				/* if not, then make it */
 	    newnod=maknod(epnt->node2a,epnt->node2b,epnt->node2c,epnt->node2d);
/* 	    setnodlst(newnod=maknod(
	       epnt->node2a,epnt->node2b,epnt->node2c,epnt->node2d),epnt); /* */
	    epnt->nodp2 = newnod;		/* set node pointer */
    }

/* ncfprintf (stderr,"checknod el %d nod1 %d nod2 %d\n",
			epnt->elnum,epnt->node1a,epnt->node2a); /* */

/* ncfprintf (stderr,"checknod elem %d  n1 %d n2 %d\n",
		epnt->elnum, epnt->nodp1,epnt->nodp2); /* */

#endif		/* XSTIMM */
}

/*------------------------------------------------------*/

void checkelemnode (elem *epnt)

/* 
   Find the nodes that connect to an element;
   for each node, set up the node's element pointer list,
   and set the element's 2 node pointers.
   If node doesn't exist, create a new node.
*/

{
   node  *npnt,*newnod;
   int found;

#ifndef XSTIMM
    if (!epnt) return;
    found = 0;
    if (epnt->node1a != NULLVAL) {
      if (epnt->node1d != NULLVAL) {
         if ((npnt= findnode(epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d, (char *)NULL))) {
	      epnt->nodp1 = npnt;			/* set node pointer */
	      setnodlst(npnt,epnt); 		/* set elem pointer */
	      found = 1;
         }
      } else {					// ignore node1d
         if ((npnt= findnode(epnt->node1a,epnt->node1b,epnt->node1c))) {
	      epnt->nodp1 = npnt;			/* set node pointer */
	      setnodlst(npnt,epnt); 		/* set elem pointer */
	      found = 1;
         }
      }
      if (! found) {				/* if not, then make it */
 	    setnodlst(newnod=maknod(
	     epnt->node1a,epnt->node1b,epnt->node1c,epnt->node1d),epnt);
	    epnt->nodp1 = newnod;		/* set node pointer */
      }
    }	/* if (epnt->node1a != NULLVAL) */

    found = 0;
    if (epnt->node2a == NULLVAL) return;	/* ignore if no second node */
      if (epnt->node2d != NULLVAL) {
         if ((npnt= findnode(epnt->node2a,epnt->node2b,epnt->node2c,epnt->node2d, (char*)NULL))) {
	      epnt->nodp2 = npnt;			/* set node pointer */
	      setnodlst(npnt,epnt);	 	/* set elem pointer */
	      found = 1;
         }
    } else {					// ignore node2d
         if ((npnt= findnode(epnt->node2a,epnt->node2b,epnt->node2c))) {
	      epnt->nodp2 = npnt;			/* set node pointer */
	      setnodlst(npnt,epnt);	 	/* set elem pointer */
	      found = 1;
         }
    }
    if (! found) {				/* if not, then make it */
 	    setnodlst(newnod=maknod(
	      epnt->node2a,epnt->node2b,epnt->node2c,epnt->node2d),epnt);
	    epnt->nodp2 = newnod;		/* set node pointer */
    }

/* ncfprintf (stderr,"checknod el %d nod1 %d nod2 %d\n",
			epnt->elnum,epnt->node1a,epnt->node2a); /* */

/* ncfprintf (stderr,"checknod elem %d  n1 %d n2 %d\n",
		epnt->elnum, epnt->nodp1,epnt->nodp2); /* */

#endif		/* XSTIMM */
}

/*------------------------------------------------------*/

void setnodlst(node *npnt, elem *epnt)
              
/* Add an element to a node's element list.  */

{
   int found;
   conlst *lpnt;

#ifndef XSTIMM
  if (!npnt) {
	ncfprintf (stderr,"setnodlst: invalid node pointer\n");
        return;
  }
  found = 0;
  for (lpnt=npnt->elemlst; lpnt; lpnt=lpnt->next) {
    if (lpnt->conpnt == (conn*)epnt) found = 1;
  }

  if (! found) 
    maklst(&npnt->elemlst,(conn *)epnt);
#endif
}

/*------------------------------------------------------*/

void unsetnodlst(node *npnt, elem *epnt)
              
/* Delete pointer to an element from a node's element list */

{
  conlst *lpnt, *olpnt, *tlpnt;

#ifndef XSTIMM
  if (!npnt) {
	ncfprintf (stderr,"unsetnodlst: invalid node pointer\n");
        return;
  }
  olpnt = (conlst *)NULL;
  for (lpnt=npnt->elemlst; lpnt;)  {
     if (lpnt->conpnt== (conn *)epnt) {		/* delete this one */
        if (olpnt) olpnt->next = lpnt->next;
        else     npnt->elemlst = lpnt->next;
        tlpnt = lpnt;
        lpnt = lpnt->next;
        efree (tlpnt);
     }
     else {
        olpnt = lpnt;
        lpnt = lpnt->next;			/* skip to next */
     }
  }
#endif
}

/*------------------------------------------------*/
    
void modrun(void)
{
  datum d1;
  Symbol *mtyp;

#ifdef DEBUG
if (debug & NCBASE) ncfprintf (stderr,"modrun start\n");
#endif

  mtyp = (Symbol *)*pc++;
  if (interp) varcopy();

  if (setdebug) debug  = setdebug;	/* "-y n" overrides "debug=n" */
  if (setdebgz) debugz = setdebgz;	/* "-z n" overrides "debugz=n" */

  switch (mtyp->type)           /* do pops */
    {
     case RUN:  break;
     case STEP:  checknum(d1=popm());
                break;
    }

#ifdef XMOD
  if (!runyet) {
      initchan();			/* set up sequential-state chans */
  }
  if (!runyet || elempnt) {
	findconnect(0);
  }
  switch (mtyp->type)
    {

     case RUN:  runsim(0.0,1);
		runyet = 0;
		break;

     case STEP: runsim(d1.val,0);
		break;
    
    }
#else
   switch (mtyp->type)           /* do pops */
     {
      case RUN:  
      case STEP:  
    		simtime = d1.val + simtime;  /* if XSTIM, set the time, at least */
    		timeptr->val = simtime;
                break;
     }
#endif   
}

/*------------------------------------------*/

char *prnode (int n1, int n2, int n3, int n4) 

{
   static char nbuf[40] = {0};
   static char undefd[] = "undef.hed";
   int test;

   nbuf[0] = 0;
   test = 0;
   if (n1 != NULLVAL) test |= 1;
   if (n2 != NULLVAL) test |= 2;
   if (n3 != NULLVAL) test |= 4;
   if (n4 != NULLVAL) test |= 8;

   switch (test) {

     default:
     case 0: 	sprintf (nbuf,"%s",undefd);
		break;
     case 1: 	sprintf (nbuf,"[%d]",n1);
		break;
     case 3: 	sprintf (nbuf,"[%d][%d]",n1,n2);
		break;
     case 7: 	sprintf (nbuf,"[%d][%d][%d]",n1,n2,n3);
		break;
     case 15: 	sprintf (nbuf,"[%d][%d][%d][%d]",n1,n2,n3,n4);
		break;
   }

   return (nbuf);
}

/*------------------------------------------*/

char *prnode (node *npnt) 

{
  if (npnt==NULL) return NULL;
  else return (prnode(npnt->nodenm1, npnt->nodenm2,  npnt->nodenm3, npnt->nodenm4));
}

/*------------------------------------------*/

