/* stimcode.c */

/* Old Version!  Use only for historical reference */

#include "nc.h"
#include "y.tab.h"
#include <stdio.h>
#include "ncelem.h"
#include "ncomp.h"
#include "stim.h"
#include "control.h"

int reccum = 0;
int recepnode = 0;                      /* receptor node from conn() */
int recepn2 = 0;
int recepn3 = 0;

/* interpreter innards: */

#define push(d) *stackp++ = (d)         /*  */
#define popm()  *--stackp               /* function pop still needed */

/* #define popm() pop()                 /*  */

#define NSTACK 256

extern char *progname;
extern char *runfile;
extern int lineno;
extern int stimhlp;
extern int stopping;
extern int continuing;
extern int breaking;
extern int returning;

extern  datum   *stackp;

extern  Inst    *progp;         /* next free spot for code generation */
extern  Inst    *pc;            /* program counter during execution */
extern  Inst    *progbase;      /* start of curent subprogram */

int curelnum=0;                 /* current elem being read in */
elem *elpnt;                    /* current elem pointer */
node *nodepnt;           	/* pointer to node list */
extern elem *elempnt;           /* pointer to elem list */

Symbol *timeptr=0;              /* set by init() */

int stimfflg = 0;

double cri=0;                           /* current RI, RM, RG */
double crm=0;
double crg=0;

int *larr[NLABELS]={0};                   /* label array */

FILE *stimout=stdout;           /* stimulus output file */

double getval();

node *maknod();
elem *makbr();

/* Remember to always change stimcode whenever
   modcode is changed.  The stimcode subroutines
   should be called usually in exactly the same way
   as their modcode counterparts (exception: stimulus
   subroutines) and they should always leave the 
   stack (popm) and program counter (pc) correct.
*/ 

/*------------------------------------------*/

initpl()        /* in "ncplot.c" */
{} 		   /* */

/*------------------------------------------*/

initrec()       /* in "ncmak.c" */
{}

initmak()       /* in "ncmak.c" */

/* make ready for new branch data to be read in
*/

{
reccum = 0;                     /* reset number of receptors */
initc(2,512);                       /* make two 256 * 256 arrays */
}

branch()

/* read in the elem number, and remember it
   in curelnum. */

{
  datum d1;

  d1 = popm();
}

double xelem()

/* return an element's number */

{
  return ((double)1);
}

xmod()

/* modify an element */

{
  datum d1;
  Symbol *param;
  double *val;
  short *typ;

  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case MODIFY: d1 = popm();
		 break;
    case ENAME: getvar(&val,&typ);	
		 break;
   }
}


xcable()

{
  datum d1;
  Symbol *param;

  param  = (Symbol *)*pc++;
  d1 = popm();
}


sphere()
{
  datum d1;
  Symbol *param;

  param = (Symbol *)*pc++;
  d1 = popm();
}

xchan()
{
  datum d1;
  Symbol *param;


  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case NA:
    case K:     break;
    case TYPE:  
    case VVREV:  
    case THRESH:
    case DENSITY: 
    case MAXCOND: d1 = popm();
		break;
   }
}

noise()
{
  datum d1;
  Symbol *param;
  static attrib *apnt; 

  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case VESNOISE:				 /* quantal vesicle noise */
    case CCHNOISE:				 /* quantal channel noise */
    		 d1 = popm();
		 break;

   }
}

synapse()
{
  datum d1;
  Symbol *param;

  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case VVREV:  
    case THRESH:
    case NFILT1: 
    case TIMEC1: 
    case NFILT2: 
    case TIMEC2: 
    case IGAIN: 
    case MAXCOND:
    case EXPON:
    case KD:    d1 = popm();
                break;

    case LINEAR:
    case OPEN: 
    case CLOSE:
                break;
   }
/* fprintf (stderr,"synapse %s %g\n",param->name,d1.val); /* */ 
}

gj()
{
  datum d1;

  d1 = popm();
}

rload()
{
  datum d1;

  d1 = popm();
}

rcap()
{
  datum d1;

  d1 = popm();
}

gcap()
{
  datum d1;

  d1 = popm();
}

resistor()
{
  datum d1;

  d1 = popm();
}

rbatt()
{
  datum d1;

  d1 = popm();
}

gbatt()
{
  datum d1;

  d1 = popm();
}

rg()
{
  datum d1;

  d1 = popm();
}

ri()
{
  datum d1;

  d1 = popm();
}

rm()
{
  datum d1;

  d1 = popm();
}

mcap()
{
  datum d1;

  d1 = popm();
}

mvrev()
{
  datum d1;

  d1 = popm();
}

mvrest()
{
  datum d1;

  d1 = popm();
}

vbuf()
{
}

conn1()
{
  datum d1,d2,d3;
  int narg;

  narg = (int)*pc++;
  if (narg > 2) d3 = popm();
  else d3.val = 0.0;
  if (narg > 1) d2 = popm();
  else d2.val = 0.0;
  d1 = popm();
  recepnode = d1.val;
  recepn2 = d2.val;
  recepn3 = d3.val;
}

conn1l()

{
  datum d1,d2,d3;
  int narg;
  node *npnt;

  narg = (int)*pc++;
  if (narg > 1) {
    if (narg > 2) 
         d3 = popm();
    else d3.val = 0.0;
    d2 = popm();	/* get location of node */
  }
  else d2.val = 0.0;
  d1 = popm();
  conn1();
}

conn2s()
{
  datum d1,d2,d3;
  int narg;

  narg = (int)*pc++;
  if (narg > 2) d3 = popm();
  else d3.val = 0.0;
  if (narg > 1) d2 = popm();
  else d2.val = 0.0;
  d1 = popm();
  recepnode = d1.val;
  recepn2 = d2.val;
  recepn3 = d3.val;
}

conn2sl()
{
  datum d1,d2,d3;
  int narg;
  node *npnt;

  narg = (int)*pc++;
  if (narg > 1) {
    if (narg > 2) 
         d3 = popm();
    else d3.val = 0.0;
    d2 = popm();	/* get location of node */
  }
  else d2.val = 0.0;
  d1 = popm();
  conn2s();
}

conn2d()
{
  datum d1,d2,d3;
  int narg;

  narg = (int)*pc++;
  if (narg > 2) d3 = popm();
  else d3.val = 0.0;
  if (narg > 1) d2 = popm();
  else d2.val = 0.0;
  d1 = popm();
  recepnode = d1.val;
  recepn2 = d2.val;
  recepn3 = d3.val;
}

conn2dl()
{
  datum d1,d2,d3;
  int narg;
  node *npnt;

  narg = (int)*pc++;
  if (narg > 1) {
    if (narg > 2) 
         d3 = popm();
    else d3.val = 0.0;
    d2 = popm();	/* get location of node */
  }
  else d2.val = 0.0;
  d1 = popm();
  conn2d();
}

gloc()

/* return x, y or z position to user */

{
    datum d1,d2,d3;
    Symbol *param;
    int nod1, nod2, nod3, recmod, narg;
    double record();

 param  = (Symbol *)*pc++;
 narg = (int) *pc++;
 if (narg > 2) d3 = popm();				/* get node */
 else d3.val = 0.0;
 if (narg > 1) d2 = popm();
 else d2.val = 0.0;
 d1 = popm();
 d1.val = 0.0;
 push(d1);

}

ndist()

/* return distance between two nodes */

{
    datum d1a,d1b,d1c,d2a,d2b,d2c;
    int narg;

 narg = (int) *pc++;
 if (narg > 2) d1c = popm();				/* get first node */
 else d1c.val = 0.0;
 if (narg > 1) d1b = popm();
 else d1b.val = 0.0;
 d1a = popm();

 narg = (int) *pc++;
 if (narg > 2) d2c = popm();				/* get second node */
 else d2c.val = 0.0;
 if (narg > 1) d2b = popm();
 else d2b.val = 0.0;
 d2a = popm();
 d1a.val = 0.0;
 push(d1a);
}

recept()

/* read in data from rod or cone statement:

                1       x position
                2       y position  (optional)
*/

{
  datum d1,d2;
  int i,narg;
  Symbol *param;
  recnod *maksnod(),*npnt;

  param = (Symbol *)*pc++;                      /* receptor type */
  narg  = (int)*pc++;
  d2.val = 0;
  if (narg > 1) {
     d2 = popm();                               /* read in y loc */
  }
  d1 = popm();                                  /* x loc */

  npnt = maksnod();
  npnt->recnm1 = recepnode;             /* recepnode is defined by conn() */
  npnt->recnm2 = recepn2;
  npnt->recnm3 = recepn3;
  npnt->xpos = (int)d1.val;
  npnt->ypos = (int)d2.val;
}

recparm()

/* get parameters for rods and cones */

{
  datum d1,d2,d3;
  Symbol *param;

  param  = (Symbol *)*pc++;
  switch (param->type)
   {
    case DIA: 
    case MAXCOND:
    case PIGM:
    case PATHL:
    case ATTF:
    case FILT:
		d1 = popm(); 
    case PHOTNOISE:
                break;
   }

}


/* variables for xstim below: */

static double start=0.0;
static double dur=0.1;
static double wavel=1.0;
static double inten=100.0;
static double backgr=10.0;
static double blurrad=11.0;
static double speriod = 20;
static double contrast = 0.5;
static double tfreq = 2.0;

static int stype,rnum1,rnum2,rnum3;
static int xstm=0, ystm=0;
static int xcent=0, ycent=0;		/* center coords for stim array */
static int stimflag=0;                  /* start stimulus */
static int blurfl=0;                    /* blur array done */

extern int stimdia;                     /* defined in "stimsub.c" */
extern int stimydist;                   /* defined in "stimsub.c" */

xstim()

/* read in data from stimulus statement:
*/

{
  datum d1,d2,d3;
  Symbol *param;
  int stimarr, backarr;
  int narg, spec;
  FILE *ftemp,*fopen();
  char stimfile[80];

  param  = (Symbol *)*pc++;
  
  switch (param->type)
   {
    case BAR: 
    case SPOT:  d1 = popm();
                stimdia   = d1.val;
                stype = param->type; 
                break;
    case SINE:  d1 = popm();
                speriod   = d1.val;
                stype = param->type; 
                break;
    case RECT:
		narg = (int)*pc++;
  		if (narg > 1) d2 = popm();
		else d2.val = 1.0;
    		d1 = popm();
		stimdia = d1.val;
		stimydist = d2.val;
		stype = param->type; 
		break;
    case NODE:
                narg = (int)*pc++;
                if (narg > 2) d3 = popm();
                else d3.val = 0.0;
                if (narg > 1) d2 = popm();
                else d2.val = 0.0;
                d1 = popm();
                rnum1 = d1.val;                 /* get node number */
                rnum2 = d2.val;
                rnum3 = d3.val;
                break;
    case ROD: 
    case CONE:  narg = (int) *pc++;
                if (narg > 2) d3 = popm();
                else d3.val = 0.0;
                if (narg > 1) d2 = popm();
                else d2.val = 0.0;
                d1 = popm();
                rnum1     = d1.val;
                rnum2     = d2.val;
                rnum3     = d3.val;
                stype = param->type; 
                break;
    case LOC:   narg  = (int)*pc++;
                if (narg > 0) {
                  if (narg > 1) {
                    d2 = popm();                /* read in y loc */
                    ystm = d2.val;
                  }
                d1 = popm();
                xstm = d1.val;
                }
                break;
    case CENTER:narg  = (int)*pc++;
                if (narg > 0) {
                  if (narg > 1) {
                    d2 = popm();                /* read in y center */
                    ycent = d2.val;
                  }
                d1 = popm();
                xcent = d1.val;
                }
                break;
    case START: d1 = popm();
                start = d1.val;
                stimflag = 1;
                break;
    case DUR:   d1 = popm();
                dur   = d1.val;
                break;
    case CONTRAST:  d1 = popm();
                contrast   = d1.val;
                break;
    case TFREQ:   d1 = popm();
                tfreq = d1.val;
                break;
    case WAVEL: spec = (int)*pc++;
		switch (spec) {
		  case 0: 
		    d1 = popm();
		    wavel = d1.val;
		    break;
		  case SUN:
                    wavel = 0;
                    break;
		  case XENON:
                    wavel = 1;
                    break;
		}
		break;
    case INTEN: d1 = popm();
                inten = d1.val;
                break;

    case BACKGR:d1 = popm();
                backgr = d1.val;
                backarr = 1;
                stype = param->type; 
                break;

    case VCLAMP:
    case CCLAMP:d1 = popm();
                inten = d1.val;
                stype = param->type; 
                break;

    case BLUR:  d1 = popm();
                blurrad = d1.val/2;
                break;
    case SFILE: if (*stimfile) fclose (stimout);
                *stimfile = 0;
                strcpy (stimfile,(char *)*pc++);  /* copy filename */
                if (stimhlp) fprintf (stderr,"stimfil '%s'\n",stimfile);
                if (*stimfile) {
                        stimfflg = 1; 
                        if ((ftemp=fopen(stimfile,"w")) == NULL) {
                          fprintf (stderr,"Can't open file '%s'\n",stimfile);
                          break;
                        }
                        else stimout = ftemp;
                }
                else {
                        stimfflg = 0;
                        stimout = stdout;
                }
                break;

    case STIM:  
        stimarr = 0;
        backarr = 1;
        
        switch (stype) {

        case 0:
        default: break;


        case VCLAMP:
                vclist(start,rnum1,rnum2,rnum3,inten,wavel,"vcl");
                vclist(start+dur,rnum1,rnum2,rnum3,inten,wavel,"eff");
                break;

        case CCLAMP:
                vclist(start,rnum1,rnum2,rnum3,inten,wavel,"ccl");
                vclist(start+dur,rnum1,rnum2,rnum3,inten,wavel,"iof");
                break;


        case BACKGR:
                  recback(stimarr,0.0,wavel);   /* zero recept stim inten */
                  recback(backarr,backgr,wavel); /* set recept backgr inten */
                  abslist(0.0, start);          /* make action list */
                  break;
        case ROD:
        case CONE: 				/* flash one receptor */
                flashrec(rnum1,rnum2,rnum3, inten,wavel,start);
                flashrec(rnum1,rnum2,rnum3,-inten,wavel,start+dur); 
                break;

        case BAR:
                if (!blurfl) {
                        makblur(blurrad);       /* make a blurring function */
                        blurfl = 1;
                }
		makrect(stimarr, stimdia, CONVSIZE,xstm-xcent,ystm-ycent,
							inten,wavel);
                if (stimflag) {
                  stimflag = 0;
                  recpat(stimarr,xcent,ycent);  /* set recept inten */
                  stimlist(1.0, start);         /* make an action list */
                  stimlist(-1.0, start+dur);
                  recback(stimarr,0.0,1.0);   /* zero recept stim inten */
                  recback(backarr,backgr,1.0);/* set recept backgr inten */
                }
                break;

        case RECT:
                if (!blurfl) {
                        makblur(blurrad);       /* make a blurring function */
                        blurfl = 1;
                }
		makrect(stimarr, stimdia, stimydist,xstm-xcent,ystm-ycent,
							inten,wavel);
                if (stimflag) {
                  stimflag = 0;
                  recpat(stimarr,xcent,ycent);  /* set recept inten */
                  stimlist(1.0, start);         /* make an action list */
                  stimlist(-1.0, start+dur);
                  recback(stimarr,0.0,1.0);   /* zero recept stim inten */
                  recback(backarr,backgr,1.0);/* set recept backgr inten */
                }
                break;

        case SPOT:
                if (!blurfl) {
                        makblur(blurrad);       /* make a blurring function */
                        blurfl = 1;
                }
                makspot(stimarr, stimdia,xstm-xcent,ystm-ycent,inten,wavel);
                if (stimflag) {
                  stimflag = 0;
                  recpat(stimarr,xcent,ycent);  /* set recept inten */
                  stimlist(1.0, start); 	/* make an action list */
                  stimlist(-1.0, start+dur);
                  recback(stimarr,0.0,1.0);     /* zero recept stim inten */
                  recback(backarr,backgr,1.0);  /* set recept backgr inten */
                }
                break;

        case SINE:
          if (!blurfl) {
                makblur(blurrad);       /* make a blurring function */
                blurfl = 1;
          }
					/* xmin is xstm; xmax is ystm. */
          if (speriod <= 0.0) speriod = 1;   /* spatial period from user */
          if (tfreq) {		/* drifting grating */
                 double timres,tperiod,tincr,stime;

            timres = .01;		/* time resolution: fractional */
            tperiod = 1.0 / tfreq;	/* temporal period */
            tincr = timres * tperiod;     /* time incr = time res * t period */
            for (stime=start; stime<(start+dur-tincr); stime += tincr) {
              maksine(stimarr,speriod,xstm+(stime-start)/tincr,ystm,
					-xcent,-ycent,inten,wavel); 
              recpat(stimarr,xcent,ycent);  /* set recept inten */
              stimlist(1.0, stime); 	    /* make an action list */
              stimlist(-1.0, stime+tincr);
              recback(stimarr,0.0,1.0);     /* zero recept stim inten */
              recback(backarr,backgr,1.0);  /* set recept backgr inten */

            }   /* for (stime= ; ; ) */
          }   /* if (drate) */
          else {		/* static grating */
              maksine(stimarr,speriod,xstm,ystm,-xcent,-ycent,inten,wavel); 

              if (stimflag) {			/* if start val is given */
                 stimflag = 0;
                 recpat(stimarr,xcent,ycent);  /* set recept inten */
                 stimlist(1.0, start); 	/* make an action list */
                 stimlist(-1.0, start+dur);
                 recback(stimarr,0.0,1.0);     /* zero recept stim inten */
                 recback(backarr,backgr,1.0);  /* set recept backgr inten */
              }
          } /* else */
	  break;

        }       /* switch (stype) */

        break;

  }     /* switch */
}

zplot()

/* read in data from plot statement
*/

{
  datum d1;
  int i,narg;

  narg  = (int)*pc++;
  for (i=0; i<narg; i++) 
    d1 = popm();
}

vplot()

/* read in data from "plot []" statement with no "max min"
*/

{
  datum d1,d2,d3;
  int narg;
  Symbol *param;

 param  = (Symbol *)*pc++;
 narg = (int)*pc++;

 if (narg > 2) d3 = popm();
 else d3.val = 0.0;
 if (narg > 1) d2 = popm();
 else d2.val = 0.0;
 d1 = popm();
 
 switch (param->type) {
   case V:
        break;
   case I:
        break;
   case L:
        break;
   case FA0:
   case FA1:
   case FA2:
   case FA3:
   case FA4:
   case FB0:
   case FB1:
   case FB2:
   case FB3:
   case FB4:
        break;
   case S:
        break;
 }
}

vplotm()

/* read in data from "plot [] max min" statement
*/

{
  datum d1,d2,d3,d4,d5;
  int narg;
  Symbol *param;

 param  = (Symbol *)*pc++;
 narg = (int)*pc++;

 if (narg > 2) d3 = popm();
 else d3.val = 0.0;
 if (narg > 1) d2 = popm();
 else d2.val = 0.0;
 d1 = popm();
 
 switch (param->type) {
   case V:
        break;
   case I:
        break;
   case L:
        break;
   case FA0:
   case FA1:
   case FA2:
   case FA3:
   case FA4:
   case FB0:
   case FB1:
   case FB2:
   case FB3:
   case FB4:
        break;
   case S:
        break;
 }
}

xrecord()

/* allow user to record from nodes directly */

{
    datum d1,d2,d3;
    Symbol *param;
    int nod, recmod, narg;
    double record();

 param  = (Symbol *)*pc++;
  narg = (int)*pc++;

 switch (param->type) {
   case V:
        recmod = VREC;
        break;
   case I:
        recmod = IREC;
        break;
   case IM:
        recmod = MREC;
        break;
   case L:
        recmod = LREC;
        break;
   case FA0:
   case FA1:
   case FA2:
   case FA3:
   case FA4:
   case FB0:
   case FB1:
   case FB2:
   case FB3:
   case FB4:
        recmod = NRECA0 + (param->type - FA0);
        break;
 }
   if (narg > 2) d3 = popm();
   else d3.val = 0.0;
   if (narg > 1) d2 = popm();
   else d2.val = 0.0;
   d1 = popm();                         /* get node */
   nod = d1.val;
/*   d1.val = record (nod1, nod2, nod3, recmod); */
   push(d1);
}

grph()

{
  datum d1,d2;
  Symbol *param;
  static int narg,i,n,np;
  static int initfl=0, plotfl=0;

 param  = (Symbol *)*pc++;

 if (param == 0) {                              /* graph a point */
        narg = (int) *pc++;
        n = narg - 1;
        if (n > PLOTNODSIZ) n=PLOTNODSIZ;
        for (i=n-1; i>=0; i-- ) {
           d2 = popm();
        }       
        d1 = popm();
 }
 else
 switch (param->type) {
   case X:                      /* bug: can't use for indiv separate graphs */
        d2 = popm();
        d1 = popm();
        break;

   case Y:
        d2 = popm();
        d1 = popm();
	np++;
        break;

  case INIT:
        initfl = 1;
        break;

  case RESTART:
        break;

  case PEN:
	if (narg > PLOTNODSIZ) n=PLOTNODSIZ;
	if (narg > np) narg = np;
	for (i=narg; i>0; i-- ) {
	   d1 = popm();
	}	

	break;

  }
}

dispnod()

{
  datum d1,d2,d3;
  Symbol *param,nulparm;
  int narg;

  param = (Symbol *)*pc++;
  nulparm.type = CONNECT;
  if (param==0) param = &nulparm;

  switch (param->type) {
  case MATCHING:
  case RANGE:
  case CONNECT:
  case ELEMENT:
 
  narg = (int)*pc++;
  if (narg) {
    if (narg > 2) d3 = popm();
    if (narg > 1) d2 = popm();
    d1 = popm();
  }

  narg = (int)*pc++;
  if (narg) {
    if (narg > 2) d3 = popm();
    if (narg > 1) d2 = popm();
    d1 = popm();
  }
  break;


  case CABLE:
  case SPHERE:
  case SYNAPSE:
  case CHAN:
  case ROD:
  case CONE:
  case GJ:
  case LOAD:
  case RESISTOR:
  case CAP:
  case GNDCAP:
  case BATT:
  case GNDBATT:
  case BUF:
    break;

  case EXCEPT: break;

  case XROT:				/* get rotation for display */
    d1 = popm();
    break;

  case YROT:
    d1 = popm();
    break;

  case ZROT:
    d1 = popm();
    break;

  case SIZE:				/* get size of display window */
    d1 = popm();
    break;

  case DSCALE:				/* get display scale */
    d1 = popm();
    break;

  case COLOR:				/* get color */
    d1 = popm();
    break;

  case CALIBLIN:			/* get calib line */
    d1 = popm();
    break;

  case HIDE:				/* make picture with hidden lines */
    break;

  case CENTER:
    narg = (int)*pc++;
    if (narg) {
      if (narg > 2) d3 = popm();
      if (narg > 1) d2 = popm();
      d1 = popm();
    }
    break;
 
  case DISPLAY:
				/* first, translate, rotate, and scale: */
   break;	/* case DISPLAY */
  } 		/* switch (param) */
}


#define NUMDIM 3			/* number of node dimensions */

foreacode()
{
	datum d;
	Inst *savepc = pc;
	int i,narg,match;
	int arg[NUMDIM];
	int val[NUMDIM];
        double *varp[NUMDIM];
	short *typ;
        node *npnt;

    narg = (int)*pc++;			/* number of snode args to expect */
    if (narg>NUMDIM) narg = NUMDIM;
    for (i=0; i<NUMDIM; i++) {
       varp[i] = NULL;
       val[i]  = NULLVAL;
       arg[i]  = (int)*pc++;
    }
    execute (savepc+6);			/* evaluate node descriptors */

    for (i=narg-1; i>=0; i--) {
         if (arg[i]) getvar(&varp[i],&typ); /* var pointer for node dimension */
	 else {
	   d = popm();			 /* else get value of node dimension */
	   val[i] = d.val;
         }
    }

   for (npnt=nodepnt; npnt; npnt=npnt->next) {	/* search all nodes */
      for (match=1,i=0; i<narg; i++) {
	if (!arg[i]) 			/* if dimension was given a value */
	  switch (i) {
	   case 0: if (val[i] != npnt->nodenm1); match = 0; break; 
	   case 1: if (val[i] != npnt->nodenm2); match = 0; break; 
	   case 2: if (val[i] != npnt->nodenm3); match = 0; break; 
	  }	
     }
     if (!match) continue;

     for (i=0; i<narg; i++) {
	if (arg[i])			/* if dimension was given blank var */
	  switch (i) {
	  case 0: *varp[i] = npnt->nodenm1; break; 
	  case 1: *varp[i] = npnt->nodenm2; break; 
	  case 2: *varp[i] = npnt->nodenm3; break; 
	  }	
     }
	execute(*((Inst **)(savepc+4)));	/* body */
	if (stopping) {
		if (returning)
			break;
		stopping = 0;
		if (breaking) {
			breaking = 0;
			break;
		}
		continuing = 0;
	}
     }
     if (!returning)
		pc = *((Inst **)(savepc+5));	/* next stmt */
}


findconnect()

/* find all nodes that connect to a branch;
   set up each node's branch entries in the node list.
   If node doesn't exist, create a new node.
 */

{
}

setnodlst(npnt,bpnt)
   node *npnt;
   conn *bpnt;

/* add to a node's branch list */

{
}
modrun()

{
  datum d1,d2;
  Symbol *mtyp;
 
  mtyp = (Symbol *)*pc++;
  varcopy();
  findconnect();
  switch (mtyp->type)
    {

     case RUN:  
/*              actrun(); */
                break;
     case STEP:  
                d1 = popm();
/*              actrun(d1.val); */
                break;
    
    }

}

/*------------------------------------*/

