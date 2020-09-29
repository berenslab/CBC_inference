/* segment ncstim in program nc */

/* Routines to make and access stimuli. */

extern "C" {

#include <stdio.h>
#include <math.h>

#ifdef CPML
#include <cpml.h>
#endif

}


#include "nc.h"
#include "y.tab.h"
#include "ndef.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "control.h"
#include "ncio.h"

extern recstim *recspnt;
extern recstim *recsend;
extern int cumrecst;

static int sinitfl=0;

extern int stimfflg;			/* =1 -> get stimulus from file */
extern FILE *stimin;			/* stimulus file */


#define DEBUG	/* */

#define RTIMESTEP 0.0001                /* time step for reading stimuli */
#define SBINSIZ 100                     /* number of bins for stimulus table */
#define BINLISTLEN 100                  /* length of bins for stimulus table */

struct recsptr {
	recstim *next;
	recstim *last;
	};

extern double specdens[TOTREC][PIGMSIZ];
extern double filts[FILTS][PIGMSIZ];
extern double lights[LIGHTS][PIGMSIZ];
extern double lumlight[LUMINS][LIGHTS];
extern double lightsens[TOTREC][LIGHTS];

int sbinsiz = SBINSIZ;
recsptr bintab[SBINSIZ][2];
recsptr *stimbintab = &bintab[0][0];  		/* stim bin table */

static double timerange = 0;
static double maxstimtime=0, minstimtime=0;

char *emalloc(unsigned int n);
void wsav (int pigm, double pathl, double wavel, int filt, double sens);
void execerror (const char *s, const char *t);
void efree(void *ptr);
recstim *makrstim(double start, int recnm1, int recnm2,
	int recnm3,int recnm4,double inten, double wavel, 
	double mask, int stimchan, int bgend, const char *action);
void delrstim(recstim *rpnt);
int setstimchan(int stimchan);

/*---------------------------------------------------*/

void setstimbins()

{
 sinitfl =  0;			/* reset binning because new stim added */
}

/*---------------------------------------------------*/

int getstimbin (double stime)

/* Find the bin that a stimulus is in, according to 
   its time.  This speeds up the stimulus routine by 
   reducing the number of stimuli to be checked.
*/

{
   int b;
   double bf;

   bf = (stime - minstimtime) / timerange;
   if (bf < 0) return -10;
   else {
	  if ((b=int(bf*sbinsiz)) > sbinsiz) return -1;
          else				     return  b;
   }
}

/*---------------------------------------------------*/

void sinstall (recstim *newstim)
              
/* Install a stimulus in bin table */
/*  Uses "direct chaining" with field "rnext". */

{
    int b;
    recstim *rpnt,*rlast;
    double stime;

   newstim->rnext = (recstim *)NULL;
   newstim->rlast = (recstim *)NULL;
   stime = newstim->time;

   b = getstimbin(stime);

   if (b >= 0) {			     /* if stimulus is within range */
     if (!(rpnt=stimbintab[b].next)) {      /* if table entry empty */
        stimbintab[b].next = newstim;          /* install directly in table */
        stimbintab[b].last = newstim;
     }
     else {                                  /* otherwise, add to list */
        rlast = stimbintab[b].last;
        rlast->rnext = newstim;              /* put new stim at end of list */
        newstim->rlast = rlast;
        stimbintab[b].last = newstim;
    }
  }
  else {
    execerror ("sinstall: Illegal stimulus, "," stopping... ");
  }
}

/*---------------------------------------------------*/

void makestimbins ()

/* Sort stimuli into a "table of bins" for faster access." */

{
  int i,bins,nrecs;
  static int osbinsiz=SBINSIZ;
  recstim *rspnt;
  double rtimestep,stime;
  static double maxstim, minstim;

  if (!sinitfl) {			/* initialize table once at start */
    sinitfl = 1;
    maxstim = -1e30;
    minstim =  1e30;
    for (nrecs=0,rspnt=recspnt; rspnt; nrecs++,
				rspnt=rspnt->next) { /* find time range */
      stime = rspnt->time;
      if (stime > maxstim) maxstim = stime;
      if (stime < minstim) minstim = stime;
    }
    maxstim =  (floor(maxstim * 10000)+1) / 10000.0;
    minstim =  (floor(minstim * 10000)-1) / 10000.0;

			/* if new bin range required */

    if ((maxstim != maxstimtime) || (minstim != minstimtime)) {
      maxstimtime = maxstim;
      minstimtime = minstim;
    }

    timerange = maxstimtime - minstimtime;
    if (timerange <= 0) timerange = .001;
    rtimestep = max(timinc,RTIMESTEP);		/* find stim timestep */
    bins = int (timerange / rtimestep) - 1;	/* make bins coarse */
    if (bins <= 0) bins = 1;

    sbinsiz = nrecs / BINLISTLEN;		/* max length of bins = 100 */
    sbinsiz = min (bins,sbinsiz);		/* make bins fit into array */

    if (sbinsiz > osbinsiz) {
        stimbintab = (recsptr *)emalloc(sizeof(recstim*)*(sbinsiz+1)*2);
        osbinsiz = sbinsiz;
    }

     for (i=0; i<=sbinsiz; i++) {		/* zero the bin table */
        stimbintab[i].next = (recstim *)NULL;
        stimbintab[i].last = (recstim *)NULL;
     }
    for (rspnt=recspnt; rspnt; rspnt=rspnt->next) {  /* bin all stimuli */
      sinstall (rspnt);
    }
  }
}

/*---------------------------------------------------*/

recstim *getstimp1 (double stime)

/* Get a pointer to the bin that a stimulus must be in.   */

{
   int b;
   recstim *rpt;

  makestimbins();			/* redo bins if new stims added */
  b = getstimbin(stime);	 	/* find bin pointer */
  if (b >= 0) { 
    rpt= stimbintab[b].next;
    return (rpt);
  }
  else return (recstim *)NULL;
}

/*---------------------------------------------------*/

recstim *getstimp2 (double stime)

/* Get a pointer to the bin after the one that a stimulus must be in.   */

{
   int b;
   recstim *rpt;

  makestimbins();			/* redo bins if new stims added */
  b = getstimbin(stime)+1;	 	/* find bin pointer */
  if (b >= 0) {
    b = min(b,sbinsiz);
    rpt= stimbintab[b].next;
    return (rpt);
  }
  else return (recstim *)NULL;
}

/*---------------------------------------------------*/

recstim *makvstim(double start, int nodnm1, int nodnm2, 
		int nodnm3, int nodnm4, double value, int bgend, const char *action)
                        
/* Make a new clamp stimulus and link it to the list. 
   A clamp stimulus looks just like a photrec stimulus, but
   has different meaning when running the simulation. 
   Check runstim() in ncsub.cc for details on meaning.
*/

{
    recstim *rspnt;
 
  setstimbins();			/* reset bins */

  if ((rspnt=(recstim *)emalloc(sizeof(recstim))) == NULL) {
     ncfprintf (stderr,"no space left for clamp stim %d\n", cumrecst+1);
     return ((recstim *)NULL);  
  }
  rspnt->next = (recstim *)NULL;
  if (!recspnt) recspnt = rspnt; 	/* save head if first synap */
  rspnt->last = recsend;
  if (recsend)
    recsend->next = rspnt;
  recsend = rspnt;

  // fprintf (stderr,"makvstim %g %g %s\n",start,value,action);
  rspnt->ctype = *action;
  rspnt->time = start;
  rspnt->bgend  = bgend;
  rspnt->recnm1 = nodnm1;
  rspnt->recnm2 = nodnm2;
  rspnt->recnm3 = nodnm3;
  rspnt->recnm4 = nodnm4;
  rspnt->val = value;
  rspnt->mask = 0;
  rspnt->stimchan = 0;
  cumrecst++; 			/* increment total */
  //ncfprintf(stderr,"makvstim start %g node %d value %g\n",start,nodnm3,value); //new for debugging
  return (rspnt); 
}

/*---------------------------------------------------*/

void makrstim(double start, double dur, photrec *rpnt, double inten, double wavel, 
	double mask, int stimchan, const char *action)
                              
/* Make a new photrec stimulus and link it to the list. */
/* Do both stim on and stim off, no bgend arg. */

{
    recstim *rspnt;

#ifdef DEBUG
   if (debug & 4 && debugz & 1) ncfprintf(stderr,"makrstim start %g %d %d %g '%c'\n",
				start,rpnt->recnm1,rpnt->recnm2,inten,*action);
#endif

  setstimbins();			/* reset bins */

  /* stimulus on */

  if ((rspnt=(recstim *)emalloc(sizeof(recstim))) == NULL) {
     ncfprintf (stderr,"no space left for recstim %d\n", cumrecst+1);
     return;  
  }
  rspnt->next = (recstim *)NULL;
  if (!recspnt) recspnt = rspnt; 	/* save head if first synap */
  rspnt->last = recsend;
  if (recsend)
    recsend->next = rspnt;
  recsend = rspnt;

  rspnt->ctype = *action;
  rspnt->time = start;
  rspnt->bgend = 0;
  rspnt->val = inten;
  rspnt->stimchan = setstimchan(stimchan);
  rspnt->mask = mask;
  rspnt->wavel = wavel;
  cumrecst++; 			/* increment total */

  rspnt->recnm1 = rpnt->recnm1;
  rspnt->recnm2 = rpnt->recnm2;
  rspnt->recnm3 = rpnt->recnm3;
  rspnt->recnm4 = rpnt->recnm4;

  /* - - - - - - - - - - */
  /* stimulus off */

  if ((rspnt=(recstim *)emalloc(sizeof(recstim))) == NULL) {
     ncfprintf (stderr,"no space left for recstim %d\n", cumrecst+1);
     return;  
  }
  rspnt->next = (recstim *)NULL;
  if (!recspnt) recspnt = rspnt; 	/* save head if first synap */
  rspnt->last = recsend;
  if (recsend)
    recsend->next = rspnt;
  recsend = rspnt;

  rspnt->ctype = *action;
  rspnt->time = start+dur;
  rspnt->bgend = 1;
  rspnt->val = -inten;
  rspnt->stimchan = setstimchan(stimchan);
  rspnt->mask = -mask;
  rspnt->wavel = wavel;
  cumrecst++; 			/* increment total */

  rspnt->recnm1 = rpnt->recnm1;
  rspnt->recnm2 = rpnt->recnm2;
  rspnt->recnm3 = rpnt->recnm3;
  rspnt->recnm4 = rpnt->recnm4;
}

/*---------------------------------------------------*/

void makrstim(double start, double dur, photorec *ept, double inten, double wavel, 
	double mask, int stimchan, const char *action)
                              
/* Make a new photorec stimulus and link it to the list. */
/* Do both stim on and stim off, no bgend arg. */

{
    recstim *rspnt;

#ifdef DEBUG
   if (debug & 4 && debugz & 1) ncfprintf(stderr,"makrstim start %g %d %d %g '%c'\n",
				start,ept->node1a,ept->node1b,inten,*action);
#endif

  setstimbins();			/* reset bins */

  /* stimulus on */

  if ((rspnt=(recstim *)emalloc(sizeof(recstim))) == NULL) {
     ncfprintf (stderr,"no space left for recstim %d\n", cumrecst+1);
     return;  
  }
  rspnt->next = (recstim *)NULL;
  if (!recspnt) recspnt = rspnt; 	/* save head if first synap */
  rspnt->last = recsend;
  if (recsend)
    recsend->next = rspnt;
  recsend = rspnt;

  rspnt->ctype = *action;
  rspnt->time = start;
  rspnt->bgend = 0;
  rspnt->val = inten;
  rspnt->stimchan = setstimchan(stimchan);
  rspnt->mask = mask;
  rspnt->wavel = wavel;
  cumrecst++; 			/* increment total */

  rspnt->recnm1 = ept->node1a;
  rspnt->recnm2 = ept->node1b;
  rspnt->recnm3 = ept->node1c;
  rspnt->recnm4 = ept->node1d;

  /* - - - - - - - - - - - - - */
  /* stimulus off */

  if ((rspnt=(recstim *)emalloc(sizeof(recstim))) == NULL) {
     ncfprintf (stderr,"no space left for recstim %d\n", cumrecst+1);
     return;  
  }
  rspnt->next = (recstim *)NULL;
  if (!recspnt) recspnt = rspnt; 	/* save head if first synap */
  rspnt->last = recsend;
  if (recsend)
    recsend->next = rspnt;
  recsend = rspnt;

  rspnt->ctype = *action;
  rspnt->time = start+dur;
  rspnt->bgend = 1;
  rspnt->val = -inten;
  rspnt->stimchan = setstimchan(stimchan);
  rspnt->mask = -mask;
  rspnt->wavel = wavel;
  cumrecst++; 			/* increment total */

  rspnt->recnm1 = ept->node1a;
  rspnt->recnm2 = ept->node1b;
  rspnt->recnm3 = ept->node1c;
  rspnt->recnm4 = ept->node1d;
}

/*---------------------------------------------------*/

recstim *makrstim(double start, int recnm1, int recnm2, 
	int recnm3, int recnm4, double inten, double wavel, 
	double mask, int stimchan, int bgend, const char *action)
                              
/* make a new photrec stimulus and link it to the list. */

{
    recstim *rspnt;

#ifdef DEBUG
   if (debug & 4 && debugz & 1) ncfprintf(stderr,"makrstim start %g %d %d %g '%c'\n",
				start,recnm1,recnm2,inten,*action);
#endif

  setstimbins();			/* reset bins */

  if ((rspnt=(recstim *)emalloc(sizeof(recstim))) == NULL) {
     ncfprintf (stderr,"no space left for recstim %d\n", cumrecst+1);
     return ((recstim *)NULL);  
  }
  rspnt->next = (recstim *)NULL;
  if (!recspnt) recspnt = rspnt; 	/* save head if first synap */
  rspnt->last = recsend;
  if (recsend)
    recsend->next = rspnt;
  recsend = rspnt;

  rspnt->ctype = *action;
  rspnt->time = start;
  rspnt->bgend = bgend;
  rspnt->recnm1 = recnm1;
  rspnt->recnm2 = recnm2;
  rspnt->recnm3 = recnm3;
  rspnt->recnm4 = recnm4;
  rspnt->val = inten;
  rspnt->stimchan = setstimchan(stimchan);
  rspnt->mask = mask;
  rspnt->wavel = wavel;
  cumrecst++; 			/* increment total */
  return (rspnt); 
}

/*---------------------------------------------------*/

void delrstim(recstim *rpnt)

/* delete a photrec stimulus. */

{
  recstim *rpt, *rlast;
  int b,found,n;
 
  if (!rpnt)  return;

  makestimbins();			/* redo bins if new stims added */
  b = getstimbin(rpnt->time);

  if (b < 0) {
    execerror ("delrstim: Illegal stimulus, "," stopping... ");
  }

       /* Patch pointers in bin list */

  if (rpnt->rlast) rpnt->rlast->rnext = rpnt->rnext; 
  if (rpnt->rnext) rpnt->rnext->rlast = rpnt->rlast;
  if (stimbintab[b].next == rpnt) stimbintab[b].next = rpnt->rnext;
  if (stimbintab[b].last == rpnt) stimbintab[b].last = rpnt->rlast;

        /* Patch pointers for main list */

  if (rpnt->last) rpnt->last->next = rpnt->next;     /* patch last pointer */
  if (rpnt->next) rpnt->next->last = rpnt->last;
  if (recspnt==rpnt) recspnt = rpnt->next;
  if (recsend==rpnt) recsend = rpnt->last;
  
  efree (rpnt);
  cumrecst--; 				/* decrement total */

  return; 
}

/*---------------------------------------------------*/

int readstim(double endtime)
                   

/* Read the current stimulus from a disk file, until
   the stimulus time is equal or greater than endtime.
   The stimulus is read once every 1 msec (simulated time).
   Erase stimuli earlier than startime.
*/

{
   static int recnm1, recnm2, recnm3, recnm4, seq, x;
   static int stimchan, needline = 1, wait = 0;
   static double inten,wavel,mask;
   static double stime=-1e6, oldstime = -1e6;
   static char linebuf[80],action[10];
   recstim *rpnt;
   int error;
   short int intaction;

  if (! stimfflg) return(1);

#ifdef DEBUG
   if (debug & 4 && debugz & 1) ncfprintf(stderr,"readstim start %g\n",endtime);
#endif
   
   //ncfprintf(stderr,"readstim endtime %g stime %g oldstime %g\n",endtime,stime,oldstime); //new for debugging
				/* then, check for new stimuli */
  error = 0;		
  do {				/* get new stimuli */
    if (endtime < oldstime) {
	 needline = 1;
         wait = 0;
    }
    if (needline) {
      if ((fgets (linebuf,80,stimin)) == NULL) { break; }  /* not an error */
      if (*linebuf == '#') continue;		/* ignore comment lines */
      oldstime = stime;

			/* can ignore seq here since it's for stim sort */

      if ((x=sscanf (linebuf,"%lf %d %d %d %d %lf %lf %lf %d %s %d\n",
		&stime,&recnm1,&recnm2,&recnm3,&recnm4,
		&inten,&wavel,&mask,&stimchan,action,&seq))< 11) {

       stimchan = 0;
       if ((x=sscanf (linebuf,"%lf %d %d %d %d %lf %lf %lf %s %d\n",
		&stime,&recnm1,&recnm2,&recnm3,&recnm4,
		&inten,&wavel,&mask,action,&seq))< 10) {

        mask   = 0;
        if ((sscanf (linebuf,"%lf %d %d %d %d %lf %lf %s %d\n",
		&stime,&recnm1,&recnm2,&recnm3,&recnm4,
		&inten,&wavel,action,&seq))< 9) {

		/* If line doesn't read correctly,
		   try old format stim file (no recnm4): */

         recnm4 = -1;
         if ((sscanf (linebuf,"%lf %d %d %d %lf %lf %s\n",
		&stime,&recnm1,&recnm2,&recnm3,
		&inten,&wavel,action))< 7) {
			 error = 1; break; }
         }
       }
      }

    }     /* if needline */

  /*   ncfprintf (stderr,"%g %g %g %d %d %d %d %g %g %g '%s' %d %d %d %d\n",
	endtime,oldstime,stime,recnm1,recnm2,recnm3,recnm4,
		inten,wavel,mask,action,seq,needline,wait,x); /* */

    if (stime < oldstime) {
       needline = 0;
       wait = 1;
       break;
    }
    if (stime > endtime) {	/* have read too far, wait for next time step */
    	needline = 0;
	break;
    }

    if (recnm4<0) recnm4 = NULLVAL;
    if (recnm3<0) recnm3 = NULLVAL;
    if (recnm2<0) recnm2 = NULLVAL;
    if (recnm1<0) recnm1 = NULLVAL;

    intaction = *action;

    switch (intaction) {

      case 'e':

	if ((rpnt=makvstim(stime,recnm1,recnm2,recnm3,recnm4,
				inten,int(wavel),action))==NULL) {
	  ncfprintf (stderr,"readstim: can't make vclamp stim %d\n",cumrecst);
	  error = 1;
	  break;
        }

      case 'g':
	//ncfprintf(stderr,"making cclamp stim from file\n");
	if ((rpnt=makvstim(stime,recnm1,recnm2,recnm3,recnm4,
				inten,int(wavel),action))==NULL) {
	  ncfprintf (stderr,"readstim: can't make cclamp stim %d\n",cumrecst);
	  error = 1;
	  break;
	}

      default:

    	if ((rpnt=makrstim(stime,recnm1,recnm2,recnm3,recnm4,
				inten,wavel,mask,stimchan,0,action))==NULL) {
	  ncfprintf (stderr,"readstim: can't make rec stim %d\n",cumrecst);
	  error = 1;
	  break;
        }
    }
    needline = 1;
  } while ((stime < endtime) && (oldstime <= stime));

  if (error)	return (0);
  else 		return (1);
}

/*------------------------------------*/

/* pigment action spectra tables */

			/* wave.h contains pigment and light tables */
// #include "wave.h"

#define RNDUPX 0.000001

double dlookup(int pigm, double wavel, double *arr)

/* calculate the specific optical density
   for given pigment type and wavelength of
   incident light.  Use quadratic interpolation.

    (Interpolation taken from Conte and de Boor (1980), 
     "Elementary Numerical Analysis", p 39. 
*/

{
   double dens,logdens;
   int w0,w1,w2;
   int dx1,dx2,dx3;
   double windx;
   double l0,l1,l2;
   double r0,r1,r2;
   double v0,v1,v2;


  windx = (wavel-MINWAV) / WAVINC;
  w1 = int(windx + RNDUPX);
  if (w1 >= (PIGMSIZ-1)) w1 = PIGMSIZ-2;
  else if (w1 < 0) 	 w1 = 0;

  if ((windx-w1) > 0.5) {			/* x val is closer to next */ 	
  //if (0) {			/* x val is closer to next */ 	
    w0 = w1;					/* x values */
    w1 = w0+1;
    w2 = w0+2;
  }
  else {					/* x val is closer to this */
    w0 = w1-1;					/* x values */
    w2 = w1+1;
  }

  v1 = *(arr+pigm*PIGMSIZ+w1);

  if (w0 < 0) v0 = v1 + v1 - *(arr+pigm*PIGMSIZ+w2);
  else	      v0 = *(arr+pigm*PIGMSIZ+w0);		/* y values */


  if (w2 >= PIGMSIZ) v2 = v1 + v1 - v0;
  else	      v2 = *(arr+pigm*PIGMSIZ+w2);		/* y values */

  r0 = windx - w0; 				/* remainder values */
  r1 = windx - w1;
  r2 = windx - w2;


  l0 = r1 * r2 * 0.5;  
  l1 = -r0 * r2;
  l2 = r0 * r1 * 0.5;

  logdens = v0 * l0 + v1 * l1 + v2 * l2;

//ncfprintf (stderr,"v0 %g v1 %g v2 %g\n",v0,v1,v2);
//ncfprintf (stderr,"l0 %g l1 %g l2 %g\n",l0,l1,l2);

//ncfprintf (stderr,"xwavel %g logdens %g\n",wavel,logdens);

  dens = exp(logdens*LN10);
  return (dens);
}

/*------------------------------------*/

struct wparm {
	int pigm;
	int filt;
	double pathl;
	double wavel;
	double sens;
	int age;
	};

#define OLDPSIZ 8
struct wparm oldparm[OLDPSIZ] = {0};	/* array to hold all combin. of */
					/*  pigm, pathl, wavel, filt calcs. */

int wfind (int pigm, double pathl, double wavel, int filt, double *oldsens)
                  
/* Look for a pigment sensitivity previously calculated
   with a specific combination of pigment, path length,
   light, and filter parameters.  
*/

{
  int i, found;
  static int first=1,lru=0,lruaddr=0;

  if (first) {
      first = 0;
      for (i=0; i<OLDPSIZ; i++) {
    	oldparm[i].pathl = 0;		/* initialize all entries first time */
    	oldparm[i].age   = 0;	
      }
  }
	 
  found = 0;
  for (i=0; !found && i<OLDPSIZ; i++) {	/* try to find same calc done before */

    if (oldparm[i].pathl == 0) break;	/* stop when list is done */

    oldparm[i].age++;		/* increment age of parameter */

    if (oldparm[i].pigm  != pigm)  continue;
    if (oldparm[i].pathl != pathl) continue;
    if (oldparm[i].wavel != wavel) continue;
    if (oldparm[i].filt  != filt)  continue;

    found = 1;
    *oldsens = oldparm[i].sens; 	/* old calc found; use it */
    oldparm[i].age = 0;			/* if found, reset age of parameter */
  }
  return (found);
}

/*------------------------------------*/

void wsav (int pigm, double pathl, double wavel, int filt, double sens)

/* Save a pigment sensitivity done with a specific combination
    of pigment, path length, light, and filter parameters
    for future use.  If not found, add this set
   of parameters, and if necesssary, remove an old entry to do so.
*/

{
  int i, found, age, lru, lruaddr;

  found = lru = lruaddr = 0;
  for (i=0; !found && i<OLDPSIZ; i++) {	/* try to find same calc done before */

    if (oldparm[i].pathl == 0) break;	/* stop when list is done */

    age = oldparm[i].age;		/* find age of parameter */
    if (age > lru) {			/*  If this is the correct one, */
	  lru = age;			/*  then no need to save it below */
	  lruaddr = i;
    }
    if (oldparm[i].pigm  != pigm)  continue;
    if (oldparm[i].pathl != pathl) continue;
    if (oldparm[i].wavel != wavel) continue;
    if (oldparm[i].filt  != filt)  continue;

    found = 1;				/* stop if correct one found */
  }

  if (!found) {
    if (i<OLDPSIZ) lruaddr = i;		/* use new location if possible */
    oldparm[lruaddr].pigm  = pigm;	/* save new parameters and sens */
    oldparm[lruaddr].pathl = pathl;
    oldparm[lruaddr].wavel = wavel;
    oldparm[lruaddr].filt  = filt;
    oldparm[lruaddr].sens  = sens;
    oldparm[lruaddr].age   = 0;		/* reset age of parameter */
  }
}

/*------------------------------------*/

double rsens(photrec *rpnt, double wavel, int filt)
               
/* Calculate receptor sensitivity given pigment, pathlength,
   wavelength, and filter.
   First, look to see if calculation has been done before;
   if so, then use old calculation.
   Second, if light is monochromatic, do the calculation directly.
   Otherwise, light must be a standard light source;
   check to see if path length and filter are standard; if so, then
   use standard lookup table () to get sensitivity.
   Lastly, if non-standard path length or filter,
   then calculate sensitivity over full wavelength range, 
   and save the calculation for future use.
*/

{
    int j,pigm,wav,lum;
    double od, pathl, logdens, specd, sens, totsens;
    double oldsens;

pigm = rpnt->chtyp->pigm;
lum = (pigm? 1 : 0);			/* scotopic or photopic calibration */
pathl = rpnt->pathl;
if (wfind (pigm,pathl,wavel,filt,&oldsens))
     sens = oldsens;		/* check to see if this has been done before */
else {
  if (filt >= FILTS) filt = FILTS-1;
  wav = (int)wavel;
  if (wav < 0) wav = 0;
  if (wav >= LIGHTS) {			/* monochromatic light */
     specd = dlookup(pigm,wavel,specdens[0]);
     od = specd * pathl;
     sens = (1.0 - exp(-od * LN10)) * dlookup(filt,wavel,filts[0]);
  }
  else {			/* (wav < LIGHTS): standard light source */
    if (pathl == rpnt->chtyp->pathl && filt == 0) {	/* if pre-computed */
       sens = lightsens[pigm][wav] / lumlight[lum][wav];
    }
    else {			/* re-calculate pigm sens to light func */
      totsens = 0.0;
      for (j=0; j<PIGMSIZ; j++) {
        logdens = specdens[pigm][j];	/* no interpolation */
        specd = exp(logdens*LN10);	/* density table is log(spec od) */
        od = specd * pathl;
        sens = 1.0 - exp(-od * LN10);
	totsens += sens * exp(-filts[filt][j]*LN10) * lights[wav][j];
      } 
      sens = totsens/(PIGMSIZ*lumlight[lum][wav]); /* divide by calib sens */
    }
   }
   wsav (pigm,pathl,wavel,filt,sens); /* save calculation */
  }

/* ncfprintf (stderr,"sens %g pigm %d wavel %g filt %d path %g\n",
			sens,pigm,wavel,filt,pathl);	/* */
  return sens;
 } 
 
