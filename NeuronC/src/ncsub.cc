/* Module ncsub in Program nc */

/* Models neuronal circuits */

/* Simulates a neuron with a compartmental model.
It assumes a dendrite is made up of small isopotential
compartments which are modeled as electrical nodes.
Each node is connected to its neighbors through axial 
resistances, and has conduction pathways to ground
through membrane resistance and capacitance.
Additional conduction pathways are allowed for any node. 

Unless otherwise noted, all length values are calibrated
in meters, all voltages are in volts, all currents are
in amperes, and capacitance is in farads.

        Oct 87                  R.G. Smith

*/

extern "C" {

#include <stdio.h>
#include <math.h>
#include "gr.h"

#ifdef CPML
#include <cpml.h>
#endif
}

#include "scheduler.h"

#include "nc.h"
#include "y.tab.h"
#include "ndef.h"
#include "nconst.h"
#include "control.h"
#include "ncsub.h"
#include "ncomp.h"
#include "ncplot.h"
#include "ncelem.h"
#include "ncio.h"
#include "drand.h"
#include "ncio.h"

#define DEBUG

#ifdef DEBUG
#include "ncdebug.h"
#endif

/*---------------------------------------*/

/* From "nclist.h": */

extern int cumelem;
extern int cumattr;
extern int cumnattr;
extern int cumcattr;
extern int cumnode;
extern int cumcomp;
extern int cumconn;
extern int cumchan;
extern int cumsynap;
extern int cumload;
extern int cumphotrec;
extern int cumrecst;
extern int cumcacomp;
extern int cumcycacomp;
extern int cumcycgcomp;
extern int cumclst;
extern int cumdlst;
extern int cumlpfilt;
extern int cumrand;
extern int cumntcomp;
extern int cumgj;
extern int reccum;
extern int delelmt;
extern int ndomains;

extern comp *compnt;
extern comp *compend;
extern conn *connpnt;
extern conn *connend;
extern synap *synpnt;
extern synap *synend;
extern chan *chanpnt;
extern chan *chanend;
extern load *loadpnt;
extern load *loadend;
extern photrec *recpnt;
extern photrec *recend;
extern recstim *recspnt;
extern recstim *recsend;

void efree(void *p);
char *emalloc(unsigned int n);

int ztotsiz=0;

/*---------------------------------------*/

extern int disp_ray;
extern double setlamcrit;
extern int pinitfl;			/* =1 -> plotinit() run */
extern synap *synpnt;
extern comp *compnt;
extern photrec *recpnt;
extern recstim *recspnt;
extern recstim *recsend;
extern node *nodepnt;
extern elem *elempnt;
extern elem *elemend;
extern elem *oelemend;

extern Symbol *timeptr;
extern Symbol *plotiptr;
extern Symbol *ncompptr;

extern int iflag;                       /* interrupt flag set by onintr()  */

plotfr plotnod[PLOTNODSIZ] = {0};       /* node numbers, vmin, vmax voltages */

int plotnum = -1;                       /* number of plots displayed */
int runyet=0;
int nocond = 0;				/* no condense */

scheduler sched;                        /* scheduler */

double record(int cnod1, int cnod2, int cnod3, int cnod4, int pmod, int pval);
double getcurv(synap *spnt, double val);
double rsens(photrec *rpnt, double wavel, int filt);
double ncabs(double x);
void initcomp(void);
void condense();
void prcomp();
void chksmall();
void setdomains();
void compcon();
void runsim(double timestep, int run);
void plotinit(int plotnum);
void runstim(double time, int nsteps, int bgend, int endfl);
void runsyn(int nsteps);
void runrecs(int t); 
double getrecval(int i);
void mplot(double val, double time, int plotnum, int i, int pflag);
void runcomp();
node *findnode(nodeint num1, nodeint num2, 
			nodeint num3, nodeint num4, const char *str);
int findphotrec (int cnod1, int cnod2, int cnod3, int cnod4,
			 photrec **rpnt, const char *s);
void execerror (const char *s, const char *t);
void getarg(int narg, int ptype, pdatum *p);   
void rchan(synap *sdyad, int stime, double chanligand);
int poisdev(double m);
double gasdev(void);
double gamdev(double a);
int readstim(double rtime);
void delrstim(recstim *rspnt);
conn *printconns(comp *pnt, int nprint);
comp *othercomp (comp *pnt,conn *cpnt);
char *prnode (int n1, int n2, int n3, int n4); 
char *prnode (node *npnt); 
elem *findelem (int num);
double gjalpha (double v, gj *gjpnt);
double gjbeta  (double v, gj *gjpnt);
double getcyca (gj *gjpnt);
double getcycg (gj *gjpnt);
double dochani (chan *cpnt, double critc);
void equilchans (void);
recstim *getstimp1 (double stime);
recstim *getstimp2 (double stime);
double sqrt(double val);
double pow(double x, double expon);
double callproc(Symbol *procp, int npar, double par1, double par2);
double getnt (comp *cpnt, int type);
double addnt(comp *cpnt, int nttype, double val);
int isnt (comp *cpnt, int type);
double getmesg(gj *gjpnt, int mesg);
char *findsym(int num);
chantype *getchantype(int ctype, int cnum);
double ccavoff(chan *cpnt);
double vext(conn *cpnt);
double ghkv (chantype *chtyp,double cao, double cai);
int smaxmin(double val);
int elemsiz(elem *epnt);
double callfunc(Symbol *funcp, int npar, double par1, double par2,
					 double par3, double par4);
void plotpen(int pen, int nplot);
void initrec();
void setvclamp (comp *cpnt, double val, int bgend);
double get_ioni (chan *ch, int ion);
double get_iono (chan *ch, int ion);
double get_ionp (chan *ch, int ion);
double setnt(comp *cpnt, int nttype, double val);
double runrodi (photreci *rpnti, double iflux);
double runrodi_orig(photreci *rpnti, double iflux);
double runchr (photrecc *rpntc, double iflux, double v);

/*------------------------------------*/

void initsim (void)
{
    int i;
    static int ocumrecst = 0;
    static double totspace,elemspace;
    elem *epnt;
    FILE *fpout;

#ifdef DEBUG
  if (debug & NCBASE) ncfprintf (stderr,"initsim\n");
  if (debug & NOCOND) nocond = 1;
#endif

    if (elemend != oelemend) {	/* if any new neural elements, also set by resetelem() */
      if (!runyet) {

	sched.init();		/* initialize scheduler */

        if (!vidmode && !disp_ray && info>0) {
		int i;
						  /* Calculate space used */
						  /* by elements */
        for (elemspace=i=0,epnt=elempnt; epnt; epnt=epnt->next,i++) {  
          /* ncfprintf (stderr,"%d type %s elemsiz %d\n", 
			i,findsym(epnt->ctype),elemsiz(epnt)); */
          elemspace += elemsiz(epnt);
        }
        totspace = elemspace + cumnode * sizeof(node) + 
		cumclst * sizeof(conlst) + cumattr * sizeof(chattrib) + 
		cumcattr * sizeof (cattrib) + cumnattr * sizeof(nattrib); 
	totspace /= MB;
        if (info >= 1) {
	  if     (cumelem> 1) ncfprintf (stdout,
			"#c %d neural elements (%.3g MB)",
					cumelem,totspace);
	  else if (cumelem==1) ncfprintf (stdout,
			"#c %d neural element (%.3g MB)",
					cumelem,totspace);
          ncfprintf (stdout,
			" converted to compartments.\n");
	  if (cumrecst>0) ncfprintf (stdout,
			"#c Stimuli %.3g MB.\n",
					cumrecst*sizeof(recstim)/MB);
        }
       }
      }
      initcomp(); 
      if (!nocond) condense();
      chksmall();		/* check for small compartments */
      setdomains();		/* set region numbers */
      equilchans();		/* equilibrate all compartment membr channels */
      compcon();		/* set up time interval in compartment */
      ncompptr->val = cumcomp;	/* set number of comps for user */
      if (!runyet) {
	fpout = stdout;
        //if (!vidmode && !disp_ray && info>0) {
        if (!disp_ray && info>0) {
	     int remelem;

	  remelem = cumelem - delelmt;
	  if      (remelem> 1) ncfprintf (stdout,
			"#c %d neural elements saved.\n",remelem);
          else if (remelem==1) ncfprintf (stdout,
			"#c %d neural element saved.\n",remelem);
	  elemspace = 0;
	  if (remelem) {
             for (i=0,epnt=elempnt; epnt; epnt=epnt->next,i++) {  
                elemspace += elemsiz(epnt);
                /* ncfprintf (stderr,"%d type %s elemsiz %d elemspace %g\n", 
			i,findsym(epnt->ctype),elemsiz(epnt),elemspace); /* */
             }
	  }
          totspace =  ((double)cumcomp*sizeof(comp)) + cumnode * sizeof(node) + 
		      cumattr * sizeof(chattrib) + cumcattr * sizeof (cattrib) + 
		      cumnattr * sizeof(nattrib) +
		      cumconn * sizeof(conn) + cumrand * sizeof(randstate) +
		      cumchan * (sizeof(chan) + 4*sizeof(stconc))  +
		      cumntcomp * sizeof(ntcomp) + cumload * sizeof(load) + 
		      cumphotrec * sizeof(photrec) + cumrecst * sizeof(recstim) +
		      cumclst * sizeof(conlst) + cumdlst * sizeof(dyadlst) +
		      cumsynap * sizeof(synap) + 
		      cumlpfilt * (sizeof(lpfilt)+2*5*sizeof(double)) +
		      cumcacomp * (sizeof(cacomp)+ 11*sizeof(double)) +
		      cumgj * sizeof(gj) + /* */
		      elemspace ;
	  totspace /= MB;
          if (info >= 1) {
	    if (vidmode) fpout = stderr;	// with vidmode, print number of compartments to stderr
	    if (info >= 2)
	       ncfprintf (fpout,
			"#c Total memory space used %.3g MB.\n",totspace);
	    if      (cumcomp> 1) ncfprintf (fpout,
			"#c %d comps",cumcomp);
            else if (cumcomp==1) ncfprintf (fpout,
			"#c %d comp",cumcomp);
	    if      (ndomains> 1) ncfprintf (fpout,
			", %d cell domains",ndomains);
            else if (ndomains==1) ncfprintf (fpout,
			", %d cell domain",ndomains);
	    if      (cumnode> 1) ncfprintf (fpout,
			", %d nodes",cumnode);
            else if (cumnode==1) ncfprintf (fpout,
			", %d node",cumnode);
	    if      (cumconn> 1) ncfprintf (fpout,
			", %d connections",cumconn);
            else if (cumconn==1) ncfprintf (fpout,", %d connection",cumconn);
	    if      (cumchan> 1) ncfprintf (fpout,", %d channels",cumchan);
            else if (cumchan==1) ncfprintf (fpout,", %d channel",cumchan);
	    if      (cumsynap> 1) ncfprintf (fpout,", %d synapses",cumsynap);
            else if (cumsynap==1) ncfprintf (fpout,", %d synapse",cumsynap);
            if      (cumgj>1)     ncfprintf (fpout,", %d gjs",cumgj);
	    else if (cumgj==1)    ncfprintf (fpout,", %d gj",cumgj);
	    if      (cumphotrec> 1) ncfprintf (fpout,", %d photorcprs",cumphotrec);
            else if (cumphotrec==1) ncfprintf (fpout,", %d photorcpr",cumphotrec);
            ncfprintf (fpout,".\n");
	    fflush (fpout);
          }
        }
       }	/* runyet = 0 */

      if (prmap & (PCOND|PCOMP)) prcomp();
    }
  else if (cumrecst!=ocumrecst) {
      if (info >= 3) {
        ncfprintf (stdout,"#c Stimulus memory space   %.3g MB\n",
				cumrecst*sizeof(recstim)/1048576.0);
        ncfprintf (stdout,"#c Total memory space used %.3g MB\n",
				totspace+cumrecst*sizeof(recstim)/1048576.0);
      }

  }
  oelemend = elemend;
  ocumrecst = cumrecst;
}
 
/*------------------------------------*/

void chksmall(void)

/* Check for possible errors in constructing compartments. */

{
   comp *pnt,*pnt2;
   node *npnt,*onpnt;
   int found;
   conlst *lpnt;
   conn *cpnt;

  for (pnt=compnt; pnt; pnt=pnt->next) {
    if (pnt->rm > (LARGERES * .9) || pnt->cap < SMALLCAP * 1.1) {
       ncfprintf(stderr,"# nc: Warning: possible compartment build error: %u.\n",
		pnt->num);
       ncfprintf(stdout,"# Warning: possible error:\n");
       ncfprintf(stdout,"# Compartment # %u is missing: ",pnt->num);
       ncfprintf(stdout," gm %6.3g, cap %6.3g",pnt->rm,pnt->cap);
       for (found=0,lpnt=pnt->nodlst; lpnt; lpnt=lpnt->next) {
            if (npnt=(node*)lpnt->conpnt) {
                ncfprintf (stdout," at node %s.\n",
			prnode(npnt->nodenm1,npnt->nodenm2,
				npnt->nodenm3,npnt->nodenm4));
                 found = 1;
            }
       }
       if (!found) 
         for (found=0,npnt=nodepnt; npnt; npnt=npnt->next) {    /* find node */
           if (npnt->comptr == pnt) { 
                ncfprintf (stdout," at node %s.\n",
			prnode(npnt->nodenm1,npnt->nodenm2,
				npnt->nodenm3,npnt->nodenm4));
	       found = 1;
               break; 
           }
       }
      if (found) 
         onpnt = npnt;
      else {
         ncfprintf (stdout,"# Small comp without node. ");
         onpnt = (node *)NULL;
      }
      ncfprintf (stdout,"# Check: \"");
      if (cpnt=printconns(pnt,1)) {
        pnt2 = othercomp (pnt,cpnt);		/* get other compartment */
        if (pnt2) {				/*  find its nodes */
          ncfprintf (stdout,"\" which connects");
           for (found=0,lpnt=pnt2->nodlst; lpnt; lpnt=lpnt->next) {
            if (npnt=(node*)lpnt->conpnt) {
                ncfprintf (stdout," to node %s.\n",
			prnode(npnt->nodenm1,npnt->nodenm2,
				npnt->nodenm3,npnt->nodenm4));
                 found = 1;
            }
         }
        }  /* if pnt2 */ } else ncfprintf (stdout,"\"\n");
       if (onpnt) ncfprintf (stdout,"# Node %s",
			prnode(onpnt->nodenm1,onpnt->nodenm2,
			onpnt->nodenm3,onpnt->nodenm4));
       else ncfprintf (stdout,"# It"); 
       ncfprintf (stdout," should also connect to a cable, sphere, or load.\n");
       fflush (stdout);

    }  /*  if (pnt-rm...) */
  }  /* for (pnt=compnt;; ) */

}

/*------------------------------------*/

static Symbol *onplpnt    = (Symbol *)NULL;
static void (*onplotpnt)()= NULL;

void __runonplot(void)
{
   double d=0;

   if (onplpnt) {
       timeptr->val = simtime;             /* save correct time for user */
       callproc (onplpnt,0,d,d);
       simtime = timeptr->val; 		   /* reset time */
       ploti = plotiptr->val;              /* reset ploti (in case changed) */
   }
   if (onplotpnt) onplotpnt(); 
}

/*------------------------------------*/

void setonplot(void (*plpnt)()) 

{
  onplotpnt = plpnt;
}


/*------------------------------------*/

void setonpl(void) 

/* look for interpreter defined onplot routine */

{
    Symbol *f;

  if (f=lookup("onplot")) onplpnt = f;
  else                    onplpnt = (Symbol*)NULL;
}
/*------------------------------------*/

static Symbol *run_on_st_pnt = (Symbol *)NULL;
static void (*run_on_step_pnt)()= NULL;

void __runonstep(void)

/* call either interpreted or compiled runtime routine */

{
     double d=0;

   if (run_on_st_pnt) { callproc (run_on_st_pnt,0,d,d); }
   if (run_on_step_pnt) run_on_step_pnt();
}

/*------------------------------------*/

void set_run_on_step(void (*runonstpnt)())

/* set pointer to compiled runtime routine */

{
   run_on_step_pnt = runonstpnt;
}


/*------------------------------------*/

void set_run_on_st(void)

/* look for interpreter defined runtime routine */

{
    Symbol *f;

  if (f=lookup("run_on_step")) run_on_st_pnt = f;
  else                         run_on_st_pnt = (Symbol*)NULL;
}

/*------------------------------------*/

void set_run_on_ex(void);   /* defined in ncm.cc */ 

void runsim(double timestep, int runfl)
                    
/* Run the simulation; 
   schedule the photrecs every msec,
   schedule the synapses to run every 100 microsec,
   and run the compartments every "timestep".
   Plot the nodes specified using "record".
*/

{
  static double ptime, oldtime=-1e38;
  static double stoptime;
  static double time5, otime5, time4, otime4;
  static int nsteps;
  extern double endexp;
  int i,newtim;

#ifdef DEBUG
 if (debug & NCBASE) ncfprintf (stderr,"runsim step %g runfl %d runyet %d\n",
			timestep,runfl,runyet);
#endif

#define ROUNDUP (stiminc * 0.001)
#define STIME5 (-1e10)

 initsim();				/* print out simulation usage  */
 simtime = timeptr->val;                /* get correct time from user */
 newtim = (simtime!=oldtime);
 oldtime = simtime;		/* initial time setting from user */
 setonpl();			/* look for user-defined plot routine */
 set_run_on_ex();		/* look for user-defined run_on_exit routine */
 set_run_on_st();		/* look for user-defined run_on_step routine */
 if (plotnum>=0 && !pinitfl)  plotinit(plotnum+1);
 if (newtim) {
      ptime = ploti;
      otime4 = otime5 = STIME5;
 }
 if (runfl) {          /* run mode: run until endexp */
   stoptime = endexp + 1e-10;
 }
 else {                          /* step mode */
   stoptime = oldtime + timestep - 1e-10;	/* stop just before end */
/*   if (stoptime > endexp) stoptime = endexp; */
 }
 
 if ((timinc / stiminc) > MAXLONG) {
	timinc = (MAXLONG - 1000) * stiminc;
	ncfprintf (stdout,
		"# runsim: 'timinc' is too large, limited to %9.5g\n",timinc);
 }

 for (simtime=oldtime; simtime<stoptime; simtime+=timinc, ptime+=timinc)
  {
    if (debug & NCBASE) ncfprintf (stderr,"runsim: simtime %g\n",simtime);

    time5 = floor (simtime / stiminc + ROUNDUP);
  /* ncfprintf (stderr,"time5 %g otime5 %g simtime %20.18g\n",time5,otime5,simtime); /* */
    if (time5 > otime5) {
  /* ncfprintf (stderr,"time5 %g otime5 %g simtime %20.18g\n",time5,otime5,simtime); /* */
        nsteps = (int)(time5 - otime5);	/* calculate steps actually taken */
        otime5 = time5;
	nsteps = 1;
        runstim(simtime,nsteps,0,0);      	/* run stimuli */
        runsyn(nsteps);			/* run synapses  */
        runrecs(nsteps);                /* run photrecs */
        time4 = floor (time5 * stiminc / MSEC);
/* ncfprintf (stderr,"time4 %g otime4 %g\n",time4,otime4); /* */
        if (time4 > otime4) {
                otime4 = time4;
                //runrecs(1);            /* run photrecs */
        }
    }
    runcomp();                           /* run compartments */
    __runonstep();			 /* run user-defined runtime subroutine */

    if (ptime >= ploti)
     {
      ptime = 1e-8; 
      __runonplot();			/* run user-defined plot */
      for (i=0; i<=plotnum; i++) {      /* possibly change colors */
		Symbol *vpen;
		double (*vpenc)(int,double,double);
		int pen;
		double yval;

	 yval = getrecval(i);
         plotnod[i].yval = yval;
         if (vpen=(Symbol*)plotnod[i].vpen) {
	    pen = int(callfunc (vpen,3,i,simtime,yval,0.0));
	    if (pen != plotnod[i].ppen) plotpen (pen,i);
	 }
         if (vpenc=(double(*)(int,double,double)) plotnod[i].vpenc) {
	    pen = (int)vpenc(i,simtime,yval);
	    if (pen != plotnod[i].ppen) plotpen (pen,i);
	 }
      }
      for (i=0; i<=plotnum; i++) {      /* display plots */

         if (plotnod[i].pmod != GREC)
		mplot (plotnod[i].yval,simtime,plotnum+1,i,1);
      }
    }

    if (iflag)
     {
       //iflag = F;
       break;                           /* stop on interrupt */
     }
     sched.step();			/* run scheduler */

  }             /* for (time;;) */

  if (stoptime + 2e-10 >= endexp) {
    ptime = 1e-8;
    __runonplot();		 	 /* run user-defined plot */
    for (i=0; i<=plotnum; i++)           /* display plots one last time */
      if (plotnod[i].pmod != GREC)
          mplot (getrecval(i),simtime,plotnum+1,i,1);
  }
  runstim(simtime,nsteps,1,0);		/* run stimuli one more time to run stim end */
  oldtime = simtime;			/* save time for next runsim() */
  if (runfl && (simtime+ROUNDUP >= endexp)) {
     runyet=0;				/* reset run flag for init */
     //runstim(1e20,1,1,1);      		/* run stimuli one last time */
  }
  else runyet=1;
  timeptr->val = simtime;               /* save correct time for user */
  timeptr->vtype = NUMBER; 

#ifdef DEBUG
  if (debug & NCBASE) ncfprintf (stderr,"runsim end\n");
#endif
}


/*------------------------------------*/

double record (int cnod1, int cnod2, int cnod3, int cnod4, int pmod, int pval)

{
  node *npnt;
  photrec *rpnt;
  double rval, recsynapse(int snum, int nf), recchan(int cnum, int nf, int pval);
  comp *cpnt;

#ifdef DEBUG
  if (debug & NCRECORD) ncfprintf (stderr,"record node %d %d %d %d mode %d %d ",
					cnod1,cnod2,cnod3,cnod4,pmod,pval);
#endif

   rval = 0.0;
   switch (pmod) {

   case VREC: if (!(npnt=findnode(cnod1,cnod2,cnod3,cnod4,
					"record V"))) return(0.0);
	      if (npnt->comptr) rval = npnt->comptr->v;
              break;

   case IREC: if (!(npnt=findnode(cnod1,cnod2,cnod3,cnod4,
					"record I"))) return(0.0);
              if ((cpnt=npnt->comptr)) {
     	         if (cpnt->miscfl & VEXT)             /* comp is v-clamped */
                     rval = cpnt->extvi;
	      else //  if (cpnt->miscfl & IEXT)           /* comp is i-clamped */
                     rval = cpnt->exti;
	      }
              break;

   case MREC: if (!(npnt=findnode(cnod1,cnod2,cnod3,cnod4,
					"record M"))) return(0.0);
              if ((cpnt=npnt->comptr)) {
                     rval = (cpnt->v - cpnt->vrev) * cpnt->rm;  // cpnt->rm = conductance
	      }
              break;

   case LREC: if (findphotrec (cnod1, cnod2, cnod3, cnod4, &rpnt, "record L"))
                    rval = rpnt->iflux / stiminc;
              break;

   case WREC: if (!(npnt=findnode(cnod1,cnod2,cnod3,cnod4,
					"record Vm"))) return(0.0);
	      if ((cpnt=npnt->comptr)) {
	         if (cpnt->pvext!=NULL)
			 rval = cpnt->v - cpnt->pvext->v; // Vm = comp voltage - vext
		 else
			 rval = cpnt->v; 		  // Vm = comp voltage
	      }
              break;

   case CREC_GLU:  case CREC_AMPA: case CREC_KAINATE: case CREC_NMDA: case CREC_CNQX:
   case CREC_GABA: case CREC_BIC:  case CREC_PTX:  case CREC_GLY:
   case CREC_STRY:  case CREC_CAMP: case CREC_CGMP: case CREC_PH: case CREC_ATP:

         if (!(npnt=findnode(cnod1,cnod2,cnod3,cnod4,
					"record nt"))) return(0.0);
	      if (npnt->comptr) {
		  int typ;

		typ = pmod - CREC_GLU + GLU;
		rval = getnt (npnt->comptr,typ);
	      }
              break;
 
   case CREC_CACONC:
   case CREC_CAS: case CREC_CICR:
   case CREC_CAV: case CREC_CAIT: case CREC_CAIP: 
   case CREC_CAIE: case CREC_CAIPE: case CREC_CABUF: case CREC_CABUFB:

         if (!(npnt=findnode(cnod1,cnod2,cnod3,cnod4,
					"record ca"))) return(0.0);
	      if (npnt->comptr) {
		    cacomp *capnt;

		if (capnt=npnt->comptr->capnt) {
		   switch (pmod) {

		case CREC_CABUF: if (capnt->cab) {
		    		  if (pval==0) rval = capnt->cao;
				  else if (pval>(capnt->cashell-1)) pval = capnt->cashell;
				  rval = capnt->cab[pval-1]; 
				 }
				break;
		case CREC_CABUFB: if (capnt->cab) {	/* bound buffer (indicator) */
		    		  if (pval==0) rval = capnt->cao;
				  else if (pval>(capnt->cashell-1)) pval = capnt->cashell;
				  if (pval==1) rval = capnt->cabti - capnt->cab[pval-1]; 
				  else         rval = capnt->cabt  - capnt->cab[pval-1]; 
				 }
				break;
		case CREC_CACONC: 
		    		if (pval<=0) {
				   int pvaln;
				   if (pval==0) pvaln = 1;
				   else         pvaln = -pval;
				   if (capnt->caos!=NULL) {
				     if (pvaln > (capnt->caoshell-1))
				          rval = capnt->caos[capnt->caoshell-1];
				     else rval = capnt->caos[pvaln-1];
				   }
				   else {
				          rval = capnt->cao;
				   }
				}
				else {
				     if (pval>(capnt->cashell-1))
				     rval = capnt->cais[capnt->cashell-1];
				     else rval = capnt->cais[pval-1];
				} break;
		    case CREC_CAV:
				rval = capnt->vrev; break;
		    case CREC_CAIT: 
				rval = -capnt->ica;   break;
		    case CREC_CAIPE: 
				rval = -(capnt->ipump+capnt->iexch); break;
		    case CREC_CAIE: 
				rval = -capnt->iexch; break;
		    case CREC_CAIP: 
				rval = -capnt->ipump; break;
		    case CREC_CAS: 
				rval = -capnt->cas; break; /* CICR Ca store */
		    case CREC_CICR: 
				rval = -capnt->cicrflux; break; /* flux from CICR store */
		   }
		}
		else {
		    static int ca_err_printed = 0;
		 if (!ca_err_printed) {
	  ncfprintf (stderr,"# nc: error, node %s does not have a Ca compartment.\n", prnode(npnt));
	           ca_err_printed = 1;
	         }
		}
	      } 
	      break; 
   case NRECA0: case NRECA1: case NRECA2: case NRECA3: case NRECA4: 
   case NRECA8: case NRECA9: 
   case NRECH0: case NRECH1: case NRECH2: case NRECH3: case NRECH4:
   case NRECB0: case NRECB1: case NRECB2: case NRECB3: case NRECB4:
   case NRECC0: case NRECC1: case NRECC2: case NRECC3: case NRECC4:
   case NRECC9: case NRECCA: 
	  rval = recsynapse(cnod1,pmod);
          break;

   case NRECG: case NRECGI: case NRECGV: case NRECGC: case NRECGM: case NRECGH:
   case NRECG0: case NRECG1: case NRECG2: case NRECG3:
	 rval = recchan(cnod1,pmod,pval);
         break;

   }    /* switch */

#ifdef DEBUG
  if (debug & NCRECORD) ncfprintf (stderr," rval %g ...end\n",rval);
#endif
 
 return (rval);
}

/*------------------------------------*/

double recsynap(synap *spnt, int nf, int snum)
           
/* Record from a synapse.
    Partly copied from modsynap() in ncmak.c.
*/

{
    elem *elpnt;
    cacomp *capnt;

    if (spnt==NULL) {
       //ncfprintf (stderr,"recsyn: can't find synapse for element %d\n",snum);
       return 0.0;
    }
    if (spnt->ctype != SYNAPSE) {
       ncfprintf (stderr,"recsyn: element %d is not a synapse\n",snum);
       return (0.0);  
    }
    if (spnt->spre && nf < NRECB0) {  /* if dyad, use other syanpse */
           spnt = spnt->spre->sdyad;/* for presyn filter A */
           if (!spnt) return 0.0;
    }
    if (nf >= NRECA0 && nf < NRECA8) {	/* look at filters before nt release */
      nf = nf - NRECA0;			/* look at any of the NUMFILT filts */
      if (spnt->filt1) {
/* ncfprintf (stderr,"f0 %d snum %d nf %d\n",spnt->filt1->lfilt[nf],snum,nf); */
      if      (nf<0) nf = 0;
      else if (nf>=spnt->filt1->nfilt) nf = spnt->filt1->nfilt-1;

        return (double)spnt->filt1->lfilt[nf]; /* voltage before filter */
      }
      else return 0.0;
   }
   else if (nf == NRECA8) {		/* look at readily releasible pool */
        return spnt->rrpool;
   }
   else if (nf == NRECA9) {		/* look at vesicle rate */
	double trvolts,fltvolts,transrate,vsize;
        int filt;

    switch (spnt->sens) {

    case V:
      if (spnt->filt1) {			/* first check for filter 1 */
            filt = spnt->filt1->nfilt-1;	/* number of filters */
            transrate = (double)spnt->filt1->lfilt[filt];  /* last filter */
      }
      else {					/* otherwise recalc tr */
	  if (oldsynapz)
	    transrate = spnt->comp1->v - spnt->comp1->vph - spnt->comp1->vext - spnt->thresh;
	  else {
	    trvolts = spnt->comp1->v - spnt->comp1->vph - spnt->comp1->vext - spnt->thresh;
	    transrate = getcurv(spnt,trvolts);
	  }
      }
      if (spnt->filt1h) {
	   double hfrac, lfrac;	

	fltvolts = transrate;
        hfrac = spnt->filt1h->tfall; 
	if (hfrac < 0) {
	   hfrac = -hfrac;
           lfrac = 1 - hfrac;
           if (lfrac < 0) lfrac = 0;
           fltvolts *= lfrac;
	}
        filt = spnt->filt1h->nfilt-1;
        fltvolts += spnt->filt1h->lfilt[filt] * hfrac + spnt->filt1h->offset;
	if (!oldsynapz) {
	  if (fltvolts < 0) fltvolts = 0;
	}
	transrate = fltvolts;
      }
      break;

    case CA:
      if (capnt=spnt->comp1->capnt) {
	      double caval;
	caval = capnt->cais[0];		 /* take calcium in first shell */
        //if (spnt->caegain != 1) fprintf (stderr,"caegain %g\n", spnt->caegain);
        transrate = exp(log(caval*dscavg)*spnt->caegain) * spnt->vgain; /* scale by vgain */
	//fprintf (stderr,"ca-trans %g\n",transrate);
      }
      else transrate = 0;
      break;

    } /* switch */

      if (oldsynapz) transrate = getcurv(spnt,transrate);

      if (spnt->maxsrate) 
        transrate *= spnt->rrpoolg + (1.0 - spnt->rrpoolg) * spnt->rrpool / spnt->maxrrpool;

      if ((vsize=spnt->vsize)<=0) {
	  if ((vsize=dvsz)<=0) vsize = 1;
      }
      transrate /= (vsize * 10.0 * HMICSEC);
//ncfprintf (stderr,"transrate %g rrpool %g maxrrpool %g\n",
//			transrate,spnt->rrpool,spnt->maxrrpool);
      return transrate;
   }
   else if (nf >= NRECH0 && nf <= NRECH4) {      /* high pass filters */
     nf = nf - NRECH0;			/* look at any of the NUMFILT filters */
     if (spnt->filt1h) {
/* ncfprintf (stderr,"f0h %d snum %d nf %d\n",spnt->filt1h->lfilt[nf],snum,nf); */
       if      (nf<0) nf = 0;
       else if (nf>=spnt->filt1h->nfilt) nf = spnt->filt1h->nfilt-1;

       return spnt->filt1h->lfilt[nf];	/* voltage before filter */
     }
     else return 0.0;
   }
   else if (nf < NRECC0) {		/* look at filters after nt release */
      nf = nf - NRECB0;			/* look at any of the NUMFILT filters */
      if (spnt->filt2) {
/* ncfprintf (stderr,"f0 %d snum %d nf %d\n",spnt->filt2->lfilt[nf],snum,nf); */
        if      (nf<0) nf = 0;
        else if (nf>=spnt->filt2->nfilt) nf = spnt->filt2->nfilt-1;
        return spnt->filt2->lfilt[nf];	/* neurotrans from filter */
      }
      else return 0.0;
   }
   else if (nf < NRECC9) { 		/* look at filters after nt sat. */
      nf = nf - NRECC0;			/* look at any of the NUMFILT filters */
      if (spnt->filt3) {
/* ncfprintf (stderr,"f0 %d snum %d nf %d\n",spnt->filt3->lfilt[nf],snum,nf); */
        if      (nf<0) nf = 0;
        else if (nf>=spnt->filt3->nfilt) nf = spnt->filt3->nfilt-1;
        return (double)spnt->filt3->lfilt[nf];
      }
      else {		/* return normalized conductance if nfilt3 == 0 */
           double c;
        if (nf==0) {
	  if (spnt->resp1)          return (spnt->conduct);
	  else if (spnt->resp2)     return (spnt->conduct);
	  else if (c=spnt->maxcond) return (spnt->conduct/c);
	  else return  spnt->conduct;
        }
        else return 0.0;
      }
   }
   else if (nf==NRECC9) {		/* return 2nd mesg [cGMP] */

   /* Problem here is that [cGMP] is not saved in synapse structure.
      We compute it from either the third filter or from the output
      conductance, using the appropriate inverse function. */

        double c,cycg,kk,p;

      if (spnt->filt3) {
        nf = spnt->filt3->nfilt-1;
        c = spnt->filt3->lfilt[nf];
	if (spnt->ntact==OPEN) {
          cycg = c;
	} else {		/* CLOSE */
          cycg =  spnt->coff - c * spnt->cgain;
	  if (cycg < 0) cycg = 0;
	}
        return cycg;
      } 
				/* if chan and secmsg is not used */
     else if (spnt->resp1 && (!spnt->secmsg || !spnt->resp1->chtyp->cgbnd)) {
         return spnt->conduct;

     } else {		    /* secmsg is used */
         c = spnt->conduct;
         if (!spnt->resp1) c /= spnt->maxcond;	/* reg syn, normalize to 1 */
	if (spnt->ntact==OPEN) {
          cycg = c;
	} else {		/* CLOSE, compute inverse sat func. */
         if ((p=spnt->chc) <=0) p = 1;   /* see hillcalcn() */
         kk = pow(spnt->ckd,p);
         cycg = pow(c * kk / (kk - c + 1),1.0/p);
	 if (cycg < 0) cycg = 0;
	}
        return cycg;
     }
  } else if (nf==NRECCA) {		/* return presyn [Ca]i */
      if (capnt=spnt->comp1->capnt) 
	return capnt->cais[0];		 /* take calcium in first shell */
  }
}

/*------------------------------------*/

double recsynapse(int snum, int nf)
           
/* Record from a synapse.
    Partly copied from modsynap() in ncmak.c.
*/

{
    elem *elpnt;
    synap *spnt;

    if (!(elpnt=findelem(snum))) {
       ncfprintf (stderr,"recsyn: can't find element %d\n",snum);
       execerror ("Missing element: "," stopping... ");
       return (0.0);
    }

    spnt = (synap *)elpnt->lptr;	/* pointer to synapse from elem */
    return recsynap(spnt,nf,snum);
}

/*------------------------------------*/


double rechan(conn *cpnt, int nf, int pval, int cnum)

/* Record from a channel given a pointer to it */
/*   Partly copied from recsynap() above. */

{
    int secmsg;
    double v;
    chan *chpnt;
    ntchan *ntchpnt;
    attrib *apnt;
    cacomp *capnt;
    int ncond,nstate;
    double cval=0;
    synap *spnt;
    photrec *rpnt;
    phvars *vpnt;

 if (cpnt==NULL) return (0);
 switch (cpnt->ctype) {

  case SYNAPSE:			/* record post-synaptic conductance */

   spnt = (synap *)cpnt;	/* get pointer to synapse from elem */
   if (spnt==NULL) {
       //ncfprintf (stderr,"rechan: can't find synapse for element %d\n",cnum);
       return 0.0;
   }
  secmsg = spnt->secmsg;
  switch (spnt->ctype) {	/* check that it is really a synapse */

   default:
           ncfprintf (stderr,"recsyn: element %d is not a synapse\n",cnum);
           return (0.0);  
	   break;

   case SYNAPSE:
    if (secmsg) ntchpnt = spnt->resp2;
    else	ntchpnt = spnt->resp1;
    if (ntchpnt==NULL) {
     if (nf==NRECG0 || (nf==NRECG && pval==0)) { /* total channel conductance */
	if (spnt->resp1 && !secmsg) cval = spnt->resp1->conduct;
	else if (spnt->resp2) 	    cval = spnt->resp2->conduct;
	else   		            cval = spnt->conduct;
       break;
     }
     if (nf==NRECGI) {			/* total ionic current */
	 v = spnt->vrev - spnt->comp2->v;
	if (spnt->resp1 && !secmsg) cval = spnt->resp1->conduct * -v;
	else if (spnt->resp2) 	    cval = spnt->resp2->conduct * -v;
	else            	    cval = spnt->conduct * -v;
	break;
     }
     if (nf==NRECGV) {
	 cval = spnt->vrev; 
         break;
     }
    }
    else {
    switch (ntchpnt->ctype) {	
     case AMPA:		/* Markov post-synaptic chans */
     case KAINATE:	
     case CGMP:	
     case GABA:
     case GLY:
     case SYN2:
     case NMDA:
     default:
//      ncond = nf - NRECG0;	/* conductance number */
//      if (ncond > NGREC) ncond = NGREC;
//      else if (ncond < 1) ncond = 1;

      if (pval < 0) pval = 0;

      switch (nf) {

	default: cval = 0; break;
	case NRECGV: cval = ntchpnt->vrev; break;
	case NRECG: 
	  if (pval==0) cval = ntchpnt->conduct;
	  else {
  	    nstate = ntchpnt->numstate;
	    if (pval > nstate) pval = nstate;
	    if (ntchpnt->nchan <= 0)
	       cval = ntchpnt->conc[pval-1].cest;
	    else
	       cval = ((double)ntchpnt->conc[pval-1].nchan) / ntchpnt->nchan;
	  }
	  break;

	  case NRECGI:			/* total ionic current */
	      v = ntchpnt->vrev - (ntchpnt->comp1->v - vext(ntchpnt));
	      cval = -ntchpnt->conduct * v;
	     break;
	  case NRECGC: 		/* calcium fraction of total conductance */
	      cval = ntchpnt->conduct * get_ionp(ntchpnt,PCA);
	     break; 

       }  /* switch (ncond) */
       break;

      } /* switch (ntchpnt->ctype) */
    }
    break;

   } /* switch (spnt->ctype==SYNAPSE) */
    break;

 case GJ: 				/* record gj conductance */
   if (nf==NRECG0) {			/* total channel conductance */
     spnt = (synap *)cpnt;		/* get pointer to synapse from elem */
     cval = spnt->conduct;
   }
    break;

 case ROD:
 case CONE:
 case CHR:
    rpnt = (photrec *)cpnt;		/* get pointer to photorec from elem */
    switch (nf) {
      default: cval = 0; break;
      case NRECG0: 
      case NRECG:
      if (pval > 11) pval = 11;
      vpnt = &rpnt->vars;

      if (rpnt->chtyp->stype==PINVG) {	/* mouse rod from invergo et al 2014 */
	      photreci *rpnti;
	      phvarsi *vpnti;

	 rpnti = (photreci *)rpnt;
         vpnti = &rpnti->vars;
	 if (nf==NRECG0) cval=rpnt->conduct; break;
         switch (pval) {

           case 0: cval = rpnt->conduct;   break;
           case 1: cval = vpnti->Rn[0];    break;
           case 2: cval = vpnti->Rn_Gt[0]; break;
           case 3: cval = vpnti->Rn_RKpost[1]; break;
           case 4: cval = vpnti->Gagtp;    break;
           case 5: cval = vpnti->PDEa_Gagtp; break;
           case 6: cval = vpnti->RGS_PDEa_Gagtp;  break;
           case 7: cval = vpnti->cGMP;     break;
           case 8: cval = vpnti->Recr_Ca;  break;
           case 9: cval = vpnti->Cafree;   break;
           case 10:cval = vpnti->Cabuff;   break;
	
         }
      } else if (rpnt->chtyp->stype==ChR2)  {	/* channel-rhodopsin */
	      photrecc *rpntc;
	      phvarsc *vpntc;
	      chan *cpnt;

	if ((cpnt=rpnt->phchan)==NULL) {	/* ChR2 in chr2.cc */
	    rpntc = (photrecc *)rpnt;
            vpntc = &rpntc->vars;
            switch (pval) {
              case 0: cval = rpntc->conduct;	    break;
              case 1: cval = vpntc->SC1;	    break;
              case 2: cval = vpntc->SC2;	    break;
              case 3: cval = vpntc->SO1;	    break;
              case 4: cval = vpntc->SO2;	    break;
              case 5: cval = vpntc->p;	    	    break;
              case 6: cval = vpntc->cond;	    break;
	    }
	} else {		// channel structure, not photrec
	     chpnt = (chrchan *)cpnt;
	     switch (nf) {
	        default: 
      		case NRECG0: cval = chpnt->conduct;   break;
		case NRECG:
		  nstate = chpnt->numstate;
		  if (pval > nstate) pval = nstate;
		  if (chpnt->nchan <= 0)
		    cval = chpnt->conc[pval-1].cest;
		  else
		    cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
		  break;
	     }
	}
      } else {				/* orig rod, cone transduction type */

	 if (nf==NRECG0) cval=rpnt->conduct; break;
         switch (pval) {

           case 0: cval = rpnt->conduct;  break;
           case 1: cval = vpnt->rhod;  break;
           case 2: cval = vpnt->rhod2; break;
           case 3: cval = vpnt->rstar; break;
           case 4: cval = vpnt->gpr;   break;
           case 5: cval = vpnt->pde;   break;
           case 6: cval = vpnt->gcyc;  break;
           case 7: cval = vpnt->cycg;  break;
           case 8: cval = vpnt->ca;    break;
           case 9: cval = vpnt->cax;   break;
           case 10:cval = vpnt->cab;   break;
           case 11:cval = vpnt->rhodk; break;
         }
    }		/* end of case NRECG */
    break;
   }
   break;

 default:
 case CACOMP:
 case CHAN:			/* record channel conductances */
 	
  if ((chpnt=(chan *)cpnt)==NULL) { /* get pnt to channel from elem */
       ncfprintf (stderr,"rechan: can't find channel for element %d\n",cnum);
       return 0.0;
  }

  /* ncfprintf(stderr,"chan type %d stype %d nf %d\n",
		chpnt->ctype,chpnt->stype,nf); /* */

  switch (chpnt->ctype) {	/* stop short on wrong element */
     case AMPA:			/* synaptic Markov chans */
     case KAINATE:
     case GABA:
     case GLY:
     case SYN2:
     case NMDA:
     case CGMP:

     case NA:			/* membrane ion chans */
     case K:
     case ClCa:
     case KCa:
     case CA: 
     case CACOMP: break;
     default:
       ncfprintf (stderr,"rechan: element %d is not a channel (%d)\n",
					cnum,chpnt->ctype);
       return (0.0);  
  }

  if (nf==NRECG0) {		/* total channel conductance */
     cval = chpnt->conduct;
  }
  else {    

//    ncond = nf - NRECG0;	/* conductance number */
//    if (ncond > NGREC) ncond = NGREC;
//    else if (ncond < 1) ncond = 1;

    switch (nf) {
        double nernstv(double val);

		/* params that are the same for all chans */

    case NRECGI: 
		 switch (chpnt->ctype) {
		   case CA: v = chpnt->gvrev; break;
		   default: v = chpnt->vrev; break;
		 }
		 v -= chpnt->comp1->v - vext(chpnt)-ccavoff(chpnt);
		 cval = chpnt->conduct * -v;
		 break;
    case NRECGV: switch (chpnt->ctype) { 
			chantype *chtyp;

		  case NMDA:
		  case CGMP:
		  case CA:
			if (chpnt->setvrev) cval = chpnt->vrev;
			else { 
			  chtyp = getchantype(chpnt->ctype,chpnt->stype);	
			  if ((capnt=chpnt->comp1->capnt) !=NULL) {
			     if (capnt->caos !=NULL) 
			        cval = ghkv(chtyp,capnt->caos[0],capnt->cais[0]);
			     else
			        cval = ghkv(chtyp,capnt->cao,capnt->cais[0]);
                           }
			  else cval = ghkv(chtyp,dcao,dcai);
			}
			break; 
		  default: cval = chpnt->vrev; break;
		}
		break;
    case NRECGC: cval = chpnt->conduct * get_ionp(chpnt,PCA); break;

    case NRECG:  if (pval==0) {
		 cval = chpnt->conduct;
		 break; 
		 }
		/* else if (pval > 0), continue on to next case: */

		/* params that are different for channel subtypes */
    default:
    switch (chpnt->ctype) {
     case NA:
	  switch (chpnt->stype) {
	    case 0:			/* m3h type Na channel */
	 	switch (nf) {
		    case NRECG1:
		    case NRECGM: cval = ((hhchan *)chpnt)->m; break;
		    case NRECG2:
		    case NRECGH: cval = ((hhchan *)chpnt)->h; break;
    		    case NRECG3: 
    		    case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
		                  cval = -chpnt->conduct * v; break;
		    default:     cval = 0.0; break;
		}
		break;
	    case 1:
	    case 2:
	    case 3:
	    default:
	      switch (nf) {
	        default: cval = 0; break;
		case NRECG:
		nstate = chpnt->numstate;
		if (pval > nstate) pval = nstate;
		if (chpnt->nchan <= 0)
		  cval = chpnt->conc[pval-1].cest;
		else
		  cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
		break;
		}
	  }
 	  break;
     case K:
	  switch (chpnt->stype) {
	    case 0:			/* n4 type K channel */
	 	switch (nf) {
		    case NRECG1:
		    case NRECGM: cval = ((hhchan *)chpnt)->m; break; /* n */
    		    case NRECG2: 
    		    case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
		                  cval = -chpnt->conduct * v; break;
		    default:     cval = 0.0; break;
		}
		break;
	    case 1:
	    case 3:			/* KA, markov */
	    case 4:			/* Ih */
	    case 5:
	    case 6:
	    case 7:
	    case 8:			/* Ih */
	    case 9:			/* Ih */
	    default:
	      switch (nf) {
	        default: cval = 0; break;
		case NRECG:
		nstate = chpnt->numstate;
		if (pval > nstate) pval = nstate;
		if (chpnt->nchan <= 0)
		  cval = chpnt->conc[pval-1].cest;
		else
		  cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
		break;
	      }
	      break;
	    case 2:			/* KA, HH */
	 	switch (nf) {
		    case NRECG1:
		    case NRECGM: cval = ((hhchan *)chpnt)->m; break; /* n */
		    case NRECG2:
		    case NRECGH: cval = ((hhchan *)chpnt)->h; break; /* h */
    		    case NRECG3: 
    		    case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
		                  cval = -chpnt->conduct * v; break;
		    default:    cval = 0.0; break;
		}
		break;
	  }
 	  break;

     case ClCa:
	      switch (nf) {
	        default: cval = 0; break;
		case NRECG:
		nstate = chpnt->numstate;
		if (pval > nstate) pval = nstate;
		if (chpnt->nchan <= 0)
		  cval = chpnt->conc[pval-1].cest;
		else
		  cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
		break;
		case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
		            cval = -chpnt->conduct * v; break;
	      }
	      break;

     case KCa:
	  switch (chpnt->stype) {
	    case 0:			/*  sKCa channel */
	    case 2:			/*  bKCa channel */
	 	switch (nf) {
		    case NRECG1:
		    case NRECGM: cval = ((hhchan *)chpnt)->m; break; /* n */
    		    case NRECG2: 
    		    case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
		                  cval = -chpnt->conduct * v; break;
		    default:     cval = 0.0; break;
		}
		break;
	    case 1:			/* sKCa */
	    case 3:			/* bKCa */
	    case 4:			/* sKCa */
	    case 5:			/* sKCa */
	    case 6:			/* bKCa */
	    default:
	      switch (nf) {
	        default: cval = 0; break;
		case NRECG:
		nstate = chpnt->numstate;
		if (pval > nstate) pval = nstate;
		if (chpnt->nchan <= 0)
		  cval = chpnt->conc[pval-1].cest;
		else
		  cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
		break;
	      }
	      break;
	  }
 	  break;
     case CA:
	capnt = chpnt->comp1->capnt;
	switch (chpnt->stype) {
        case 0:			/* m3 type Ca channel */
        case 1:			/* m3 type Ca channel */
 	switch (nf) {
	    case NRECG1:
	    case NRECGM: cval = ((hhchan *)chpnt)->m; break; /* m */
	    case NRECGH: cval = ((hhchan *)chpnt)->h; break; /* h */
    	    case NRECG2: 
    	    case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
	                  cval = -chpnt->conduct * v; break;
	    default:    cval = 0.0; break;
		}
		break;
	 case 2:
	 default:
	  switch (nf) {
	    default: cval = 0; break;
	    case NRECG:
	    nstate = chpnt->numstate;
	    if (pval > nstate) pval = nstate;
 	    if (chpnt->nchan <= 0)
	       cval = chpnt->conc[pval-1].cest;
	    else
	       cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
	    break;
    	    case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
	          cval = -chpnt->conduct * v; break;
	    }
	  break;
	}
        break;

     default:			/* any other channel type */
     case CGMP:
	  switch (chpnt->stype) {
	  case 1:
	  default:
	    switch (nf) {
	      default: cval = 0; break;
	      case NRECG:
	      nstate = chpnt->numstate;
	      if (pval > nstate) pval = nstate;
	      if (chpnt->nchan <= 0)
		  cval = chpnt->conc[pval-1].cest;
	      else
		  cval = ((double)chpnt->conc[pval-1].nchan) / chpnt->nchan;
              break;
    	     case NRECGI: v = chpnt->vrev - (chpnt->comp1->v-vext(chpnt)-ccavoff(chpnt));
	          cval = -chpnt->conduct * v; break;
	     }
 	  break;
    }  /* switch (chpnt->stype) */
    break;
    
    }  /* switch (chpnt->ctype) */
   }  /* switch (nf)    */
  } /* else nf > NRECG0 */

  break; 	/* case CHAN */

 }  /* switch (cpnt->ctype) */

 return cval;
}

/*------------------------------------*/

double recchan(elem *elpnt, int nf, int pval)

/* Record from a channel, given element pointer.  */

{
    int found;
    conn *cpnt;
    attrib *apnt;
    node *npnt;
    comp *pnt;
    conlst *lpnt;

 if (elpnt->ctype==CHAN) {
   if ((apnt=elpnt->attpnt)==NULL) {
    ncfprintf (stderr,"recchan: can't find channel attrib for element %d\n",elpnt->elnum);
       return 0.0;
   }
   cpnt = (conn *)apnt->lptr;
 }
 else cpnt = (conn *)elpnt->lptr;	/* get pointer to channel from elem */
 return rechan(cpnt,nf,pval,elpnt->elnum);
}

/*------------------------------------*/

double recchan(int cnum, int nf, int pval)

                /* cnum = element number of channel */

/* Record from a channel.  */

{
    int found;
    elem *elpnt;
    conn *cpnt;
    attrib *apnt;
    node *npnt;
    comp *pnt;
    conlst *lpnt;

 if (!(elpnt=findelem(cnum))) {
      ncfprintf (stderr,"recchan: can't find element %d\n",cnum);
      execerror ("Missing element: "," stopping... ");
      return 0.0;
 }
 return recchan(elpnt,nf,pval);
}

/*------------------------------------*/

double reccable(int elemnum, double fracdist)

/* record from a cable with element number "elemnum",
   at a fractional distance "fracdist" from the first
   end (the end that contacts node 1, the node given first).
*/

{
    int i,n,cnum;
    double r, v, vold;
    node *np1;
    elem *epnt;
    conn *cpnt,*cpnt1,*cpnt2,*tcpnt;
    comp *pnt,*pnt1,*pnt2;
    conlst *lpnt;


 if (!(epnt=findelem(elemnum))) {
	ncfprintf (stderr,"reccable: can't find element %d\n",elemnum);
	execerror ("Missing element: "," stopping... ");
	return 0.0;  
 }

 np1 = epnt->nodp1;
 if (!np1) {
         ncfprintf (stderr,"reccable: can't find node for element %d\n",elemnum);
         return 0.0;  
 }
 cpnt = (conn *)epnt->lptr;
 if (!cpnt) {
         ncfprintf (stderr,"reccable: can't find low-lev conn for element %d\n",
							elemnum);
  	 return 0.0;
 } 
 		/* find first node's compartment */

 pnt = np1->comptr;
 if (!pnt) {
      ncfprintf (stderr,"reccable: can't find comp for node %d\n",np1->nodenm1);
      return 0.0;
 } 
 		/* Find connection from first comp to second comp. */ 
		/* Must be same as connection saved in cable. */

 cpnt1 = (conn*)NULL;
 for (lpnt=pnt->clst; lpnt && (cpnt1 != cpnt); lpnt=lpnt->next) {
    if (!(cpnt1=lpnt->conpnt)) {
	ncfprintf (stderr,"reccable: missing conn pt\n");
	return 0.0;
    }
 }
 tcpnt = cpnt1;		/* remember the conn */

		/*  (Code below taken from dcomp() in "prcomp.c") */

		/* Make sure conn points back to orig comp. */

 if 	 ((pnt1=cpnt1->comp1)==pnt) pnt1 = cpnt1->comp2;
 else if ((pnt1=cpnt1->comp2)==pnt) pnt1 = cpnt1->comp1;
 else return 0.0;

 if (cpnt1->ctype!=AXIALRES) return 0.0;
 
 for (n=0; pnt1->nodlst==NULL; n++) {	       /* traverse cable, count comps */

    lpnt = pnt1->clst;
    if (!lpnt) break;
    if (!(cpnt2=lpnt->conpnt)) break;		 /* find conn */
    if (cpnt1==cpnt2) {
	lpnt  = lpnt->next;
	if (lpnt) cpnt2 = lpnt->conpnt; 	 /* wrong conn? */
    }
    if (!cpnt2) break;
    while (cpnt2->ctype!=AXIALRES) {
      lpnt  = lpnt->next;
      if (!lpnt) break;
      cpnt2 = lpnt->conpnt;
      if (!cpnt2) break;
    }

    if ((pnt2=cpnt2->comp1)==pnt1) pnt2 = cpnt2->comp2;
    pnt1=pnt2;
    cpnt1=cpnt2;

   }   /* for (n=0;;) */

   n++;					/* final comp at second node */
   cnum = (int)(fracdist * n);		/* integer number of requested comp */
  					/* end comps (0,n) are 1/2 size */
 
			/* Do it again; this time stop at comp "cnum" */
   cpnt1 = tcpnt;
   if 	   ((pnt1=cpnt1->comp1)==pnt) pnt1 = cpnt1->comp2;
   else if ((pnt1=cpnt1->comp2)==pnt) pnt1 = cpnt1->comp1;
   else return 0.0;

   if (cpnt1->ctype!=AXIALRES) return 0.0;

   vold = pnt->v;
   v = pnt1->v;
   for (i=1; i<=cnum; i++) {	       /* traverse cable, count comps */

      lpnt = pnt1->clst;
      if (!(cpnt2=lpnt->conpnt)) break;		 /* find conn */
      if (cpnt1==cpnt2) {
	  lpnt  = lpnt->next;
	  if (lpnt) cpnt2 = lpnt->conpnt; 	 /* wrong conn? */
      }
      if (!cpnt2) break;
      while (cpnt2->ctype!=AXIALRES) {
        lpnt  = lpnt->next;
        if (!lpnt) break;
        cpnt2 = lpnt->conpnt;
        if (!cpnt2) break;
      }

      if ((pnt2=cpnt2->comp1)==pnt1) pnt2 = cpnt2->comp2;
      vold = v;
      v = pnt2->v;
      pnt1=pnt2;
      cpnt1=cpnt2;

     }   /* for (i=0;;) */

   r = fracdist * n - cnum;		/* interpolate across axial resistors */
					/* between adjacent compartments */

/*  ncfprintf (stderr,"n %d fracdist %g cnum %d v %g vold %g r %g\n",
			n,fracdist,cnum,v,vold,r);	/* */
    
   return (vold * (1.0-r) + v * r);
}

/*------------------------------------*/

double getrecval (int plotnum)
{
  int cnod1,cnod2,cnod3,cnod4,pmod,pval;
  Symbol *spnt, *func;
  double (*funcc)(double pl, double t);
  int narg;
  pdatum p={0};
  double rval, *var, *garrayval(Symbol *s);
  double recsynapse(int snum, int nf),reccable(int elemnum, double fracdist);

#ifdef DEBUG
  if (debug & NCRECORD) ncfprintf (stderr,"getrecval plotnum %d\n",plotnum);
#endif

   cnod1 = smaxmin(plotnod[plotnum].cnod1);
   cnod2 = smaxmin(plotnod[plotnum].cnod2);
   cnod3 = smaxmin(plotnod[plotnum].cnod3);
   cnod4 = smaxmin(plotnod[plotnum].cnod4);

   pmod  = plotnod[plotnum].pmod;
   pval  = plotnod[plotnum].pval;

   switch (pmod) {

   case VREC:
   case IREC:
   case MREC:
   case LREC: 
   case WREC: 
   
   case CREC_GLU:  case CREC_AMPA: case CREC_KAINATE: case CREC_NMDA: case CREC_CNQX:
   case CREC_GABA: case CREC_BIC:  case CREC_PTX:  case CREC_GLY:
   case CREC_STRY:  case CREC_CAMP: case CREC_CGMP: case CREC_PH: case CREC_ATP:

   case CREC_CACONC:
   case CREC_CAS: case CREC_CICR:
   case CREC_CAV: case CREC_CAIT: case CREC_CAIP: case CREC_CAIE: 
   case CREC_CAIPE: case CREC_CABUF: case CREC_CABUFB:

             rval = record (cnod1,cnod2,cnod3,cnod4,pmod,pval);
             break;
   case GREC:
	     break;
   case FREC:			/* run user-defined plot func */ 

             if (func=(Symbol*)plotnod[plotnum].func) {
		  double plotvl;

	       if (plotnod[plotnum].plotval!=LARGENUM) 
		  plotvl = plotnod[plotnum].plotval;
	       else
		  plotvl = plotnum;
	       rval = callfunc (func,2,plotvl,simtime,0.0,0.0);
	     }
	     else if (funcc=(double(*)(double p, double t))plotnod[plotnum].funcc) {
		  double plotvl,plotvl2,plotvl3,plotvl4,plotvl5,plotvl6,plotvl7;
	          plotvl=plotnod[plotnum].plotval;
	          plotvl2=plotnod[plotnum].plotval2;
	          plotvl3=plotnod[plotnum].plotval3;
	          plotvl4=plotnod[plotnum].plotval4;
	          plotvl5=plotnod[plotnum].plotval5;
	          plotvl6=plotnod[plotnum].plotval6;
	          plotvl7=plotnod[plotnum].plotval7;

	       if (plotvl7!=LARGENUM) 
	          rval = ((double (*)(double,double,double,double,double,double,double,double))funcc) 
			  (plotvl,plotvl2,plotvl3,plotvl4,plotvl5,plotvl6,plotvl7,simtime);
	       else if (plotvl6!=LARGENUM) 
	          rval = ((double (*)(double,double,double,double,double,double,double))funcc) 
			  (plotvl,plotvl2,plotvl3,plotvl4,plotvl5,plotvl6,simtime);
	       else if (plotvl5!=LARGENUM) 
	          rval = ((double (*)(double,double,double,double,double,double))funcc) 
			  (plotvl,plotvl2,plotvl3,plotvl4,plotvl5,simtime);
	       else if (plotvl4!=LARGENUM) 
	          rval = ((double (*)(double,double,double,double,double))funcc) 
			  (plotvl,plotvl2,plotvl3,plotvl4,simtime);
	       else if (plotvl3!=LARGENUM) 
	          rval = ((double (*)(double,double,double,double))funcc) 
			  (plotvl,plotvl2,plotvl3,simtime);
	       else if (plotvl2!=LARGENUM) 
	          rval = ((double (*)(double,double,double))funcc) 
			  (plotvl,plotvl2,simtime);
	       else if (plotvl!=LARGENUM) 
	          rval = ((double (*)(double,double))funcc) (plotvl,simtime);
	       else 
	          rval = ((double (*)(double,double))funcc) (plotnum,simtime);
	     }
	     else rval = 0;
	     break;
   case SREC: 
        if (spnt = (Symbol *)plotnod[plotnum].spnt) {  /* code is from getvar() */
          narg = (long int)spnt;
          if (narg > 0 && narg < 20) {    /* if var is arg to a function */
                getarg(narg,0,&p);   
		rval = *p.val;
          }
          else {
            switch (spnt->type) {
            case UNDEF:
            case CONST:
            case VAR:    rval = spnt->val;
                      break;
            case ARRAY:  
                      rval = *garrayval(spnt); 
                      break;
            default: execerror("attempt to eval non-var", spnt->name);
                      break;
            } /* switch */
          } /* else */
	} else if (var=(double *)plotnod[plotnum].var) { // C++ version
	    rval = *var;
	}
        break;
			/* record from synapse */

case NRECA0: case NRECA1: case NRECA2: case NRECA3: case NRECA4:
case NRECA8: case NRECA9:
case NRECH0: case NRECH1: case NRECH2: case NRECH3: case NRECH4:
case NRECB0: case NRECB1: case NRECB2: case NRECB3: case NRECB4:
case NRECC0: case NRECC1: case NRECC2: case NRECC3: case NRECC4: 
case NRECC9:

			/* record conductance from channel */

case NRECG0: case NRECG1: case NRECG2: case NRECG3:
case NRECG: case NRECGI: case NRECGV: case NRECGC: case NRECGM: case NRECGH: 

			/*  cnod1 = element number, not node number */ 
		rval = record(cnod1,cnod2,cnod3,cnod4,pmod,pval);
        	break;
case CABLE:
     		rval = reccable(cnod1,((double) cnod2)/CABLFRAC);
		break;

   }    /* switch */
 
#ifdef DEBUG
  if (debug & NCRECORD) ncfprintf (stderr,"getrecvalend\n");
#endif

 return (rval);
}

/*------------------------------------*/

double getcurv(synap *spnt, double val)

/* Given a synapse (pointer) and a voltage (volts), return
instantaneous value for amount of neurotransmitter released per
time step (100 usec), in numbers of half-saturating vesicles
(i.e. vsize=1). */

{
   double rval=0;

  switch (spnt->curve) {

   case LINEAR:	if (val < 0) val = 0;		/* linear, sharp shoulder */
		rval = val * 10.0 * spnt->ngain;/* val=0.1V -> rval = 1 */
		break;

   case EXPON:	rval =  0.025 * exp(val*VTOMV/spnt->ngain); /* exponential */
		if (rval < 0.00001) rval = 0.0;	/* val is in volts */ 
		break;
  }
/*  ncfprintf (stderr,"curv %d, ngain %g val %g, rval %g cgain %g vgain %g\n",
			spnt->curve,spnt->ngain,val,rval,spnt->cgain,
			spnt->vgain); /* */
 return (rval * spnt->vgain);
}

/*------------------------------------*/

double rpresyn(synap *spnt, int stime)

/* run a synapse's vesicle release */

{
     int i,nsites;
     double trvolts,fltvolts,vesrel,transrate,p;
     double vsize,refr,meanves,pmean,nvtime,sd;       
     double chanligand, maxsrate, caval;
     int binomdev(double pp, int n);
     double synfilt(lpfilt *lpnt, double transrate);
     cacomp *capnt;

  switch (spnt->sens) {
   case V:
   if (stime) {
     trvolts = (spnt->comp1->v - spnt->comp1->vph - spnt->comp1->vext - spnt->thresh);  /* calib in volts */
     if (!oldsynapz) 
          transrate = getcurv(spnt,trvolts);         /* static tranfer func */
     else transrate = trvolts;

     for (i=0; i<stime; i++) {                   /* one or more steps */
       if (spnt->filt1)
         fltvolts = synfilt(spnt->filt1,transrate); /* first lowp filt */
       else fltvolts = transrate;
       if (spnt->filt1h) {		/* highpass filt */
	   double hfrac, lfrac;

	 hfrac = spnt->filt1h->tfall;	/* highpass gain */
	 if (hfrac < 0) {		/* if hp gain < 0, scale lp gain to keep initial edge same */
	    hfrac = -hfrac;
	    lfrac = 1 - hfrac;
	    if (lfrac < 0) lfrac = 0;
            fltvolts *= lfrac; 
	 }
         fltvolts += synfilt(spnt->filt1h,transrate) * hfrac + spnt->filt1h->offset; 
         if (!oldsynapz) { if (fltvolts < 0) fltvolts = 0; }
       }
       transrate = fltvolts;
    }
   }
   else {
     trvolts = (spnt->comp1->vest - spnt->comp1->vph - spnt->comp1->vext - spnt->thresh); /* calib in volts */
     if (!oldsynapz) transrate = getcurv(spnt,trvolts);         /* static tranfer func */
     else fltvolts = trvolts;
   }
   if (oldsynapz) transrate = getcurv(spnt,fltvolts);         /* static tranfer func */
   break;

   case CA:				/* calcium sensitivity */
    if (capnt=spnt->comp1->capnt) {
	    caval = capnt->cais[0];		 /* take calcium in first shell */
	    transrate = exp(log(caval*dscavg)*spnt->caegain) * spnt->vgain; /* scale by vgain */
	    // fprintf (stderr,"caval %g dscavg %g\n",caval,dscavg);
    }
    else transrate = 0;
    break;
  } /* switch (spnt->sens) */

  if (stime) {

    if ((maxsrate=spnt->maxsrate) > 0) { /* "readily releasible pool" is in use */
		 // transrate *= spnt->rrpool / spnt->maxrrpool;
		 transrate *= spnt->rrpoolg + (1.0 - spnt->rrpoolg) * spnt->rrpool / spnt->maxrrpool;
    }
    for (i=0; i<stime; i++) {                   /* one or more steps */

      if (nsites=spnt->vsites && spnt->vsize) {  /* if noisy release */
        vsize=spnt->vsize * 10.0 * HMICSEC / stiminc;   /* ves rate*timeinc */

		/* for binomial, stdev = sqrt(np(1-p)) */

        meanves = transrate / vsize; 
        if (maxsrate > 0) { 			/* "readily releasible pool" is in use */
           if (meanves > spnt->rrpool) meanves = spnt->rrpool;
	}
	sd = meanves * spnt->cov; 
	refr = spnt->vrefr;			/* user set refractory period */
	if (refr < 0.0) refr = 0.0;
	if ((refr > 1) && (meanves > 1)) {
	    meanves = 1; 
	}
        if (meanves > 1) {
          if (spnt->cov > 0.3) {
            p = 1e-2;
	  }
	  else {
	    p = 1.0 - (sd * sd / meanves);
	  }
	  nsites = (int) (meanves / p);
	  if (nsites < 1) nsites = 1;
	  p = meanves / nsites;
	  if (p < 0) p = 0;
	  else if (p > 1.0) p = 1.0;
	  if (spnt->vstate) drand_setstate (spnt->vstate);
          vesrel = binomdev(p,nsites);    
	}
	else {		/* Less than 1 ves/time step, */
			/*  use refractory period: */
	   double gorder;

	/* Mark van Rossum's way for calculationg refr period: */

/*	   refr = 1.0 / meanves - sd / (meanves*sqrt(meanves)); */

	   if (spnt->cov==0) gorder = 1e6;
	   else gorder = 1.0 / (spnt->cov*spnt->cov);

	   if (meanves<MINVES) {	/* If very low rate, */
	      meanves = MINVES;		/*  just reset */
	      pmean = 1.0 / meanves;
	      if (spnt->vtime > 0) spnt->vtime = 0;
	      if (spnt->vstate) drand_setstate (spnt->vstate);
	      spnt->vtint = gamdev(gorder)/gorder;
	      spnt->vflag = 0;
	      vesrel = 0; 
	   }
	   else {	/* meanves is non-zero */
		double x;
	     x = 1.0 / meanves;			/* x ~= pmean = time until release */
	     if (spnt->vflag) x -= refr;
	     if (x >= 0) { 	     /* asymptote to x for pos pmean */
		  if (x > 100) pmean = x;
	          else         pmean = exp(-x/2) + x;
	     }
	     else {		     /* asymptote to 0 for neg pmean*/ 
		if (x < -100) pmean=0; 
	        else          pmean = exp(x/2);
	     }
	     if (pmean > MAXVINT) pmean = MAXVINT; 	

/* ncfprintf (stderr,"vtime %g vtint %g pmean %g meanves %g gorder %g refr %g\n",
  spnt->vtime, spnt->vtint, pmean, meanves,gorder, refr);  /* */

			/* if rate increases, reset interval */

						/* vesicle event */
	   if (spnt->vtime >= spnt->vtint) {	
	         if (spnt->vstate) drand_setstate (spnt->vstate);
	         spnt->vtint = gamdev(gorder)/gorder;	/* vtint = ves. thresh interval for reset*/
	         spnt->vtime = -refr / pmean;		/* vtime = ves. release time */
	         spnt->vflag = 1;			/* ves release flag, after refr done */
	         vesrel = 1;
	   }
	   else {		/* no event, just count interval */
	      vesrel = 0;
	    }

	  } /* else meanves is non-zero */

	  nvtime = spnt->vtime * pmean;
	  nvtime++;
      	  if (nvtime > MAXVINT) {
	      spnt->vtime = -refr / pmean;	
	  }
	  spnt->vtime = nvtime / pmean;

/* ncfprintf (stderr, "vtime %g vtint %g meanves %g pmean %g vesrel %g\n",
		spnt->vtime,spnt->vtint,meanves,pmean,vesrel); /* */

	}  /* less than 1 ves/step */

 /* ncfprintf (stderr, "time %g vtime %g vtint %g meanves %g vesrel %g\n",
		time,spnt->vtime,spnt->vtint,meanves,vesrel); /* */

        if (spnt->vcov > 0) {
	     double r,voldev;

	   if (spnt->gstate) drand_setstate (spnt->gstate);
	   r = gasdev() * spnt->vcov + 1.0;
	   restorstate(); 
	   if (r < 0) r = 0;
	   voldev = r * r * r; 
	   transrate = vesrel * vsize * voldev; 
	}
	else {				/* no ampl. variation */
	   transrate = vesrel * vsize; 
	}
			/* The orig Poisson way: */

//	/*  if (spnt->vstate) drand_setstate (spnt->vstate);
//	    vesrel = poisdev (transrate/vsize); 
//	    transrate = vesrel * vsize; 
//	*/
//
//			/* the old binom way: */
//
//  /*  #define TMAXP 5.0     			     /* max prob at 10.0 */ 
//    /*		   transrate /= TMAXP;
//		   if ((p=transrate/(vsize*nsites))<=1) {
//		     if (vsize<1) {		   
//		       p = transrate/nsites;	    
//		       nsites = (int)(nsites / vsize);
//		     }
//		     nsites = (int)(nsites * TMAXP);
//		     if (spnt->vstate) drand_setstate (spnt->vstate);
//	             vesrel = binomdev(p,nsites);    
//		     transrate = vesrel * vsize; 
//		   }
//		   else {
//		     transrate *= TMAXP;
//		   }
//	*/

/* ncfprintf(stderr,"vesrel %g transrate %g p %g size %g nsites %d\n",
			vesrel,transrate,p,vsize,nsites); 	/* */
	 } /* if (vsites && vsize) */

	 else {			/* no noise in release */
           vsize=  dvsz * 10.0 * HMICSEC/stiminc;   /* ves rate*timeinc */ // Added stiminc
           meanves = transrate / vsize; 
	   vesrel = meanves;
	   if (vesrel < 0.0001) vesrel = 0;
	}

/* ncfprintf(stderr,"transrate %g vesrel %g\n", transrate,vesrel); 	/* */

	if (maxsrate > 0) { /* Set "readily releasible pool" */
	      double rrpool;

	  rrpool = spnt->rrpool;
          if (nsites=spnt->vsites && spnt->vsize) {  /* if noisy release */
	     rrpool += poisdev(maxsrate * stiminc) * (spnt->maxrrpool - rrpool) / spnt->maxrrpool;
	  }
	  else 
	     rrpool += maxsrate * stiminc * (spnt->maxrrpool - rrpool) / spnt->maxrrpool;

// fprintf (stderr,"rrpool %g maxsrate %g vesrel %g\n",rrpool,maxsrate,vesrel);

	  if (vesrel>0) {
	     rrpool -= vesrel;
	     if (rrpool < 0) rrpool = 0;
	  }
	  if (rrpool > spnt->maxrrpool) rrpool = spnt->maxrrpool;
	  spnt->rrpool = rrpool;
	} // if (maxsrate>0)

     } /* for (i=0;;) */
     restorstate();   /* restore random state from binomdev and gamdev */
  }
  spnt->transrate = transrate;
  return transrate;
}

/*------------------------------------*/

void runsyn(int stime)
  			/* number of 100 microsec steps */

/* Run all the synapses, depending on their
   type and history.  Normally, call once every 100 microsec.
   For static runs, call with stime=0 and only
   the amplitude transfer function is produced, not the
   time functions.
*/

{
   int i, npre;
   double chanligand, chanligand2, transrate;
   double synfilt(lpfilt *lpnt, double transrate);
   synap *spnt;
   dyadlst *lpnt;
   gj *gjpnt;
   double vj, cyc, normcond, gnv;
   double alpha, beta, n, rpar;
   double v;
   dbuf *dlpnt;
   comp *cpnt;

#ifdef DEBUG
  if (debug & NCBASE) ncfprintf (stderr,"runsynap");
#endif

 if (synaptau > 0) stime /= synaptau;		/* run synapse faster when synaptau is small */

 for (spnt=synpnt; spnt; spnt= (synap *)spnt->next)
   {
     switch (spnt->ctype)
      {
       case SYNAPSE:
         if (!spnt->spre) {
	    transrate = rpresyn(spnt,stime) * spnt->trconc;
            npre = 1;				/* number of pre-syns summed */
	 }
	 else {
	    transrate = 0; 
	    npre = 0;
         }
	 if (spnt->spre) {
	   for (lpnt=spnt->spre; lpnt; lpnt=lpnt->next) {
             if (lpnt->sdyad) {
               transrate += lpnt->sdyad->transrate; /* trans from other presyn */
	       npre++;
	     }
           }
	   transrate /= npre;
	 }

	 if (stime) {
           for (i=0; i<stime; i++) {            /* one or more steps */
	    if (spnt->filt2)
	      chanligand = synfilt(spnt->filt2,transrate); /* 2nd lowp filt*/
	    else chanligand = transrate;
           }
	 }
	 else chanligand = transrate;

         rchan(spnt,stime,chanligand);		/* run the channels */

         break;   /* case SYNAP */

       case GJ:                                 /* GJ modulated by voltage */
       case PNX:                                 /* PNX modulated by voltage */
	gjpnt = (gj *)spnt;
	gnv = gjpnt->gnv;
	if (stime) {
         for (i=0; i<stime; i++) {		/* one or more steps */

	cyc = getmesg(gjpnt,gjpnt->modtyp);

	if (gjpnt->sign) cyc = 1.0 - cyc;

			/* For details, see Harris, Spray, Bennett, */
			/* J.  Gen Physiol. 77:95-117, 1981.        */

#define GJLAMBDA 0.0013
#define Aa 0.077
#define Ab 0.14

          vj = gjpnt->comp1->v - gjpnt->comp2->v;
	  alpha = gjalpha(vj, gjpnt);
	  beta =  gjbeta (vj, gjpnt);
	  n = gjpnt->n;
	  n += -beta * n + alpha * (1.0 - n);
	  if (n < 0.0) n = 0.0;
	  if (n > 1.0) n = 1.0;
	  gjpnt->n = n;
	  if (gjpnt->rect != 0) {
	     if ((gjpnt->comp1->v - gjpnt->comp2->v) > 0)
		  gjpnt->conduct = 1;
	     else
		  gjpnt->conduct = 0;
	     if (gjpnt->rect < 0) gjpnt->conduct = !gjpnt->conduct;
	     gjpnt->conduct *= gjpnt->maxcond;
	  } 
	  else gjpnt->conduct = gjpnt->maxcond * ( gnv + n * (1.0-gnv)) * (1.0-cyc);

 	  }   /* for (i;;) */ 
	 }  /* if (stime;;) */

	 else { 	/* stime==0,  static gj, update compartment implf */

	   cyc = getmesg(gjpnt,gjpnt->modtyp);
	   if (gjpnt->sign) cyc = 1.0 - cyc;

           vj = ncabs(gjpnt->comp1->v - gjpnt->comp2->v);
	   alpha = gjalpha(vj, gjpnt);
	   beta =  gjbeta (vj, gjpnt);
	   n = alpha / (alpha + beta);
	   gjpnt->n = n;
	   if (gjpnt->rect != 0) {
	     if ((gjpnt->comp1->v - gjpnt->comp2->v) > 0) normcond = 1;
	     else					  normcond = 0;
	     if (gjpnt->rect < 0) normcond = !normcond;
	   } 
	   else normcond = (gnv + n * (1.0-gnv)) * (1.0-cyc);

	   cpnt = gjpnt->comp1;		/* check first comp */
	   if (cpnt) {
	     if (cpnt->implf) {	/* if implf has been calculated yet... */
	      rpar = ((1. / cpnt->implf) - 1.) / cpnt->k*implfk - gjpnt->conduct;
	      rpar += normcond * gjpnt->maxcond; 
	      cpnt->implf = 1. / (1. + rpar * cpnt->k*implfk);  /* implicit factor*/
	    }
	  }

	   cpnt = gjpnt->comp2;		/* check second comp */
	   if (cpnt) {
	     if (cpnt->implf) {	/* if implf has been calculated yet... */
	      rpar = ((1. / cpnt->implf) - 1.) / cpnt->k*implfk - gjpnt->conduct;
	      rpar += normcond * gjpnt->maxcond; 
	      cpnt->implf = 1. / (1. + rpar * cpnt->k*implfk);  /* implicit factor*/
	    }
	  }
          gjpnt->conduct = gjpnt->maxcond * normcond;
	 } /* else static */
	 
	 if (spnt->ctype==PNX) {		/* set ATP flux to outside */
		 double pnxcurr, atp_level, ph;
		 pnx *pnxpnt;

	     pnxpnt = (pnx*)gjpnt;
	     if (pnxpnt->pnxfunc==NULL) {
                vj = pnxpnt->comp2->v - pnxpnt->comp1->v;
	        pnxcurr = pnxpnt->conduct * vj * 1e14;   
	        atp_level = addnt(pnxpnt->comp2,ATP,pnxcurr);
	        atp_level *= pnxpnt->atp_decr;
	        setnt(pnxpnt->comp2,ATP,atp_level);
	     } else {
		ph = pnxpnt->pnxfunc(pnxpnt,0);     /* run user-specified function */
		// setnt(pnxpnt->comp2,PH,ph);	    // now done in user func, phbuf_inc.cc
		// atp_level = pnxpnt->atp;
	     }
	 }
         break;

       case BUF:                                /* delay buffer */
	 dlpnt = (dbuf *)spnt;
         for (i=0; i<stime; i++) {                   /* one or more steps */
	   v = (dlpnt->comp1->v - dlpnt->offset) * dlpnt->gain;
	   if (dlpnt->filt) v = synfilt(dlpnt->filt,v);
	   if (((dbuf *)spnt)->delay > 0) {	/* if delay */
	      if (dlpnt->delbuf) {
	        *dlpnt->delpnt = v;
	        if (++dlpnt->delpnt >= dlpnt->delbuf+ (int)(dlpnt->delay))
		   dlpnt->delpnt = dlpnt->delbuf;  /* incr. pntr, make circ. */
					/* leave pointing to next output val */
	      }  /* if (dlpnt->delbuf) */
              v = *dlpnt->delpnt;  /* output value */
	   }  /* if (dlpnt->delay>0) */
         }  /* for (stime;;) */
	 dlpnt->v = v;
         break;  /* case BUF */

      } /* switch */
   }        /* for */
#ifdef DEBUG
  if (debug & NCBASE) ncfprintf (stderr,"...end\n");
#endif
}

/*------------------------------------*/

double hillcalc (double ligand, double kd, double p)

/* Calculate receptor binding and normalized conductance from
   kd and Hill coefficient.
*/
 
{
    double kk, ll;

  kk = pow(kd,p);
  ll = pow(ligand,p);
  return (ll / (ll + kk));  	/* fraction bound */
}

/*------------------------------------*/

double hillcalcn (double ligand, double kd, double p)

/* Calculate second-messenger receptor binding and normalized conductance 
   from kd and Hill coefficient.  

   Problem is that input range is restricted to 0-1 so standard
   saturation with Hill coeff doesn't give a reasonable amount
   of saturation when kd=1, and when kd=0.1 the standard
   saturation gives a gain of 0.1.

   Therefore, normalize so input of 1 gives output of 1. For best
   results should have kd of .1 (to make saturation). 
*/

 {  
    double kk, ll;

  kk = pow(kd,p);
  ll = pow(ligand,p);
  return (ll / (ll + kk) * (1 + kk));  	/* fraction bound */
}

/*------------------------------------*/

double setbind(synap *spnt, double chligand)

/* calc receptor binding and norm. conductance for standard synapse */

{
   double bound;

 bound = hillcalc (chligand,spnt->nkd,spnt->npow);
 return bound;
}

/*------------------------------------*/

double setcond(synap *spnt, double bound)

/* Calc conductance for standard synapse */
/* If synaptic 2nd mesg with cGMP channel,
    return [cGMP] as the conductance. 
   Otherwise, if 2nd mesg with "generic" synapse,
   assume instantaneous binding, and return bound cGMP
   as the conductance.
*/

{
   double active;

   switch (spnt->ntact) {		/* calc active signal */
     case OPEN: 
	  active = bound;
	  break;

	  /* Here we simulate the action of the cGMP cascade in an
     	     on-bipolar.  See Shiells and Falk, 1994. */

     case CLOSE:
         active = spnt->coff - bound * spnt->cgain;
         if (active < 0.0) active = 0.0;
	 break;
   }
   return active;
}

/*------------------------------------*/

double setbind2 (synap *spnt, double active)

	/* Here we simulate binding and saturation of second messenger */
{
   double active2;
   ntchan *resp;
   
  if (spnt->secmsg) {		/* if second messenger, calc binding */
				/* otherwise pass linear signal */ 

     if (resp=spnt->resp2) {	/* if markov channel */
	if (resp->chtyp->cgbnd)	/* if artif cCMP binding is needed */
 	   active2 = hillcalcn (active,spnt->ckd,spnt->chc);
        else active2 = active;
     }
     else {
	if (spnt->maxcond > 0)
	  active2 = hillcalcn (active,spnt->ckd,spnt->chc);
        else			/* no channel, must use mesg out */
	  active2 = active; 
     }
  }
  else active2 = active;
  return active2;
}

/*------------------------------------*/

void rchan(synap *spnt, int stime, double chligand)

/* Run the channel conductance part of a synapse.  */

/*  AMPA, KAINATE, CGMP, GABA, GLY, SYN2, and NMDA channels are run as state machines 
    instead of having Hill coefficient binding and temporal filters. 
    Their conductances are found and assigned to the synaptic conductance 
    just like the standard synapse.  These channels are run during the 
    compartment iteration so are not part of this procedure.
*/

{
   int i,nchan;
   double bound,active,active2,cycg,cdur,rpar;
   double synfilt(lpfilt *lpnt, double transrate);
   double a,x,k,kk,ll,pp,cycgmax;
   comp *cpnt;


  if (spnt->resp1) {			/* if postsynaptic markov channel, */
    switch (spnt->resp1->ctype){	/* calc "active" signal sent to chan */

	case AMPA:			/* ligand-binding channels */ 
	case KAINATE: 
	case NMDA: 
	case GLY: 
	case GABA: active2 = active = chligand; 
		   break;
    }
  }
				/* Since we don't have a Markov version of mGluR6, */
				/*  we approximate using hill eqn and use kd = trconc  */

  if (spnt->resp2) {			/* if 2nd messenger markov channel */
    switch (spnt->resp2->ctype){	/* calc "active" signal sent to chan */
	case CGMP:			/* "second-messenger" chans */
	case SYN2:			/* Need "generic" post-syn receptor */
	  if (!spnt->resp1) {
   	    bound = setbind (spnt,chligand/dstr);  /*  calc Glu binding, normalized */
   	    if (stime && spnt->filt3) bound = synfilt(spnt->filt3,bound);
   	    active  = setcond(spnt,bound);     /* calc norm cond, 2nd mesg  */
   	    active2 = setbind2(spnt,active);   /* calc 2nd mesg binding */
          }
	  break;
    }
  }
  if (!(spnt->resp1 || spnt->resp2)) {	/* else standard postsynaptic mechanism or mesg*/

   bound = setbind (spnt,chligand/dstr);    /*  calc postsyn rec binding, normalized */
   if (stime && spnt->filt3) bound = synfilt(spnt->filt3,bound);
   active = setcond(spnt,bound);     /* calc norm cond, possib w/2nd mesg */
   active2 = setbind2(spnt,active);   /* calc 2nd mesg binding */

   if (!stime) { /* stime==0,  static syn, update compartment implf */
     cpnt = spnt->comp2;
     if (cpnt) {
       if (cpnt->implf) {	/* if implf has been calculated yet... */
         rpar = ((1. / cpnt->implf) - 1.) / cpnt->k*implfk - spnt->conduct;
         rpar += active2 * spnt->maxcond; 
         cpnt->implf = 1. / (1. + rpar * cpnt->k*implfk);  /* implicit factor*/
       }
     }
     else ncfprintf (stderr,"rchan: panic, no second compartment\n");
   }
  } 

  if (spnt->mesg1){	/* neurotransmitter output signal */
        double a, x;
     a = (active - spnt->oldc1);
     x = addnt(spnt->comp2,spnt->mesg1, a);
     spnt->oldc1 = active;
  }
  if (spnt->mesg2){	/* second messenger output signal */
			/* mesgconc is set from chantype->trconc */
			/* which can be different than spnt->trconc */
     a = (active - spnt->oldc2)*spnt->mesgconc;
     x = addnt(spnt->comp2,spnt->mesg2, a);
     spnt->oldc2 = active;
  }
  if (spnt->resp1)
     spnt->conduct = active2; /* save norm. conductance of chan */
   else
     spnt->conduct = active2 * spnt->maxcond;

/* ncfprintf (stderr,"maxcond %g active2 %g chlig %g bound %g nkd %g ntact %d\n",
			spnt->maxcond,active2,chligand,bound,spnt->nkd,spnt->ntact);  /* */
}

/*------------------------------------*/

void runrec(photrec *rpnt, int dtimestep) 

/* run a photoreceptor for one or more timesteps */

{
   int i,tspeed;
   double mask;
   recpar *cpt;
   photreci *rpnti;
   comp *cpnt;
   static double flux,iflux,pflux,rflux,cag,cak,dcab,cond;
   phvars *vpnt;
  
    flux = 0; 
    for (i=0; i<numstimchan; i++) {  /* weight masked flux by chan value */
       if (rpnt->stimchan>0 && rpnt->stimchan != i) continue;  /* skip if wrong channel */
       mask = rpnt->mask[i];
       if (mask > 1) mask = 1;
       if (mask < 0) mask = 0;
       if (rpnt->mflux[i] <= unmaskthr) mask = 0;	    /* intensity masking */
       flux += rpnt->chanw[i] * (rpnt->aflux[i]*(1-mask) + rpnt->mflux[i]*mask); /* masking flux */
    }

    vpnt = &rpnt->vars;				/* pointer to state variables */
    switch (rpnt->ctype) {

      case VTRANSDUCER:
	 rpnt->iflux = flux; 			/* save iflux for record */
 	 if (cpnt=rpnt->comp1) {
	    if (stimonh > stimonl) {		// turn on vclamp within stimonl, stimonh window

 	      if (flux > stimonl && flux < stimonh) setvclamp(cpnt,flux,0);
	      else                                  setvclamp(cpnt,flux,1);
	    }
	    else {				// turn on vclamp outside stimonl, stimonh window

 	      if (flux > stimonl || flux < stimonh) setvclamp(cpnt,flux,0);
	      else                                  setvclamp(cpnt,flux,1);
	    }
	 }
	 break;
      case ITRANSDUCER:				/* continuous cclamp, reset each time step */
	 rpnt->iflux = flux; 			/* save iflux for record */
 	 if (cpnt=rpnt->comp1) {
              cpnt->exti -= cpnt->exti_old;	/* stop previous current clamp */
              cpnt->exti += flux;		/* set up current clamp */
              cpnt->exti_old = flux;		/* save in previous current clamp */
	 }
	 break;

      case ROD:
      case CONE:
      case CHR:
         if (!(cpt=rpnt->chtyp)) {		/* get pointer to const */
	    break;
	 }
	 else {
	 if (tempcel != cpt->oldtempcel) initrec(); /* recalc if temp changes*/
	 tspeed = (int) (rpnt->tspeed * dtimestep);
	 rflux = flux;
				/* calib for invergo rod, DellOrco Rect */

	 if (cpt->pigm==RODI) rflux *= ((recpari*)cpt)->gain / 
			    	       ((recpari*)cpt)->gain_spr_mul; 

         if (rpnt->pnois) {
              if (rpnt->pstate) drand_setstate (rpnt->pstate);
	      pflux = poisdev(rflux);
	      restorstate();
	      iflux = pflux;
	 }
	 else iflux = rflux; 
	 rpnt->iflux = iflux;  		/* save iflux for record */

	 switch (cpt->pigm) {		/* transduction without cascade */
		double isens;
	  case 7:
	  case 8:
	  case 9:
	    isens = vpnt->ksens * stiminc;
            vpnt->cond = isens / (iflux + isens);
	    cond = vpnt->cond;
 	    break;

	 default:
	  iflux *= rpnt->intenadj;		/* adjust inten for timec */
          for (i=0; i<tspeed; i++) {             /* maybe do several steps */
		double t1, t2, t3;

            vpnt->rhod  += iflux * cpt->lgain * (1.0 - vpnt->resp_b);         /* */

	    					/* bleaching from van hateren cone */
	    vpnt->resp_b  += cpt->cn*vpnt->rhod - cpt->rk_b/cpt->tau_b * 
					vpnt->resp_b/(vpnt->resp_b + cpt->rk_b);

	    t1 = vpnt->rhod  * cpt->rgain1;
	    vpnt->rhod  -= t1;
            vpnt->rhod2 += t1;

            t2 = vpnt->rhod2 * cpt->rgain2;
            vpnt->rhod2 -= t2;
            vpnt->rstar += t2;

            vpnt->gpr   +=  vpnt->rstar * cpt->rgain3;

            t3 = vpnt->gpr * cpt->ggain;
	    vpnt->gpr -= t3; 
	    vpnt->pde += t3;
	    if (rpnt->dnois>0) {
                 if (rpnt->dstate) drand_setstate (rpnt->dstate);
		 vpnt->pde  += (drand() - 0.5) * .0004856 * rpnt->dnois;
 		 restorstate();  	/* restore default random state */
	    }

	      /* Calcium feedback */

            // cag = pow(rpnt->ca,cpt->powcag);
            // cak = pow(rpnt->ca,cpt->powcak);
            cag = pow(vpnt->ca,cpt->powcag);
            cak = pow(vpnt->cax,cpt->powcak);	/* rhodk has extra delay */

            vpnt->gcyc  = 1.0 - (cag / (cag + cpt->kdgcyc));
	    vpnt->rhodk = 1.0 - (cak / (cak + cpt->kdrhodk));

	    // fprintf (stderr,"cak %g rhodk %g\n",cak,vpnt->rhodk);

            vpnt->cycg -= vpnt->pde * vpnt->cycg * cpt->pdegain;
            vpnt->cycg += vpnt->gcyc * cpt->gcygain;
            if (vpnt->cycg < 0) vpnt->cycg = 0;

            vpnt->cond += pow(vpnt->cycg,cpt->powcycg) * 
			 (cpt->totcond - vpnt->cond) / cpt->totcond * 
			  cpt->condgain;

	    // This implements Russ Hamer's (2005, VN 22:417) calcium buffer. 
	    // We had been using cax instead.
	    // The effect is similar, to delay the effect of calcium.
	    //
	    dcab = vpnt->ca * cpt->gcab * (cpt->cat - vpnt->cab) - vpnt->cab * cpt->gcabr;
            vpnt->cab  += dcab; 
            vpnt->ca   += vpnt->cond * cpt->gca - dcab; 
            // vpnt->ca   += vpnt->cond * cpt->gca; 
            vpnt->cax  += vpnt->ca * cpt->cxgain;

            vpnt->rstar *= cpt->decrstar - pow(vpnt->rhodk,1) * cpt->rkgain;	

            vpnt->pde  = cpt->pdedec * (vpnt->pde- cpt->pdebase) + cpt->pdebase;
            vpnt->ca   -= cpt->capump * vpnt->ca / (vpnt->ca + 1.0);
            vpnt->cond *= cpt->deccond;
            vpnt->cax  *= cpt->deccax;
 	  }	/* for */
	  cond = vpnt->cond;
          break;

	  case 19:
	  case 20:
	  case 21:		/* van hateren cone */

				/* van Hateren, JH, Snippe HP. (2007). */
				/*  Simulating human cones from mid-mesopic up to high-photopic luminances.*/
	 			/*  Journal of Vision, 7(4):1, 1-11. */
				/* http://journalofvision.org/7/4/1/, doi:10.1167/7.4.1. */
						 
	     double alpha, beta, beta_e;
	     double resp_im;

	    iflux *= rpnt->intenadj;		/* adjust inten for timec */

            for (i=0; i<tspeed; i++) {             /* maybe do several steps */

		vpnt->resp_r  += (iflux * cpt->lgain* (1-vpnt->resp_b-cpt->cn*vpnt->resp_r) - 
					vpnt->resp_r) / cpt->tau_r;
		vpnt->resp_b  += cpt->cn*vpnt->resp_r/cpt->tau_r - cpt->rk_b/cpt->tau_b * 
					vpnt->resp_b/(vpnt->resp_b + cpt->rk_b);
		vpnt->resp_e  += (vpnt->resp_r-vpnt->resp_e)/cpt->tau_e;

		beta     = cpt->c_beta + cpt->rk_beta*vpnt->resp_e;
		beta_e   = beta / (1 + beta/cpt->beta_e_max);
		alpha    = 1 / (1 + pow(cpt->a_c*vpnt->resp_c,cpt->rnc));

		vpnt->resp_x  += alpha - beta_e*vpnt->resp_x;
		vpnt->resp_os  = pow(vpnt->resp_x,cpt->rnx);

		vpnt->resp_c  += (vpnt->resp_os - vpnt->resp_c)/cpt->tau_c;

		vpnt->resp_is += (vpnt->resp_os/vpnt->atten_i - vpnt->resp_is)/cpt->tau_vc;
		resp_im  = pow(vpnt->resp_is,cpt->gamma_is);

		vpnt->atten_i += (cpt->a_is*resp_im - vpnt->atten_i)/cpt->tau_is;
		 
		vpnt->cond = vpnt->resp_is / cpt->dark_resp_is;

		// vpnt->cond = (vpnt->resp_os - cpt->dark_resp_os) / cpt->dark_resp_os;
		//vpnt->cond = vpnt->resp_os / cpt->dark_resp_os;
 	     }	/* for */
	     cond = vpnt->cond;
          break;

	  case RODI:				/* invergo et al (2014) mouse rod */
	      double k, d;			/*  reactions described in supplemental data */

	   rpnti = (photreci *)rpnt;
	   iflux *= rpnti->intenadj;		/* adjust inten for timec */
           for (i=0; i<tspeed; i++) {             /* maybe do several steps */
	        rpnti->vars.cond = runrodi(rpnti, iflux);	/* in rodi.cc */
	      //  rpnti->vars.cond = runrodi_orig(rpnti, iflux);	/* in rodi.cc */
	   }
	   cond = rpnti->vars.cond;

          break;

	  case ChR2:				/* channel-rhodopsin, Williams et al (2013) */
	       double vg,g;
	       photrecc *rpntc;
	   rpntc = (photrecc *)rpnt;
	   iflux *= rpntc->intenadj;		/* adjust inten for timec */
	   if (rpntc->phchan!=NULL) {		/* use channel defined in chanchr2.cc */
	      ((chrchan*)(rpntc->phchan))->iflux = iflux;  /* set inten for markov CHRC chan */
	      cond = ((chrchan*)(rpntc->phchan))->conduct;  /* get cond for markov CHRC chan */

	   } else {				/* use ChR2 in chr2.cc */
              for (i=0; i<tspeed; i++) {             /* maybe do several steps */
	        rpntc->vars.cond = runchr(rpntc, iflux, rpntc->comp1->v);	/* in chr2.cc */
	      }
	      cond = rpntc->vars.cond;
	   }
	   vg = rpnt->comp1->v * VTOMV;
	   if (vg > -10) vg = -10;
	   cond *= (10.6408 - 14.6408 * exp(-vg/42.7671)) / vg;
	// ncfprintf (stderr,"cond %g\n", cond);                /* */
	   break;

         }    /* switch (cpt->pigm) */
        }    /* else */
// fprintf (stderr,"rpnt->cond %g resp %g lgain %g\n",rpnt->cond, rpnt->resp_os, cpt->lgain);

         rpnt->conduct = rpnt->maxcond * cond;   /* chan current */

// ncfprintf (stderr,
//        "tspeed %d iflux %g cond %g conduct %g 2\n", tspeed, iflux, cond, rpnt->conduct);                /* */

if (rpnt->ctype==ROD || rpnt->ctype==CONE) {

  if (debug & NCPHOT && debugz & 16)
    ncfprintf (stderr,
        "rec %d num %d tspeed %d rpntinten %g chanw %g, flux %g rflux %g pflux %g iflux %g aflux %g mflux %g pnois %d rh %g pd %g gc %g cg %g c %g cond %g\n",
                rpnt->recnm1, rpnt->recnm2, tspeed, rpnt->intenadj, rpnt->chanw[0], flux, rflux, pflux, 
		iflux, rpnt->aflux[0], rpnt->mflux[0], rpnt->pnois, vpnt->rhod, 
		vpnt->pde, vpnt->gcyc,
                 vpnt->cycg, vpnt->cond, rpnt->conduct);                /* */
} 

         break;
    }       /* switch */
}

/*------------------------------------*/

void runrecs(int dtimestep)

/* Run all the photrecs. 
 depending on their type and history.
 Call at multiples of tenths of milliseconds only. */

{
   photrec *rpnt;

#ifdef DEBUG
  if (debug & NCPHOT)  ncfprintf (stderr,"runrecs\n");
#endif

  for (rpnt = recpnt; rpnt; rpnt=(photrec *)rpnt->next) {
     runrec(rpnt,dtimestep);
  }     /* for (rpnt; ; )  */

#ifdef DEBUG
  if (debug & NCPHOT) ncfprintf (stderr,"...end\n");
#endif
}

/*------------------------------------*/

void setvclamp (comp *cpnt, double val, int bgend)

/* set a voltage clamp for a compartment */

{
	int nvclamps;

 if (cpnt) {
  if (!(cpnt->miscfl & VBUF)) { /* not buf already */
    if (bgend==0) {  		/* turn on vclamp */
          cpnt->extv = val;     /* set voltage clamp */
          cpnt->v    = val;     /* set comp voltage (for plot)*/
	  nvclamps = cpnt->miscfl & VEXT;
	  nvclamps++;		/* incr. num of vclamps */
	  if (nvclamps > VEXT) nvclamps=VEXT;
	  cpnt->miscfl &= ~VEXT;
	  cpnt->miscfl |= nvclamps;
     } else {  			/* turn off vclamp */
          // cpnt->extv = val;   
	  nvclamps = cpnt->miscfl & VEXT;
	  if (nvclamps>0) nvclamps--;	/* incr. num of vclamps */
	  cpnt->miscfl &= ~VEXT;
	  cpnt->miscfl |= nvclamps;
     }
    }  
  }  /* if (cpnt) */ 
  else {
      ncfprintf (stderr,"runstim, setvclamp: can't find compartment for node.\n");
  }
}

/*------------------------------------*/

void runstim(double rtime, int dmsec, int bgend, int endfl)

/* Advance all the stimuli one timestep.
 Run all light stimuli, voltage and current clamps.
 Call at multiples of tenths of milliseconds only. */
/* bgend = 1 => run ending of delta stimulus only */

{
#define RECTIM 1			/* how many 0.1 msec intervals */
#define READRECTIM 10			/* how often to read stimuli */
// #define RTIMERND  0.00001		/* down-round value */
// #define RTIMESTEP 0.0001		/* time step for reading stimuli, now srtimestep */
// #define RTIMERND  0.0000001		/* down-round value */
// #define RTIMESTEP 0.000001		/* time step for reading stimuli */

   int i,c,nt;
   recstim *rspnt,*next,*rbinpnt;
   photrec *rpnt;
   node *npnt;
   comp *cpnt;
   static int readmsec=READRECTIM;
   static double stop,sens,rtimestep;

#ifdef DEBUG
  if (debug & NCSTIM && debugz & 1)  ncfprintf (stderr,"runstim time %g\n",rtime);
#endif

// rtimestep = max(timinc,RTIMESTEP);
 rtimestep = max(timinc,srtimestep);

/* The input from the file is sorted by time before being read in. 
   This is done as a final step by "stim".
*/

 if(readmsec++ >= READRECTIM) {
   if (!readstim(rtime+readmsec*rtimestep)) {  /* read new photrec inputs */
         ncfprintf (stderr,"runstim: invalid stimulus at time %g\n",rtime);
	 execerror ("Missing stimulus: "," stopping... ");
   }   
   readmsec = 0;                
 }

 if(dmsec >= RECTIM) {        /* calculate photrecs for 0.1 msec intervals: */

	/* Problem here is that "rtime" comes from a continuing
	 * sum which has some roundup error, causing it to
	 * gradually increase to a greater value than the correct
	 * stimulus run time.  Therefore, we offset "rtime" and
	 * the "stop" time by a small amount so a stimulus won't be missed. 
	 */
           
  // rtime -= RTIMERND;
  rtime -= rtimestep * 0.1;
  stop  = rtime + rtimestep;

  for (i=0; i<2; i++) {
    if (i==0) rbinpnt = getstimp1(rtime);
    else      rbinpnt = getstimp2(rtime);
    if (endfl) {
              rbinpnt = recspnt;
	      stop = rtime;
	      rtime = simtime;
    }

    for (rspnt=rbinpnt; rspnt; ) {

  /*   ncfprintf(stderr,
 "stim num %d %d type '%c' val %10g wavel %5g time %15.15g rtime %15.15g stop %15.15g\n",
 rspnt->recnm1,rspnt->recnm2,rspnt->ctype,rspnt->val,rspnt->wavel,rspnt->time,rtime,stop);/**/
 
   if (debug & NCSTIM)  ncfprintf (stderr,"rspnt->time %-10.5g start %g stop %g val %6.3g c '%c' mask %g bgend %d %d\n",
			  rspnt->time,rtime, stop, rspnt->val,rspnt->ctype,rspnt->mask,bgend,rspnt->bgend);

   if ((rspnt->time >= rtime) && (rspnt->time <stop) && (bgend==0 || bgend==rspnt->bgend)) {
     c = rspnt->stimchan;
     switch (rspnt->ctype) {

       case 'a':                /* absolute additive light intensity */
          if (findphotrec (rspnt->recnm1, rspnt->recnm2, 
			 rspnt->recnm3, rspnt->recnm4,
					&rpnt, "runstim")) {
	   if (rpnt->ctype==VTRANSDUCER || rpnt->ctype==ITRANSDUCER) {
	      rpnt->aflux[c] = rspnt->val;
	   }
	   else {	/* rod or cone */
              sens = rsens(rpnt,rspnt->wavel,rpnt->filt);
              rpnt->aflux[c] = rpnt->area * rspnt->val * sens * rpnt->attf * dqeff * stiminc;
              if (rpnt->aflux[c] > MAXPHOT) rpnt->aflux[c] = MAXPHOT;
              if (rpnt->aflux[c] < 1e-10)   rpnt->aflux[c] = 0.0;
	   }
	   rpnt->chanw[c] = 1.0;
	 }
#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"stim abs %d %g %g aflux %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel,rpnt->aflux[c]);
#endif
         break;

       case 'b':                /* delta masking light intensity */
          if (findphotrec (rspnt->recnm1, rspnt->recnm2, 
			   rspnt->recnm3, rspnt->recnm4,
					&rpnt, "runstim")) {
	    if (rpnt->ctype==VTRANSDUCER || rpnt->ctype==ITRANSDUCER) {
	       rpnt->mflux[c] += rspnt->val;
	    }
	    else {
               sens = rsens(rpnt,rspnt->wavel,rpnt->filt);
               rpnt->mflux[c] += rpnt->area * rspnt->val * sens * rpnt->attf * dqeff * stiminc;
	    }
	    rpnt->mask[c] = rspnt->mask;
	    rpnt->chanw[c] = 1.0;
	 }
/* ncfprintf (stderr,"xwavel %g mflux %g\n", rspnt->wavel,rpnt->mflux[c]);  /* */
/* ncfprintf (stderr,"absorb %g\n", rpnt->area * sens * * rpnt->attf * dqeff); /* */

#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"mask delta %d %g %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel);
#endif
         break;

       case 'c':                /* absolute masking light intensity */
          if (findphotrec (rspnt->recnm1, rspnt->recnm2, 
			   rspnt->recnm3, rspnt->recnm4,
					&rpnt, "runstim")) {
	    if (rpnt->ctype==VTRANSDUCER || rpnt->ctype==ITRANSDUCER) {
	       rpnt->mflux[c] = rspnt->val;
	    }
	    else {
               sens = rsens(rpnt,rspnt->wavel,rpnt->filt);
               rpnt->mflux[c] = rpnt->area * rspnt->val * sens * rpnt->attf * dqeff * stiminc;
	    }
	    rpnt->mask[c] = rspnt->mask;
	    rpnt->chanw[c] = 1.0;
	  }
/* ncfprintf (stderr,"xwavel %g mflux %g\n", rspnt->wavel,rpnt->mflux[c]);  */
/* ncfprintf (stderr,"absorb %g\n", rpnt->area * sens * rpnt->attf * dqeff); */

#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"mask abs  %d %g %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel);
#endif
         break;

       case 'd':                /* delta additive light intensity */
          if (findphotrec (rspnt->recnm1, rspnt->recnm2, 
			   rspnt->recnm3, rspnt->recnm4,
					&rpnt, "runstim")) {
	    if (rpnt->ctype==VTRANSDUCER || rpnt->ctype==ITRANSDUCER) {
	       rpnt->aflux[c] += rspnt->val;
	       rpnt->chanw[c] = 1.0;
	    }
	    else {	/* rod or cone */
	       sens = rsens(rpnt,rspnt->wavel,rpnt->filt);
               rpnt->aflux[c] += rpnt->area * rspnt->val * sens * rpnt->attf * dqeff * stiminc;
               if (rpnt->aflux[c] > MAXPHOT) rpnt->aflux[c] = MAXPHOT;
               if (rpnt->aflux[c] < 1e-10)   rpnt->aflux[c] = 0.0;
	    }
	    rpnt->chanw[c] = 1.0;
	  }
// ncfprintf (stderr,"inten %g\n", rspnt->val); /*  */
// ncfprintf (stderr,"xwavel %g aflux %g\n", rspnt->wavel,rpnt->aflux[c]); /*  */
// ncfprintf (stderr,"absorb %g sens %g\n", rpnt->area * sens * rpnt->attf * dqeff,sens); /* */
// ncfprintf (stderr,"area %g sens %g attf %g dqeff %g\n", rpnt->area, sens, rpnt->attf, dqeff); /* */
 /* ncfprintf(stderr,
 "stim num %d %d type '%c' val %10g wavel %5g time %15.15g rtime %15.15g stop %15.15g\n",
 rspnt->recnm1,rspnt->recnm2,rspnt->ctype,rspnt->val,rspnt->wavel,rspnt->time,rtime,stop);/**/

#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"stim delta %d %g %g\n",
                                        rspnt->recnm1,rspnt->val,rspnt->wavel);
#endif
         break;

       case 'e':                /* vclamp */
          if (npnt=findnode (rspnt->recnm1, rspnt->recnm2, 
			     rspnt->recnm3, rspnt->recnm4,
					 "runstim vclamp")) {
            if (npnt->comptr) setvclamp(npnt->comptr,rspnt->val,rspnt->bgend);
          }
#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"stim vclamp %d %d: %g %d\n",
                                        rspnt->recnm1,rspnt->recnm2,rspnt->val,rspnt->bgend);
#endif
         break;

       case 'g':                /* cclamp */
          if (npnt=findnode (rspnt->recnm1, rspnt->recnm2, 
			     rspnt->recnm3, rspnt->recnm4,
				 	"runstim cclamp")) {
	    if (npnt->comptr) {
              npnt->comptr->exti += rspnt->val;    /* set up current clamp */
	    }
	    else {
	      ncfprintf (stderr,"runstim: can't find compartment for node %s.\n",prnode(npnt));
	    }
          }

#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"stim cclamp %d %g\n",
                                        rspnt->recnm1,rspnt->val);
#endif
         break;

#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"stim cclamp %d %g\n",rspnt->recnm1,-rspnt->val);
#endif
         break;

       case 'i':                /* nt puff */
       case 'j':                
       case 'k':                
       case 'l':                
       case 'm':                
       case 'n':                
       case 'o':                
       case 'p':                
       case 'q':                
       case 'r':                
       case 's':                
       case 't':                
	   nt = GLU + rspnt->ctype - 'i';
          if (npnt=findnode (rspnt->recnm1, rspnt->recnm2,
			     rspnt->recnm3, rspnt->recnm4,
				 	"runstim puff")) {
	    if (!addnt(npnt->comptr,nt,rspnt->val))
	      ncfprintf (stderr,"runstim: can't find compartment for node %s.\n",prnode(npnt));
          }
#ifdef DEBUG
          if (debug & NCSTIM)  ncfprintf (stderr,"puff %s %d %g\n",
                                        findsym(nt),rspnt->recnm1,rspnt->val);
#endif
         break;

       case 'h':                /* Ca puff */
	 nt = CA;
         if (npnt=findnode (rspnt->recnm1, rspnt->recnm2,
			     rspnt->recnm3, rspnt->recnm4,
				 	"runstim puff")) {
		int cashell;
		double caval;
		cacomp *capnt;

	  if ((cpnt=npnt->comptr)==NULL) {
	      ncfprintf (stderr,"runstim: can't find compartment for node %s.\n",prnode(npnt));
	  } 
	  else if ((capnt=cpnt->capnt)==NULL) {
	      ncfprintf (stderr,"runstim: can't find Ca comp for node %s.\n",prnode(npnt));
	  }
	  else {			/* give puff of [Ca]i */
		caval = rspnt->val;
		cashell = capnt->cashell;
		capnt->cai += caval;
		for (i=0; i<cashell; i++) {
		  capnt->cais[i] += caval;
		}
	     }
         }
         break;

        default:
         break;

       }         /* switch */

	if (endfl) next = rspnt->next;	/* next in stim list */ 	 
        else       next = rspnt->rnext;	/* next in bin */
        delrstim(rspnt);                /* delete stimulus after use */
        rspnt = next;
     }        /* if */

    else {   /* no stimulus here */
      if (endfl) rspnt = rspnt->next;
      else       rspnt = rspnt->rnext;
    }
   }       /* for (rspnt=;;) */
  }       /* for (i=;;) */
 }      /* if (dmsec>=RECTIM) */
#ifdef DEBUG
  if (debug & NCSTIM && debugz & 1)  ncfprintf (stderr,"runstim end\n");
#endif
}

 
