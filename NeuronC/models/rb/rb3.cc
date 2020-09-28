/* file rb3.cc */
/* Simulation of rod and rod bipolar circuit */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ncio.h"
#include "ncsub.h"
#include "ncomp.h"
#include "ncfuncs.h"
#include "retcolors.h"
#include "gprim.h"
// #include "ncinit.h"
// #include "setexpt.h"

#include "scheduler.h"

extern scheduler sched;

 /* Identity numbers for neurons */ 

int xrod = 1;
int rbp  = 2;
int hbat = 3;
int soma    = 1000;		/* node numbers for parts of the cell */
int axon    = 1001;		/* axon of rb */
int axontip = 1002;

/* array parameters */

#define RODDENS   450000
#define MAXROD    500
#define MAXRODARRSIZ  5

//int rodarrsiz = 1;		/* width of rod array */
int rodarrsiz = MAXRODARRSIZ;		/* width of rod array */
double rodspac;			/* rod spacing, um */
int roddens = RODDENS;		/* number of cells / mm2 */
int maxrod  = maxrod;		/* maximum number of rods */

int rbzpos   = -8;         	/* z location of rod bipolar soma */
int hbzpos   = -2;         	/* z location of HB at terminal */

int disp_phot = 1;		/* display photons in flash during exp */

int makerods = 1;		/* = 1 -> make rods */
int makerbp	 = 1;		/* = 1 -> make rod bipolar cells */

/* biophysical parameters */

double rodvrev	= -.1056;	/* battery for rod Rm */
double rodrest	= -.040;	/* rod initial rest   */
double rbvrev	= -.050;	/* battery for RB Rm  */
double rbrest	= -.05;		/* RB initial rest    */
double rbprm    = 10000; 	/* RB Rm */

double rodrm   = 5000;
double rodtimec= .7;		

double hbvrev   = -0.05; 	/* HB vrev */
double hbrest   = -0.05; 	/* HB vrest */
double hbrm     = 5000; 	/* HB Rm */

/* Synaptic parameters */

double cexpon  = 0.73;//9;		/* gain of rod -> RB synapse: mv/efold change */
double sthresh = -0.0434;	/* -0.0434 */
double scov	  = 1;		/* Coeff of variation of vesicle timing =1 -> Poisson */
double svsiz	  = 300;	/* vesicle size */
double svgain	  = 2;		/* linear gain factor for vesicle release rate */
double skd	  = 1;		/* Kd for Glut binding at synapse */
//double scgain    = 1.5;	/* gain for cGMP cascade */
double scgain    = 1.2;		/* gain for cGMP cascade */
double condrb    = 650e-12;	/* forward conductance  rod -> rb */
double condhb    = 250e-12;	/* forward conductance  rod -> hb */
double scaperm   = .05;		/* Ca partial permeability through cGMP chan */
double scakd     = 100e-9;	/* Kd for Ca binding to cGMP channel */
double staua     = 1;
double staub     = 1;
double stauc     = 1;
double staud     = 1;

double bapta     = 0;
double bapta_kd  = 1e-6;
double scai      = 10e-9;
int    outt      = 0;
int    Gchan	 = 1; 

double rbtipdia  = .3;		/* diameter of rod bipolar dendritic tips */

double fbgain = 1e9;		/* gain of electrical feedback */
double hbsynfilt = 4;		/* time constant for hb PSPs */

double mxrot    = 120;		/* X rotation */

/* experiment parameters */

int only_report=0;

// const char *recmode    = "cclamp";	/* recording mode from RB */
const char *recmode    = "vclamp";

// const char *expt = "disp_model";	/* display model */
//const char *expt = "disp_flash";	/* display flash photons in rods */
//const char *expt = "flash_ir";	/* flash intensity - response run */
//const char *expt = "step_ir";	/* step intensity - response run */
//const char *expt = "flash_recov";	/* recovery, 2-flash */
const char *expt = "flash_8";		/* 8 flashes at different frequencies */


int ninfo = 0;			/* don't print anything except run */
//int ninfo = 1;		/* print just a little, stop */
//int ninfo = 3;		/* print info about rods connecting to rb */
//int ninfo = 1;

// int ca_mod = 0;
int ca_mod = 1;		/* =1 -> Ca through cGMP chan in RB dendrite */
			/*  also activates Ca plot */

int nonoise = 1;	/* =1 -> turn all noise off */

double dnoise;	/* dark continuous noise in rod */
int pnoise;	/* photon noise */
int vnoise;	/* vesicle noise */
int cnoise;	/* postsyn channel noise */

double srange;	/* synaptic range in terms of 1 photon response */

int rsyn[MAXROD];		/* array to hold enames for rod->rb synapse */
int rchan[MAXROD];		/* array to hold enames for rb channels */
int cabf[MAXROD];		/* array to hold enames for cacomps */
int hbsyn[MAXROD];		/* array to hold enames for rod->rb synapse */

int drod=0;
   
double spotdia;
int nflash;
int nrepeat;
double flash_inten;
double fbase;		/* flash step factor */
double flashdur;
double stimstart; 	/* pre-stimulus int */
double stimint; 	/* ISI */
double stimrecover; 	/* recovery time */
int use_inten_list;	/* use inten_list[] */
double inten_scale;	/* scale for inten_list*/
double ifstart;		/* inter-flash start factor */
double predur;

/**********************************************/

void setptrs(void)

/* must include here all variables to be set from command line */

{
   setptrn("ninfo",    &ninfo);
   setptrn("expt",     &expt);
   setptrn("recmode",  &recmode);

   setptrn("nonoise",  &nonoise);

   setptr("dnoise",  &dnoise);
   setptr("pnoise",  &pnoise);
   setptr("vnoise",  &vnoise);
   setptr("cnoise",  &cnoise);

   setptr("srange",  &srange);

   setptrn("rodvrev",  &rodvrev);
   setptrn("rodrest",  &rodrest);
   setptrn("rbvrev",    &rbvrev);
   setptrn("rbrest",   &rbrest);
   setptrn("rbprm",    &rbprm);
   setptrn("rodrm",    &rodrm);
   setptrn("rodtimec", &rodtimec);

   setptrn("fbgain",	&fbgain);
   setptrn("drod",	&drod);

   setptrn("hbrm",      &hbrm);
   setptrn("hbvrev",    &hbvrev);
   setptrn("hbrest",    &rbrest);

   setptrn("disp_phot", &disp_phot);
   setptrn("outt", &outt);
   setptrn("Gchan", &Gchan);

   setptrn("rodarrsiz", &rodarrsiz);

   setptr("spotdia",     &spotdia);
   setptr("nflash",      &nflash);
   setptr("nrepeat",     &nrepeat);
   setptr("flash_inten", &flash_inten);
   setptr("fbase",       &fbase);		/* flash step factor */
   setptr("flashdur",    &flashdur);
   setptr("stimstart",   &stimstart); 		/* pre-stimulus int */
   setptr("stimint",     &stimint); 		/* ISI */
   setptr("stimrecover", &stimrecover); 	/* recovery time */
   setptr("use_inten_list",&use_inten_list);	/* use inten_list[] */
   setptr("inten_scale", &inten_scale);		/* scale for inten_list*/
   setptr("ifstart",	&ifstart);		/* inter-flash int start factor */
   setptr("predur",	&predur);		/* interval before start of plot */
}

/*--------------------------------------------*/

int mid(int siz)

/* function to return center element of array */

{
   int m;

  if (int(siz/2)*2==siz) m = (siz+1) * siz / 2;
  else                   m = (siz*siz-1) / 2;
  return m;
}


/*--------------------------------------------*/

void makrod (double xpos, double ypos, int n) 

 /* procedure to make simple rod: one compartment */

{
   double dia;
   node *np;
   sphere *s;
   photorec *p;

   np = loc(nd(xrod,n), xpos,ypos,0);
   s=make_sphere(np,dia=8); s->Rm=rodrm; s->vrest=rodrest, s->vrev=rodvrev;
   p = make_rod(np); p->xpos=xpos; p->ypos=ypos; p->timec1=rodtimec;
       p->maxcond=35e-12; p->photnoise=pnoise; p->darknoise=dnoise;
}

/*-----------------------------------------*/

void make_rbp (int nrod, double xpos, double ypos) 

/* Make rod bipolar and connect rods to it. */

{
    int r, cn, nfilt;
    double dia;
    double catau, dtcai;

    elem *e; 
    cattrib *ca;
    node *np;
    cable *c;
    sphere *s;
    synapse *sp;
    nattrib *napnt;
    chattrib *cpnt;

  /* first, make rb skeleton */

   cn = 0;
   s = make_sphere(loc(nd(rbp, cn, soma),xpos,ypos,rbzpos), dia=5, rbprm);
   s->vrest = rbrest;
   s->vrev  = rbvrev;
   c = make_cable (ndn(rbp,cn,soma), loc(nd(rbp,cn,axon),xpos,ypos,rbzpos-8));    
	c->Rm = rbprm; c->dia=0.8; c->vrest=rbrest; c->vrev=rbvrev;
   c = make_cable (ndn(rbp,cn,axon), loc(nd(rbp,cn,axontip),xpos,ypos,rbzpos-16)); 
	c->Rm = rbprm; c->dia=1.2; c->vrest=rbrest; c->vrev=rbvrev;

  /* next, connect rods */

  for (r=0; r<nrod; r++) {
				/* make rb dendrite */
   np = nd(xrod,r);
   c = make_cable(loc(nd(rbp,cn,r),np->xloc,np->yloc,np->zloc), ndn(rbp,cn,soma));
	c->Rm = rbprm; c->dia=rbtipdia; c->cplam=0.05; c->vrest=rbrest; c->vrev=rbvrev;
   
   //s = make_sphere(ndn(rbp,cn,r),dia=1, rbprm); s->vrest=rbrest; s->vrev=rbrev;  /* add some volume */

   catau = 200e-3;			/* sets vmax for capump */
   dtcai = scai;			/* [Ca]i threshold */

   if (ca_mod) {
	e = at(ndn(rbp, cn, r), CACOMP);                /* Ca comp at soma */
	ca = (cattrib *)make_chan(e, CACOMP, 0);
	ca->pump=1; ca->vmax=2e-6/catau; ca->pkm=10e-6;
	ca->cai=scai;
	ca->cashell=3;
	ename(e,&cabf[r]); 
	if (bapta) {
	   ca->cabuf=1; ca->bmax=1e8; ca->bkd=bapta_kd; ca->btot=1e4; ca->btoti=1e4;
	}
    }

				/* connect rod to rb with synapse */	

     sp = make_synapse(nd(xrod,r), nd(rbp,cn,r));
          sp->ntact  = CLOSE;
          sp->curve  = EXPON;
          sp->ngain  = cexpon;
          sp->thresh = sthresh;
          sp->vgain  = svgain;
          sp->maxcond= 0;
          sp->trconc = 8e-6;
          sp->nkd = skd;
          sp->npow = 2;
	  nfilt = 2;
          sp->nfilt2 = nfilt;
          sp->timec2 = makfiltarr(nfilt,0,sp->timec2,4);
          sp->nfilt3 = nfilt;
          sp->timec3 = makfiltarr(nfilt,0,sp->timec3,25);
          sp->secmsg = 1;
          sp->cgain=scgain;
          sp->mesg2 = CGMP;
          if (vnoise) {
            napnt=make_vesnoise(sp);
            napnt->vsize = svsiz;
            napnt->cov  = scov;
          }
          ename(sp,&rsyn[r]);
   
          e = at(nd(rbp,cn,r),CHAN);
            cpnt = make_chan(e, CGMP, 11);
                cpnt->maxcond = condrb; 
                cpnt->vrev    = 0;
                cpnt->caperm  = scaperm;
                cpnt->taua    = staua; cpnt->taub  = staub;
                cpnt->tauc    = stauc; cpnt->taud  = staud;
                if (cnoise>0) {
                  napnt=make_chnoise(e);
                  napnt->unitary=20e-12;
                }
                ename (e,&rchan[r]);

        if (ninfo>=3) printf ("#  conn %2d %3d to %2d %2d %2d\n",
				xrod,r,rbp,cn,r); /* */

				/* connect rod to HBat */

     s  = make_sphere(loc(nd(hbat, cn, r),xpos,ypos,hbzpos), dia=3, hbrm);
     s->vrest = hbrest;
     s->vrev  = hbvrev;
     sp = make_synapse(ndn(xrod,r), nd(hbat,cn,r));
          sp->ntact  = OPEN;
          sp->curve  = EXPON;
	  sp->dyadelem = rsyn[r];
          sp->ngain  = cexpon;
          sp->thresh = sthresh;
          sp->vgain  = svgain;
          sp->maxcond= condhb;
	  nfilt = 2;
          sp->nfilt2 = nfilt;
          sp->timec2 = makfiltarr(nfilt,0,sp->timec2,hbsynfilt);
          ename(sp,&hbsyn[r]);
   
 } /* for (r=0; r<nrod; ) */

if (ninfo>=2) {
   printf ("#\n");
   printf ("# make_rbp: %d rod connections\n",nrod);;
 }
} /* make_rbp */

/* - - - - - - - - - - - - - - - - - - - - - - */

void rod_phot_dr (int type, int pigm, int color, double dscale, double dia, 
		double foreshortening, int hide) 

/* Procedure to draw rods with color */
/*  Runs automatically when displaying rods */
/*  Set up with "set_phot_dr()", */
/*  see "phot_dr" in manual */

{
      int num;
      char cbuf[20];

    num = color;
    gpen (color);
    gcirc (dia*0.4,1);
    if (strcmp(expt,"disp_flash")==0) {
      gframe("rodphot");
      grotate(-45);
      gpen (0);
      gmove (-0.12,-0.12);
      gcwidth (0.4);
      sprintf (cbuf,"%d",num);
      gtext (cbuf);
      gpen (6);
      gframe("..");
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - */

int rodnum;
int midrod;

double fdur = 0.0001;

double phtot (double val, double xtime)

/* Count all the photons captured in the rod array. */

{
     int i;
     double ptot;

   for (ptot=i=0; i<rodnum; i++) {
      ptot += l(ndn(xrod,i))*fdur;
   }
   return ptot;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

#define CCMAPSIZ 100

int ccmap[CCMAPSIZ];

void init_phot(double x, double y, double siz)

/* procedure to initialize photon display */

{
     int i;
     double s, size;
     double xrot, yrot, zrot;
     double xcent, ycent, zcent;

  gframe ("phot");
  gorigin (x,y);
  gsize (siz);

  size = rodarrsiz*rodspac*2.5;
  set_disp_rot (xrot=0,yrot=0,zrot=0,xcent=0.8,ycent=0.8,zcent=rbzpos+2,0,0,0,size);
  for (i=0; i<CCMAPSIZ; i++) ccmap[i] = i;
  setcmap(ccmap,CCMAPSIZ);
  gpen (1);
  s=.5;
  gmove  (.2,.21);
  grdraw (0,s);
  grdraw (s,0);
  grdraw (0,-s);
  grdraw (-s,0);
  gframe ("..");
};

/* - - - - - - - - - - - - - - - - - - - - - - */

void print_phot (int flashn, double photn)

{
    char cbuf[50];

  if (photn==1)
  sprintf(cbuf,"Flash %d, %3.3g photon.",flashn,photn);
  else
  sprintf(cbuf,"Flash %d, %3.3g photons.",flashn,photn);
  gtext(cbuf);
}


/* - - - - - - - - - - - - - - - - - - - - - - */

void run_phot(int i)

/* procedure to run photon display */

/* Note that several special colors allow spatial displays: */

/*  vcolor = voltage
    lcolor = light
    cacolor = [Ca]i
*/

{
    double ph;
    double x,y;
    int color, cmap;
    double dscal,vmax,vmin;
    static double oldph;
 
  if (!outt) {
   if (i==0) oldph=0;
   gframe ("phot");
   display (ROD,MATCHING,ndt(xrod,-1,-1),0,color=LCOLOR,vmax=CCMAPSIZ,vmin=0,cmap=10,dscal=0.4);
   gframe ("..");
   ph = phtot(0,0);
   gframe ("phot_count");
   gcwidth(.02);
   x = .05;
   y = .96;
   gmove (x,y);
   gpen (0);
   print_phot(i,oldph);
   gmove (x,y);
   gpen (i+1);
   print_phot(i+1,ph);
   gpurge();
   oldph = ph;
   gframe ("..");
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - */

void stimrods (double sinten, double sstart, double sdur) 

/* stimulate rods with single-photon responses, without photon noise */

{
    int i, x, nphot;
    double eff; 
    int phot_given;		/* number of photons given in stimulus so far */
    static int stimrod[MAXROD];
    double inten, start, dur;

   eff = .8889;			/* rod area * sens * filt * qeff */
				/* 2.2      * .67  * .9   * .67  */ 

   nphot = int(sinten*sdur*rodarrsiz*rodarrsiz*eff);
   for (i=0; i<rodnum; i++) {
       stimrod[i]=0;
   }
   for (x=phot_given=0; phot_given < nphot; x) {
     for (i=0; i<rodnum; i++) {
       stimrod[i] ++;			/* Make array that contains number */
       if (++phot_given>=nphot) break;  /* of photons each rod should receive */
     }
   }
   for (i=0; i<rodnum; i++) {
       stim_rod(ndn(xrod,i), inten=stimrod[i]/(eff*sdur),start=sstart,dur=sdur);
   }
}

/* - - - - - - - - - - - - - - - - - - - - - - */

void disp_plots (void)

{
   int j, nstates;
   int cn;
   double glur_siz;
   double plg;			/* gain of trace */
   double offb;			/* zero for trace */
   double offtr;		/* position of trace within graph */
   int pen,plnum; 
   double plsize;
   double ft, max, min;
   elem *epnt;
   chan *cpnt;
   char chanstr[20];

  if (disp_phot) disp = 33;	/* set display mode to save elems */
  cn = 0;

  if (strcmp(recmode,"vclamp")==0) {	 /* plot current from RB voltage clamp */
   plg   =  30e-12;
   offb  = -30e-12;
   offtr = -0.1;
   plg *= (double)rodnum / (1+ca_mod);
   offb *= (double)rodnum / (1+ca_mod);
   plot (I, ndn(rbp,cn,soma), max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
   				plot_param( "Irb", pen=2, plnum=1, plsize=1.2);

  } else {			/* plot voltage in current clamp mode */
   plg   =  .03;
   offb  = -.06;
   offtr = 0;
   plot (V, ndn(rbp,cn,soma), max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
				plot_param("Vrb",pen=2, plnum=2, plsize=0.8);
  }
 
  if (ca_mod) {				/* Ca postsyn to center rod */
    if(1) {	
    plg = 10e-6;
    offb  = 0;
    offtr = 0;
    plot (CA, 1, ndn(rbp,cn,drod), max=(0-offtr)*plg+offb, min=(1-offtr)*plg+offb);
				plot_param("CaRB", pen=5, plnum=3, plsize=0.5);
    }  
  }

 if (Gchan) {

  plg = 1;
  offtr = 0;
  offb = 0;

  if ((nstates=get_chan_nstate(rchan[drod]))>0) {
    for (j=0; j< nstates; j++) {
      sprintf (chanstr,"Gchan%d",j);
      plot (G, j+1, rchan[drod],  max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
				 plot_param(chanstr,pen=j+1,plnum=5,plsize=0.7);
    }  /* for (j;;) */
  }
 }

  /* Plot cGMP concentration */
  if (0) {
  plg = 2;
  offtr = 0;
  offb = 0;
/* plot FC9 rsyn[drod] max (1-offtr)*plg+offb 
		        min (0-offtr)*plg+offb pen 9 plnum 6 plname "cGMP";  */

  plg *= 20e-6;		/* cGMP in center RB dend */  
  plot (CGMP,ndn(rbp,cn,drod), max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb); 
			plot_param("cGMP",pen=13,plnum=6, plsize=1);
  }

/* plot rod->rb vesicle release rate, ves/s  */
	plg = 200;
	offtr = 0;
	offb = 0;
	plot (FA9,rsyn[drod],max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
		 plot_param("VesRate", pen=11, plnum=8,plsize=0.8);
  if (0) {
    for(j=0;j<(rodarrsiz*rodarrsiz-1);j++) {
	plot (FA9,rsyn[j],max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
		 plot_param("VesRate", pen=(j+1), plnum=8,plsize=0.8);

    }
  }  /* if (0) */

	glur_siz = 0.5;

	plg = 1;
	offtr = 0;
	offb = 0;
	plot (FC0,rsyn[drod],max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
			plot_param("GluRbnd", pen=4, plnum=7,plsize=glur_siz);
	plg = 10;
	offtr = 0;
	offb = 0;
	plot (FB2,rsyn[drod],max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
			plot_param("Ves", pen=5, plnum=7,plsize=glur_siz);
	plg = 1;
	offtr = 0;
	offb = 0;
	plot (FC2,rsyn[drod],max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
			plot_param("GbndFilt", pen=2,plnum=7,plsize=glur_siz);
	plg = 20e-6;
	offtr = 0;
	offb = 0;
	plot (CGMP,ndn(rbp,cn,drod),max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
			plot_param("CycG", pen=1,plnum=7,plsize=glur_siz);

  if (0) {

  ft = .05;				/* filter time const */
  plg   = .04;
  offtr = .5;
  offb  = 0;
  double timec[]={ft,ft,ft,ft,ft}; 
  plotfilt (5,timec);    				/* total photons */ 
  plot_func (phtot,   0,  max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb);
		          plot_param("PhTot", pen=1,plnum=9,plsize=1);
  }

  plg = .005;
  offtr = .6;
  offb = -.04;
  plot(V,ndn(xrod,drod), max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb); 
		       plot_param("RodV", pen=1,plnum=10,plsize=1);
  if (0) {
  for(j=0;j<24;j++) {

//  plot (V, ndn(xrod,midrod), max=(1-offtr)*plg+offb,	/* voltage in center rod */ 
//		       min=(0-offtr)*plg+offb,
//		       "RodV", pen=6, plnum=10, plsize=1); 
							/* voltage in first rod */ 
  plot(V,ndn(xrod,j), max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb); 
		       plot_param("RodV", pen=(j+1),plnum=10,plsize=1);
  }
  } /* if (0) */

  /* plot the HB current and voltage */

  plg = condhb * 0.05;
  offtr = 0.9;
  offb =  0;
  plot(G,I,hbsyn[drod], max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb); 
		       plot_param("Ihb", pen=9,plnum=12,plsize=0.5);

  plg   = 0.05;
  offtr = 0.1;
  offb =  -0.05;
  plot(V,ndn(hbat,cn,drod), max=(1-offtr)*plg+offb, min=(0-offtr)*plg+offb); 
		       plot_param("Vhb", pen=10,plnum=13,plsize=0.5);

  //plot S time pen 0 plnum 15; /* dummy plot to allow space for photon disp */

} /* end proc disp_plots() */

/**********************************************/

synap *dyadhb[MAXROD];
synap *synrodrb[MAXROD];
chan *synrodrbchan[MAXROD];


void init_elec_feedback(void)

/* initialize electrical feedback */

{
   elem *epnt;
   int i;

  for (i=0; i<rodnum; i++) {

     /* find rod -> rb synapse, save pointer */
    if ((epnt=findelem(rsyn[drod]))==NULL) {
       fprintf (stderr,"rb: can't find synapse %d\n",rsyn[drod]);
    }
    if ((synrodrb[i]=(synap*)epnt->lptr)==NULL) {
       fprintf (stderr,"rb: rb synapse %d not created.\n",rsyn[drod]);
    }

     /* find rod -> rb chan, save pointer */
    if ((epnt=findelem(rchan[drod]))==NULL) {
       fprintf (stderr,"rb: can't find chan %d\n",rchan[drod]);
    }
    if ((synrodrbchan[i]=(chan*)epnt->lptr)==NULL) {
       fprintf (stderr,"rb: rb synapse %d not created.\n",rchan[drod]);
    }

     /* find rod -> hb synapse, save pointer */
    if ((epnt=findelem(hbsyn[drod]))==NULL) {
       fprintf (stderr,"rb: can't find synapse %d\n",hbsyn[drod]);
    }
    if ((dyadhb[i]=(synap*)epnt->lptr)==NULL) {
       fprintf (stderr,"rb: hb synapse %d not created.\n",hbsyn[drod]);
    }
  }
}

/*--------------------------------------------*/

void elec_feedback(void)

{
    int i;
    static int runyet=0;
    double v;
    double dyadhbi, rbi, svoffset;


  // fprintf (stderr,"elec_feedback: %g\n",simtime);

  if (!runyet) init_elec_feedback();
  runyet = 1;
  for (i=0; i<rodnum; i++) {
    svoffset = synrodrb[i]->thresh - sthresh;
    v = dyadhb[i]->vrev + svoffset - dyadhb[i]->comp2->v;
    dyadhbi = dyadhb[i]->conduct * -v;
    v = synrodrb[i]->vrev +svoffset - synrodrb[i]->comp2->v;
    rbi = synrodrbchan[i]->conduct * -v;
    // rbi = 0;
    svoffset = (dyadhbi + rbi) * fbgain;
    synrodrb[i]->thresh = sthresh + svoffset;
    if (ninfo >= 2) fprintf (stderr,"elec_feedback: n %d rbi %g\n",i, rbi);
    if (ninfo >=2) fprintf (stderr,"elec_feedback: n %d voff %g\n",i, svoffset);
  }
}

/**********************************************/

main(int argc, char **argv)
{
   int i, j, n, ct;
   double rb_version;;
   double scal, nscal;
   FILE *filout;
   char *chseed;
   char cbuf[20];
   int rwid, seed;
   double inten_list[] = {2,6,18,82,230,606,1900,4410,10000};
   int nf = sizeof(inten_list)/sizeof(double);

 ncinit(argc,argv);
 filout = stderr;
 rb_version = 3.001;

 timinc   = 1e-4;
 ploti    = 1e-3;
 crit     = 1e-8; 
 implicit = 0;
 plsep = 1;		/* =1 -> separate plots */

// chseed = system ("date +%m%d%H%M%S");		/* random seed */
// sscanf (chseed,"%d",&seed);	
                               /* Note, can set setrand(-1) for similar effect. */
 setrand(15);

/*-------------------------------------------------------*/

 setvar();		/* set variables from command line */

if (notinit(nonoise)) nonoise = 1;	/* =1 -> turn all noise off */

if (nonoise) {
 if (notinit(dnoise)) dnoise = 0;	/* dark continuous noise in rod */
 if (notinit(pnoise)) pnoise = 0;	/* photon noise */
 if (notinit(vnoise)) vnoise = 0;	/* vesicle noise */
 if (notinit(cnoise)) cnoise = 0;	/* postsyn channel noise */
}
else {
 if (notinit(dnoise)) dnoise = 0.2;	/* dark continuous noise in rod */
 if (notinit(pnoise)) pnoise = 1;	/* photon noise */
 if (notinit(vnoise)) vnoise = 1;	/* vesicle noise */
 if (notinit(cnoise)) cnoise = 0;	/* postsyn channel noise */
}
pnoise = 1;

if (notinit(srange)) srange=1;

fprintf (stderr,"rodtimec %g srange %g\n",rodtimec,srange);

/* synaptic parms for different dynamic ranges */

/* To calibrate, set big plot on GluRbnd and set resting
   level to .9, and peak to .1 by varying cexpon, sthresh. */

if (rodtimec>=0.9) {

  if (srange==.5) {		/* 0.5 photon dynamic range */
     cexpon  = 0.44;		/* harder saturation */
     sthresh = -0.0418; 
  }
  else if (srange==1) {		/* 1 photon dynamic range */
     cexpon  = 0.59;		/* hard saturation */
     sthresh = -0.0424; 
  }
  else if (srange==1.5) {		/* 1.5 photon dynamic range */
     cexpon  = 0.84;
     sthresh = -0.0434; 
  }
  else if (srange==2) {		/* 2 photon dynamic range */
     cexpon  = 1.12;
     sthresh = -0.0446;
  }
  else if (srange==3) {		/* 3 photon dynamic range */
     cexpon  = 1.68;
     sthresh = -0.0468;
  }
  else if (srange==4) {		/* 4 photon dynamic range */
     cexpon  = 2.2;
     sthresh = -0.0490;
  }
}
else if (rodtimec>=0.8){

  if (srange==.5) {		/* 0.5 photon dynamic range */
     cexpon  = 0.42;		/* harder saturation */
     sthresh = -0.0423; 
  }
  else if (srange==1) {		/* 1 photon dynamic range */	
     cexpon  = 0.52;		/* hard saturation */
     sthresh = -0.0427; 
  }
  else if (srange==1.5) {		/* 1.5 photon dynamic range */
     cexpon  = 0.80;
     sthresh = -0.0438; 
  }
  else if (srange==2) {		/* 2 photon dynamic range */
     cexpon  = 1.12;
     sthresh = -0.0451;
  }
}
else if (rodtimec>=0.6){

  if (srange==.5) {		/* 0.5 photon dynamic range */
     cexpon  = 0.42;		/* harder saturation */
     sthresh = -0.0423; 
  }
  else if (srange==1) {		/* 1 photon dynamic range */	
     cexpon  = 0.55;		/* hard saturation */
     sthresh = -0.04365; 
  }
  else if (srange==1.5) {		/* 1.5 photon dynamic range */
     cexpon  = 0.83;
     sthresh = -0.0448;
  }
  else if (srange==2) {		/* 2 photon dynamic range */
     cexpon  = 1.12;
     sthresh = -0.0460;
  }
  else if (srange==3) {		/* 2 photon dynamic range */
     cexpon  = 1.66;
     sthresh = -0.0460;
  }
}

/* for single rod -> rb sim */

//  sthresh = -0.0439; 
   
if (notinit(ca_mod)) ca_mod = 0;  /* allow user to override default */
if(info>0) printf("# vsize= %g , dnoise= %g  CoV= %g \n", svsiz, dnoise, scov);

if (only_report==1) exit(0);

/*--------------------------------------------*/

 set_phot_dr(rod_phot_dr);

//erase model;

/* Make array of rods */

if (makerods) {

  rodspac   = 1 / sqrt (roddens* 1e-6);         /* rod spacing */
  if (notinit(rodarrsiz)) rodarrsiz=1;

  if (info>=2) {
    printf ("# rod array size: %g\n",rodarrsiz);
  }
  rodnum = 0;
  rwid = rodarrsiz / 2;
    for (j= -rwid; j<=rwid; j++) 
    for (i= -rwid; i<=rwid; i++) {
      makrod (i*rodspac, j*rodspac, rodnum++);
    }
  if (info>=2)
    printf ("# Number of rods: %d ", rodnum);
  midrod = mid(rodarrsiz);
  if (info>=2) printf ("# midrod=%d\n", midrod);

}  /* make rods */

else {
  midrod = 0;
}

/*-------------------------------------------------*/

   /* Make array of rod bipolars */

if (makerbp) {
  if (makerods) {
     double rbxloc,rbyloc;
    make_rbp (rodnum, rbxloc=0, rbyloc=0);
  }  /* if makerods */
  if (info>=2) printf ("# rod bipolar cell done.\n");
} /* if makerbp */

/*-------------------------------------------------*/

/* Display */

if (strcmp(expt,"disp_model")==0) {
     double len;
     double dscale;
     int dcolor;
     double xrot, yrot, zrot;
 
  disp=1;					/* allow display of anatomy */
  display_center (0,0,rbzpos-2);
 // display_size (rodarrsiz*rodspac*2);
  display_size (30);
  display_calib (len=5.0,red);
  display_rot (xrot=mxrot,yrot=-5,zrot=0);

  display (SPHERE, MATCHING, ndt(xrod,-1), dcolor=4, dscale= 0.15);

  if (makerbp) {
      display (SPHERE,  MATCHING, ndt(rbp,-1,-1), dcolor=2);
      display (CABLE,   MATCHING, ndt(rbp,-1,-1), dcolor=2);
      display (SYNAPSE, MATCHING, ndt(rbp,-1,-1), dcolor=2);
  } 
  exit(0);
}

/*-------------------------------------------------*/

 gmove (0.5,0.95);
 gpen (1);
 sprintf (cbuf,"fbgain %g",fbgain);
 gtext(cbuf);

/*-------------------------------------------------*/

  sched.addTask(1e-4, elec_feedback);		/* set electrical feedback at runtime */

/*-------------------------------------------------*/

/* Experiments	*/

if (strcmp(expt,"flash_ir")==0) {	/* flash intensity-response, overlaid */

    int r,cn;
    double start, dur;
    double runlen;
    double fl_inten, inten;
   cattrib *apnt;
   elem *epnt;

  cn = 0;
  if (strcmp(recmode,"vclamp")==0) {
     vclamp (ndn(rbp,cn,soma), rbrest,start=0,dur=100);
  }

  disp_plots();				/* make plots */

   if (notinit(spotdia))      spotdia     = 100;
   if (notinit(nflash))       nflash      = nf;
   if (notinit(flash_inten))  flash_inten = 1e3;
   if (notinit(fbase))        fbase       = 2;		/* flash step factor */
   if (notinit(flashdur))     flashdur    = .0001;
   if (notinit(stimstart))    stimstart   = 0.1; 	/* pre-stimulus int */
   if (notinit(stimint))      stimint     = 2; 		/* ISI */
   if (notinit(stimrecover))  stimrecover = .01; 	/* recovery time */
   if (notinit(use_inten_list)) use_inten_list = 1;	/* use inten_list[] */
   if (notinit(inten_scale))  inten_scale = 20;	/* scale for inten_list*/
   runlen = stimstart+stimint;
   endexp = runlen;				/* end of expt */

   if (notinit(disp_phot)) disp_phot=0;

   if (disp_phot) init_phot(.78,.8,.25);	/* initialize photon display */

   if (notinit(predur))	predur = 0.5; 		/* equilibrium time */
   simtime = 0 - predur;
   setxmin = 0;
   step (predur);
   for (i=0; i<nflash; i++){

      	simtime = 0;
	if (use_inten_list) fl_inten = inten_list[i]*inten_scale;
	else                fl_inten = flash_inten * pow(fbase,i);

      	if (!pnoise) stimrods(fl_inten, stimstart, flashdur);
     	else stim_spot (spotdia,0,0,inten=fl_inten,start=stimstart,dur=flashdur);

     	step (stimstart+flashdur);
     	if (disp_phot) run_phot(i);   /* display photons in flash during expt */
	fprintf (stderr,"expt %s: Flash %d: inten %7.3g\n",expt, i+1, fl_inten*flashdur);

	step (stimint-stimrecover);

     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = 1000;
        apnt->taub = 0.001;
        apnt->tauc = 1000;
        apnt->taud = 0.001;
        apnt->caflg = 1;
        apnt->cai = scai;
     }
     step (stimrecover);
     for (r=0; r<rodnum; r++) {		/* modify chans for normal recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = staua;
        apnt->taub = staub;
        apnt->tauc = stauc;
        apnt->taud = staud;
     }
     if (outt) printf ("\n");
   }

 }  /* expt==flash_ir */

else

/*-------------------------------------------------*/

if (strcmp(expt,"step_ir")==0) {			/* step intensity-response, overlaid */

   int i, r,cn;
   double start, dur;
   double runlen;
   double fl_inten, inten;
   cattrib *apnt;
   elem *epnt;

  cn = 0;
  if (strcmp(recmode,"vclamp")==0) {
     vclamp (ndn(rbp,cn,soma), rbrest, start=0, dur=100);
  }

  disp_plots();				/* make plots */

   if (notinit(spotdia))      spotdia     = 100;
   if (notinit(nflash))       nflash      = nf;
   if (notinit(flash_inten))  flash_inten = 5;
   if (notinit(fbase))        fbase       = 2;		/* flash step factor */
   if (notinit(flashdur))     flashdur    = 1500e-3;
   if (notinit(stimstart))    stimstart   = 0.2; 	/* pre-stimulus int */
   if (notinit(stimint))      stimint     = 3000e-3; 	/* ISI */
   if (notinit(stimrecover))  stimrecover = .05; 	/* recovery time */
   if (notinit(use_inten_list)) use_inten_list = 1;	/* use inten_list[] */
   if (notinit(inten_scale))  inten_scale = 0.2;		/* scale for inten_list*/

   //endexp = nflash*stimint+stimstart;		/* end of expt */

   runlen = stimstart+stimint;
   endexp = runlen;

   disp_phot = 0;
   if (notinit(disp_phot)) disp_phot=0;
   if (disp_phot) init_phot(.78,.8,.25);	/* initialize photon display */

   if (notinit(predur)) predur = .5;		/* equilibrium time */
   simtime = 0 - predur;
   setxmin = 0;
   step (predur);
   for (i=0; i<nflash; i++){

     simtime = 0;
     if (use_inten_list) fl_inten = inten_list[i]*inten_scale;
     else                fl_inten = flash_inten * pow(fbase,i);
     fprintf (stderr,"expt %s: Flash %d: inten %5.3g\n",expt, i+1, fl_inten);

     if (!pnoise) stimrods(fl_inten, stimstart, flashdur);
     else stim_spot (spotdia, 0, 0, inten=fl_inten,start=stimstart,dur=flashdur);

     step (stimstart+0.0002);
     if (disp_phot) run_phot(i);    /* display photons in flash during expt */
     step (stimint-stimrecover);

     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = 1000;
        apnt->taub = 0.001;
        apnt->tauc = 1000;
        apnt->taud = 0.001;
        apnt->caflg = 1;
        apnt->cai = scai;
     }
     step (stimrecover);
     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = staua;
        apnt->taub = staub;
        apnt->tauc = stauc;
        apnt->taud = staud;
     }
     if (outt) printf ("\n");
   }

 }  /* expt==step_ir */

else

/*-------------------------------------------------*/

if (strcmp(expt,"flash_recov")==0) {			/* recovery, 2-flash */

   int i,r,cn;
   double start, dur;
   double stimstart2, stimint2;
   double inten, runlen;
   cattrib *apnt;
   elem *epnt;

  cn=0;
  if (strcmp(recmode,"vclamp")==0) {
     vclamp (ndn(rbp,cn,soma), rbrest, start=0, dur=100);
  }

  disp_plots();				/* make plots */

   if (notinit(spotdia))      spotdia     = 100;
   if (notinit(nrepeat))      nrepeat     = 6;
   if (notinit(flash_inten))  flash_inten = 4e4;
   if (notinit(fbase))        fbase       = 1;	 /* flash step factor */
   if (notinit(flashdur))     flashdur    = .0001;
   if (notinit(stimstart))    stimstart   = .1;  /* pre-stimulus int */
   if (notinit(stimint))      stimint     = 2;   /* repeat interval */
   if (notinit(ifstart))      ifstart     = .2;  /* inter-flash start factor */
   if (notinit(stimrecover))  stimrecover = .01; /* recovery time */
 
   runlen = stimstart+stimint;
   endexp = runlen;				/* end of expt */

   if (notinit(disp_phot)) disp_phot=0;

   if (disp_phot) init_phot(.78,.8,.25);	/* initialize photon display */

   if (notinit(predur))	      predur	  = .5; 	/* equilibrium time */
   simtime = 0 - predur;
   setxmin = 0;
   step (predur);
   for (i=0; i<nrepeat; i++){
        simtime = 0;

	fprintf (stderr,"expt %s: Flash %d: interval %5.3g\n",
				expt, i+1, ifstart*(i+1));

        /* first flash */

      	if (!pnoise) stimrods(flash_inten * pow(fbase,i), stimstart, flashdur);
      	else stim_spot (spotdia,0,0, inten=flash_inten*pow(fbase,i), stimstart, flashdur);

        step (stimstart+flashdur);
     	if (disp_phot) run_phot(i*2);  /* disp photons in flash during expt */

        /* second flash */
        stimstart2 = stimstart + ifstart * (i+1);

      	if (!pnoise) stimrods(flash_inten * pow(fbase,i), stimstart2, flashdur);
      	else stim_spot (spotdia,0,0,inten=flash_inten * pow(fbase,i), stimstart2, flashdur); 

     	stimint2 = stimstart2-stimstart;
     	step (stimint2);
     	if (disp_phot) run_phot(i*2+1); /* disp photons in flash during expt */

	step (stimint-stimstart-stimint2-stimrecover);

     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = 1000;
        apnt->taub = 0.001;
        apnt->tauc = 1000;
        apnt->taud = 0.001;
        apnt->caflg = 1;
        apnt->cai = scai;
     }
     step (stimrecover);
     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = staua;
        apnt->taub = staub;
        apnt->tauc = stauc;
        apnt->taud = staud;
     }
     if (outt) printf ("\n");
   };

}  /* expt==recov */

else

/*-------------------------------------------------*/

if (strcmp(expt,"flash_8")==0) {			/*8-flash sequence */

#define NREPEAT 6

   int i, r, cn;
   double start, dur;
   double inten, runlen;
   double frame_period, flashstart;
   double flashintervals[NREPEAT] = {25,20,10,5,2,1};
   cattrib *apnt;
   elem *epnt;

  if (strcmp(recmode,"vclamp")==0) {
     vclamp (ndn(rbp,cn=0,soma), rbrest, start=0, dur=100);
  }

  disp_plots();				/* make plots */

   if (notinit(spotdia))      spotdia     = 100;
   if (notinit(nrepeat))      nrepeat     = NREPEAT;
   if (notinit(nflash))       nflash      = 8;
   if (notinit(flash_inten))  flash_inten = 4e4;
   if (notinit(fbase))        fbase       = 1;	 /* flash step factor */
   if (notinit(flashdur))     flashdur    = .0001;
   if (notinit(stimstart))    stimstart   = .2;  /* pre-stimulus int */
   if (notinit(stimint))      stimint     = 3;   /* repeat interval */
   if (notinit(ifstart))      ifstart     = .2;  /* inter-flash start factor */
   if (notinit(stimrecover))  stimrecover = .01; /* recovery time */

   frame_period = 0.01333333;

   for (i=0; i<NREPEAT; i++) {
     flashintervals[i]= flashintervals[i]* frame_period;
   }

   runlen = stimstart+stimint;
   endexp = runlen;				/* end of expt */

   if (notinit(disp_phot)) disp_phot=0;

   if (disp_phot) init_phot(.78,.8,.25);	/* initialize photon display */

   if (notinit(predur))	predur = .5; 		/* equilibrium time */
   simtime = 0 - predur;
   setxmin = 0;
   step (predur);
   for (i=0; i<nrepeat; i++){
        simtime = 0;
	fprintf (stderr,"expt %s: Flash %d: interval %5.3g\n",
				expt,i+1, flashintervals[i]);
	step (stimstart);
        /* 8 flashes */
	for(j=0;j<nflash;j++) {
        	flashstart = simtime;
 		//  fprintf(stderr,"%g %g\n", simtime,flashstart);

	      	if (!pnoise) stimrods(flash_inten * pow(fbase,i), flashstart, flashdur);
      		else stim_spot (spotdia, 0, 0, inten=flash_inten * pow(fbase,i), flashstart, flashdur);
		step (flashdur);
    	        if (disp_phot) run_phot(j);   /* display photons in flash */
		step (flashintervals[i]-flashdur);
 
	} /* flash sequence */


	step (stimint-stimstart-8*flashintervals[i]-stimrecover);

     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = 1000;
        apnt->taub = 0.001;
        apnt->tauc = 1000;
        apnt->taud = 0.001;
        apnt->caflg = 1;
        apnt->cai = scai;
     };
     step (stimrecover);
     for (r=0; r<rodnum; r++) {		/* modify chans for fast recovery */
        epnt = modify (rchan[r]);
        apnt = (cattrib*)chanattr(epnt, CGMP);
        apnt->taua = staua;
        apnt->taub = staub;
        apnt->tauc = stauc;
        apnt->taud = staud;
     }
     if (outt) printf ("\n");
   }

}  /* expt== flash_8 */

else

/*-------------------------------------------------*/

if (strcmp(expt,"disp_flash")==0) {		/* display flash */

     int i, t, t2;
     double start, dur;
     double runlen;
     double fl_inten, inten;

   timinc   = 2e-5;
   disp=1;				/* set display mode to save elems */

   if (notinit(spotdia))      spotdia     = 100;
   if (notinit(nflash))       nflash      = 10;
   if (notinit(flash_inten))  flash_inten = 1e3;
   if (notinit(fbase))        fbase       = 1.5;	/* flash step factor */
   if (notinit(flashdur))     flashdur    = 0.0001;
   if (notinit(stimstart))    stimstart   = 0.001; 	/* pre-stimulus int */
   if (notinit(stimint))      stimint     = 0.01; 		/* ISI */
   if (notinit(stimrecover))  stimrecover = .01;        /* recovery time */
   if (notinit(use_inten_list)) use_inten_list = 0;	/* use inten_list[] */
   if (notinit(inten_scale))  inten_scale = 1;		/* scale for inten_list*/

   //endexp = nflash*stimint+stimstart;		/* end of expt */

   runlen = stimstart+stimint;

   if (disp_phot) init_phot(0,0,.9);

   for (i=0; i<nflash; i++){

     simtime = 0;
     if (use_inten_list) fl_inten = inten_list[i]*inten_scale;
     else                fl_inten = flash_inten * pow(fbase,i);
     fprintf (stderr,"expt %s: Flash %d: inten %5.3g\n",expt, i+1, fl_inten);

     step (stimstart);
     if (!pnoise) stimrods(fl_inten, simtime, flashdur);
     else stim_spot (spotdia,0,0, inten=fl_inten, start=simtime,dur=flashdur); 

     step (flashdur);
     if (disp_phot) run_phot(i);	/* display photons */
     step (stimint-stimrecover);
     simwait(2.0);
   }
 
}  /* expt==disp_flash */

else

/*-------------------------------------------------*/

if (strcmp(expt,"single_ph_trace")==0) {	/* single photon in rod and RBC */

    int cn=0;
    double start, dur;
    double inten, flash_inten;

   disp=33;				/* set display mode to save elems */

   if (strcmp(recmode,"vclamp")==0) {
      vclamp (ndn(rbp,cn=0,soma), rbrest, start=0, dur=100);
   }

   if (notinit(spotdia))      spotdia     = 3;
   if (notinit(nflash))       nflash      = 1;
   if (notinit(flash_inten))  flash_inten = 1e3;
   if (notinit(fbase))        fbase       = 1.5;	/* flash step factor */
   if (notinit(flashdur))     flashdur    = .0001;
   if (notinit(stimstart))   stimstart    = .75; 	/* pre-stimulus int */
   if (notinit(stimint))      stimint     = 1.5; 		/* ISI */
   if (notinit(stimrecover))  stimrecover = .01;        /* recovery time */

   flash_inten = 1e3;

   endexp = 3;		/* end of expt */

   disp_plots();				/* make plots */

   for (i=0; i<nflash; i++){

     step (stimstart);
     stim_rod (ndn(xrod,0), inten=flash_inten*0.7,start=simtime,dur=flashdur);
     step (stimstart);
     stim_rod (ndn(xrod,0), inten=flash_inten*1.5,start=simtime,dur=flashdur);

     step (flashdur);
     step (stimint-stimrecover);
   }
 
 }  /* expt==single_ph_trace */
}  /* end of runmod() */

