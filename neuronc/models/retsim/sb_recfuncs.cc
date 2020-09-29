/* module sb_recfuncs.n */
/* Functions for setting up recording sites, plotting, in SB amac cell */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

#define NSBNODES 2000
#define NSB      20
#define NDSGC    10
#define NDSGCNODES 3300
#define SBDSCONNS  5 
#define NDENDREC   6			/* max number of recording sites */
#define MAXNCELREC 6			/* max number of cells to record from */
#define NUMDEN     30			/* number of dendrites, set to max number */
#define BIGNEGVAL -1e10
#define BIGPOSVAL  1e4

#define CT       0
#define CN       1
#define NOD      2
#define VPEAK    3
#define VDIP     4
#define T_VPEAK  5
#define T_VDIP   6
#define IPEAK    7
#define IDIP     8
#define T_IPEAK  9
#define T_IDIP   10
#define CAPEAK   11
#define CADIP    12
#define NRECPARS 13

#define DENFOUND  0
#define RMIN      1
#define RMAX      2
#define MAXNOD    3
#define NDENPARMS 4

#define PEAK  0
#define VOLT  1
#define CALC  2
#define CURR  3
#define SYNIN 4
#define DIP   1

#define NORANGE 0		/* =0 -> no range, proc must find one */

extern const char *expt;
extern double stfreq;

static int sbnodes=0;
static int xnsb=0;
static int xndsgc=0;
static int sb_dsconns;

int dsconninfo=0;
double dsconn_gen; 		/* 0=least, 1=most general conn */
int ndendrec;
int addca;
int DSloc;
int plotsep=0;
int somaclamp;
int sbreccel;
int sbrecnod;
double cadisptime=0;		/* time for [Ca]i display    */

int cbpsb_synsa [NSB][NSBNODES];	/* cbp syns per sb node */
int cbpsb_cbp_cn[NSB][NSBNODES];	/* celnum of presyn cbp, per sb node */
int cbpdsgc_cbp_cn[NDSGC][NDSGCNODES]; 	/* celnum of presyn cbp, per dsgc node */

 /* for each dsgc cel, per nodnum save connected type,cel,nod numbers */

int dsgc_in_syns[NDSGC][NDSGCNODES][3];
int sb_out_syns[NSB][SBDSCONNS];

 /* array holds dir of ds, cellnum of sb, nodenum on sb, dist to gccable */
 /* make last dim one larger than required, to allow overflowing later */

double sb_dsconn_dist[NSB][SBDSCONNS] = {1e6};
int sb_dsconn_nodenrs[NSB][SBDSCONNS] = {-1};
int dsgc_sbconn_nodenrs[NDSGC][SBDSCONNS] = {-1};

double dsgc_pref[NDSGC] = {0};

int recpoints[NCELTYPES][MAXNCELREC][NDENDREC] = {-1}; /* recording sites */
int recpointsR[NCELTYPES][MAXNCELREC][NDENDREC] = {-1};/* sites for computeR*/

double Vresponses[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* max V resp for stim */
double Caresponses[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* max Ca responses */
double iresponses[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGPOSVAL}; /* max current responses */
double Vmaxresp[2][NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* max V resp / direction */
double camaxresp[2][NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* max Ca responses */
double imaxresp[2][NCELTYPES][MAXNCELREC][NDENDREC] = {BIGPOSVAL}; /* max current responses */
double zerovals[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGPOSVAL}; /*  V value at time=0; */
double Vpeaks[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* peaks for all rec points */
double Vdips[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGPOSVAL};
double capeaks[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* [Ca] peaks for all rec pnts */
double cadips[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGPOSVAL};
double ipeaks[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* peaks and dips in current I */
double idips[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGPOSVAL};
double inputpeaks[NCELTYPES][MAXNCELREC][NDENDREC]= {BIGNEGVAL}; /* peaks in synaptic input */
double peakVtimes[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* time at peak voltage */
double peakinputtimes[NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL}; /* time at peak syn input */
double peakVtimesCFCP[2][NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL};
double peakinputtimesCFCP[2][NCELTYPES][MAXNCELREC][NDENDREC] = {BIGNEGVAL};

/* array to hold recording sites + associated values of measured params */

int recpts[NDENDREC][NRECPARS] = {{-1}};

//int R [NUMDEN][NDENDREC];	/* array to hold all R-values computed */

double denparms [NDENPARMS];	/* holds denfound,rmin,rmax,xdmaxnod  */
				/* (see func searchdend) */


/*---------------------------------------------------------------------*/

void sb_recinit(void)

{
 //sbnodes = celnode[sb][1]*2 + 1;	/* rough estimate of sb nodes */
 sbnodes = NSBNODES;			/* estimate of max sb nodes */

 if (!notinit(nsbac))   xnsb = nsbac;
 else                   xnsb = 1;
 if (!notinit(ndsgc)) xndsgc = ndsgc;
 else                 xndsgc = 1;

 sb_dsconns = 5;         /* max nr of sb-dsgc synapses   */
 setptrn("dsconninfo",&dsconninfo);
 setptr ("dsconn_gen",&dsconn_gen);
 setptr ("addca",&addca);
 setptr ("DSloc",&DSloc);
 setptrn("plotsep",&plotsep);
 setptr ("somaclamp",&somaclamp);
 setptr ("sbreccel",&sbreccel);
 setptr ("sbrecnod",&sbrecnod);
 setptrn("cadisptime",&cadisptime);
}

/*---------------------------------------------------------------------*/

void conn_sbdsgc (int prect, int precn, int postct, int postcn)

{
   int n, nfilt;
   int sbnod, dsgcnod;
   int ctindex, cnindex, nindex;
   int rcs;
   int sb_dsconns = SBDSCONNS;
   synapse *spnt;
   nattrib *napnt;

  rcs  = getconn(prect,postct);         /* synapse type prect -> postct */
  if (dsconninfo>=1) {
    fprintf(stderr,"sb_dsconn_nodenrs[%d]:\n", precn);
    for(n=0;n<sb_dsconns;n++)
      fprintf(stderr,"  %d",sb_dsconn_nodenrs[precn][n]);
    fprintf(stderr,"\n");
  };

  for (n=0;n<sb_dsconns;n++) {
    ctindex =0;
    cnindex =1;
    nindex  =2;

    sbnod   = sb_dsconn_nodenrs[precn][n];
    dsgcnod = dsgc_sbconn_nodenrs[postcn][n];

    if(sbnod!=-1 && dsgcnod!=-1) {

      spnt = make_synapse (nd(prect,precn,sbnod),nd(postct,postcn,dsgcnod));
      spnt->ntact  = OPEN;
      spnt->ngain  = gs(SGAIN);
      spnt->maxcond=gs(SCOND);
      spnt->thresh =gs(STHRESH);
      nfilt=(int)gs(SNFILTH);
      spnt->nfilt1h=(short int)nfilt;
      spnt->timec1h=makfiltarr(nfilt,0,spnt->timec1h,gs(SDURH));
      spnt->filt1hg=gs(SHGAIN);
      spnt->vrev   =gs(SVREV);
      spnt->vgain  =1;
      nfilt =(int)gs(SNFILT);
      spnt->nfilt2 =(short int)nfilt;
      spnt->timec2 = makfiltarr(nfilt,0,spnt->timec2,gs(SDUR));
      if (gs(SCNOISE)>0) {
        napnt=make_chnoise(spnt);
        napnt->unitary=22e-12;
      }
      if (gs(SVNOISE)>0) {
        napnt=make_vesnoise(spnt);
        napnt->vsize  =gs(SVSIZ);
      }
      sb_out_syns[precn][n]=spnt->elnum;    /* save unique identifier */
 
      dsgc_in_syns[postcn][dsgcnod][ctindex]=prect;
      dsgc_in_syns[postcn][dsgcnod][cnindex]=precn;
      dsgc_in_syns[postcn][dsgcnod][nindex] =sbnod;

      if (dsconninfo>=1)
         fprintf(stderr,
        "connections made: [%s][%d][%d] to [%d][%d][%d]; n=%d\n",
                cname[prect],precn,sbnod, postct,postcn,dsgcnod,n);
    }
  }  /* for (n;n<sb_dsconns;) */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int sbdsgc_conn(int prect, int precn, int postct, int postcn)

/* Call for 2 cells of known type, connect them according to specificity param. */
/* Prect is starburst, postct is dsgc. */
 
{
    int c, c2, dd, nn;
    int i,n, count;
    double xsoma,ysoma;
    double x1,y1,dx,dy;
    double dsgcpref_deg, sb_denpref, lowbound, upbound;
    double dist, pn;
    char *elemtyp;
    node *npnt;
    elem *epnt,*epnt2;

 if (notinit(dsconn_gen)) dsconn_gen = .5;   /* 0=least, 1=most general conn */

 npnt = ndn(prect,precn,soma);
 xsoma = npnt->xloc;
 ysoma = npnt->yloc;

 for (epnt=elempnt; epnt=foreach(epnt,CABLE,prect,precn,-1,-1,NULL,NULL,&nn,NULL);
		     epnt=epnt->next) {
  c = epnt->elnum;
  elemtyp = get_elabel(epnt,ELABL);
  if(nsbac==0) nn=0;    /* if no starbursts, set nn to 0 */

  /* if sb dendrite is distal, and node has cbp associated with it */

  if (!strcmp(elemtyp,"dend_dist") && !notinit(cbpsb_synsa[precn][nn])) {

      npnt = ndn(prect,precn,nn);
      x1 = npnt->xloc;
      y1 = npnt->yloc;

      dx = x1 - xsoma;
      dy = y1 - ysoma;
      sb_denpref = atanx(dx,dy)*DEG;    /* get sb node pref angle to soma */
      while (sb_denpref>180)            /* limit to range -180 to 180 deg */
        sb_denpref -= 360;

      dsgcpref_deg = dsgc_pref[postcn]; /* get dsgc preferred direction */
      while (dsgcpref_deg>180)          /* limit to range -180 to 180 deg */
        dsgcpref_deg-=360;

      //if (dsgcpref_deg<=45 && dsgcpref_deg>-45) dsgcdir = EAST
      //else if (dsgcpref_deg<=135 && dsgcpref_deg>45) dsgcdir = NORTH
      //else if (dsgcpref_deg<=-135 || dsgcpref_deg>135) dsgcdir = WEST
      //else if (dsgcpref_deg<=-45 && dsgcpref_deg>-135) dsgcdir = SOUTH;

      lowbound = dsgcpref_deg - dsconn_gen*180; /* compute angle range */
      upbound  = dsgcpref_deg + dsconn_gen*180;

     if (lowbound<=sb_denpref && sb_denpref<=upbound) {
       for (epnt2=elempnt; epnt2=foreach(epnt2,CABLE,postct,postcn,-1,-1,
						NULL,NULL,&dd,NULL,
						10,dist2d,ndt(prect,precn,nn));
			   epnt2=epnt2->next) {
	epnt2->elnum = c2;
        dist=endist2d(c2, ndn(prect,precn,nn));

        if (dsconninfo>=2)
 	  fprintf(stderr,
	"  sbdsgc_conn: node [%s][%d][%d]: dist=%g sbpref=%g dsgcpref=%g postcn=%d\n",
                        cname[prect],precn,nn,dist,sb_denpref,dsgcpref_deg,postcn);

    /* keep the closest distances in an ordered array: */
    /* shortest distance,closest node first */

        for (n=0;n<sb_dsconns;n++) {
          if (dist<sb_dsconn_dist[precn][n]) {
            for (i=sb_dsconns;i>n;i--) {          /* shift everything to right */
              sb_dsconn_dist[precn][i]=sb_dsconn_dist[precn][i-1];
              sb_dsconn_nodenrs[precn][i]=sb_dsconn_nodenrs[precn][i-1];
              dsgc_sbconn_nodenrs[postcn][i]=dsgc_sbconn_nodenrs[postcn][i-1];
            };
                                                      /* add in new values: */
            sb_dsconn_dist[precn][n]=dist;              /* dist */
            sb_dsconn_nodenrs[precn][n]=nn;             /* sb nodenr */
            dsgc_sbconn_nodenrs[postcn][n]=dd;  /* dsgc nodenr */

            count=0;
            for (i=0;i<sb_dsconns;i++) {                /* check for doubles */
              if (sb_dsconn_nodenrs[precn][i]==nn)
                count++;
                if (count>1) {   /* if found, erase by shifting array to left */

                 sb_dsconn_dist[precn][i]=sb_dsconn_dist[precn][i+1];
                 sb_dsconn_nodenrs[precn][i]=sb_dsconn_nodenrs[precn][i+1];
                 dsgc_sbconn_nodenrs[postcn][i]=dsgc_sbconn_nodenrs[postcn][i+1];
                };
             };

             if (dsconninfo>=2) {                       /* print out arrays */
                fprintf(stderr, "  array sb_dsconn_nodenrs:  \n");
                for (n=0;n<sb_dsconns;n++)
                  fprintf(stderr, "  %d", sb_dsconn_nodenrs[precn][n]);
                fprintf(stderr, "\n  array sb_dsconn_dist:  \n");
                for (n=0;n<sb_dsconns;n++)
                  fprintf(stderr, "  %-4.3g", sb_dsconn_dist[precn][n]);
                fprintf(stderr, "\n    sbdsgc_conn: sb[%d][%d] and dsgc[%d][%d]; dist=%g\n",
                                     precn,nn,postcn,dd,dist);
              };
              break;
            }
          }
        }
      }
    }
  }
  return 1;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void save_cbp_syns (int prect, int precn, int postct, int postcn, int hnod, int synout)

{
 if (postct==sbac) {
   if(nsbac==0) hnod=0;		     /* for each sb node, save cbp synapse num */
   if (postcn > NSB) postcn = NSB;
   if (hnod > NSBNODES) hnod = NSBNODES;
   cbpsb_synsa[postcn][hnod] = synout; /* save cbp celnr per sb node */
   cbpsb_cbp_cn[postcn][hnod] = precn; /* save cbp celnr per dsgc node */
 }
 else if (postct==dsgc) {
   cbpdsgc_cbp_cn[postcn][hnod] = precn; 
 };

 /* fprintf(stderr,"cbpdsgc_cbp_cn[%g][%g] = %g \n",postcn,hnod,precn); /* */
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void print_sb_out_syns()

{
    int i,n;

  fprintf(stderr,"retsim1: sb_out_syns:\n");
  for (i=1; i<=nsbac; i++) {
    for (n=0;n<sb_dsconns;n++)
      fprintf(stderr," %d",sb_out_syns[i][n]);
    fprintf(stderr,"\n");
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void display_sb_out_syns(int i)

{
    int n;
    int color;
    double dscale, z1, z2;
 
  for (n=0;n<sb_dsconns;n++) {
    if(!notinit(sb_out_syns[i][n])) {
      if (info_disp>=1)
        fprintf(stderr,
            "retsim1: make_sb cel %d: sb_dsconns=%d n=%d synapse=%d\n",
                                    i, sb_dsconns,n,sb_out_syns[i][n]);

       display (ELEMENT, sb_out_syns[i][n], color=white, dscale=4);

    }
    else
      break;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int searchdend(int ct, int cn, double anglrad, double dt)
		/* dt (searchrange), anglrad both in RAD-units! */
{
    int denfound1,cnod,cbrnum,xdmaxnod;
    double xsoma,ysoma;
    double xden,yden,dx,dy;
    double r,rmax,rmin;
    double anglmax,anglmin,dangl,d1,d2,d3;
    elem *epnt;
    node *npnt;

//fprintf(stderr,"searchdend %g\n", dt);
  npnt = ndn(ct,cn,soma);
  xsoma = npnt->xloc;		/* xsoma is the x-loc of node 0 (soma) */
  ysoma = npnt->yloc;

//fprintf(stderr,"searchdend ct%g cn%g\n", ct,cn);
  /* First look at nodes within most exact range of angles possible: */
  rmax = -1e10;
  rmin =  1e10;

  anglmax = anglrad + dt;	//in radians
  anglmin = anglrad - dt;

  for (epnt=elempnt; epnt=foreach(epnt,CABLE,ct,cn,-1,-1,NULL,NULL,&cnod,NULL);
			epnt=epnt->next) {
    cbrnum = epnt->elnum;
    npnt = ndn(ct,cn,cnod);
    xden  = npnt->xloc;
    yden  = npnt->yloc;
    dx = xden - xsoma; 		/* dist between den-node and soma */
    dy = yden - ysoma; 
    dangl = atanx(dx,dy);

    /* atanx gives angle in RAD: dy/dx is tan(angle), */
    /* so atan(dy/dx) gives that angle */

    r = sqrt(dx*dx+dy*dy);
    d1=dangl-anglmax;
    d2=dangl-anglmin;

    if (inrange(anglmin,anglmax,dangl)) {  /* dend must be in angle range */
      if (rmin > r)
        rmin = r; 		/* if r is smaller than rmin, change rmin */
      if (rmax < r) {
        rmax = r;		/* if r is bigger than rmax, change rmax */
        xdmaxnod = cnod;	/* that node then becomes the farthest node */
      }
    }
  }
  denfound1 = (rmax != -1);	/* denfound=1 when a dendrite is found */
  denparms[DENFOUND]=denfound1;
  denparms[RMIN]=rmin;
  denparms[RMAX]=rmax;
  denparms[MAXNOD]=xdmaxnod;
  return denfound1;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int findrecpts (int ct, int cn, double angl, int searchrange, int ndendrc)

 /* Finds recording points and puts nodnrs in array "recpts".
    If it finds a SBAC dendrite within a certain range of angles,
       1) Get first node within radius range,
       2) Try to get node that receives connection from bipolar
           (for bp synapse recording).
        3) Assumes soma is node 0, and that its primary dendrites
            start at node 1.
 */
{
    int i, n, q, color, nborders;
    int cbrnum, cnod, nod1, nod2, nconn, nelem;
    int denfound, xdmaxnod, dend;
    double r, rmax, rmin, dscale;
    double xden, yden, xsoma, ysoma;
    double dx,dy,dt,dangl;
    double anglrad,anglmax, anglmin,dtexact;
    double recsitewidth, whigh, wlow, temp_angle;
    int nfound[NDENDREC]={0};	/* if recsite x found, nfound[x]=1 */
    node *npnt;
    elem *epnt;
    electrode *rp;

  anglrad = angl/DEG;
  dtexact = 2*PI/30;

  ndendrec = ndendrc;			/* set global variable */

  if (ndendrec > NDENDREC) ndendrec = NDENDREC;

  denfound=searchdend(ct,cn,anglrad,dtexact);/* look for dend at exact angle */

  //fprintf(stderr,"denfound: %d ct=%d (%s)\n", denfound, ct, cname[ct]);

  if (!denfound) {	/* if no dendrite found in the exact search */
    if(searchrange == NORANGE)	{	/* if no specific searchrange entered */
      denfound=searchdend(ct,cn,anglrad,2*PI/10);  /* expand search to 2PI/10 */
      if (!denfound) denfound=searchdend(ct,cn,anglrad,2*PI);/* expand to 2PI */
    }
    else {		/* if specific range entered, search only that range */
      dt = .5 * searchrange;
      denfound=searchdend(ct,cn,anglrad,dt);
    }
  }

  /* If a dendrite was found, follow through parent nodes to soma. */
  /* First initialise bounds for placing the recording sites, */

  if (denfound) {
    rmin=denparms[RMIN];
    rmax=denparms[RMAX];
    xdmaxnod=int(denparms[MAXNOD]);
    npnt = ndn(ct,cn,0);
    xsoma = npnt->xloc;
    ysoma = npnt->yloc;
    npnt = ndn(ct,cn,xdmaxnod);
    xden  = npnt->xloc;
    yden  = npnt->yloc;
    dx = xden - xsoma;
    dy = yden - ysoma;

    /* Interval defined as:  space between sites / num of intervals needed */
    /* Because soma is always recsite, make max ndendrec-1 recsites on dend */

    recsitewidth= 1 / (ndendrec-1) * 0.3;/* width of recsite as fraction of dend length*/

    /* fprintf(stderr,"recsitewidth=%g\n", recsitewidth); /* */
    /* fprintf(stderr,"denfound=%g; rmin=%g; rmax=%g; xdmaxnod=%g\n", 
		denfound,rmin,rmax,xdmaxnod); /* */

    /* find angle (in deg) of farthest node and remember it in global var */

    temp_angle=atanx(dx,dy)*DEG;

    recpoints[ct][cn][0]=soma;		/* make soma first rec-point! */
    nfound[0]=1;			/* so first rec-point is always found */

    recpoints[ct][cn][ndendrec-1]=xdmaxnod;	/* make tip last rec-point! */
    nfound[ndendrec-1]=1;

    for (n=xdmaxnod; n != soma; n) {
      npnt = ndn(ct,cn,n);
      xden  = npnt->xloc;
      yden  = npnt->yloc;
      dx = xden - xsoma;
      dy = yden - ysoma;
      r = sqrt(dx*dx+dy*dy);

      for (q=1;q<ndendrec;q++) {
        whigh = (q/ndendrec + recsitewidth)*rmax;
        wlow  = (q/ndendrec - recsitewidth)*rmax;
        if (ct==sbac) {

  if (!notinit(cbpsb_cbp_cn[cn][n]))
 /* fprintf(stderr, "cbpsb_cbp_cn[%g][%g]= %g lowbound=%g,r=%g,hibound=%g,nfound[%g]=%g\n",
       cn,n,cbpsb_cbp_cn[cn][n],wlow,r,whigh,q,nfound[q]);
/* */
	  if (!notinit(cbpsb_cbp_cn[cn][n]) && r>wlow && 
			r<=whigh && !nfound[q]) {

/*  fprintf(stderr, "cbpsb_cbp_cn[%g][%g]= %g lowbound=%g,r=%g,hibound=%g,nfound[%g]=%g\n",
       cn,n,cbpsb_cbp_cn[cn][n],wlow,r,whigh,q,nfound[q]);
/* */
	    recpoints[ct][cn][q] = n;
	    nfound[q] = 1;
	  };
	}
	else if (ct==dsgc) {
	  if (!notinit(cbpdsgc_cbp_cn[cn][n])) {

/* fprintf(stderr, "cbpdsgc_cbp_cn[%g][%g]= %g lowbound=%g,r=%g,hibound=%g,nfound[%g]=%g\n",
       cn,n,cbpdsgc_cbp_cn[cn][n],wlow,r,whigh,q,nfound[q]);
/* */

	    if (r>wlow && r<=whigh+0.001 && !nfound[q]) {
	      recpoints[ct][cn][q] = n;
	      nfound[q] = 1;

	   /* fprintf(stderr, "nfound[%d] = %d: recpoints[dsgc][%d][%d]=%d \n",
                  q, nfound[q],cn,q,recpoints[ct][cn][q]); /* */

            }
	  }
	}
      }

      nconn = int(get_nfield(ndn(ct,cn,n),NUMCONN));
      //fprintf (stderr,"node = %d; nconn=%d\n",n,nconn);

      for (i=1; i<=nconn; i++) {  /* Find path back to soma */
        nelem = int(get_nfield(ndn(ct,cn,n),ELNUM,i));
        epnt = get_elempnt(nelem);

        if (epnt->ctype != CABLE) continue;

	if (ninfo>2)
	  printf ("tracing from node %d conn %d\n",n,i);

	nod1=epnt->node1c;	//nr of first node
        nod2=epnt->node2c;	//nr of second node
	//fprintf (stderr,"nod1 %d nod2 %d n %d r %g\n",nod1,nod2,n,r);
        if (nod2 == n)
	  continue;        	/* descendent node */
        if (nod1 == n) {
	  n = nod2;        	/* follow parent node */
 	  break;
        }
      }  /* for (i;;) */

    //fprintf (stderr,"updated node = %g\n",n);
    } /* for (n;;)   */
  } /* if(denfound) */

  if(ninfo >= 2) {
    for (q=0;q<ndendrec;q++) {
        if (recpoints[ct][cn][q] >= 0)
           rp = (electrode *)at(nd(ct,cn,recpoints[ct][cn][q]), ELECTRODE); rp->dia=10;
	display (ELECTRODE,MATCHING,ndt(ct,cn,recpoints[ct][cn][q]), color=q+1,dscale=1);
    };
    //if (disp) exit(0);
  }
  if (ninfo >= 3) {
    printf ("# sb_dend %d recording nodes:\n",dend=1); // fix this
    for (q=0;q<ndendrec;q++)
      printf ("%d; ", recpoints[ct][cn][q]);
  }

  if(notinit(addca)) addca = 0;
  addca = int(getn(sbac,BIOPHYS));

  //if (addca==2) //add ca-channel at one node only
  //  add_oneCachan();

  return denfound;

} /* findrecpts */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int sb_distnode(int cn)

/* Finds distal node in certain sb, using array recpoints,
   which, if made for sb, is ordered so that last node in array is most distal.
   Alternatively, user can specify node using DSloc (0=soma,1=next one,etc.).
*/

{
    int recnod,cellnr;

  if (notinit(DSloc))	{	/* unless node specified, get outermost node */
     for (recnod=ndendrec-1;recnod>=0;recnod--)
        {  if (recpoints[sbac][cn][recnod]!=-1)  /* take the outermost rec-node */
   	   { if (addca)		   /* only take node w/Ca-compartment */
	       { if (get_nfield(ndn(sbac,cn,recpoints[sbac][cn][recnod]),CACOMP))
                 break;
	       }
	     else
	       break;
	   }
        }
  }
  else if (addca) {	/* if addca, only accept DSloc node with Ca-comp */
    if (get_nfield(ndn(sbac,cn,recpoints[sbac][cn][DSloc]),CACOMP))
       recnod=DSloc;
     else {
	fprintf(stderr,
	   "# Func findnode error; for DSloc please select node with Ca-comp.");
        exit(2);
       }
  }
  else if (DSloc>=ndendrec || DSloc<0 || recpoints[sbac][cn][DSloc]== -1) {
    fprintf (stderr,"# The node entered is not a valid rec-point.\n");
    fprintf (stderr,"# Please check func 'findnode()'.\n");
    fprintf (stderr,"# With DSloc enter range of 0 (=soma) - %d (distal)\n",
		ndendrec-1);
    exit(2);
  }
  else
    recnod=DSloc;

  return recnod;
}

/*------------------------ plot -----------------------------*/

void plot_cbpsb_v(int cn, int n, double vmin, double vmax, double siz)

/* display bipolar cell voltages presynaptic to sb number cn */

{
  int first;
  double plg, offtr, offb, pmax, pmin, plsize;
  int rec, pen, plnum;
  char *plname;
  
  plg   = vmax-vmin;	/* gain of trace */
  offtr = .0;		/* position of trace within graph */
  offb  = vmin;		/* trace offset base */
  if (siz<0) siz = 0.5;

  if (make_dbp1_sbac && script==0) {
    if (n==-1) {	/* if node = -1, show all bipolars connected to recpoints */
      for (rec=0;rec<ndendrec;rec++) {
        if (!notinit(cbpsb_synsa[cn][recpoints[sbac][cn][rec]])) {
          sprintf (plname,"Vcbp_%d",recpoints[sbac][cn][rec]);
	  plot (V,ndn(dbp1,cbpsb_cbp_cn[cn][recpoints[sbac][cn][rec]], soma), 
			pmax=(1-offtr)*plg+offb,
	    	        pmin=(0-offtr)*plg+offb);
	  plot_param(plname, pen=rec+1, plnum=ndendrec+3, plsize=1);
        }
      }  /* for (rec;;) */
    }
    else if (!notinit(cbpsb_cbp_cn[cn][n])) {

      /* make string with input bipolar and receiving sb node */

      sprintf(plname,"Vcbp%d->sb%d_%d", cbpsb_cbp_cn[cn][n], cn, n);

      plot (V, ndn(dbp1,cbpsb_cbp_cn[cn][n],soma), 
			pmax=(1-offtr)*plg+offb,
		        pmin=(0-offtr)*plg+offb);
      plot_param(plname, plnum=ndendrec+3, plsize=1);
    }
  }  /* if (make_cbpsb) */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_cbpsb_rate(int cn, int n, double maxvesrate, double siz)

/* display bipolar cell vesicle release rate at synapses presynaptic to sb number cn */

{
  int first, rec;
  int plnum, pen;
  double plg, offtr, offb;
  char *plname;
  double plsize;
  double pmax, pmin;
 
  plg   = maxvesrate;		/* gain of trace */
  offtr = .0;		/* position of trace within graph */
  offb  = 0;		/* trace offset base */
  if (siz<0) siz = 0.5;

  if (make_dbp1_sbac && script==0) {
    if (n==-1) {	/* if node = -1, show all bipolars connected to recpoints */
      for (rec=0;rec<ndendrec;rec++) {
        if (!notinit(cbpsb_synsa[cn][recpoints[sbac][cn][rec]])) {
          sprintf (plname,"Rcbp_%d",recpoints[sbac][cn][rec]);
	  plot (FA9,cbpsb_synsa[cn][recpoints[sbac][cn][rec]],
			pmax=(1-offtr)*plg+offb,
	    	        pmin=(0-offtr)*plg+offb);
	  plot_param(plname, pen=rec+1, plnum=ndendrec+2,plsize=siz);
        }
      }  /* for (rec;;) */
    }
    else if (!notinit(cbpsb_synsa[cn][n])) {

      /* make string with input bipolar and receiving sb node */

      sprintf(plname,"cbp%d->sb%d_%d", cbpsb_cbp_cn[cn][n], cn, n);

      plot (FA9,cbpsb_synsa[cn][n],
			pmax=(1-offtr)*plg+offb,
		        pmin=(0-offtr)*plg+offb);
      plot_param(plname,plnum=ndendrec+2,plsize=1);
    }
  }  /* if (make_cbpsb) */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_v_recpnt(int ct, int cn, double vmin, double vmax, double size)

/* Display voltage at dendritic sites in array recpoints */

{
   static int plotnum=0;
   int i;
   int recpt, pen, plnum;
   double plg, offtr, offb, pmax, pmin;
   double plsize;
   char *plname;

  if (script==0) {
    if(ct==dsgc || ct==gca || ct==gcb || ct==gcaoff || ct==gcboff) /* plotnum sets where to put the plot */
				/* (lowest num = bottom) */
      plotnum+=1;		/* put dsgc plot at bottom of display */
    else if(ct==sbac)
      plotnum+=ndendrec+1;
    else if(ct==dbp1 || ct==dbp2 || ct==rbp)
      plotnum+=ndendrec+3;
    else if (ct==xcone || ct ==xrod)
      plotnum+=ndendrec+4;
    else plotnum+=ndendrec+5;
    if (size<0) size = 1;

    for (i=0;i<ndendrec;i++) {
      if (recpoints[ct][cn][i]>=0) {
        recpt=recpoints[ct][cn][i];
        sprintf(plname,"V%s_%d_%d",cname[ct],cn,recpt);
        plot (V, ndn(ct,cn,recpt), pmax=vmax, pmin=vmin);
	plot_param(plname, pen=i+1, plnum=plotnum,plsize=size);
	if(plotsep) plotnum++;	/* change number so each node has own plot */
      };
    };
    plg = 200;		/* gain of trace */
    offtr = .0;		/* position of trace within graph */
    offb  = 0;		/* trace offset base */
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


void defaultplots(double vmin, double vmax, double maxvesrate, double maxCa)

/* Sets the standard plots for most experiments */

{
   int ct,cn,n;
   double siz;

  if (script==0 && !disp) {     //don't plot if running higher level script or -d
     plot_cbpsb_rate(1,-1,maxvesrate,0.3);
     if (!notinit(somaclamp)) plot_i_nod(ct=sbac,cn=1,soma,-1e-10,1e-10,magenta,"",-1,siz=1);
     plot_ca_nod(ct=sbac,cn=1,n=46,maxCa,cyan,"",-1,siz=0.2);
     plot_v_recpnt(sbac,1,-0.08,-0.02,siz=0.5);
     //plot_v_nod(ct=sb,cn=1,soma,-0.08,-0.02,red,"Vsb",-1,siz=1);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void getresponses(int ct,int cn, int fugorpet, char *stimtype)
{
  int n,recnod;

  for (n=0;n<ndendrec;n++) {
    recnod=recpoints[ct][cn][n];
    if (recnod>=0) {
      if (!strcmp(stimtype,"sine") || !strcmp(stimtype,"sineann")) {
          if (notinit(somaclamp))
            Vresponses[ct][cn][n] = Vpeaks[ct][cn][n] - Vdips[ct][cn][n];
          else
            iresponses[ct][cn][n] = ipeaks[ct][cn][n] - idips[ct][cn][n];
      }
      else
        Vresponses[ct][cn][n] = Vpeaks[ct][cn][n]-zerovals[ct][cn][n];

    // fprintf(stderr, "Vresponses[%g][%g][%g] = %g\n",ct,cn,n,Vresponses[ct][cn][n]);

      if (addca && int(get_nfield(ndn(ct,cn,recnod),CACOMP))) {
        if (!strcmp(stimtype,"sine") || !strcmp(stimtype,"sineann"))
           Caresponses[ct][cn][n] = capeaks[ct][cn][n]-cadips[ct][cn][n];
        else
            Caresponses[ct][cn][n] = capeaks[ct][cn][n] - dtcai;
      }
      else
        Caresponses[ct][cn][n] = -1;

      Vmaxresp[fugorpet][ct][cn][n] = Vresponses[ct][cn][n];
      peakVtimesCFCP[fugorpet][ct][cn][n] = peakVtimes[ct][cn][n];
      peakinputtimesCFCP[fugorpet][ct][cn][n] = peakinputtimes[ct][cn][n];
      camaxresp[fugorpet][ct][cn][n] = Caresponses[ct][cn][n];
      imaxresp[fugorpet][ct][cn][n] = iresponses[ct][cn][n];
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* find sb cel and node to record from that connects to dsgc */

void find_sb_connect(void)

{
   int n;

  if (make_dsgc && make_sbac) {
    if (notinit(sbreccel)) {
      for (n=1;n<=nsbac;n++)
        if (sb_dsconn_nodenrs[n][0]>=0) { /* if elem 0 is neg, array is empty */
          sbreccel = n;                   /* if arr nonempty, use this cell */
          break;
        }
    }
    if (notinit(sbrecnod)) {
      if (sb_dsconn_nodenrs[sbreccel][0]>=0)	 /* if elem 0 is neg, array is empty */
        sbrecnod=sb_dsconn_nodenrs[sbreccel][0]; /* just take elem 0 */
      else
       fprintf(stderr,"# retsim: sbreccel has no nodes connecting to the dsgc");
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void savecones(void)

/* save all cone voltages for later */

{
   elem *epnt;
   photorec *ppnt;

  for (epnt=elempnt; epnt=foreach(epnt,CONE); epnt=epnt->next) {
   if (!(epnt->modif))
     ppnt = (photorec*)makelem(CONE,modify (epnt->elnum));
     ppnt->save = 1;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void restorcones(void)

/* restore all cone voltages */

{
   elem *epnt;
   photorec *ppnt;

  for (epnt=elempnt; epnt=foreach(epnt,CONE); epnt=epnt->next) {
   if (!(epnt->modif))
     ppnt = (photorec*)makelem(CONE,modify (epnt->elnum));
     ppnt->restore = 1;
  }
}

/*-------------------------------------------------------------------*/

void getextr(int ct, int cn, int ptype, double t1, double t2, int maxormin)

/* Find max or min V or [Ca] between two times at a node. */

{
    int recnode, n;
    node *npnt;
    static double farval[NDENDREC]     ={BIGNEGVAL};
    static double farcaval[NDENDREC]   ={0};
    static double farival[NDENDREC]    ={BIGNEGVAL};
    static double farsyninval[NDENDREC]={BIGNEGVAL};

  if (simtime<=0) {	/* (re)set all values before an experiment begins */
    for (n=0; n<NDENDREC; n++) {
      capeaks[ct][cn][n] = BIGNEGVAL;
      cadips[ct][cn][n] = BIGPOSVAL;
      Vpeaks[ct][cn][n] = BIGNEGVAL;
      Vdips[ct][cn][n]  = BIGPOSVAL;
      ipeaks[ct][cn][n] = BIGNEGVAL;
      idips[ct][cn][n]  = BIGPOSVAL;
      inputpeaks[ct][cn][n] = BIGNEGVAL;
      peakVtimes[ct][cn][n] = BIGNEGVAL;
      peakinputtimes[ct][cn][n] = BIGNEGVAL;
    };
  };

  if (simtime>0 && simtime<=ploti) {     /* store V at time=0, for every rec node */
    for (n=0; n<NDENDREC; n++)
      if (recpoints[ct][cn][n]>=0) 
		zerovals[ct][cn][n] = v(ndn(ct,cn,recpoints[ct][cn][n]));
  }

  if (simtime>t1 && simtime<t2) {
    for (n=0; n<ndendrec; n++)	{ /* get voltage extrema for all rec points */
      if (recpoints[ct][cn][n]>=0) {
        npnt = ndn(ct,cn,recpoints[ct][cn][n]);
        if (ptype==VOLT) {	  	 /* assign which plot farval looks at */
	  farval[n] = v(npnt);

	  if (maxormin==PEAK) {  /*  find either voltage peak or dip: */
	    if (farval[n] > Vpeaks[ct][cn][n]) {
	      Vpeaks[ct][cn][n] = farval[n];
              peakVtimes[ct][cn][n] = simtime;
	    };
	  }
	  else if (maxormin==DIP) {
	    if (farval[n] < Vdips[ct][cn][n]) Vdips[ct][cn][n] = farval[n];
	  }
	  else
	    fprintf (stderr,"#6th arg of func getextr must be PEAK or DIP\n");
        }
        else if (ptype==CALC) {		/* check for ca compart. */
          if (get_nfield(npnt,CACOMP)) { 
            farcaval[n] = record_ca(npnt,1);

	    if (maxormin==PEAK) {		/* find [Ca] peak or dip */
	      if (farcaval[n] > capeaks[ct][cn][n]) 
		capeaks[ct][cn][n] = farcaval[n];
            }
	    else if (maxormin==DIP) {
	      if (farcaval[n]<cadips[ct][cn][n]) cadips[ct][cn][n]=farcaval[n];
	    }
	    else
	      fprintf (stderr,"#6th arg of func getextr must be PEAK or DIP\n");

	  }	/* if Ca-compartment */
        }

        else if (ptype==CURR) {		/* get extrema for currents */
	  farival[n] = i(npnt);

	  if (maxormin==PEAK) {		/*  find voltage peak or dip */
	    if (farival[n] > ipeaks[ct][cn][n]) ipeaks[ct][cn][n] = farival[n];
	  }
	  else if (maxormin==DIP) {
	    if (farival[n] < idips[ct][cn][n]) idips[ct][cn][n] = farival[n];
 	  }
	  else
	    fprintf (stderr,"#6th arg of func getextr must be PEAK or DIP\n");
        };

        if (ptype==SYNIN) {	/* get extrema for synaptic inputs */
          if (!notinit(cbpsb_synsa[cn][recpoints[ct][cn][n]])) {
            farsyninval[n] = record_synapse(cbpsb_synsa[cn][recpoints[ct][cn][n]],FA9);

	    if (maxormin==PEAK) {
	      if (farsyninval[n] > inputpeaks[ct][cn][n]) {
                inputpeaks[ct][cn][n] = farsyninval[n];
                peakinputtimes[ct][cn][n] = simtime;
	      }
	    }
	  } /* if cbpsb_synsa is initialised */
        } /* ptype SYNIN */
      } /* if recpoints exist */
    } /* forloop all recpoints */
  } /* if within time window */
} /* proc getextr */

/*----------------------------------------------------*/

void sb_onplot(void)

{
     int ct,cn,recnod, n;
     double t1, t2, peakCaval, resetdur;

   ct = sbac; cn = 1;
   peakCaval = BIGNEGVAL;	/* [Ca] peak at certain time (cadisptime); */
   if (addca && disp && simtime >  cadisptime-0.0001 && simtime <= cadisptime)
   {
      //recnod=findnode();

      /* get [Ca] at cadisptime */

      peakCaval=5e-7;  //Ca(1)[ct][cn][recpoints[ct][cn][recnod]];
      //disp_ca(peakCaval);	/* show anatomy with [Ca]i colored */
   };

    /* times in in which to look at plot: */
   t1 = .015; 			/* start looking a bit after start of simulation */
   t2 = endexp;

    /* for sine stims, find min and max responses */
    if (!strcmp(expt,"sb_cc_sine") || !strcmp(expt,"sb_cc_sineann")) {

     t1 = endexp - 1/stfreq;    	/* set t1 to look at last cycle only */

     getextr(ct,cn,VOLT,t1,t2,PEAK);
     getextr(ct,cn,VOLT,t1,t2,DIP);
     getextr(ct,cn,SYNIN,t1,t2,PEAK);

     if (!notinit(somaclamp)) {
       getextr(ct,cn,CURR,t1,t2,PEAK);
       getextr(ct,cn,CURR,t1,t2,DIP);
     }
     if (addca) {
       getextr(ct,cn,CALC,t1,t2,PEAK);
       getextr(ct,cn,CALC,t1,t2,DIP);
     }
   }
   else /* for all other stims, find only peak responses */
      getextr(ct,cn,VOLT,t1,t2,PEAK);

   if (addca)					/* find peak Ca level */
      getextr(ct,cn,CALC,t1,t2,PEAK);

//fprintf(stderr,"time=%g\n",time);

}

