/* module sb_makfuncs.cc */

/* functions for generating amacrine cells */
/* for use with retsim1.n */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

double sbrangles;
double ang_brpt1;
double ang_brpt2;
double ang_brpt3;
double ang_brpt4;
double ang_brpt5;

double branchdist;
double branchdist1;
double branchdist2;
double branchdist3;
double branchdist4;
double branchdist5;

int nbranchpoints;

double xradius;
double varicos;
double sbac_dend_dia;

double am_seglen;
double am_dia_prox_rad;
double am_dia_prox_factor;
double amdia;

int am_den_seg;
int n_amseg;
int makesbfile;
int sbac_info;
int rand_amth_sb;
int sbac_dend;
double sbac_dend_cplam;

double denddiataper;
int nonrand;

#define NBRPT 5
double arr_sbrangles[NBRPT] = {0};

#define NBRDIST 5
double arr_branchdist[NBRDIST] = {0};

void mak_real_cel(int ct, int cn, double xpos, double ypos, double thetax, 
		double thetay, double thetaz, int flip, 
		int dend1, int dend2, int dend3, int dend4, int dend5, double cplam, int taper);

/* ---------------------------------------------------------------*/

void sb_init(void)

{
  int i;
  static int runyet = 0;

  if (runyet) return;
  runyet = 1;

  setptr("ang_brpt1",&ang_brpt1);
  setptr("ang_brpt2",&ang_brpt2);
  setptr("ang_brpt3",&ang_brpt3);
  setptr("ang_brpt4",&ang_brpt4);
  setptr("ang_brpt5",&ang_brpt5);

  setptr("branchdist1",&branchdist1);
  setptr("branchdist2",&branchdist2);
  setptr("branchdist3",&branchdist3);
  setptr("branchdist4",&branchdist4);
  setptr("branchdist5",&branchdist5);

  setptr("xradius",&xradius);
  setptr("varicos",&varicos);
  setptr("sbac_dend_dia",&sbac_dend_dia);
  setptr("denddiattaper",&denddiataper);
  setptr("am_dia_prox_rad",&am_dia_prox_rad);
  setptr("am_dia_prox_factor",&am_dia_prox_factor);
  setptr("am_den_seg",&am_den_seg);
  setptr("am_seglen", &am_seglen);
  setptr("n_amseg",&n_amseg);
  setptr("amdia",&amdia);
  setptr("makesbfile",&makesbfile);
  setptr("sbac_info",&sbac_info);
  setptr("rand_amth_sb",&rand_amth_sb);
  setptr("sbac_dend_cplam",&sbac_dend_cplam);

  if (notinit(sbac_dend_dia))  sbac_dend_dia = 0.5; /* dendritic cable diameter */
  if (notinit(xradius)) xradius=150;	      /* set dendrite radius */
  if (notinit(varicos)) varicos = 0; 	      /* default no varicosities */

  //branching-angles for each branchpoint: set one value for all
  if (notinit(sbrangles)) sbrangles = 15;

  //or set different angle for each branchpoint
  if (notinit(ang_brpt1)) ang_brpt1=sbrangles;
  if (notinit(ang_brpt2)) ang_brpt2=sbrangles;
  if (notinit(ang_brpt3)) ang_brpt3=sbrangles;
  if (notinit(ang_brpt4)) ang_brpt4=sbrangles;
  if (notinit(ang_brpt5)) ang_brpt5=sbrangles;

  i=0;
  arr_sbrangles[i++] = ang_brpt1; 
  arr_sbrangles[i++] = ang_brpt2; 
  arr_sbrangles[i++] = ang_brpt3; 
  arr_sbrangles[i++] = ang_brpt4; 
  arr_sbrangles[i++] = ang_brpt5; 

 //set distance between branchpoints [NOT distance of branchpoints from soma]

  if (notinit(branchdist)) branchdist = xradius/2.5; //-> 2brpts per den

  //allow setting of individual branchdistances by hand
  if (!notinit(branchdist1)) arr_branchdist[0] = branchdist1;
  if (!notinit(branchdist2)) arr_branchdist[1] = branchdist2;
  if (!notinit(branchdist3)) arr_branchdist[2] = branchdist3;
  if (!notinit(branchdist4)) arr_branchdist[3] = branchdist4;
  if (!notinit(branchdist5)) arr_branchdist[4] = branchdist5;

  if (notinit(nbranchpoints)) nbranchpoints = 2;
  if (nbranchpoints>=NBRDIST) nbranchpoints = NBRDIST-1;
  arr_branchdist[nbranchpoints]=xradius; //make last branchdistance very long
}

/*--------------------------------------------------------*/

int extend_branch(int brnum, int obrnum, int ctype,int cellnum,
			double xsoma, double ysoma, double sdia,
			double nonbrlen, double rad, double seglen,
                        double branchthresh, double growthresh) 

  /* Extend an existing branch on a starburst cell dendrite. */
  /* Includes competition between branches of same cell type. */

{
    int i, nconn,nelem;
    int cnod, cbrnum;
    int n1a, n1b, n1c, n1d;
    int n2a, n2b, n2c, n2d, csn;
    int done, pvals;
    double dx,dy;
    double leng, dist, nsdia;
    double mdist, mdist1, mdist2, mdist3;
    double slen, sangl;
    double angl1, angl2, angl3;
    int nsn, nsn2, nsn3;
    double xm1, ym1, zm, xm2, ym2, xm3, ym3;
    double xden, yden, zden, xden2,yden2;
    double bt,bt2,gt;
    double radfrac,radtaper,currad,dangl;
    double newrad, newrad2, newrad3;
    double pvaricos, slope, taper;
    double cplam;
    elem *epnt;
    node *npnt,*npnt2;

  if (notinit(denddiataper)) denddiataper=1; //default do dend diam tapering
  cplam = -1;

  pvals = 0;
  if (ninfo>3)
    pvals = 1;
  if (ninfo>4)
    pvals = 2;
  if (pvals)
    printf ("entering extend_branch %d %d l %g\n", brnum,obrnum,nonbrlen);

  /* find branches at distal end of segment */

  epnt = get_elempnt(brnum);
  n1a = epnt->node1a;		/* descendent node */
  n1b = epnt->node1b;
  n1c = epnt->node1c;
  n1d = epnt->node1d;

  n2a = epnt->node2a;		/* parent node */
  n2b = epnt->node2b;
  n2b = epnt->node2c;
  n2b = epnt->node2d;

  if (n2a < 0) {		/* must be sphere, copy second node */
    n2a = n1a;
    n2b = n1b;
    n2c = n1c;
    n2d = n1d;
  };
  npnt  = ndn(n1a,n1b,n1c,n1d);
  npnt2 = ndn(n2a,n2b,n2c,n2d);
  
  nconn = int(get_nfield(npnt,NUMCONN));
  xden  = get_nfield(npnt,XLOC);
  yden  = get_nfield(npnt,YLOC);
  zden  = get_nfield(npnt,ZLOC);

  xden2 = get_nfield(npnt2,XLOC);   /* location of parent node */
  yden2 = get_nfield(npnt2,YLOC);

  dx = xden - xden2;
  dy = yden - yden2;

  //sangl = atanx(dx,dy) * DEG;
  if (dx == 0) dx = 1e-6;
  slope = dy/dx;
  sangl = atan(slope); 
  if (slope >= 0) {
        if(dx < 0) sangl -= PI;
  }
  else {
        if(dx < 0) sangl += PI;
  };
  if (sangl < 0) sangl += 2 * PI;
  sangl *= DEG;

  /* if (pvals>1) printf ("xden %d %g %g %d %g %g\n",
                        n1c,xden2,yden2,n2c,xden2,yden2); /* */
  csn = n1c;

  if (pvals>1)
    printf ("tracing from node %d %d %d\n", n1a,n1b,n1c);
  if (pvals>1)
    printf ("nconn %d\n",nconn);
  if (pvals>1)
    printf ("sangl %g\n",sangl);

  /* follow each branch to its end */

  if (nconn>1) {
    done = 1;
    if ((nconn>=3) && (brnum!=obrnum))
      nonbrlen = 0;	/* reset dist from br pt, but skip soma */
    for (i=1; i<=nconn; i++) {
      nelem = int(get_nfield(ndn(n1a,n1b,n1c),ELNUM,i));
      epnt = get_elempnt(nelem);
      if (pvals>1)
        printf ("conn %d: elem %d type %s\n", i,nelem,get_elabel(epnt,TYPE));
      if (nelem==brnum) continue;
      if (nelem==obrnum) continue;
      if (epnt->ctype != CABLE) continue;
      if (pvals>1)
        printf ("following conn %d: elem %d\n\n", i,nelem);

      leng = get_efield(epnt,LENGTH);
      done &= extend_branch(nelem,brnum,ctype,cellnum,xsoma,ysoma,
		     sdia,nonbrlen+leng,rad,seglen,branchthresh,growthresh);
    }
    if (pvals)
      printf ("extend_branch endm done %d\n",done);
    return done;
  }
  else {       /* If this is the end, add to it */
    if (pvals>1)
      printf ("found the end at node %d %d %d\n", n1a,n1b,n1c);

    slen  = rrange (0.8,1.2) * seglen;
    dx = xsoma - xden;
    dy = ysoma - yden;
    currad = sqrt (dx*dx + dy*dy);
    if (currad > rad*rrange(.1,1.9)) { //build in randomness in radlength
      if (pvals)
        printf ("extend_branch end done %d\n",1);
      return 1;
    }
    radfrac = currad/rad;		/* fraction of final radius */
    radtaper = (1-radfrac)*.6 + .4; 	/* 1 - fraction of final radius */

    /* try 3 locations for new tip */

    dangl = 30;
    angl1 = sangl + rrange (-1,1);
    angl2 = angl1 + dangl;
    angl3 = angl1 - dangl;

    if (angl1 < 0) angl1 += 360;
    if (angl2 < 0) angl2 += 360;
    if (angl3 < 0) angl3 += 360;

    xm1 = cos(angl1/DEG) * slen + xden;
    ym1 = sin(angl1/DEG) * slen + yden;
    zm = zden;
    xm2 = cos(angl2/DEG) * slen + xden;
    ym2 = sin(angl2/DEG) * slen + yden;
    xm3 = cos(angl3/DEG) * slen + xden;
    ym3 = sin(angl3/DEG) * slen + yden;

    /* check to make sure branch doesn't curve back towards soma */
    dx = xsoma - xm2;
    dy = ysoma - ym2;
    newrad2 = sqrt(dx*dx+dy*dy);
    if (currad > newrad2) {
      xm2 = xm1;
      ym2 = ym1;
    };
    dx = xsoma - xm3;
    dy = ysoma - ym3;
    newrad3 = sqrt(dx*dx+dy*dy);
    if (currad > newrad3) {
      xm3 = xm1;
      ym3 = ym1;
    };
    dx = xsoma - xm1;
    dy = ysoma - ym1;
    newrad = sqrt(dx*dx+dy*dy);
    if (currad > newrad) {
      if (pvals)printf ("extend_branch end done %d\n",1);
      return 1;
    };

    nsn = ++celnode[ctype][cellnum];
    nsn2 = nsn + 1;
    nsn3 = nsn + 2;
    loc(nd(ctype,cellnum,nsn), xm1,ym1,zm);
    loc(nd(ctype,cellnum,nsn2),xm2,ym2,zm);
    loc(nd(ctype,cellnum,nsn3),xm3,ym3,zm);

    /* Check proximity of other branches of same cell. */

    mdist1 = mdist2 = mdist3 = 1e10; /* now find closest dendrite */
    for (epnt=elempnt; epnt=foreach(epnt,CABLE,ctype,cellnum,-1,-1,
						NULL,NULL,&cnod,NULL);
						epnt = epnt->next) {
      cbrnum = epnt->elnum;
      if (cbrnum==brnum) continue;             /* ignore this cable */
      if (cbrnum==obrnum) continue;            /* ignore this cable */

      if ((dist=endist3d(epnt,ndn(ctype,cellnum,nsn))) < mdist1) {
        mdist1 = dist;                     /* remember this dist */
      }
      if ((dist=endist3d(epnt,ndn(ctype,cellnum,nsn2))) < mdist2) {
        mdist2 = dist;
      }
      if ((dist=endist3d(epnt,ndn(ctype,cellnum,nsn3))) < mdist3) {
        mdist3 = dist;
      }

      //if (pvals>1) printf ("cnod %d n1c %d dist %g\n", cnod,n1c,dist);
    }   /* foreach cable, remember dist. */

    /* find which direction is farthest away from nearest branches */

    erase (ndn(ctype,cellnum,nsn2));
    erase (ndn(ctype,cellnum,nsn3));
    if (mdist1 >= mdist2 && mdist1 >= mdist3) {
      mdist = mdist1;
      loc(nd(ctype,cellnum,nsn), xm1,ym1,zm);
    }
    else if (mdist2 >= mdist1 && mdist2 >= mdist3) {
      mdist = mdist2;
      loc(nd(ctype,cellnum,nsn), xm2,ym2,zm);
    }
    else if (mdist3 > mdist1 && mdist3 > mdist2) {
      mdist = mdist3;
      loc(nd(ctype,cellnum,nsn), xm3,ym3,zm);
    }

    /* Extend this branch if other branches are far enough away. */
    if (pvals>1)
      printf ("closest branch %g\n", mdist);

    gt  = growthresh * radtaper * rrange(.8,1.2);
    bt  = growthresh * 2 * radtaper * rrange(.8,1.2);
    bt2 = branchthresh * rrange(.8,1.2);

    if (varicos) {			//default 0
      pvaricos=radfrac*rrange(0,1);
      if (pvaricos>.5)
        nsdia = sdia*4;			//make varicosities
      else
        nsdia = sdia;
    }
    else if (denddiataper)		//default 1
      nsdia = sdia * radtaper;
    else
      nsdia = sdia;

    if (pvals)
      printf ("mdist %g gt %g bt %g bt2 %g nonbrlen %g\n",
                                   mdist, gt, bt, bt2, nonbrlen);
    if (mdist >= gt) {               	//extend dendrite
      if (mdist >= bt && nonbrlen > bt2) {   /* branch into 2 */
        if (pvals>1)
	  printf ("splitting branch into 2\n");
        erase (ndn(ctype,cellnum,nsn));
	radtaper = (1-radfrac)*.5 + .5;     /* 1 - fraction of final radius */
        dangl = 55 * radtaper * rrange (.7,1.3);
        angl2 = sangl + dangl;
        angl3 = sangl - dangl;
        if (angl2 < 0)
	  angl2 += 360;
        if (angl3 < 0)
	  angl3 += 360;
        zm  = zden;

	xm2 = cos(angl2/DEG) * slen + xden;
        ym2 = sin(angl2/DEG) * slen + yden;

	xm3 = cos(angl3/DEG) * slen + xden;
        ym3 = sin(angl3/DEG) * slen + yden;

	nsn2 = nsn;
        nsn3 = ++celnode[ctype][cellnum];

        loc(nd(ctype,cellnum,nsn2),xm2,ym2,zm);
        loc(nd(ctype,cellnum,nsn3),xm3,ym3,zm);
        if (pvals>1)
	  printf ("dividing branch from %d to %d and %d\n", csn, nsn2,nsn3);
        make_celseg (ctype, cellnum, nsn2, csn, nsdia, taper=0, DEND,cplam);
        make_celseg (ctype, cellnum, nsn3, csn, nsdia, taper=0, DEND,cplam);
      }
      else {                                 /* straight */
	if (pvals>1)
	  printf ("extending branch from %d to %d\n", csn, nsn);
        make_celseg (ctype, cellnum, nsn, csn, nsdia, taper=0, DEND,cplam);
      };
    }
    else {
      erase (ndn(ctype,cellnum,nsn));
      if (pvals)
        printf ("branch too close, stopping\n");
      if (pvals)
        printf ("extend_branch end done %d\n",0);
      return 1;
    };
  };   /* add onto the end */

  if (pvals)
    printf ("extend_branch end done %d\n",0);
  return 0;
}

/*----------------------------------------------------------------------*/

int extend_branch2(int brnum, int obrnum, int ctype, int cellnum, double xsoma, 
	double ysoma, double sdia, double nonbrlen, double rad, double seglen,
	double branchthresh, double growthresh, int nbranchpts)  {

    /* Extend an existing branch on a starburst cell dendrite. */
    /* Includes competition between branches of same cell type. */

    int i, nconn,nelem;
    int cnod, cbrnum;
    int n1a, n1b, n1c, n1d;
    int n2a, n2b, n2c, n2d, csn;
    int nsn, nsn2, nsn3;
    int pvals, done;
    double dx,dy;
    double leng, dist;
    double mdist,dangl,taper;
    double nsdia, slen, sangl;
    double angl1, angl2, angl3;
    double xm1, ym1, zm, xm2, ym2, xm3, ym3;
    double xden, yden, zden, xden2, yden2;
    double bt,bt2,gt;
    double radfrac,radtaper,currad;
    double newrad, newrad2, newrad3;
    double cplam;

    elem *epnt;
    node *npnt,*npnt2;

  if (notinit(denddiataper)) denddiataper=0; //default no dend diam tapering
  cplam = -1;

  pvals = 0;
  if (ninfo>3)
    pvals = 1;
  if (ninfo>4)
    pvals = 2;
  if (pvals)
    printf ("entering extend_branch2 %d %d l %g %d\n",
		          brnum,obrnum,nonbrlen,nbranchpts);

  /* find branches at distal end of segment */

  epnt = get_elempnt(brnum);
  npnt  = get_elemnode1(epnt);
  npnt2 = get_elemnode1(epnt);

  n1a = int(get_efield(epnt,NODE1A));	/* descendent node */
  n1b = int(get_efield(epnt,NODE1B));
  n1c = int(get_efield(epnt,NODE1C));
  n1d = int(get_efield(epnt,NODE1D));

  n2a = int(get_efield(epnt,NODE2A));	/* parent node */
  n2b = int(get_efield(epnt,NODE2B));
  n2c = int(get_efield(epnt,NODE2C));
  n2d = int(get_efield(epnt,NODE2D));

  if (n2a < 0) {		/* must be sphere, copy second node */
    n2a = n1a;
    n2b = n1b;
    n2c = n1c;
    n2d = n1d;
  };

  nconn = int(get_nfield(npnt,NUMCONN));	//nr of conns from n1
  xden  = get_nfield(npnt,XLOC);
  yden  = get_nfield(npnt,YLOC);
  zden  = get_nfield(npnt,ZLOC);

  xden2 = get_nfield(npnt2,XLOC);   /* location of parent node */
  yden2 = get_nfield(npnt2,YLOC);

  dx = xden - xden2;		//x-distance from n1 to n2
  dy = yden - yden2;		//y-distance from n1 to n2
  sangl = atanx(dx,dy) * DEG;	//angle with respect to x-axis between n1 and n2

  csn = n1c;

  if (pvals>1) printf ("tracing from node %d %d %d\n",
                        n1a,n1b,n1c);
  if (pvals>1) printf ("nconn %d\n",nconn);
  if (pvals>1) printf ("sangl %g\n",sangl);

  /* follow each branch to its end */
  if (nconn>1) {	//if n1 has >1 connection (ie >1 parents)
    done = 1;
    if ((nconn>=3) && (brnum!=obrnum)){
      nonbrlen = 0;/* reset dist from br pt */
      nbranchpts++;
    };
						     /* but skip soma */
    for (i=1; i<=nconn; i++) {
       elem *e;
      nelem = int(get_nfield(npnt,ELNUM,i));
      e = get_elempnt(nelem); 
      if (pvals>1) printf ("conn %d: elem %d type %s\n",
                         i, nelem,  get_elabel(e,TYPE));
      if (nelem==brnum)
        continue;
      if (nelem==obrnum)
        continue;
      if (get_efield(nelem,NTYPE)==CABLE)
        continue;
      if (pvals>1)
        printf ("following conn %d: elem %d\n\n", i,nelem);

      leng = get_efield(e,LENGTH);
      done &= extend_branch2(nelem,brnum,ctype,cellnum,xsoma,ysoma,sdia,
      			     nonbrlen+leng,rad,seglen,branchthresh,
			     growthresh,nbranchpts);
    }
    if (pvals)printf ("extend_branch2 endm done %d\n",done);
      return done;
  }
  else {       /* If this is the end, add to it */
    if (pvals>1)
      printf ("found the end at node %d %d %d\n", n1a,n1b,n1c);

    slen  =  seglen;
    dx = xsoma - xden;		     //x-distance of n1 to soma
    dy = ysoma - yden;		     //y-distance of n1 to soma
    currad = sqrt (dx*dx + dy*dy);   //distance of n1 to soma (current radius)
    if (currad > rad)  {	//if desired final radius < current radius, stop
      if (pvals)printf ("extend_branch2 end done %d\n",1);
      return 1;
    }

    radfrac = currad/rad;	     //current fraction of desired final radius
    radtaper = (1-radfrac)*.6 + .4;  //1 - fraction of final radius

    /* make new tip */
    angl1 = sangl;
    if (angl1 < 0)
      angl1 += 360;

    xm1 = cos(angl1/DEG) * slen + xden;
    ym1 = sin(angl1/DEG) * slen + yden;
    zm = zden;
    nsn = ++celnode[ctype][cellnum];
    loc(nd(ctype,cellnum,nsn),xm1,ym1,zm);

    /* Extend this branch if other branches are far enough away. */

    gt  = growthresh * radtaper ;
    bt  = growthresh * 2 * radtaper ;
    bt2 = arr_branchdist[nbranchpts] ;
    if (denddiataper)
      nsdia = sdia * radtaper;
    else
      nsdia = sdia;

    //for now, mdist is set to some large value so always mdist>gt,bt
    mdist = 1000;

    // if (pvals>1) printf ("closest branch %g\n", mdist);

    if (mdist >= gt) {               	      /* extend dendrite */
      if (mdist >= bt && nonbrlen > bt2) {   /* branch into 2 */
        if (pvals>1)
	  printf ("splitting branch into 2\n");
        erase (ndn(ctype,cellnum,nsn));
        radtaper = (1-radfrac)*.5 + .5;     /* 1 - fraction of final radius */
        dangl = arr_sbrangles[nbranchpts]; /* arr hold angles for every brchpt */
        angl2 = sangl + dangl;
        angl3 = sangl - dangl;

        if (angl2 < 0) angl2 += 360;
        if (angl3 < 0) angl3 += 360;

        zm = zden;
        xm2 = cos(angl2/DEG) * slen + xden;
        ym2 = sin(angl2/DEG) * slen + yden;

        xm3 = cos(angl3/DEG) * slen + xden;
        ym3 = sin(angl3/DEG) * slen + yden;
        nsn2 = nsn;
        nsn3 = ++celnode[ctype][cellnum];

        loc(nd(ctype,cellnum,nsn2),xm2,ym2,zm);
        loc(nd(ctype,cellnum,nsn3),xm3,ym3,zm);

        if (pvals>1) printf ("dividing branch from %d to %d and %d\n",
                                       csn, nsn2,nsn3);
        make_celseg (ctype, cellnum, nsn2, csn, nsdia, taper=0, DEND,cplam);
        make_celseg (ctype, cellnum, nsn3, csn, nsdia, taper=0, DEND,cplam);
      }
      else {                                 /* straight */
        if (pvals>1)
	  printf ("extending branch from %d to %d\n", csn, nsn);
        make_celseg (ctype, cellnum, nsn, csn, nsdia, taper=0, DEND,cplam);
      };
    }
    else {
      erase (nd(ctype,cellnum,nsn));
      if (pvals) printf ("branch too close, stopping\n");
      if (pvals) printf ("extend_branch2 end done %d\n",0);
      return 1;
    };

  };	/* add onto the end */
  if (pvals)
    printf ("extend_branch2 end done %d\n",0);
  return 0;
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void maksbacs(int ct, int cellnum, double x, double y, double t, double scale, int nbr) 

    /* Make one amacrine cell, given position and rotation. */
    /* Make straight, non-branching dendrites. */
{

    int d, i, amregion, dbr;
    double dt, r, u, z;
    double amdia, amdia_fac, rmrange;
    double taper;
    double cplam;

  z = getn(ct,SOMAZ);
  loc(nd(ct,cellnum,soma),x,y,z);
  make_celseg(ct, cellnum, soma, soma, amdia=8, taper=0, amregion=SOMA,cplam);

  dt = 2 * PI / nbr;

  for (d=1; d<=am_den_seg; d++) {
    r = d * am_seglen;       /* extend branches one varicosity at a time */
    for (i=0; i<nbr; i++) {           /* make "nbr" branches */
      dbr = int(d + i*am_den_seg);
      u = t + i*dt;
      loc(nd(ct,cellnum,dbr), x+r*cos(u),y+r*sin(u),z);
      //make_celseg(ct, cellnum, dbr, dbr, amdia=2, taper=0, amregion=VARICOS,cplam);

      if (r<=am_dia_prox_rad)
        amdia_fac= am_dia_prox_factor;
      else
        amdia_fac=1;

      if (d==1)
        make_celseg(ct, cellnum, dbr, soma, amdia=sbac_dend_dia*amdia_fac,
	   				    taper=0, amregion=DEND,cplam);
      else if (r<xradius*0.66)
	make_celseg(ct, cellnum, dbr, dbr-1, amdia=sbac_dend_dia*amdia_fac,
	   			            taper=0, amregion=DEND,cplam);
      else
	make_celseg(ct, cellnum, dbr, dbr-1, amdia=sbac_dend_dia*amdia_fac,
	   			            taper=0, amregion=DEND_DIST,cplam);
    };
  };
  n_amseg = am_den_seg * nbr;

};    /* proc makamacs() */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void maksbacx(int ct, int cellnum,double x,double y, double t, double scale, int nbr) {

    /* Make one amacrine cell, given position and rotation. */
    /* Split branch whenever growing tip is farther than a threshold distance
      from nearest node in the cell.  Tip grows in direction farthest away
      from nearest node in the cell. */

    int d, dbr, i, done;
    int soma_elnum, amregion;
    double dt, r, u, ut, ta, z;
    double amdia, amsomdia;
    double nonbrlen, rad, sdia, taper;
    double branchthresh, growthresh, seglen;
    double varbt;
    double *randa;
    double cplam;

  randa = (double *)emalloc(nbr*sizeof(double));
  if (notinit(denddiataper))
    denddiataper=1; 	//default do dend diam tapering
  amsomdia = 10;	//soma diameter
  cplam = -1;

  z = getn(ct,SOMAZ);
  loc(nd(ct,cellnum,soma),x,y,z);
  make_celseg(ct, cellnum, soma, soma, amdia=amsomdia, taper=0, amregion=SOMA, cplam);

  /* start the dendritic branches */
  dt = 2 * PI / nbr;
  r = am_seglen;
  dbr = 1;			/* dendrite node number */
  for (ut=i=0; i<nbr; i++) {    /* get "nbr" random numbers */
    u = rrange(.6,1.4);
    randa[i] = u;
    ut += u;
  }
  for (i=0; i<nbr; i++)        	/* normalize to total of 1 */
    randa[i] *= nbr/ut;
  for (ta=t,i=0; i<nbr; i++) {       /* start "nbr" branches */
    ta += dt * randa[i];	//total angle
    loc(nd(ct,cellnum,dbr),x+r*cos(ta),y+r*sin(ta),z);
    make_celseg(ct, cellnum, dbr, soma, amdia=sbac_dend_dia*2, taper=0, amregion=DEND, cplam);
    dbr++;
    celnode[ct][cellnum]++;
  }

  /* grow the dendritic branches */
  soma_elnum = int(get_nfield(ndn(ct,cellnum,soma),ELNUM,1));

  for (done=i=0; i<35 && !done; i++) {
    varbt = (xradius/4.5) - pow(1.05,i);
    done=extend_branch(soma_elnum,soma_elnum,ct,cellnum,x,y,sdia=sbac_dend_dia,
	       nonbrlen=25,rad=xradius,seglen=4,branchthresh=varbt,growthresh=6);
  }

  if (notinit(makesbfile)) makesbfile=0;
  if (makesbfile==1) {
    makanatfile (sbac,cellnum);
    return;
  }
  free (randa);
}    /* proc makamacx() */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void maksbacc(int ct, int cellnum, double x,double y, double t, double scale, int nbr) {

 /* Make one amacrine cell, given position and rotation. */
 /* Make nonrandom symmetric cell, according to user specifications */

    int dbr, i, done;
    int amregion, soma_elnum, nbranchpts;
    double dt, r, u, ut, ta, z;
    double amdia;
    double nonbrlen, rad, sdia, taper;
    double branchthresh, growthresh, seglen;
    double *randa;
    double cplam;

  cplam = -1;
  randa = (double *)emalloc(nbr*sizeof(double));
  z = getn(ct,SOMAZ);
  loc(nd(ct,cellnum,soma),x,y,getn(ct,SOMAZ));

  if (notinit(nonrand)) nonrand = 1;	     //make default nonrandom
  if (notinit(denddiataper)) denddiataper=0; //default no tapering

  make_celseg(ct, cellnum, soma, soma, amdia=8, taper=0, amregion=SOMA,cplam);

  /* start the dendritic branches */
  if (nbr==1 && nbranchpoints==1 && branchdist==0){
  /* Special case: make 2 straight cables from soma */
    dt=sbrangles/DEG;
    arr_branchdist[0]=xradius*2; // make cable without branchpoints
    nbr=2;	                 // make 2 cables
  }
  else
    dt = 2 * PI / nbr;
  r = am_seglen;
  dbr = 1;			 /* dendrite node number */
  for (ut=i=0; i<nbr; i++) {     /* get "nbr" random numbers */
    u = rrange(.6,1.4);
    randa[i] = u;
    ut += u;
  };
  for (i=0; i<nbr; i++)          /* normalize to total of 1 */
    randa[i] *= nbr/ut;

  for (ta=t,i=0; i<nbr; i++) {   /* start "nbr" branches */
    ta += dt * randa[i];		   // total angle
    loc(nd(ct,cellnum,dbr),x+r*cos(ta),y+r*sin(ta),z);
    make_celseg(ct, cellnum, dbr, soma, amdia=sbac_dend_dia*2, taper=0, amregion=DEND,cplam);
    dbr++;
    celnode[ct][cellnum]++;
  };

  /* grow the dendritic branches */
  soma_elnum = int(get_nfield(nd(ct,cellnum,soma),ELNUM,1));

  for (done=i=0; i<35 && !done; i++)
    done=extend_branch2(soma_elnum,soma_elnum,ct,cellnum,x,y,sdia=sbac_dend_dia,
                        nonbrlen=25,rad=xradius,seglen=6,branchthresh=branchdist,
			growthresh=7,nbranchpts=0);
  free (randa);
};    /* proc makamacc() */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/* realistic amacrine morphology */

/* Here make just one of the dendrites, set by "sbac_dend"
   and rotate it with "sbtheta". Then rotate it
   randomly within a small range to give some randomness.  
   
   The dendrites are labeled by the "dendr" column in the anatomy file.  
   To set the dendrite number this way is a little arbitrary since a "dendrite"
   invariably branches. Five sbac_dend1-5 variables allow several dendrites
   to be selected */

/* orig combination for sbac3 anatomy -- sbac_dend=3,sbtheta=-40*PI/180 */

/* Another parameter, sbac_dend_cplam, added to allow other dendrites
   to be assigned a different complam than specified by the density file. */
 
/* Allow the user to change the diameter of the dendrites in a
   convenient way.  The diameter of dendrites can be specified with a
   number or with the the "amdia" variable in the anatomy file (e.g.
   sb1).  This variable sets the value of the diameter that is
   stored in the "cabldata" array for use by the script.

   The factor "am_dia_prox_factor" multiplies only the dendrites
   within the radius specified by "am_dia_prox_rad".  This is to
   allow the primary dendrites to be a little larger as originally
   described by Tauchi & Masland (1984). */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* - - - - - - - - - - - make starburst cells  - - - - - - - - - - - */

void mak_sbac(int ct, int n, int nden, double xpos, double ypos, double thetax, 
		double thetay, double thetaz, int flip, int morph) 

{
  double arbor_scale;
  double sbac_theta;

  sb_init();
  if(notinit(sbac_info)) sbac_info=0;

  if (notinit(rand_amth_sb)) rand_amth_sb=0; /* default: sb at random angle*/
  if (rand_amth_sb) {
    sbtheta = drand()*360;
  }				/* else use user-specified or preset angles */
  else if (notinit(sbtheta) && sbac_file=="morph_sb1")    sbac_theta=PI;
  else if (notinit(sbtheta) && sbac_file=="morph_sbac3b") sbac_theta=PI;
  else if (notinit(sbtheta)) sbac_theta = thetaz;  	/* rot to make dend = 0 deg */
  else                       sbac_theta = sbtheta;

  if (sbac_info >=2)
    fprintf(stderr, "#sbacfuncs, mak_sbac: sbac_theta for sb %d = %g deg\n", n, sbac_theta);

  if (notinit(sbac_dend_cplam)) sbac_dend_cplam = -1; /* sets complam for dendrites not sbac_dend  */
  celnode[ct][n] = 0;				// count nr of dendrites

  if (notinit(amdia)) amdia = 0.5; 	/* diameter of sb dendrites in file */
  if (notinit(dia_prox_rad))    dia_prox_rad    = 25; /* radius for prox dia */
  if (notinit(dia_prox_factor)) dia_prox_factor = 2; /* dia factor */
  if (notinit(taper)) taper = 1; /* make cables with taper */

  if      (morph==MORPH_REAL) {
         arbor_scale = getn(ct,ARBSCALE);
         if (arbor_scale==0) arbor_scale=1.0;
	  mak_real_cel (ct,n,xpos,ypos,thetax,thetay,sbac_theta,flip,
		     sbac_dend1,sbac_dend2,sbac_dend3,sbac_dend4,sbac_dend5,
		     sbac_dend_cplam,taper);  /* use morph file */
  }
  else if (morph==MORPH_A1) maksbacs     (ct,n,xpos,ypos,sbac_theta/180*PI,arbor_scale,nden);    /* artif neuron */
  else if (morph==MORPH_A2) maksbacx     (ct,n,xpos,ypos,sbac_theta/180*PI,arbor_scale,nden);
  else if (morph==MORPH_A3) maksbacc     (ct,n,xpos,ypos,sbac_theta/180*PI,arbor_scale,nden);
  else fprintf (stderr,"retsim: mak_sbac, unknown sb_morph, %d\n",morph);
};

