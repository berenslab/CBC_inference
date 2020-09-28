/* module celfuncs for script retsim */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ncio.h"

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "drand.h"
#include "nconst.h"

static int ptrace = 0;

int tot_ncel_ind(int ct, int cn);

/*-------------------------------------------------*/

double mod (double a,double b)

{
    double div,a1,a2;

  if (b>0)
  {
   a1=int(a/b);
   a2=a1*b;
   return (a-a2);
  }
  else return 0;
}

/*-------------------------------------------------*/

double round(double n,double p)

{
   return (int((n / p) + 0.5) * p);
};

/*-------------------------------------------------*/

double modangl(double a)

/* set angle between 0 and 2*PI, in radians */

{
 while (a<0)      { a+= 2*MPI;}
 while (a>=2*MPI) { a-= 2*MPI;}
 return a;
};

/*-------------------------------------------------*/

double sindeg (double theta)

/* sine function in degrees */

{
  return (sin(theta*MPI/180));
};
/*-------------------------------------------------*/

double cosdeg (double theta)

/* cosine function in degrees */

{
  return (cos(theta*MPI/180));
};

/*-------------------------------------------------*/

/* function to test range of angle */

int inrange(double a1, double a2, double t)

/* a1, a2 define a chord range. a1 is supposed to be less than a2. */
/* t is supposed to be 0 <= t < 360 */
/* If a1 is less than a2, and if t is in that range, return 1. */
/* However if a1 is greater than a2, if t is outside that range   */
/*  return 1, else return 0.       */

{
   int retval, sign;

  a1 = modangl(a1);
  a2 = modangl(a2);
//  t  = modangl(t);

//  if (a1>a2) {
//     a2 += 2*MPI;
//     if (t < a1) t += 2*MPI;
//  };

  if (a1 < a2) retval = 1;
  else         retval = 0;
  if (a1>t)    retval ^= 1;
  if (a2<t)    retval ^= 1;

  if (abs(a1-a2) < 1e-6) retval = 1;
//   fprintf (stderr,"a1 %20.20g a2 %20.20g t %g ret %d\n",a1,a2,t,retval);
  return retval;
};

/*-------------------------------------------------*/

double atanx (double dx, double dy)

/* similar to the atan2(y,x) function in the C library */

{
    double slope, theta;

  if (dx == 0) dx = 1e-20;
  slope = dy/dx;
  theta = atan(slope);
  if (slope >= 0) {
        if(dx < 0) theta -= MPI;
  }
  else {
        if(dx < 0) theta += MPI;
  };
   if (theta < 0) theta += 2 * MPI;
  return theta;
};

/*-------------------------------------------------*/

int mid(int siz)

/* Returns center element of array */
/* Useful to find "center cone" */
 
{
   int m;

  if (int(siz/2)*2==siz) {
/*    print "even"; */
    m = (siz+1) * siz / 2;
  }
  else {
/*    print "odd"; */
    m = (siz*siz-1) / 2;
  };
  return m;
};

/*-------------------------------------------------*/

int midrow(int siz)

/* function to return start of middle row of array */

{
  return (siz*int(siz/2));
};

/*-------------------------------------------------*/
 
int ff(void)

/* return fifty percent chance */

{
  return (drand() > 0.5);
};

/*-------------------------------------------------*/

int ffs(void)

/* return 1 or -1 with fifty percent chance */

{
  return ((drand() > 0.5) * 2 - 1);
};

/*-------------------------------------------------*/

double rrange (double l, double h)

/* return a number between L and H */

{
    double t__;

  if (h < l) { 
      t__ = l;
      l = h; 
      h = t__;
   };

  return ( l + drand() * (h - l) ); 
};

/*-------------------------------------------------*/

double gauss(double r, double rad)

{
 return exp(-(r*r)/(rad*rad));
}


/*-------------------------------------------------*/

int bptype (int ct)
{
   int retval;
   switch (ct) {
   	case RBP:
   	case DBP1:
   	case DBP2:
   	case HBP1:
   	case HBP2: retval = 1; break;
   	  default: retval = 0; break;
   }
   return retval;
}

/*--------------------------------------------------*/

double node_angle (node *npnt1, node *npnt2)

/* Find angle between two nodes */

{
    double nx, ny;

  nx = npnt1->xloc - npnt2->xloc;
  ny = npnt1->yloc - npnt2->yloc;
  return (atanx(nx,ny));
  
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - */

double node_angle (node *npnt)

/* Find angle between a node and the soma */

{
  int ct,cn;

  ct = npnt->nodenm1;
  cn = npnt->nodenm2;
  return (node_angle (npnt, ndn(ct,cn,soma)));
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - */

double node_angle (int ct, int cn, int n)

/* Find angle between a node and the soma */

{
  return (node_angle (ndn(ct,cn,n)));
}

/*--------------------------------------------------*/

void getangles(int ct, int cn, double *angles, int arrsize)

/* Function to determine orientation of dendrites */
/* Takes as arguments: celltype, cell number, empty array of doubles,  */
/* and integer size of array.  */
 
/* It stores angles of dendrites in array,  indexed by dendrite number. */
/* The dendrite number is read from the node, originating from the morphology file. */

{
     node *npnt;
     int i;
     int numdend=0;		/* number of dendrites */

#define NX 0
#define NY 1

  for(npnt=nodepnt; 
      npnt=foreach(npnt,ct,cn,-1, NULL, NULL,NULL); npnt=npnt->next) {
        if (numdend < npnt->dendn) numdend = npnt->dendn; 		  /* count dendrites */
  }
  if (numdend >= arrsize) {
      fprintf(stderr,"getangles: Error, number of dendrites %d exceeds maximum array size %d\n",
		      			numdend, arrsize);
      return;
  }

double dendrxy[numdend][2]; 		/* array for summed x, y positions of dendrites */

  for (i=0; i<numdend; i++) { 	/* initialize array values */
	dendrxy[i][NX]=0;
	dendrxy[i][NY]=0;
  }
	
  for(npnt=nodepnt; 			/* Sum vectors of node positions for each dendrite. */
      npnt=foreach(npnt,ct,cn,-1, NULL, NULL,NULL); npnt=npnt->next) { 
	   dendrxy[npnt->dendn][NX] += npnt->xloc;
	   dendrxy[npnt->dendn][NY] += npnt->yloc;
	   fprintf(stderr,"%g,%g\n",dendrxy[npnt->dendn][NX],dendrxy[npnt->dendn][NY]);
  }

  for (i=0; i<arrsize; i++) { 			/* store angles in array */
    if (i>0 && i<numdend) {			/* skip dend. number=0 */

      angles[i] = atan2(dendrxy[i][NY],dendrxy[i][NX]) * 180 / PI;
      //fprintf(stderr,"Dendrite %d: (%g,%g), angle: %g \n", i,dendrxy[i][NX],dendrxy[i][NY],angles[i]);

    } else {
      angles[i]=-361;				/* dend = 0, overflow, -> set angle to -361.*/
    }
  }
}	/* getangles */


#undef NX
#undef NY

/*--------------------------------------------------*/

double rad_dist2 (int celtyp, int n, int tnod, int parent, int stflag) 

/* Compute distance along dendrite back to cell body
   as the radial distance from the distal node to
   the soma.
*/

{
  return (dist3d (ndn(celtyp,n,tnod), ndn(celtyp,n,0)));
};

/*-------------------------------------------------*/

double rad_dist (int celtyp, int n, int tnod, int parent, int stflag) 

/* Compute distance along dendrite back to cell body.
   Assumes that dendrites with lower node numbers are 
   created first (closer to soma).
*/

{
  int i,diam,nconn,ncabl,nelem;
  int nod1, nod2, tracenode;
  double this_length, desc_length;
  elem *epnt;

// printf ("rad_dist %g %g %g %g %g\n", celtyp,n,tnod,parent,stflag);

  if (ptrace) ncfprintf (stdout,
	"rad_dist, node %d %d %d, parent %d start %d\n",celtyp,n,tnod,parent,stflag);

  nconn = (int)get_nfield(ndn(celtyp,n,tnod),NUMCONN);

  if (ptrace) ncfprintf (stdout,"rad_dist, node %d, parent %d\n",tnod,parent);

   if (tnod == 0) return 0;
   if (nconn <= 1 && !stflag) return 0;

   for (ncabl=0,i=1; i<=nconn; i++) {        /* Count the cable connections */
        epnt = get_elempnt(ndn(celtyp,n,tnod), i);
	if (epnt->ctype == CABLE) ncabl++;
   };
   if (ptrace) ncfprintf(stdout," %d cable connections.", ncabl);

   for (i=1; i<=nconn; i++) {  /* Find the path back to soma */
        epnt = get_elempnt(ndn(celtyp,n,tnod), i);
	if (epnt->ctype != CABLE) continue;
	if (ptrace) ncfprintf (stdout,"tracing conn %d\n", i);
        nod1= epnt->node1c;
        nod2= epnt->node2c;
        if (nod1 == parent || nod2 == parent) continue; 
        if (nod1 != tnod && nod2 == tnod) {
	   tracenode = nod1;
	}
	else {
	  if (nod2 != tnod && nod1 == tnod) {
	    tracenode = nod2;
	  }
	  else {
	   ncfprintf (stderr,"rad_dist: error, connection to node %d not found.\n",tnod);
	   return 0;
	  };
	};
	if (tracenode < 0) tracenode = 0;

	if (ptrace) ncfprintf (stdout,"tracenode %d\n", tracenode);

	this_length= get_length(epnt);
     if (ptrace)
        ncfprintf (stdout,"tracing element %d nodes %d %d length %g\n",
			nelem,nod1,nod2,this_length);
        if (tracenode==0) {
	   return this_length; 
	} 
	else {
	   if ((desc_length=rad_dist (celtyp, n, tracenode, tnod, 0)) != 0) {
      if (ptrace) ncfprintf (stdout,"rad_dist returning with %g\n",desc_length);
		return (this_length + desc_length); 
	   };
	};

   }; /* for (i;;) */
   return 0; 
}

/*-------------------------------------------------*/

int rad_dir (int elnum, int c1, int c2, int c3) 

/* Find which end of element is closest to center of cell,
   and return a 1 if it should be rotated. */

{
   int n1a,n1b,n1c;
   int n2a,n2b,n2c;
   elem *epnt;

  epnt = get_elempnt(elnum);
  n1a = epnt->node1a;
  n1b = epnt->node1b;
  n1c = epnt->node1c;

  n2a = epnt->node2a;
  n2b = epnt->node2b;
  n2c = epnt->node2c;

  if (epnt->ctype == CABLE) {
     if (dist3d(ndn(n1a,n1b,n1c), ndn(c1,c2,c3)) >
         dist3d(ndn(n2a,n2b,n2c), ndn(c1,c2,c3)))
       return 1;
       else return 0;
   }
  else return 0;
}

/*-------------------------------------------------*/

double rad_diam (double dist,double diaspc,double diabs)

{
  double diam;

  diam = exp (-dist/diaspc) * diabs;
  return diam;
};

/*-----------------------------------------*/

int taperden  (int ct,int cn,int n1,int n2, double cd1, double cd2,int nden) {

/* Make a tapered dendrite starting at existing nodes 
   [n1a][n1b][n1c], [n2a][n2b][n2c], adding extra cable 
   segments if nodes are far apart.  Return the minor 
   node number of the newly created node.
*/

    int i, nsegs, region;
    double x1,y1,x2,y2,z1,z2;
    double dist,newdist,nlength;
    double xincr,yincr,zincr,dincr;
    double xrange,yrange,zrange,drange;
    double cplam;
    node *npnt;

 npnt = ndn(ct,cn,n1);
 x1 = npnt->xloc;
 y1 = npnt->yloc;
 z1 = npnt->zloc;

 npnt = ndn(ct,cn,n2);
 x2 = npnt->xloc;
 y2 = npnt->yloc;
 z2 = npnt->zloc;

 xrange = (x2-x1);
 yrange = (y2-y1);
 zrange = (z2-z1);

 nlength = sqrt (xrange*xrange + yrange*yrange);
 nsegs = int (nlength / 10) + 1;
 xincr = xrange / nsegs;
 yincr = yrange / nsegs;
 zincr = zrange / nsegs;

 drange = cd1 - cd2;

 if (nsegs > 2) dincr = drange / (nsegs-1);
 else           dincr = .1;

 if (ptrace) 
ncfprintf (stdout,"taperden: cable from node %d %d %d to %d %d %d dia %g\n\n",
						ct,cn,n1, ct,cn,n2, cd1);

 cplam = -1;

 if (nsegs==1) {
    make_celseg(ct,cn,n1,n2,cd1,0,region=DEND,cplam);
 }
 else {					/* nsegs = 2 or more */

    npnt = nd(ct,cn,nden+1); npnt->xloc=x1+xincr; 
			     npnt->yloc=y1+yincr;
			     npnt->zloc=z1+zincr;
    make_celseg(ct,cn,n1,nden+1,cd1,0,region=DEND,cplam);

   for (i=1; i<nsegs-1; i++) {		/* middle segments */

     npnt = nd(ct,cn,nden+i+1); npnt->xloc= x1+(i+1)*xincr; 
			        npnt->yloc= y1+(i+1)*yincr; 
			        npnt->yloc= z1+(i+1)*zincr;

     make_celseg(ct,cn,nden+1,nden+i+1,cd1-dincr*i,0,region=DEND,cplam);

   }

    make_celseg(ct,cn,nden+i,n2,cd2,0,region=DEND,cplam);
 };

 return nden+nsegs-1;
};

/*-------------------------------------------------*/

int taperdistden(int ct,int cellnr,int nod1,int nod2,double diabs,double diaspc,int region,int nden) {

/* Make a tapered dendrite starting at existing nodes 
   [n1a][n1b][n1c], [n2a][n2b][n2c], adding extra cable 
   segments if nodes are far apart.  Return the minor 
   node number of the last newly created node.

   Diameters are calculated from distance back to the soma.
*/
    double xrange,yrange,drange;
    double x1,y1,x2,y2;
    double d1,d2;
    double dist,newdist,nlength,cplam;
    node *npnt;

 npnt = ndn(ct,cellnr,nod1);
 x1 = npnt->xloc;
 y1 = npnt->yloc;

 npnt = ndn(ct,cellnr,nod2);
 x2 = npnt->xloc;
 y2 = npnt->yloc;

 xrange = (x2-x1);
 yrange = (y2-y1);

 dist = rad_dist2 (ct, cellnr, nod1, nod2, 1);	/* find distance to soma */

 if (ptrace) ncfprintf (stdout,"taperdistden: node %d %d %d dist %g\n",ct,cellnr,nod1,dist);

 newdist = sqrt (xrange*xrange + yrange*yrange);

 d1 = rad_diam (dist,diaspc,diabs);
 d2 = rad_diam (newdist,diaspc,d1);
 cplam = -1;

/*  return (taperden(ct,celnr,nod1, ct,celnr,nod2, d1,d2,nden)); */

 make_celseg(ct,cellnr,nod1,nod2,d1,d2,region,cplam);

 return nden;
}

/*-------------------------------------------------*/

double taperdia  (int n1a,int n1b,int n1c, int n2a,int n2b,int n2c, 
					double diabs, double diaspc) 
{

/* Find diameter of far end of a tapered dendrite, given start
   diameter. Calculate from space constant. 
*/
    double xrange,yrange,drange; 
    double x1,y1,x2,y2; 
    double edia;
    double dist,newdist,nlength;
    node *npnt;

 npnt = ndn(n1a,n1b,n1c);
 x1 = npnt->xloc;
 y1 = npnt->yloc;

 npnt = ndn(n2a,n2b,n2c);
 x2 = npnt->xloc;
 y2 = npnt->yloc;

 xrange = (x2-x1);
 yrange = (y2-y1);

 dist = rad_dist2 (n1a, n1b, n1c, n2c, 1);	/* find distance to soma */

 if (ptrace) ncfprintf (stdout,"taperdistden: node %d %d %d dist %g\n", n1a,n1b,n1c,dist);

 newdist = sqrt (xrange*xrange + yrange*yrange);

 edia = rad_diam (newdist,diaspc,diabs);

 return (edia);
}

/*-------------------------------------------------*/

double sigm(double xmin,double xmax,double ymin,double yrange,double x)

/* Returns y-val on specified sigm func for pointx */

{
     double y, midpoint;

 midpoint=(xmin+xmax)/2;
 y = ymin + yrange/(1+exp(-x+midpoint));
 return y;
}

/*-------------------------------------------------*/

double comp_phase(double tfreq,double delaytime) 

{
  double phase, ncycles;

  while (delaytime < 0)	      //make delaytimes positive
    delaytime = delaytime + 1/tfreq;
  ncycles = delaytime * tfreq; //nr of cycles in delaytime
  while (ncycles>=1)	      //subtract all whole cycles
     ncycles--;
  phase = ncycles*360;	      //convert remainder into degrees

  return phase;
}

/*-------------------------------------------------*/

double sinewaves(double phase1,double phase2,double ampl1,double ampl2) 

  //add 2 sinewaves, return ampl of sum
{
  int nsteps;
  double shiftrad1, shiftrad2;
  double x, ysum, ypeak, yratio;
  double y1, y2;

  shiftrad1 = (phase1/360)*2*MPI;  //deg -> rad
  shiftrad2 = (phase2/360)*2*MPI;

  ysum = 0;
  ypeak = 0;
  yratio = 0;

  nsteps=100;

  for (x=0;x<=2*MPI;x+=(2*MPI)/nsteps) {
    y1 = ampl1*sin(x+shiftrad1);
    y2 = ampl2*sin(x+shiftrad2);
    ysum = y1+y2;
    if (ysum>ypeak) ypeak=ysum;
  };

  yratio = ypeak/(ampl1+ampl2);

  return ypeak;

}

/*-------------------------------------------------*/

int findmid(int ct, double xoffset, double yoffset)

/* Procedure to find a cell in array to record from. */
/* Find cell closest to (xoffset, yoffset) from center of array. */
{
   int i,n,cn,found,midn,midcn;
   double sumx, sumy;
   double dx, dy, xcent, ycent, dist, mindist;
   node *npnt;

  midcn = found = n = 0;
  sumx = sumy = 0;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
     sumx += npnt->xloc;		// find x,y center
     sumy += npnt->yloc;
     n++;
  }
  if (n==0) n = 1;
  xcent = sumx / n;
  ycent = sumy / n;
  mindist = 1e10;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
     dx = npnt->xloc - xcent - xoffset + cn * 1e-6;
     dy = npnt->yloc - ycent - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       midcn = cn;
       found = 1;
     }
  }
  if (!found) { 			// find first one that still exists
     for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
       midcn = cn;
       break;
     }
  }
  return midcn;
}

/*-------------------------------------------------*/

int findmida(int ct, double xoffset, double yoffset)

/* Procedure to find a cell to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   int i,n,cn,nod,found,closecn;
   double dx, dy, dist, mindist;
   node *npnt;

  cn = found = 0;
  mindist = 1e10;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
     dx = npnt->xloc - xoffset;
     dy = npnt->yloc - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       closecn = cn;
       found = 1;
     }
  }
  if (!found) closecn = 0;
  return closecn;
}

/*-------------------------------------------------*/

int findcell(int ct, double xoffset, double yoffset)

/* Procedure to find a cell soma closest to absolute (xoffset, yoffset) */

{
    int cn,found,closecell;
    double dx, dy, dist, mindist;
    node *npnt;
 
  closecell = found = 0;
  mindist = 1e10;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
     dx = npnt->xloc - xoffset;
     dy = npnt->yloc - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       closecell = cn;
       found = 1;
     }
  }
  if (!found) closecell = -1;
  return closecell;
}

/*-------------------------------------------------*/

int findcella(int ct, double roffset, double theta)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (roffset, theta) from its soma. */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findcell(ct, xoffset, yoffset);
}

/*-------------------------------------------------*/

int findnodloc(int ct, int cn, double xoffset, double yoffset)

/* Procedure to find a node in cell to record from. */
/* Find node closest to absolute (xoffset, yoffset). */
{
   int nod,found,closenod;
   double dx, dy, dist, mindist;
   node *npnt;

  closenod = found = 0;
  mindist = 1e10;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, cn, -1, NULL,NULL,&nod); npnt=npnt->next) {
     dx = npnt->xloc - xoffset;
     dy = npnt->yloc - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       closenod = nod;
       found = 1;
     }
  }
  if (!found) closenod = -1;
  return closenod;
}

/*-------------------------------------------------*/

int findnodloc(int ct, int cn, double xoffset, double yoffset, double maxdist)

/* Procedure to find a node in cell to record from. */
/* Find node closest to absolute (xoffset, yoffset), but nearer than maxdist. */
{
   int i,n,nod,found,closenod;
   double dx, dy, dist, mindist;
   node *npnt;

  closenod = found = 0;
  mindist = 1e10;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, cn, -1, NULL,NULL,&nod); npnt=npnt->next) {
     dx = npnt->xloc - xoffset;
     dy = npnt->yloc - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       closenod = nod;
       found = 1;
     }
  }
  if (mindist > maxdist) closenod = -1;
  if (!found) closenod = -1;
  return closenod;
}

/*-------------------------------------------------*/

int findnodlocr(int ct, int cn, double xoffset, double yoffset)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (xoffset, yoffset) relative to soma. */
{
   int i,n,nod,found,closenod;
   double dx, dy, dist, mindist;
   node *npnt, *somanode;

  closenod = found = 0;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (npnt=nodepnt; npnt=foreach(npnt, ct, cn, -1, NULL,NULL,&nod); npnt=npnt->next) {
     dx = npnt->xloc - xoffset;
     dy = npnt->yloc - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       closenod = nod;
       found = 1;
     }
  }
  if (!found) closenod = -1;
  return closenod;
}

/*-------------------------------------------------*/

int findnodlocrz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (xoffset, yoffset) relative to soma. */
/* Works with cn = -1 */

{
   int i,n,ncn,nod,found,closenod;
   double dx, dy, dist, mindist;
   node *npnt, *somanode;

  closenod = found = 0;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (npnt=nodepnt; npnt=foreach(npnt, ct, cn, -1, NULL,&ncn,&nod); npnt=npnt->next) {
     if ((npnt->nodenm1==ct && npnt->nodenm2==ncn) && 
        (zmin <= npnt->zloc && npnt->zloc <= zmax)) {
        dx = npnt->xloc - xoffset;
        dy = npnt->yloc - yoffset;
        dist = sqrt(dx*dx+dy*dy);
        if (mindist>dist) {
          mindist = dist;
          closenod = nod;
          found = 1;
        }
     }
  }
  if (!found) closenod = -1;
  return closenod;
}

/*-------------------------------------------------*/

int findnodlocr(int ct, int cn, double xoffset, double yoffset, double maxdist)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (xoffset, yoffset) from its soma. */
{
   int i,n,nod,found,closenod;
   double dx, dy, dist, mindist;
   node *npnt, *somanode;

  closenod = found = 0;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (npnt=nodepnt; npnt=foreach(npnt, ct, cn, -1, NULL,NULL,&nod); npnt=npnt->next) {
     dx = npnt->xloc - xoffset;
     dy = npnt->yloc - yoffset;
     dist = sqrt(dx*dx+dy*dy);
     if (mindist>dist) {
       mindist = dist;
       closenod = nod;
       found = 1;
     }
  }
  if (mindist > maxdist) closenod = -1;
  if (!found) closenod = -1;
  return closenod;
}

/*-------------------------------------------------*/

int findnodlocra(int ct, int cn, double roffset, double theta)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (roffset, theta) from its soma. */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findnodlocr(ct, cn, xoffset, yoffset);
}

/*-------------------------------------------------*/

int findnodlocraz(int ct, int cn, double roffset, double theta, double zmax, double zmin)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (roffset, theta) from its soma. */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findnodlocrz(ct, cn, xoffset, yoffset, zmax, zmin);
}

/*-------------------------------------------------*/

int findnodlocra(int ct, int cn, double roffset, double theta, double maxdist)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (roffset, theta) from its soma. */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findnodlocr(ct, cn, xoffset, yoffset,maxdist);
}

/*-------------------------------------------------*/

synapse *findsyn(int ct, int cn, int nod)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Works with cn, nod = -1 */
{
   int i,n,ncn,nnod,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, nod, NULL,&ncn,&nnod); epnt=epnt->next) {
     if (epnt->node1a==ct && epnt->node1b==ncn && epnt->node1c==nnod) {
         closesyn = (synapse *)epnt;
         found = 1;
	 break;
     }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsyn(int ct, int cn, int nod, int ct2)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Limit to ct, cn, nod, ct2, cn2. */
/* Works with cn, nod = -1 */

{
   int i,n,ncn,nnod,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, nod, NULL,&ncn,&nnod); epnt=epnt->next) {
     if (epnt->node1a==ct  && epnt->node1b==ncn && epnt->node1c==nnod && 
         epnt->node2a==ct2) {
         closesyn = (synapse *)epnt;
         found = 1;
	 // fprintf (stderr,"syn elem  %d %d %d %d %d\n",ct,ncn,nnod,ct2,closesyn->elnum);
	 break;
     }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsyn(int ct, int cn, int nod, int ct2, int cn2)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Limit to ct, cn, nod, ct2, cn2. */
/* Works with cn, nod, cn2 = -1 */

{
   int i,n,ncn,nnod,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, nod, NULL,&ncn,&nnod); epnt=epnt->next) {
     if (cn2<0) { 
        if (epnt->node1a==ct  && epnt->node1b==ncn && epnt->node1c==nnod && 
            epnt->node2a==ct2) {
            closesyn = (synapse *)epnt;
            found = 1;
	    break;
        }
     } else {
        if (epnt->node1a==ct  && epnt->node1b==ncn && epnt->node1c==nnod && 
            epnt->node2a==ct2 && epnt->node2b==cn2) {
            closesyn = (synapse *)epnt;
            found = 1;
	    break;
        }
     }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

int findsynlocx(int ct, int cn, double xoffset, double yoffset)

/* Procedure to find a node in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   int i,n,nod,found,closenod;
   double dx, dy, dist, mindist;
   elem *epnt;

  closenod = found = 0;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1, NULL,NULL,&nod); epnt=epnt->next) {
     if (epnt->node1a==ct && epnt->node1b==cn) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closenod = nod;
         found = 1;
       }
    }
  }
  if (!found) closenod = -1;
  return closenod;
}

/*-------------------------------------------------*/

synapse *findsynloc(int ct, double xoffset, double yoffset)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
{
   int i,n,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, -1, -1); epnt=epnt->next) {
     if (epnt->node1a==ct) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = (synapse *)epnt;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynloc(int ct, int cn, double xoffset, double yoffset)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset). */
{
   int i,n,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1); epnt=epnt->next) {
     if (epnt->node1a==ct && epnt->node1b==cn) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = (synapse *)epnt;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

int findsynlocn(int ct, int cn, double xoffset, double yoffset)

/* Find synapse in a cell, return presynaptic node */
{
   synapse *s;
  s = findsynloc(ct, cn, xoffset, yoffset);
  if (s==NULL) return 0;
  else return s->node1c;
}

/*-------------------------------------------------*/

synapse *findsynloc(int ct, int cn, double xoffset, double yoffset, double vrev, double maxdist)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
{
   int i,n,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1); epnt=epnt->next) {
     s = (synapse *)epnt;
     if (epnt->node1a==ct && epnt->node1b==cn && (abs(s->vrev-vrev) < 0.001)) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = s;
         found = 1;
       }
    }
  }
  if (mindist > maxdist) closesyn = NULL;
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynloc(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset). */
{
   int i,n,nod,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, ct2, cn2, 0); epnt=epnt->next) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = (synapse *)epnt;
         found = 1;
       }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynloc(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
/* limit to ct, cn, ct2, cn2 */

{
   int i,n,nod,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, ct2, cn2, 0); epnt=epnt->next) {
     s = (synapse *)epnt;
     if (abs(s->vrev-vrev) < 0.001) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = s;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynloc(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev, double maxdist)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset), but closer than maxdist. */
/* limit to ct, cn, ct2, cn2 */

{
   int i,n,nod,found;
   double dx, dy, dist, mindist;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, ct2, cn2, 0); epnt=epnt->next) {
     s = (synapse *)epnt;
     if (abs(s->vrev-vrev) < 0.001) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = s;
         found = 1;
       }
    }
  }
  if (mindist > maxdist) closesyn = NULL;
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynloca(int ct, double roffset, double theta)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;
   
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynloc(ct, xoffset, yoffset);
}

/*-------------------------------------------------*/

synapse *findsynloca(int ct, int cn, double roffset, double theta)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;
   
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynloc(ct, cn, xoffset, yoffset);
}

/*-------------------------------------------------*/

synapse *findsynloca(int ct, int cn, int ct2, int cn2, double roffset, double theta)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;
   
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynloc(ct, cn, ct2, cn2, xoffset, yoffset);
}

/*-------------------------------------------------*/

synapse *findsynloca(int ct, int ct2, int cn2, double roffset, double theta)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) from (0,0). */
{
  return findsynloca(ct, -1, ct2, cn2, roffset, theta);
}

/*-------------------------------------------------*/

synapse *findsynlocr(int ct, int cn, double xoffset, double yoffset)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to (xoffset, yoffset) relative to its soma. */
/* Works with cn = -1 */
{
   int i,n,ncn,nod,found;
   double dx, dy, dist, mindist;
   node *somanode;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1, NULL,&ncn,&nod); epnt=epnt->next) {
     if (epnt->node1a==ct && epnt->node1b==ncn) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = (synapse *)epnt;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynlocr(int ct, int cn)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find synapse closest to its soma. */
{
    double xoffset, yoffset;

    return findsynlocr(ct, cn, xoffset=0.0, yoffset=0.0);
}

/*-------------------------------------------------*/

synapse *findsynlocr(int ct, int cn, double xoffset, double yoffset, double vrev)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
/* Works with cn = -1 */
{
   int i,n,ncn,nod,found;
   double dx, dy, dist, mindist;
   node *somanode;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1, NULL,&ncn,&nod); epnt=epnt->next) {
     s = (synapse *)epnt;
     if (epnt->node1a==ct && epnt->node1b==ncn && (abs(s->vrev-vrev) < 0.001)) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = s;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynlocr(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
/* limit to ct, cn, ct2, cn2 */
/* Works with cn = -1 */

{
   int i,n,ncn,nod,found;
   double dx, dy, dist, mindist;
   node *somanode;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, ct2, cn2, 0); epnt=epnt->next) {
     s = (synapse *)epnt;
     if (abs(s->vrev-vrev) < 0.001) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = s;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynlocr(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev, double maxdist)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
/* limit to ct, cn, ct2, cn2 */
/* Works with cn = -1 */
{
   int i,n,nod,found;
   double dx, dy, dist, mindist;
   node *somanode;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, ct2, cn2, 0); epnt=epnt->next) {
     s = (synapse *)epnt;
     if (abs(s->vrev-vrev) < 0.001) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = s;
         found = 1;
       }
    }
  }
  if (mindist > maxdist) closesyn = NULL;
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynlocra(int ct, int cn, double roffset, double theta)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocr(ct,cn,xoffset,yoffset);
}

/*-------------------------------------------------*/

synapse *findsynlocra(int ct, double roffset, double theta)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocr(ct,-1,xoffset,yoffset);
}

/*-------------------------------------------------*/

synapse *findsynlocra(int ct, int cn, double roffset, double theta, double vrev)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocr(ct,cn,xoffset,yoffset,vrev);
}

/*-------------------------------------------------*/

synapse *findsynlocra(int ct, int cn, int ct2, int cn2, double roffset, double theta, double vrev)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (roffset, theta) from (0,0). */
/* limit to ct, cn, ct2, cn2 */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocr(ct,cn,ct2,cn2,xoffset,yoffset,vrev);
}
/*-------------------------------------------------*/

synapse *findsynlocra(int ct, int cn, int ct2, int cn2, double roffset, double theta, double vrev, double maxdist)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (roffset, theta) from (0,0). */
/* limit to ct, cn, ct2, cn2 */
{
   double xoffset, yoffset;

  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocr(ct,cn,ct2,cn2,xoffset,yoffset,vrev,maxdist);
}

/*-------------------------------------------------*/

synapse *findsynlocz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to absolute loc (xoffset, yoffset). */
/* Works with cn = -1 */

{
   int i,n,nod,ncn,found;
   double dx, dy, dist, mindist;
   node *npnt;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1, NULL,&ncn,&nod); epnt=epnt->next) {
     npnt = epnt->nodp1;
     if ((epnt->node1a==ct && epnt->node1b==ncn) && 
        (zmin <= npnt->zloc && npnt->zloc <= zmax)) {
       dx = npnt->xloc - xoffset;
       dy = npnt->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = (synapse *)epnt;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynlocrz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) relative to its soma. */
/* Works with cn = -1 */

{
   int i,n,nod,ncn,found;
   double dx, dy, dist, mindist;
   node *somanode, *npnt;
   elem *epnt;
   synapse *s,*closesyn;

  found = 0;
  closesyn = NULL;
  mindist = 1e10;
  if ((somanode=nde(ct,cn,soma)) != NULL) {
     xoffset += somanode->xloc; 
     yoffset += somanode->yloc;
  } 
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1, NULL,&ncn,&nod); epnt=epnt->next) {
     npnt = epnt->nodp1;
     if ((epnt->node1a==ct && epnt->node1b==ncn) && 
        (zmin <= npnt->zloc && npnt->zloc <= zmax)) {
       dx = npnt->xloc - xoffset;
       dy = npnt->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closesyn = (synapse *)epnt;
         found = 1;
       }
    }
  }
  if (!found) closesyn = NULL;
  return closesyn;
}

/*-------------------------------------------------*/

synapse *findsynlocaz(int ct, int cn, double roffset, double theta, double zmax, double zmin)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;
   
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocz(ct, cn, xoffset, yoffset, zmax, zmin);
}

/*-------------------------------------------------*/

synapse *findsynlocaz(int ct, double roffset, double theta, double zmax, double zmin)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   int cn;
   double xoffset, yoffset;
 
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocz(ct, cn=-1, xoffset, yoffset, zmax, zmin);
}

/*-------------------------------------------------*/

synapse *findsynlocarz(int ct, int cn, double roffset, double theta, double zmax, double zmin)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   double xoffset, yoffset;
   
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocrz(ct, cn, xoffset, yoffset, zmax, zmin);
}

/*-------------------------------------------------*/

synapse *findsynlocarz(int ct, double roffset, double theta, double zmax, double zmin)

/* Procedure to find a synapse in a cell with a synapse to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0). */
{
   int cn;
   double xoffset, yoffset;
   
  xoffset = roffset *  cos (theta*PI/180.0);
  yoffset = roffset * -sin (theta*PI/180.0);
  return findsynlocrz(ct, cn=-1, xoffset, yoffset, zmax, zmin);
}

/*-------------------------------------------------*/

int findnodlocz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin)

/* Procedure to find a node in cell to record from. */
/* Find node closest to (xoffset, yoffset) from (0,0) that is within zmax, zmin. */
{
   int i,n,nod,found,closenod;
   double dx, dy, dist, mindist,temp;
   node *npnt;

  if (zmax<zmin) { temp = zmax; zmax = zmin; zmin = temp; }
  closenod = found = 0;
  mindist = 1e10;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, cn, -1, NULL,NULL,&nod); npnt=npnt->next) {
     if (zmin <= npnt->zloc && npnt->zloc <= zmax) {
       dx = npnt->xloc - xoffset;
       dy = npnt->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closenod = nod;
         found = 1;
       }
    }
  }
  if (!found) closenod = 0;
  return closenod;
}

/*-------------------------------------------------*/

int findsynlcpre(int ct, int cn, int ct2, int cn2)

/* Procedure to find the presynaptic node for synapse from ct,cn, to ct2, cn2. */

{
   int prenode,found;
   elem *epnt;

  found = 0;
  prenode = 0;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1); epnt=epnt->next) {
     if ((epnt->node1a==ct  && epnt->node1b==cn) && 
         (epnt->node2a==ct2 && epnt->node2b==cn2)) {
         prenode = ((synapse *)epnt)->node2c;
         found = 1;
    }
  }
  if (!found) prenode = 0;
  return prenode;
}

/*-------------------------------------------------*/

int findsynlocp(int ct, int cn, int ct2, int cn2)

/* Procedure to find the postsynaptic node for synapse from ct,cn, to ct2, cn2. */

{
   int postnode,found;
   elem *epnt;

  found = 0;
  postnode = 0;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1); epnt=epnt->next) {
     if ((epnt->node1a==ct  && epnt->node1b==cn) && 
         (epnt->node2a==ct2 && epnt->node2b==cn2)) {
         postnode = ((synapse *)epnt)->node2c;
         found = 1;
    }
  }
  if (!found) postnode = 0;
  return postnode;
}

/*-------------------------------------------------*/

int findsynlocpa(int ct, int cn, double xoffset, double yoffset, int ct2, int cn2)

/* Procedure to find the postsynaptic node for synapse from ct,cn, to ct2, cn2. */
/* First, find node in ct,cn closest to (xoffset, yoffset) that has a synaptic output. */
/* Then return the postsynaptic node in ct2,cn2. */

{
   int found;
   double dx, dy, dist, closenod, mindist;
   elem *epnt;

  found = 0;
  closenod = 0;
  mindist = 1e10;
  for (epnt=elempnt; epnt=foreach(epnt, SYNAPSE, ct, cn, -1); epnt=epnt->next) {
     if ((epnt->node1a==ct  && epnt->node1b==cn) && 
         (epnt->node2a==ct2 && epnt->node2b==cn2)) {
       dx = epnt->nodp1->xloc - xoffset;
       dy = epnt->nodp1->yloc - yoffset;
       dist = sqrt(dx*dx+dy*dy);
       if (mindist>dist) {
         mindist = dist;
         closenod = ((synapse *)epnt)->node2c;
         found = 1;
       }
    }
  }
  if (!found) closenod = 0;
  return closenod;
}

/*-------------------------------------------------*/

int findmidc (int ct, int offset)

/* find middle cell by count */

{
   int i,n,cn,found,midcn,midn;
   node *npnt;

  found = n = 0;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
     n++;
  }
  midn = n/2 + offset;
  i = found = 0;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
    if (i==midn) {
      midcn = cn;
      found = 1;
      break;
    }
    i++;
  }
  return midcn;
}

/*-------------------------------------------------*/

int findnext (int ct, int celcnt)

/* find next cell by count */

{
   int i,cn,found;
   node *npnt;

  i = found = 0;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, NULL,&cn,NULL); npnt=npnt->next) {
    if (i==celcnt) {
      found = 1;
      break;
    }
    i++;
  }
  if (!found) cn = -1;
  return cn;
}

/*-------------------------------------------------*/

int find_gtconn(int ct, int nconn)

/* procedure to find a cell in array to record from */

{
   int i,cn,found,retcn;
   node *npnt;

  found = 0;
  for (npnt=nodepnt; npnt=foreach(npnt, ct, -1, soma, -1, NULL,&cn,NULL,NULL); npnt=npnt->next) {
     if (tot_ncel_ind(ct,cn)==nconn) {
       retcn = cn;
       found = 1;
       break;
     }
  }
  if (!found) {
     fprintf (stderr,"find_gtconn: can't find %s cell with %d inputs.\n",cname[ct], nconn);
     retcn = 0;
  }
  return retcn;
}


/*-------------------------------------------------*/


void find_maxmin(int celltype,int cellnum)

 /* Find size of cell arrays; if celltype=-1, look at all celltypes, */
 /*  if cellnum = -1, look at all cells of that type */

{
    int ct,cn,n;
    double xl, yl, zl, rad;
    elem *epnt;
    node *npnt;

    xmax = -1e6;
    xmin =  1e6;
    ymax = -1e6;
    ymin =  1e6;
    zmax = -1e6;
    zmin =  1e6;
    xmaxnode = -1;
    xminnode = -1;
    ymaxnode = -1;
    yminnode = -1;
    zmaxnode = -1;
    zminnode = -1;
    ct=cn=0;

    /* find existing dendritic tree, ignore axons, synapses, etc */

    for (epnt=elempnt; epnt=foreach(epnt, ELEMENT, celltype, cellnum, -1, -1, 
					      NULL,	NULL,	 &n, NULL); epnt=epnt->next) {  

       if (epnt->ctype  == CABLE || epnt->ctype == SPHERE ||
                      epnt->ctype  == ROD || epnt->ctype == CONE) {

 	 ct = epnt->node1a;
 	 cn = epnt->node1b;

	 if (!bptype(ct) && 
			     (
			     (epnt->region==HILLOCK) ||
             		     (epnt->region==AXON_THIN) || 
			     (epnt->region==AXON_PROX) || 
			     (epnt->region==AXON_DIST) 
			     )) continue;

        /* if ((ct==xcone) || (ct==xrod)) printf ("xcone %d, xrod %d, ct %d\n",xcone,xrod,ct); */
          npnt = epnt->nodp1;
	  xl = npnt->xloc;
	  yl = npnt->yloc;
	  zl = npnt->zloc;
          if (epnt->ctype  == ROD || epnt->ctype == CONE) {
            zl += 10;
	  }
	  if (xl < xmin) { xmin = xl; xminnode = n; };
	  if (xl > xmax) { xmax = xl; xmaxnode = n; };
	  if (yl < ymin) { ymin = yl; yminnode = n; };
	  if (yl > ymax) { ymax = yl; ymaxnode = n; };
	  if (zl < zmin) { zmin = zl; zminnode = n; };
	  if (zl > zmax) { zmax = zl; zmaxnode = n; };
	  if (epnt->ctype==SPHERE) {	// check radius 
	     rad = ((sphere*)epnt)->dia*0.5;
	      if ((xl-rad) < xmin) { xmin = xl-rad; xminnode = n; };
	      if ((xl+rad) > xmax) { xmax = xl+rad; xmaxnode = n; };
	      if ((yl-rad) < ymin) { ymin = yl-rad; yminnode = n; };
	      if ((yl+rad) > ymax) { ymax = yl+rad; ymaxnode = n; };
	      if ((zl-rad) < zmin) { zmin = zl-rad; zminnode = n; };
	      if ((zl+rad) > zmax) { zmax = zl+rad; zmaxnode = n; };
	  }
	  if (epnt->ctype==CABLE) {	// check second node
            npnt = epnt->nodp2;
	    xl = npnt->xloc;
	    yl = npnt->yloc;
	    zl = npnt->zloc;
	    if (xl < xmin) { xmin = xl; xminnode = n; };
	    if (xl > xmax) { xmax = xl; xmaxnode = n; };
	    if (yl < ymin) { ymin = yl; yminnode = n; };
	    if (yl > ymax) { ymax = yl; ymaxnode = n; };
	    if (zl < zmin) { zmin = zl; zminnode = n; };
	    if (zl > zmax) { zmax = zl; zmaxnode = n; };
	  }
       }
    }
    if (xmax < xmin)
    {  xmax = 1;
       xmin = 0;
    }
    if (ymax < ymin)
    {  ymax = 1;
       ymin = 0;
    }
}

/*-------------------------------------------------*/

double find_maxrad(int ct,int cn)

 /* proc to find size of cell arrays; if celltype=-1, look at all celltypes */

{
    int n;
    double xd, yd, maxrad;
    double xl, yl;
    double xmax, xmin;
    double ymax, ymin;
    elem *epnt;
    node *npnt;

    xmax = -1e6;
    xmin =  1e6;
    ymax = -1e6;
    ymin =  1e6;

    /* find existing dendritic tree, ignore photoreceptors and axons */

    for (epnt=elempnt; epnt=foreach(epnt, ELEMENT, ct, cn, -1, -1, 
					      NULL,NULL,&n,NULL); epnt=epnt->next) {  
        if (bptype(ct) ||  ((epnt->region!=HILLOCK) &&
             		   (epnt->region!=AXON_THIN) && 
			   (epnt->region!=AXON_PROX) && 
			   (epnt->region!=AXON_DIST)))

       {
          npnt = ndn(ct,cn,n);
	  xl = npnt->xloc;
	  yl = npnt->yloc;
	  if (xl < xmin) xmin = xl;
	  if (xl > xmax) xmax = xl;
	  if (yl < ymin) ymin = yl;
	  if (yl > ymax) ymax = yl;
       };
    };
    xd = xmax - xmin;
    yd = ymax - ymin;
    maxrad = (xd + yd) / 2;
    return maxrad;
}

/*-------------------------------------------------*/

double qfact(double q10) {
  return (exp(log(q10)*(tempcel-22)/10));
}

/*-------------------------------------------------*/

void makanatfile (int celltype,int cellnr) 

/* make an anatomy file out of an existing artificial cell morphology */

{
    int d,n;
    int node1, node2, dendcount, dendnr;
    const char *region;
    double xpos, ypos, zpos, diam ;
    elem *epnt;
    node *npnt;
    static int dendnodes[30][100];	/* hold the nodes belonging to each dend */

  ncfprintf (stdout, "#  node  parent   dia     xbio    ybio     zbio       region  dendnr\n");

  /* initialise local array dendnodes */

  for (d=0;d<30;d++)
    for (n=0;n<100;n++)
      dendnodes[d][n] = -1;

  dendcount=0; /* count nr of dendrites */

  /* go through all cables, for every node print out parameters for amacfile */

  for (epnt=elempnt; epnt=foreach (epnt, ELEMENT); epnt=epnt->next) {
   
   node1 = epnt->node1c; 	/* descendent node */
   node2 = epnt->node2c; 	/* parent node */

   /* soma has neg nodenr as parent, call the rest dend */

   if (node2<0){
   	region="SOMA";
	node2=0;
	dendnr=0;
   }
   else region="DEND";

   /* count nr of dendrites leaving soma, make array with parents for every dend nr */

   if (node2==0) {  	/* if parent is soma, give child its own dendr-nr. */
     dendnodes[dendcount][0] = node1;
     dendnr=dendcount;
     dendcount++;
   }
   else  /* go through dendnodes, find parent, use its dendnr, add node to arr */
     for(d=0;d<dendcount;d++)
       for (n=0;n<100;n++) {
         if (dendnodes[d][n]>-1 && dendnodes[d][n]==node2) {
	   dendnr=d;
	   dendnodes[d][n+1]=node1;
	   break;
	};
     };

   npnt = ndn(am,cellnr,node1);
   xpos = npnt->xloc;
   ypos = -1* (npnt->yloc);  /* -1 is needed.. */
   zpos = npnt->zloc;
   diam = ((cable *)epnt)->dia;
   zpos = -5;

   fprintf (stdout, "%7d %6d   %-7.3g %-8.4g %-8.4g %-8.4g   ",
   		node1, node2, diam, -1*xpos, -1*ypos, zpos);
   ncfprintf(stdout, "%-s  ", region);
   ncfprintf(stdout, "%6d\n", dendnr);
  }
}
