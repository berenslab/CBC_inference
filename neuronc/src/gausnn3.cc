/* Program Gausnn */

/* Makes a gaussian nearest-neighbor distance */

/*  S.S.Kumar and R.G.Smith */

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define MAXPTS 2000
#define NUMN 10

#define FRAMEX  200  /* Frame size in microns; X-dimension */
#define FRAMEY  200  /* Frame size in microns; Y-dimension */
#define BORDERSIZ 0  /* Border size in fraction of frame */

#define NUMCELLS 400

/* #define M_PI 3.14159265358979323846264  /* */

double density = ((double)NUMCELLS) / (FRAMEX*FRAMEY);
double gmean = 10.0;
double gstdev = 1.0;
double gms = 10.0;
double cmean = 0;
double cstdev = 0;
int initcells = 8;

int num_cells = NUMCELLS;
int num_bins = 0;
int num_reject=MAXPTS;
int densityfl = 0;
int meanfl = 0;
int msflag = 0;
int initfl = 0;
int numfl = 0;
int nnfl  = 0;			/* =0 -> nearest neighbor mean, stdev */
				/* =1 ->  avg.   neighbor mean,stdev */
int printfl = 0;
int xfrflag = 0;
int yfrflag = 0;
long int rseed = 1234;

double framex  = FRAMEX;
double framey  = FRAMEY;
double drand();
double nndist(int i, double *xv, double *yv, int n);
double packing(double reg);
double elastprob(int i, double *xv, double *yv, int n);

double xval[MAXPTS] = {0};
double yval[MAXPTS] = {0};

FILE *ncstdout = stdout;
FILE *ncstderr = stderr;

FILE *textout;
FILE *outfil;
FILE *fp;

/* -------------------------------------------------------------- */

main(int argc, char **argv)
{
   char *cptr;
   FILE *freopen(const char *, const char *, FILE *);
         
 outfil = NULL;
 textout = ncstdout;
 if (argc==1)                   /* if user needs help */
   run();
 else
 do                                     /* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {
     
          case 'd': 
                argv++; argc--;
                density = atof(*argv);
                if (density==0.0) density=1e-6;
                densityfl = 1;
                break;

          case 'i': 
                argv++; argc--;
                initcells = atof(*argv);
                initfl = 1;
		break;

          case 'm': 
                argv++; argc--;
                gmean = atof(*argv);
                meanfl = 1;
                break;

          case 'n': 
                argv++; argc--;
                num_cells = atof(*argv);
                numfl = 1;
                break;

          case 'p': 
                printfl = !printfl;
                break;

          case 'r': 
                argv++; argc--;
                rseed = atof(*argv);
                break;

          case 's': 
                argv++; argc--;
                gstdev = atof(*argv);
                break;

          case 't': 
                argv++; argc--;
                gms = atof(*argv);
		msflag = 1;
                break;

          case 'x': 
                argv++; argc--;
                framex = atof(*argv);
                xfrflag = 1;
                break;

          case 'y': 
                argv++; argc--;
                framey = atof(*argv);
                yfrflag = 1;
                break;

          default:
                ncfprintf (ncstderr,"gausnn: unknown switch '%s'\n",*argv);
                exit();

        }  /* switch */
      }    /* if */
     else
      {
       if((outfil=fopen(cptr,"w"))==NULL)
         {
           ncfprintf(ncstderr,"gausnn: cannot open %s\n",cptr);
           fflush (stderr);
           continue;
         }
       run();
       if (argc <= 1) break;
     }
    }
   else run();
  }
 while (argc > 0);
}

/* -------------------------------------------------------------- */

run(void)
{

   int i,j,numpts;
   long int tries,tottries,itercrit,maxiter;
   long int iter,slowcount;
   double p,okdist,mdist,mstd;
   double nnd[NUMN];
   int nn[NUMN],varn[4];
   double sum[4],sumsq[4],val,area,ms;
   double tvar=0.0, tstdev=0.0,tmean=0.0;
   double tstdevx=0.0, tmeanstd=0.0;
    
/* MAXITER is how many incorrect points per cell is not allowable.  */
/* ITERCRIT is how many incorrect points per cell is allowable,
    as long as the next cell is found quickly.  */

#define MAXITER 1.5
#define ITERCRIT .15

  setrand(rseed);		/* initialize random number generator */

  if (xfrflag && !yfrflag) {
	framey = framex;
  }
  area = framex * framey;

  if (numfl) {			/* number specified overrides everything */
    density = num_cells / area;
    if (!meanfl) gmean = sqrt (1.0/density);  
  }
  else if (densityfl) {		/* density specified overrides mean */
    num_cells = density * area;
    if (!meanfl) gmean = sqrt (1.0/density);  
  }
  if (msflag) {			/* "mean/stdev" ratio specified */
    ms = gms;
    if (!meanfl) gmean = sqrt (1.0/density) * packing(ms);  
    gstdev = gmean / ms;
  }
  else if (meanfl) {
     ms = gmean / gstdev;
     density = packing(ms)*packing(ms)/(gmean*gmean);
     num_cells = density * area;
 /*    if (ms<6.7) gstdev *= pow((10/ms-.45),1.4)*.9;		/* */
     if (ms<10) gstdev *= pow((10/ms),1.3)*.75;		/* */
  }
/*  if (!meanfl) gmean *= packing(ms);  			/* */
/*  if (!numfl) num_cells *= packing(ms)*packing(ms);	/* */

  ms = gmean / gstdev;
  cstdev = gstdev * 2.0; 		/* stdev created by this alg is low  */
  cmean = gmean + gstdev*2.0; 		/* nnd is always less than 2nnd */

  if (num_cells > MAXPTS){
	ncfprintf (ncstderr,"Gausnn: too many cells: %d\n",num_cells);
	num_cells = MAXPTS;
  }
ncfprintf (textout,"Gausnn,  Sept 1992:\n\n");

ncfprintf (textout,"Array size:       %-4g x %-4g\n",framex,framey);
ncfprintf (textout,"Area:             %-6.4g\n",area);
ncfprintf (textout,"Target num cells: %-6d\n",num_cells);
ncfprintf (textout,"Target density:   %-6.4g\n",density);
ncfprintf (textout,"Target mean dist: %-6.4g\n",gmean);
ncfprintf (textout,"Target stdev:     %-6.4g\n",gstdev);

  if (outfil==NULL) {
     fp = fopen ("gausn.dat","w");
     outfil = fp;
  }
  tottries = 0;
  numpts = 0;
  slowcount = 0;
 
  for (i=0; i<num_cells; i++) {		/* zero the x,y locations */
        xval[i] = 0;
        yval[i] = 0;
  }

  i=0;
  xval[i] = framex/2.0; 		/* make one point in center */
  yval[i] = framey/2.0; 
  sendxy (i,xval[i],yval[i]);
  numpts++;
					/*  make some initial seed points */
					/*  but not too close */
  if (!initfl) {
    if (ms>6) initcells = min(4,initcells);  /* */
    if (ms>10) initcells = 1; 		/* only one start if very regular */
  }

  okdist = 4.0*gmean;            
  for (i=numpts; i<initcells; numpts++,i++) { 
     for (mdist=0,tries=0; mdist<okdist && tries<100; tries++) { 
        tottries++;
        xval[i] = drand() * framex; 
        yval[i] = drand() * framey; 
        mdist = nndist(i,xval,yval,numpts);
     }
     sendxy(i,xval[i],yval[i]);
  }

		/* make maxiter proportional to number of cells */
		/*   and square of regularity */

  mstd = min (ms,20);
  maxiter  = MAXITER  * num_cells * ms * mstd;
  itercrit = ITERCRIT * num_cells * ms * mstd;

  for (i=numpts; i<num_cells; numpts++,i++) {
     for (iter=0; iter<maxiter; iter++) {
        tottries++;
        xval[i] = drand() * framex; 
        yval[i] = drand() * framey; 
        p = elastprob(i,xval,yval,numpts);
	/* ncfprintf (ncstderr,"n %d p %g\n",numpts,p); /* */
        if (p > drand()) break; 
     }
     if (iter >= maxiter) {
        if (printfl) 
	 ncfprintf(textout,"Number of tries (%d) too many, stopping...\n",iter);
	   break;
     }
     sendxy(i,xval[i],yval[i]);
     if (printfl) ncfprintf (textout,"iter %d\n",iter);
     if (iter >= itercrit ) {
	slowcount++;
        if (printfl) 
           ncfprintf(textout,"Incrementing slowcount: %d. \n",slowcount);
	
	if (slowcount >= 3) {
          if (printfl) 
           ncfprintf(textout,"Slowing down on last %d tries, stopping.\n",
				slowcount);
           numpts++;
	   break;
	}
     }
     else {
	if (slowcount>0) slowcount--;
     }

  }  /* for (i;;) */

  for (j=0; j<1; j++) {			/* clear arrays for stdev calc. */
     sum[j] = 0;
     sumsq[j] = 0;
  }

  for (i=0; i<numpts; i++) { 
    nndistm(i,xval,yval,numpts,nnd,nn);		/* find nearest 4 neighbs */
    for (j=0; j<2; j++) {
      val = nnd[j];
      sum[j] += val;
      sumsq[j] += val*val;
    }
  } 
  varn[0] = numpts;
  varn[1] = numpts;
/*  varn[1] += varn[0];			/* */
/*  sum[1]  += sum[0];			/* add nearest to next neighbor dist */
/*  sumsq[1] += sumsq[0];			/* */
  if (varn[j] <= 1) varn[j] = 2;
  for (j=0; j<2; j++) {
    tmean = sum[j] / varn[j];
    tvar = (sumsq[j] - (sum[j] * tmean)) / (varn[j] - 1);
    tstdev = sqrt (tvar);
    tstdevx = tstdev;
    if (tstdevx == 0.0) tstdevx = 1.0;
    tmeanstd = tmean / tstdevx;
    ncfprintf (textout,"nn #%d:  ",j+1);
    ncfprintf (textout,"n %d mean %g stdev %g m/s ratio %g\n",
		varn[j],tmean,tstdev,tmeanstd);
  }
 ncfprintf (textout,"Density %g Area %g\n", numpts/area, area);
 ncfprintf (textout,"Total tries %d\n",tottries);
}

/* -------------------------------------------------------------- */

double packing (double reg)
              

/* Find approximate packing radius.
   Weight sqrt ( 2 / sqrt(3)) == 1.07457 for triangular packing
   when regularity is high.  When packing is more
   square or random, reduce packing radius.
*/

{
   double a,rfactor;


 reg = max(reg,1);
 a = 2/reg;
 if (a<0) a = 0;
 else if (a>1) a=1;
 rfactor = (a*0.70) + ((1-a) * 1.05);
 return (rfactor);    
}

/* -------------------------------------------------------------- */

sendxy(int n, double x, double y)
             

/* print a cell's position */

{

    ncfprintf(outfil,"%g  %g\n",x,y);
    fflush (outfil);
    if (printfl) ncfprintf (textout,"N %d\n",n+1);
}

/* -------------------------------------------------------------- */

double gauss(double x, double mu, double sigma)
{
double value,r;
        r = (x-mu) / sigma;
        value = exp(-r*r);
        return(value);
}

/* ---------------------------------------------------------- */

nndistm(int i, double *xv, double *yv, int n, double *nnd, int *nn)
                    
            
               
            

/* Find nearest 4 neighbors, and their distances from point i */

{
int j;
double dist,xt,yt,md1,md2,md3,md4;
int n1=0,n2=0,n3=0,n4=0;

        for (md4=md3=md2=md1=1e10,j=0; j<n; j++){
                if (i!=j) {
                   xt = xv[i]-xv[j];
                   yt = yv[i]-yv[j];
                   dist = xt*xt + yt*yt;        
                   if (dist < md1) {
                         md4 = md3;
                         md3 = md2;
                         md2 = md1;
                         md1 = dist;
                         n4  = n3;
                         n3  = n2;
                         n2  = n1;
                         n1  = j;
                   }
                   else if (dist < md2) {
                         md4 = md3;
                         md3 = md2;
                         md2 = dist;
                         n4  = n3;
                         n3  = n2;
                         n2  = j;
                   }
                   else if (dist < md3) {
                         md4 = md3;
                         md3 = dist;
                         n4  = n3;
                         n3  = j;
                   }
                   else if (dist < md4) {
                         md4 = dist;
                         n4  = j;
                   }

                 }
        }
   nnd[3] = sqrt(md4);
   nnd[2] = sqrt(md3);
   nnd[1] = sqrt(md2);
   nnd[0] = sqrt(md1);
   nn[3]  = n4;
   nn[2]  = n3;
   nn[1]  = n2;
   nn[0]  = n1;
}

/* ---------------------------------------------------------- */

double nndist(int i, double *xv, double *yv, int n)
                    
            

/* find nearest neighbor distance */

{
int j;
double dist,xt,yt,min_dist;

        for (min_dist=1e10,j=0; j<n; j++){
                if (i!=j){
                   xt = xv[i]-xv[j];
                   yt = yv[i]-yv[j];
                   dist = xt*xt + yt*yt;        
                   if (dist < min_dist) {
                         min_dist = dist;
                   }
                }
        }
return sqrt(min_dist);
}

/* ---------------------------------------------------------- */

double elastprob(int i, double *xv, double *yv, int n)
                    
            

/* find probability of acceptance */

{

  int j;
  double prob;
  double nnd[NUMN];
  int nn[NUMN];

   nndistm(i,xv,yv,n,nnd,nn); /* find nearest 4 neighbs */
   for (prob=1.0,j=0; j<2; j++) {
         if (j>0 && nnd[j] > cmean*2+cstdev*2) continue; /* */
        /* if (j>0 && nnd[j] > 1000) continue; /* */
        prob *= gauss(nnd[j],cmean,cstdev);
  /* ncfprintf (ncstderr,"n %d j %d dist %g prob %g\n",n,j,nnd[j],prob);  /*  */
   }

/* for (prob=1.0,j=0; j<n; j++){
        if (i!=j) {
           xt = xv[i]-xv[j];
           yt = yv[i]-yv[j];
           dist = sqrt(xt*xt + yt*yt);  
           if (dist < cmean) {
                prob *= gauss(dist,cmean,cstdev);
 /* ncfprintf (ncstderr,"prob %g\n",prob);  /*  */
/*         }
        }
  }
*/

  return (prob);
}

/* ---------------------------------------------------------- */

