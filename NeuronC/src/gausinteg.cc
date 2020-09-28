
#include <stdio.h>
#include <math.h>
#include "ncinit.h"
#include "ncfuncs.h"

#define MAXINTBIN 20000
int intbin[MAXINTBIN];

void makbin (double val)
{

  if (val>MAXINTBIN) val = MAXINTBIN;
  if (val<0) val = 0;
  intbin[int(val)]++;
};

double gausfunc(double x, double m, double s)
{
  return (exp(-((x-m)*(x-m)/(s*s*2.0))));
}

main (int argc, char **argv)

{
    int n;
    double m,s,tot2=0;
    int i,tot;

   setptr("m",      &m);
   setptr("n",      &n);
   setptr("s",      &s);
   ncinit(argc,argv);
   m = 10000;
   n = 10000000;
   s = 100;
   setvar();

  //for (i=0; i<n; i++) {
  //   makbin (gasdev()*s+m);
  //}
  for (tot=i=0; i<MAXINTBIN; i++) {
    //tot += intbin[i];
    tot2 += gausfunc(i,10000,s);
    printf ("%d %g %g %g\n",i,((double)tot)/n,gausfunc(i,10000,s),tot2/(s*sqrt(2*M_PI)));
  } 
}
