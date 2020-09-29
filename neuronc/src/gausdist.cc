
#include <stdio.h>
#include "ncfuncs.h"
#include "ncinit.h"

main (int argc, char **argv)

{
    int n;
    double m, s;
    int i,c;

   setptr("m",      &m);
   setptr("n",      &n);
   setptr("s",      &s);
   setptr("c",      &c);
   ncinit(argc,argv);
   m = 100;
   n = 250000;
   n = 2497580;
   s = 10;
   c = 2;
   setvar();

  for (i=0; i<n; i++) {

     if (c==1) printf ("%g\n",gasdev()*s + m);
     else      printf ("%d %g\n",i, gasdev()*s + m);

  }
}
