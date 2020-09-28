
#include <stdio.h>
#include "ncinit.h"
#include "ncfuncs.h"

main (int argc, char **argv)

{
    int n;
    double m;
    double corr,s1,s2,s3,s4,s5;
    int i,t;

  setptr("m",      &m);
  setptr("n",      &n);
  setptr("t",      &t);
  setptr("corr", &corr);
  ncinit(argc,argv);
  m = 100;
  n = 100000;
  t = 0;
  corr = -1;
  setvar();

  if (corr > 1) corr = 1;

  for (i=0; i<n; i++) {
   if (corr>=0) {
     s1 = poisdev(m*corr); 
     s2 = poisdev(m*(1-corr));
     s3 = poisdev(m*(1-corr)); 
     s4 = poisdev(m); 
     s5 = poisdev(m); 
     s2 += s1;
     s3 += s1;
     if (t==0) printf ("%8.5g %8.5g %8.5g %8.5g\n",s2,s3,s4,s5);
     else      printf ("%d %8.5g %8.5g %8.5g %8.5g\n",i,s2,s3,s4,s5);
   }
   else {   
     if (t==0) printf ("%d\n",poisdev(m));
     else      printf ("%d %d\n",i, poisdev(m));
   }

  }
}
