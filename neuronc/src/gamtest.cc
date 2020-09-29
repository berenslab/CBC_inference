
#include <stdio.h>


/* Use these definitions for "random.cc" (from glibc2 library) */

double gamdev(int a);
double gamdev2(double a);

int cumrand;
int rseed = 1235938;

main()

{
   long int i,n,x,x1;
   double v;
   double a,b;

  n = 100000;
  a = 1e4;
  b = 1.0;
  for (i=0; i<n; i++) {
     v=gamdev2(a);
     //v=gamdev((int)a);
     printf ("%g\n",v);
  }
}

