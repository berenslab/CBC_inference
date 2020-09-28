

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

main  (int argc, int argv)
{
double a;

a = pow (2.,3.);
fprintf (stderr,"pow(2,3) = %g\n",a);
}


