
#include <stdio.h>
#include "ncinit.h"
#include "ncfuncs.h"

double data[] = { 
	.07, .13, .19, .26, .32, .38, .44, .51,
        .57, .63, .69, .76, .82, .88, .94,

	.24, .35, .43, .49, .55, .61, .66, .71,
        .75, .79, .83, .87, .90, .94, .97 };

#define DATASIZ 15

double p[] = { 1, 1, 1 };       // use any starting value except { 0,0,0 }

double fit_func(double x, double *p, int n)
{
    //if (p[0] > 5) p[0] = 5; 
    return (p[0] * x + (1 - p[0] + p[1] + p[2]) * x * x) /
        (1 + p[1] * x + p[2] * x * x);
}

int info;

main()

{
    // perform the fit:

    lmfit (fit_func,DATASIZ,&data[0],3,&p[0]);

   fprintf (stderr,"%g %g %g\n",p[0],p[1],p[2]);
}
