
#include <stdio.h>
#include <math.h>

#include "ncfuncs.h"
#include "digfilt.h"

int rseed=333221;
int cumrand=0;

double tau, tstep;

/*------------------------------------*/

main (int argc, char **argv)

{
  double t, exptdur, scontrast, meanval, val;

  exptdur = 1;

  meanval = 1;
  tau = 0.01;
  tstep = 0.001;
  scontrast = 0.3;

  init_digfilt  (meanval, tau, tstep);
  init_digfilt2 (meanval, tau, tstep);

  for (t=0; t<exptdur; t+= tstep){
        double start, dur, inten;

    val = gasdev()*scontrast + meanval;
    //inten = digfilt2(digfilt(val));
    inten = digfilt(val);
    if (inten < 0) inten = 0;
    printf ("%g %g %g\n",t,val,inten);
  }
}

