/* program mglur6.cc   simulates second-messenger cascade */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>

char *strtok(char *s, const char *delim);
int strcmp(const char *s1, const char *s2);
char *strstr(const char *haystak, const char *needle);
double atof(char *s);
double log(double val);

#ifdef __cplusplus
}
#endif

/*------------------------------------*/

double calctau(double tau, double timestep)

/* tau input in msec, output is multiplier for digital filter */
/* M.v.Rossum, 1996 */

{
   double k,taugran;
   int sign;

#define TIMERES  (.1)
#define MINTAU  (TIMERES)
#define TAUGRAN (.632120558)  /* (1 - 1/e) */

 sign = (tau >= 0);
 if (tau<0) tau = -tau;
 if (tau < timestep) tau =  timestep;
 if (tau==0.0) tau = 1.0;

/*  taugran = timestep * pow(TAUGRAN,timestep/tau);  /* the old
 *  way */
/*  k = taugran / tau; /* */

 k = 1.0 - exp(-timestep/tau);          /* fix error for small
tau*/
 if (k > 1.0) k = 1.0;                  /* limit k to reasonable
values */
 else if (k<=0.0) k = .0001;
 return ( sign ? k : -k);
}

/*------------------------------------*/

main(void)
{
  int i,j,startfl,plotinc;
  double stiminc, startval;
  double qc, qt, qcond, qrate, tempcel;
  double sfactor, rrate;
  double cabuf, cycgbuf, loopgain;
  double pdebase, kdgcyc, kdrhodk;
  double rhod, rhodk, rhod2, rstar;
  double gpr, pde, gcyc, cycg, cond, totcond;
  double t1, t2, t3;
  double ca, cax, ca2;
  double lgain, rgain1, rgain2, rgain3, rkgain;
  double ggain, gcygain, pdegain, condgain, gca, cxgain;
  double decrstar, pdedec, deccond, capump, deccax;

 
 stiminc = .0001;

 qt = 37.0;             /* base temperature for channel rate & cond */
 qc = 1.4;              /* Q10 factor for conductance shape */
 tempcel = 35;

 qrate = exp(log(qc) * (tempcel - qt) / 10.0);

 
 sfactor = .5;               /* speedup of cascade & intensity */
 sfactor = .1;               /* speedup of cascade & intensity */
 rrate = stiminc * qrate /sfactor;    /* stiminc = timinc for syn/stimuli */

 cabuf    = 1.0;
 cycgbuf  = 1.0;
 loopgain = 1.0;

 pdebase = .04;
 kdgcyc  = .020;
 kdrhodk = .002;

 rhod  = 0.0;
 rhodk = 0.0;
 rhod2 = 0.0;
 rstar = 0.0;
 gpr   = 0.0;
 pde   = pdebase;

 cycg  = 1.03805;       /* values taken from tcomp21b */
 cond  = .902444;
 totcond = 20.0;
 ca    =  0.315626;
 cax   =  0.184104;

 lgain    = rrate * .25;        /* adjust to give R* resp ~ 0.02 pA */
 rgain1   = calctau (0.0004,rrate);
 rgain2   = calctau (0.0004,rrate); /* 3 rhod steps give .5 ms tau */
 rgain3   = rrate * 800;
 rkgain   = rrate * 50;      /* adjust this for amount of adaptation */

                                          /*  for R* *  activation */
 ggain    = rrate * 100.0;
 gcygain  = rrate * 11.2 / cycgbuf * loopgain;
 pdegain  = rrate * 100.0 / cycgbuf;
 condgain = rrate * 448 * loopgain;
 gca      = rrate * 270.0 / cabuf * loopgain;
 cxgain   = rrate * 5.0 * loopgain;

 decrstar = 1 - (rrate * 20);
 pdedec   = 1 - (rrate * 1000); /* controls tail of cond */
 deccond  = 1 - (rrate * 480);
 capump   =      rrate * 1000 / cabuf;
 deccax   = 1 - (rrate * 10);

 startfl = 1;
 plotinc = 10;
 for (i=0,j=plotinc; i<1000; i++) {

  if (i == 500) mglur6 += 1000.0 * lgain;			/* */

  t1 = mglur6  * rgain1;
  rhod  -= t1;
  rhod2 += t1;

  t2 = rhod2 * rgain2;
  rhod2 -= t2;
  rstar += t2;
 
  gpr   +=  rstar * rgain3;

  t3 = gpr * ggain;
  gpr -= t3;
  pde += t3;

    /* Calcium feedback */

  ca2 = cax;
  ca2 = ca2 * ca2;

  gcyc  = 1.0 - (ca2 / (ca2 + kdgcyc));
  rhodk = 1.0 - (ca2 / (ca2 + kdrhodk));

  cycg -= pde * cycg * pdegain;
  cycg += gcyc * gcygain;
  if (cycg < 0) cycg = 0;

  cond += cycg * cycg * cycg *
           (totcond - cond) /
            totcond * condgain;
  ca   += cond * gca;
  cax  += ca * cxgain;

  rstar *= decrstar - rhodk * rkgain;

  pde  = pdedec * (pde- pdebase) + pdebase;
  ca   -= capump * ca / (ca + 1.0);
  cond *= deccond;
  cax  *= deccax;

  if (i >= 00) {
    if (startfl) {
	startfl = 0;
	startval = cond;
    } 
  if (++j >= plotinc) {
	j = 0;
	//printf ("%d %10.9g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",
	//		i,cond,cond-startval,pde,gcyc,cycg,ca,cax,ca2,rhodk);
	printf ("%d %10.9g %8.5g %8.5g %8.5g %8.5g\n",
		i,cond,pde,gcyc,ca2,rhodk);
   }
  }
 }
}


