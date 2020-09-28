/* program cone.c   simulates fifth-order cone kinetics */

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
  int i,j,startfl;
  double stiminc, startval, plotinc;
  double qc, qt, qcond, qrate, tempcel;
  double sfactor, rrate;
  double cabuf, cycgbuf, loopgain;
  double pdebase, kdgcyc, kdrhodk;
  double rhod, rhodk, rhod2, rstar;
  double gpr, pde, gcyc, cycg, cond, totcond;
  double t1, t2, t3;
  double ca, cax;
  double cab, cag, cak, cat;
  double dcab, gcab, gcabr;
  double lgain, rgain1, rgain2, rgain3, rkgain;
  double ggain, gcygain, pdegain, condgain, gca, cxgain;
  double decrstar, pdedec, deccond, capump, deccax;
  double powcag, powcak, powcycg;
   
 stiminc = .0001;

 qt = 37.0;             /* base temperature for channel rate & cond */
 qc = 1.4;              /* Q10 factor for conductance shape */
 tempcel = 35;

 qrate = exp(log(qc) * (tempcel - qt) / 10.0);

 
 sfactor = 1.1;               /* speedup of cascade & intensity */
 sfactor = .5;               /* speedup of cascade & intensity */
 rrate = stiminc * qrate /sfactor;    /* stiminc = timinc for syn/stimuli */

 cabuf    = 1.0;
 cycgbuf  = 1.0;
 loopgain = 1.0;

 powcag  = 2.3;
 powcak  = 2;
 powcycg = 3;

 pdebase = .04;
 kdgcyc  = .020;
 kdrhodk = .01;

 rhod  = 0.0;
 rhodk = 0.0;
 rhod2 = 0.0;
 rstar = 0.0;
 gpr   = 0.0;
 pde   = pdebase;

 cycg  = 1.0307;       /* values taken orig from tcomp21b */
 cond  = .902444;
 totcond = 20.0;
 ca    =  0.405;
 cax   =  0.23119;
 cat   = 250;
 cab   = 9.7458;

 lgain    = rrate * 4.0;        /* adjust to give R* resp ~ 0.02 pA */
 rgain1   = calctau (0.0004,rrate);
 rgain2   = calctau (0.0004,rrate); /* 3 rhod steps give .5 ms tau */
 rgain3   = rrate * 800;
 rkgain   = rrate * 50000;      /* adjust this for amount of adaptation */

                                          /*  for R* *  activation */
 ggain    = rrate * 40.0;
 gcygain  = rrate * 30.0 / cycgbuf * loopgain;
 pdegain  = rrate * 100.0 / cycgbuf;
 condgain = rrate * 448 * loopgain;
 gca      = rrate * 325.0 / cabuf * loopgain;
 gcab     = rrate * 20  * loopgain;
 gcabr    = rrate * 200 * loopgain;
 cxgain   = rrate * 2.0 * loopgain;

 decrstar = 1 - (rrate * 20);
 pdedec   = 1 - (rrate * 100); /* controls tail of cond */
 deccond  = 1 - (rrate * 480);
 capump   =      rrate * 1000 / cabuf;
 deccax   = 1 - (rrate * 4);
 
 startfl = 1;
 plotinc = 10;
 for (i=0,j=plotinc; i<2500; i++) {

  if (i == 500) rhod += 1000.0 * lgain;			/* */
/*  if (i > 400  && i <= 600) rhod += 100.0 * lgain;	/* */
/*  if (i <= 400) rhod += 1.0 * lgain;			/* */
/*  rhod += .1 * lgain;					/* */

  t1 = rhod  * rgain1;
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

  cag = pow(ca, powcag);
  cak = pow(cax,powcak);   /* rhodk has extra delay */

  gcyc  = 1.0 - (cag / (cag + kdgcyc));
  rhodk = 1.0 - (cak / (cak + kdrhodk));

  cycg -= pde * cycg * pdegain;
  cycg += gcyc * gcygain;
  if (cycg < 0) cycg = 0;

  cond += pow(cycg,powcycg) * (totcond - cond) / totcond * condgain;

  // This implements Russ Hamer's (2005, VN 22:417) calcium buffer.
  // It delays the effect of calcium. 
 
  dcab = ca * gcab * (cat - cab) - cab * gcabr;
 
  cab += dcab;

  ca   += cond * gca - dcab;
  cax  += ca * cxgain;

  rstar *= decrstar - pow(rhodk,2) * rkgain;

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
	printf ("%g %10.9g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",
			i/plotinc,cond,pde,gcyc,cycg,ca,cak,rhodk);
	//printf ("%d %10.9g %8.5g %8.5g %8.5g %8.5g\n",
	//	i,cond,pde,gcyc,cak,rhodk);
   }
  }
 }
}


