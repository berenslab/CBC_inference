/* program rod.c   simulates fifth-order cone kinetics */


main(void)
{
  double dec1, dec2, dec3, dec4, dec5, dec6;
  double gpr2, gpr3, rhodgain, gprgain, ggain2, ggain3;
  double pdegain, gcygain, condgain, cxgain;
  double rhod, gprt, pde, pdebase, gcyc, cycg, ca, gca, cax, cabuf;
  double ca2, capump, kdgcyc;
  double rise, fall, startval, cond, totcond, thresh, thresh6;
  double cycgbuf, loopgain, factor, factor5, factor6;
  int i,j,startfl,plotinc;

/* rise = .020;
 fall = .001;
*/

 factor = 1.34 * .5;

 rise = .020 / factor;
 fall = .001 / factor;

 rhod = 0.0;
 gprt = 0.0;
 gpr2 = 0.0;
 gpr3 = 0.0;
 pde  = 0.0;
 pdebase = .02;
 kdgcyc = .010;
 cycg = 1.03564096155416; 
 cond = 1.00009243; 
 totcond = 20.0;
 ca   =   0.340753892653351; 
 cax  =  0.156197785692584; 
 cabuf  = .2;
 cycgbuf = .5;
 loopgain = .2;

 rhodgain = rise * 46.;
 gprgain  = rise * 1.0;
 ggain2   = rise * 1.0;
 ggain3   = rise * 1.0;
 gcygain  = rise * .35 / cycgbuf * loopgain;
 pdegain  = rise * 2.0 / cycgbuf;
 condgain = rise * 7.53 * loopgain;
 gca      = rise * 3.45 / cabuf * loopgain;
 cxgain   = rise * 1.0 * loopgain;

 dec1 = 1 - (fall * 180);
 dec2 = 1 - (fall * 180);
 dec3 = 1 - (fall * 180);
 dec4 = 1 - (fall * 5);
 dec5 = 1 - (fall * 30);
 capump =    fall * 50 / cabuf;
 dec6 = 1 - (fall * 10);

 startfl = 1;
 plotinc = 20;
 for (i=0,j=plotinc; i<2500; i++) {

  if (i == 1000) rhod += 20.0 * rhodgain * factor;	/* */
/*  if (i > 400  && i <= 600) rhod += 100.0 * rhodgain;	/* */
/*  if (i <= 400) rhod += 1.0 * rhodgain;		/* */
/*  rhod += .1 * rhodgain;				/* */
  gprt += rhod * gprgain;
  gpr2 += gprt * ggain2;
  gpr3 += gpr2 * ggain3;

  pde  = gpr3 + pdebase;
  ca2 = cax * cax;
  gcyc  = 1.0 - (ca2 / (ca2 + kdgcyc));

  cycg -= pde * cycg / (cycg + 1.0) * pdegain;
  cycg += gcyc * gcygain;

  cond += cycg * cycg * (totcond - cond) / totcond * condgain;		/* */
  ca   += cond * gca;
  cax  += ca * cxgain;

  rhod *= dec1;
  gprt *= dec2;
  gpr2 *= dec3;
  gpr3 *= dec4;
  /* ca   -= capump * ca / (ca + 1.0);  */
  ca   -= capump * ca / (ca + 1.0); 
  cond *= dec5;
  cax  *= dec6;

  if (i >= 00) {
    if (startfl) {
	startfl = 0;
	startval = cond;
    } 
  if (++j >= plotinc) {
	j = 0;
	printf ("%d %10.9g %8.5g %15.15g %15.15g %15.15g %15.15g\n",
		i,cond,pde,gcyc,cycg,ca,cax);
   }
  }
 }
}


