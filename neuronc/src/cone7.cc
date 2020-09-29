/* program cone.c   simulates fifth-order cone kinetics */


main(void)
{
  double dec1, dec2, dec3, dec4, dec5, dec6;
  double gpr2, gpr3, rhodgain, gprgain, ggain2, ggain3;
  double pdegain, gcygain, condgain, cxgain;
  double rhod, gprt, pde, pdebase, gcyc, cycg, ca, gca, cax, cabuf;
  double ca2, capump, kdgcyc;
  double timinc, startval, cond, totcond, thresh, thresh6;
  double cycgbuf, loopgain;
  int i,j,startfl,plotinc;

 timinc = .001;

 rhod = 0.0;
 gprt = 0.0;
 gpr2 = 0.0;
 gpr3 = 0.0;
 pde  = 0.0;
 pdebase = .01;
 kdgcyc  = .02;
 cycg = 1.4454;
 cond = 1.0;
 totcond = 20.0;
 ca   = 0.48447;
 cax  = .12463;
 cabuf  = .1;
 cycgbuf = .02;
 loopgain = .70;

 rhodgain = timinc * 5.;
 gprgain  = timinc * 20.0;
 ggain2   = timinc * 20.0;
 ggain3   = timinc * 20.0;
 gcygain  = timinc * .6 / cycgbuf * loopgain;
 pdegain  = timinc * 40.0 / cycgbuf;
 condgain = timinc * 80.0 * loopgain;
 gca      = timinc * 26.2 / cabuf * loopgain;
 cxgain   = timinc * 8. * loopgain;

 dec1 = 1 - (timinc * 200);
 dec2 = 1 - (timinc * 200);
 dec3 = 1 - (timinc * 200);
 dec4 = 1 - (timinc * 50);
 dec5 = 1 - (timinc * 100);
 capump =    timinc * 50 / cabuf;
 dec6 = 1 - (timinc * 30);

 startfl = 1;
 plotinc = 10;
 for (i=0,j=plotinc; i<1000; i++) {

  if (i == 500) rhod += 2000.0 * rhodgain;			/* */
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
	printf ("%d %10.9g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",
		i,cond,cond-startval,pde,gcyc,cycg,ca,cax);
   }
  }
 }
}


