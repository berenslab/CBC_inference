/* program cone.c   simulates fifth-order cone kinetics */


main(void)
{
  double dec1, dec2, dec3, dec4, dec5, dec6, dec7;
  double gpr2, gpr3, rhodgain, gprgain, ggain2, ggain3;
  double pdegain, gcygain, condgain, modgain, ca4gain;
  double rhod, gprt, pde, pdebase, gcyc, cycg, ca, gca, mod, cabuf;
  double ca2, ca4, capump, kdmod, kdca4;
  double timinc, startval, cond, totcond, thresh, thresh6;
  double cycgbuf, loopgain;
  int i,j,startfl,plotinc;

 timinc = .001;

 rhod = 0.0;
 gprt = 0.0;
 gpr2 = 0.0;
 gpr3 = 0.0;
 pde  = 0.0;
 pdebase = .04;
 kdca4  = 1.;
 kdmod  = 1.;
 cycg = .1;
 gcyc = .1;
 cond = 1.0;
 totcond = 20.0;
 ca   = 0.1;
 mod  = .12463;
 cabuf   = .1;
 cycgbuf = .05;
 loopgain = 1.0;

 rhodgain = timinc * 15.;
 gprgain  = timinc * 20.0;
 ggain2   = timinc * 20.0;
 ggain3   = timinc * 20.0;
 gcygain  = timinc * 1.0 / cycgbuf * loopgain;
 pdegain  = timinc * 40.0 / cycgbuf;
 condgain = timinc * 200.0 * loopgain;
 gca      = timinc * 10. / cabuf * loopgain;
 ca4gain  = timinc * 800. * loopgain;
 modgain  = timinc * 100. * loopgain;

 dec1 = 1 - (timinc * 200);
 dec2 = 1 - (timinc * 200);
 dec3 = 1 - (timinc * 200);
 dec4 = 1 - (timinc * 20);		/* gpr3 */
 dec5 = 1 - (timinc * 100);		/* gcyc */
 dec6 = 1 - (timinc * 100);		/* cond */
 capump =    timinc *  15 / cabuf;
 dec7 = 1 - (timinc * 100);

 startfl = 1;
 plotinc = 10;
 for (i=0,j=plotinc; i<1000; i++) {

  if (i == 500) rhod += 20000.0 * rhodgain;		/* */
/*   if (i == 100) rhod += 2000.0 * rhodgain;		/* */
/*   if (i == 700) rhod += 2000.0 * rhodgain;		/* */
/*   if (i > 400  && i <= 900) rhod += 100.0 * rhodgain;	/* */
/*  if (i <= 400) rhod += 1.0 * rhodgain;		/* */
  rhod += 100.0 * rhodgain;				/* */
  gprt += rhod * gprgain;
  gpr2 += gprt * ggain2;
  gpr3 += gpr2 * ggain3;

  pde  = gpr3 + pdebase;

  gcyc  += modgain * (1 - (mod / (mod + kdmod)));

  cycg -= pde * cycg / (cycg + 1.0) * pdegain;
  cycg += gcyc * gcygain;

  cond += cycg * cycg * (totcond - cond) / totcond * condgain;		/* */
  ca   += cond * gca;
  ca2   = ca * ca;
  ca4   = ca2 * ca2;
  mod  += ca4gain * (1 - mod) * ca4 / (ca4 + kdca4);

  rhod *= dec1;
  gprt *= dec2;
  gpr2 *= dec3;
  gpr3 *= dec4;
  gcyc *= dec5;
  ca   -= capump * ca / (ca + 1.0); 
  cond *= dec6;
  mod  *= dec7;

  if (i >= 00) {
    if (startfl) {
	startfl = 0;
	startval = cond;
    } 
  if (++j >= plotinc) {
	j = 0;
	printf ("%d %10.9g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",
		i,cond,cond-startval,pde,gcyc,cycg,ca,mod);
   }
  }
 }
}


