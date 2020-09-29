/* program cone.c   simulates fifth-order cone kinetics */


main(void)
{
  double gain1, gain2, gain3, gain4, gain5, gain6, gain7, gain8;
  double dec1, dec2, dec3, dec4, dec5, dec6;
  double gpr2, gpr3, rhodgain, gprgain, ggain2, ggain3;
  double pdegain, gcygain, condgain, cxgain;
  double rhod, gprt, pde, pdebase, gcyc, cycg, ca, gca, cax, cabuf;
  double ca2, capump;
  double rise, fall, startval, cond, totcond, thresh, thresh6;
  double cycgbuf;
  int i,j,startfl;

 rise = .020;
 fall = .001;

 rhod = 0.0;
 gprt = 0.0;
 gpr2 = 0.0;
 gpr3 = 0.0;
 pde  = 0.0;
 pdebase = .1;
 cycg = 0.6;
 cond = 0.7;
 totcond = 20.0;
 ca   = 0.4;
 cax  = .4;
 cabuf  = .5;
 cycgbuf = .1;

 rhodgain = rise * .2;
 gprgain  = rise * 1.0;
 ggain2   = rise * 1.0;
 ggain3   = rise * 1.0;
 gcygain  = rise * 1.0 / cycgbuf;
 pdegain  = rise * 10.0 / cycgbuf;
 gca      = rise * 8.0 / cabuf;
 condgain = rise * 2.5;
 cxgain   = rise * 1.0;

 dec1 = 1 - (fall * 100);
 dec2 = 1 - (fall * 100);
 dec3 = 1 - (fall * 100);
 dec4 = 1 - (fall * 20);
 capump =    fall * 150 / cabuf;
 dec5 = 1 - (fall * 100);
 dec6 = 1 - (fall * 6);

 startfl = 1;
 for (i=0,j=5; i<800; i++) {

  if (i == 400) rhod += 100.0 * rhodgain;		/* */
/*  if (i > 400  && i <= 600) rhod += 100.0 * rhodgain;	/* */
/*  if (i <= 400) rhod += 1.0 * rhodgain;		/* */
  rhod += .001 * rhodgain;
  gprt += rhod * gprgain;
  gpr2 += gprt * ggain2;
  gpr3 += gpr2 * ggain3;

  pde  = gpr3 + pdebase;

  ca2 = cax * cax;
  gcyc  = 1.0 - (ca2 / (ca2 + 1.0));

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

  if (i >= 400) {
    if (startfl) {
	startfl = 0;
	startval = cond;
    } 
  if (++j >= 5) {
	j = 0;
	printf ("%d %10.9g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",
		i,cond, cond-startval,pde,gcyc,cycg,ca,cax);
   }
  }
 }
}


