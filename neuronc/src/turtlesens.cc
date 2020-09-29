/* turtlesens.cc */
 
/*  See explanation in the sensitivity function file "perlman.txt".  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ncfuncs.h"

#define turtle(row, col) (*(turtlearr+(row)*cols+(col)))

/*------------------------------------------------*/
 
int lshift = 19;	/* shift the L cone curve left by 19 points */
int sshift = 17;	/* shift the S cone curve right by 17 points */

double sl = 16700;
double  b = -28.0;

double pigmlen = 35;
double lpeaksens = .7;
double mpeaksens = .7;
double speaksens = .7;
double pigml = 35;			/* vary this to change peak width */

#define PL 0
#define PM 1
#define PS 2

#define maxwav 800
#define minwav 380
#define wavinc 5

#define wavsiz ((maxwav - minwav) / wavinc + 1)
 
double wav[wavsiz];
double pigms[3][wavsiz];

/*------------------------------------------------*/

#define STRSIZ 480

double *xfread(const char *filnam, int *nlongp, int *nwidp)

{
    int j,nlong,nwid,arrsiz;
    FILE *ftemp;
    char *cp,*np,*strp;
    double *p,*ptr;
    static char str[STRSIZ];
    static char cbuf[STRSIZ];

#define NLONG 100
#define NWID  4

 arrsiz = NLONG*NWID*sizeof(double);
 ptr = (double *)malloc (arrsiz);

 for (cp=(char*)ptr,j=0; j<arrsiz; j++) *cp++= 0;
 p = ptr;

 if (strcmp(filnam,"stdin")==0) ftemp = stdin;
 else if ((ftemp=fopen(filnam,"r")) == NULL) {
   fprintf (stderr,"fread: can't open file '%.50s'.\n",filnam);
   return NULL;
 }
 while ((strp=fgets(str,STRSIZ,ftemp)) && *str== '#'); /* get first line */
 if (!strp) {
    fclose (ftemp);
    fprintf (stderr,"fread: file '%.50s' is zero length.\n",filnam);
    return NULL;
 }
 cp=strstr(strp,"#");                   /* look for # after values */
 if (cp) {
    *cp = '\0';                         /* terminate string at # */
 }
                                /* must use buffered read for stdin */
 strncpy(cbuf,str,STRSIZ-1);

 for (cp=np=cbuf,nwid=0; cp; nwid++) {  /* how many numbers on line */
    cp = strtok(np," ,\t\n");           /* get pointer to next token */
    np = (char *)NULL;
 }
 nwid--;

 for (cp=np=str,j=0; j<nwid; j++) {                     /* scan first line */
   cp = strtok(np," ,\t\n");
   if (cp==NULL) {
      fprintf (stderr,"Error in format of file '%s' while reading array.\n",filnam);
   }
   np = (char *)NULL;
   *ptr++ = atof(cp);
 }   /* for j */

                                                /* read in lines from file */
 for (nlong=1; ; nlong++) {
   while ((strp=fgets(str,STRSIZ,ftemp)) && *str=='#');/* get next line */
   if (!strp) break;
   np = str;
   cp=strstr(strp,"#");                 /* look for # after values */
   if (cp) {
     *cp = '\0';                                /* terminate string at # */
   }

 for (j=0; j<nwid; j++) {                       /* scan lines */
   cp = strtok(np," ,\t\n");
   if (cp==NULL) {
      fprintf (stderr,"Error in format of file '%s' while reading array, line %d.\n",
                                                        filnam,nlong);
            break;
   }
   np = (char *)NULL;
  *ptr++ = atof(cp);
 }   /* for j */

 }  /* for (nlong=1;;) */
 *nlongp = nlong;
 *nwidp  = nwid;
  fclose (ftemp);
  return p;
}

/*------------------------------------------------*/

double longwav (double j) {
    return sl * 1/(minwav+j*wavinc) + b;
}

/*------------------------------------------------*/
   
main(int argc, char **argv)

{
   int i, j, p, lines, cols;
   int startwav, startplot, stoplot;
   int linelen = 100;
   double *turtlearr;

 turtlearr = xfread ("perlman.txt",&lines,&cols);
 
/* the turtle sensitivity functions */

 for (j=0; j<3; j++) {
  for (i=0; i<wavsiz; i++) {
   wav[i] = 0;
   pigms[j][i]   = 0;
  }
 }

/****************************************************/
/* Model of long-wave log(sens) function, after Lewis, 1955) */
/*  This relies on the fact that at long wavelengths the     */
/*  log(sensitivity) of visual pigments asymtotes to a       */
/*  straight line plotted against wavenumber (1/wavelength). */

/* The three pigments' log(sens) functions are supposed to superimpose */
/*  when they are plotted against wavenumber and shifted horizontally */
/*  to superimpose their peaks. The slope of the line may vary a little
/*  but is assumed here to be the same for all three functions when   */
/*  superimposed. */


/****************************************************/

/* Uncomment this part to print out the log10 of turtle cone sensitivity */
/*  from "perlman.txt". */

/*
for (i=0; i<wavsiz; i++) {
 printf ("%6.4g %6.4g %6.4g %6.4g\n",turtle(i,0),
		log10(turtle(i,1)),
		log10(turtle(i,2)),
		log10(turtle(i,3)));
}; /* */


/* Read in data, set L, M, and S curves to be log10 of sensitivity */

for (i=0; i<lines; i++) {

 wav[i] = turtle(i,0);
 pigms[PL][i]   = log10(turtle(i,1));
 pigms[PM][i]   = log10(turtle(i,2));
 pigms[PS][i]   = log10(turtle(i,3));

}; /* */

/* Extend the wavelength numbers to 800 nm */
/*  to allow extrapolation */

for (i=lines; i<wavsiz; i++) {
   wav[i] = minwav + i * wavinc;
}; /* */

/****************************************************/

/*    Shift L and S curves to overlap M curve */
/*     Do this only for test purposes */

/*
for (i=0; i<wavsiz; i++) {

 if (i+lshift < wavsiz) pigms[PL][i] = log10(turtle(i+lshift.1))
 else                   pigms[PL][i] = log10(turtle(i.1));
 pigms[PM][i]   = log10(turtle[i.2]);
 if (i-sshift >= 0)     pigms[PS][i] = log10(turtle(i-sshift.3))
 else                   pigms[PS][i] = log10(turtle(i.3));

}; /* */


/* Extrapolate the L curve beyond 680 nm */
/*  Between 680 and 730, use the S curve, shifted right */
/*  Beyond 730, use longwave extrapolation */

for (i=0; i<wavsiz; i++) {

 if (wav[i]<=680) pigms[PL][i] = log10(turtle(i,1));
 else 
 if (wav[i]<=730) pigms[PL][i] = log10(turtle(i-lshift-sshift,3));
 else             pigms[PL][i] = longwav(i-lshift);
   
}; /* */

/* Extrapolate the M curve beyond 650 nm */

for (i=0; i<wavsiz; i++) {

 if (wav[i]<=645) pigms[PM][i] = log10(turtle(i,2));
 else             pigms[PM][i] = longwav(i);
   
}; /* */

/* Extrapolate the S curve beyond 540 nm */

for (i=0; i<wavsiz; i++) {

 if (wav[i]<=535) pigms[PS][i] = log10(turtle(i,3));
 else             pigms[PS][i] = longwav(i+sshift);
   
}; /* */

/****************************************************/

/* Convert to specific optical density (for testing only). */

/*
for (i=0; i<wavsiz; i++) {

   pigms[PL][i] = -log(1-exp(pigms[PL][i]*log(10))*lpeaksens)/log(10) / pigmlen;
   pigms[PM][i] = -log(1-exp(pigms[PL][i]*log(10))*mpeaksens)/log(10) / pigmlen;
   pigms[PS][i] = -log(1-exp(pigms[PL][i]*log(10))*speaksens)/log(10) / pigmlen;

}; /* */

/****************************************************/

/* Plotting variables */

startwav = 380;
startplot = (startwav-minwav)/wavinc;
stoplot = wavsiz;

pigml = pigmlen;

for (i=startplot,p=1; i<stoplot; i++,p++) { /* */


/*printf ("%g %g %g %g \n", wav[i], pigms[PL][i], pigms[PM][i], pigms[PS][i]);

/* printf ("%g %g %g %g \n", wav[i], 1-exp(-pigms[PL][i]*pigml*log(10)), 
				  1-exp(-pigms[PM][i]*pigml*log(10)), 
				  1-exp(-pigms[PS][i]*pigml*log(10)));
*/

/*
printf ("%g %g %g %g \n", wav[i], exp(pigms[PL][i]*log(10)), 
				  exp(pigms[PM][i]*log(10)), 
				  exp(pigms[PS][i]*log(10)));
*/

//printf ("%g %g %g %g %g\n", - 1/wav[i], pigms[PL][i], 
//				pigms[PM][i], pigms[PS][i], longwav(i));
//printf ("%g %g %g\n", - 1/wav[i],  pigms[PL][i], longwav(i));
//printf ("%g %g %g\n", wav[i],  pigms[PM][i],     longwav(i));
//printf ("%g %g %g\n", wav[i-sshift],  pigms[PS][i], longwav(i));
//printf ("%g %g %g\n", wav[i+lshift],  pigms[PL][i], longwav(i));

};

linelen = 8;

printf ("/* Turtle sensitivity functions (log(sens)) */\n");
printf ("/*  from 380 to 800 nm in 5 nm increments.  */\n");
printf ("/*  Order: L, M, S cones */\n\n");
printf ("double turtlepigm[3][PIGMSIZ] = {\n ");


for (j=PL; j<3; j++) {
  for (i=startplot,p=1; i<stoplot; i++,p++) { /* */
    printf ("%7.5g,", pigms[j][i]);
    if (p>=linelen) {
	printf ("\n ");
	p = 0;
    }
  }
  printf ("\n\n ");
}
printf ("  };\n");
}

