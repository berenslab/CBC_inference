/* program stimplot */

/* Program to make stimulus and blur it by a gaussian.
 Prints array on standard output */

/* This program was never finished or debugged.    R.G.Smith */

#include <stdio.h>
#include "stim.h"

#define SPOT 1
#define BAR 2
#define BLUR 3
#define STIM 4
#define BACK 5
#define STIMBLUR 6
#define UDEF 7

int reccum=0;			/* number of receptors */
extern int stimdia;
static double blurrad;
extern recnod reclist[];
extern double blurarr[];

extern double atof();
extern double getblur();

int operation=UDEF;

FILE *stimout=stdout;
FILE *ncstdout=stdout;
FILE *ncstderr=stderr;

FILE *pictin;

/*------------------------------------*/

main(argc,argv)
   int argc;
   char **argv;

begin
   char *cptr;
   int i;
   FILE *temp,*freopen();
	 
 pictin = stdin;
 if (argc==1)			/* if user needs help */
   stimrun();
 else
 do					/* if there are any switches */
  begin
   argc--; argv++;
   cptr = *argv;
   if (argc)
    begin
     if (*cptr == '-')
      begin
       cptr++;
       switch (*cptr)
        begin

	  case 's': operation = SPOT; break;
	  case 'l': operation = BAR; break;

	  case 'w': operation = STIMBLUR; break;
	  case 'x': operation = BLUR; break;
	  case 'y': operation = STIM; break;
	  case 'z': operation = BACK; break;
	   
	  case 'b': 
		argv++; argc--;
		blurrad = atof(*argv);
		break;
	
	  case 'd': 
		argv++; argc--;
		stimdia = atof(*argv);
		break;
	
 	/*     	  case 'e': 
		argv++; argc--;
		stop = atof(*argv);
		stopfg = T;
		break;
	*/

	  default:
		ncfprintf (ncstderr,"stim: unknown switch '%s'\n",*argv);
		exit();

        end  /* switch */
      end	   /* if */
     else
      begin
       if((pictin=fopen(cptr,"r"))==NULL)
         begin
           ncfprintf(ncstderr,"stim: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         end
       stimrun();
       if (argc <= 1) break;
     end
    end
   else stimrun();
  end
 while (argc > 0);
end

/*------------------------------------*/

stimrun()

begin
    int i,j,n1,stimarr,backarr;
    int xstim,ystim;
    double stimval,backval;
    double start, dur;
    recnod *npt;

stimval = 100.0;
backval = 10.0;
stimarr = 0;
backarr = 1;
xstim = ystim = 0;		/* no offset */
start = 0.0;			/* start at time 0 */
dur = 100e-3;			/* 100 msec stimulus */
initc(2,256);			/* set up two 256 * 256 arrays */

setc(backarr, backval);		/* make a background level in array */
if (operation != BACK && operation != STIM) {
	blursize = makblur(blurrad);		/* make a blurring function */
	ncfprintf (ncstderr,"gaussian done...\n"); /* */
}
switch (operation) {

  case BACK:	printa (backarr);
		break;

  case STIM:	setc (stimarr,backval);
		makrect(stimarr, stimdia, 200, xstim, ystim, stimval);
		printa (stimarr);
		break;

  case STIMBLUR:setc (stimarr,backval);
		makrect(stimarr, stimdia, 200, xstim, ystim, stimval);
		reccum = 20; 
		for (npt=reclist,i=0; i<reccum; npt++, i++) {
		  npt->xpos = (i * 5) - 50;
		  npt->ypos = 0;
		}
		recpat(stimarr,xstim,ystim); /* set recept inten */
		recpat(backarr,xstim,ystim); /* set recept inten */
		for (npt=reclist,i=0; i<reccum; npt++, i++) {
		  ncfprintf (ncstdout,"%d %g\n",npt->xpos, npt->stim);
		}
		break;

  case BLUR:	printb ();
		break;

  case UDEF:	
	 	/* readnod(reclist);	/* read stimulus descriptor file */
		break;

  case SPOT:	setc(stimarr,backval);
		makspot(stimarr,  stimdia, xstim, ystim, stimval);
	 	/* readnod(reclist);	/* read x,y pos of nodes from file */
		recpat(stimarr,xstim,ystim); /* set recept inten */
		recpat(backarr,xstim,ystim); /* set recept inten */
		stimlist(1.0, start);	/* make an action list */
		stimlist(-1.0, start+dur);

		ncfprintf (ncstderr,"spot normal exit...\n");
		break;
 
 case BAR:	setc(stimarr,backval);
		ncfprintf (ncstderr,"beginning rect...\n");
		makrect(stimarr, 10, 100, xstim, ystim, stimval);
		ncfprintf (ncstderr,"rect done...\n");

	 	/* readnod(reclist);	/* read x,y pos of nodes from file */
		recpat(stimarr,xstim,ystim); /* set recept inten */
		recpat(backarr,xstim,ystim); /* set recept inten */
		stimlist(1.0, start);	/* make an action list */
		stimlist(-1.0, start+dur);

		ncfprintf (ncstderr,"rect normal exit...\n");
		break;

}	/* switch */

end

/*------------------------------------*/


