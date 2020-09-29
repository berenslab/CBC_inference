#
/* Filter for positioning.  R.G.Smith 
 *
 * Modified from:
 *
 * vector plot filter for Tektronix 4012
 *
 * author- Rob Clayton, Dept. of Geophysics, Stanford Univ.
 * Modified for Version 7, 12/3/79
 * Mitchell Gart, Ampex Corp.
 */


#include	<stdio.h>
#include	<string.h>
#include	<signal.h>
#include	<sys/types.h>
#include	<sys/stat.h>
  
typedef void (*sighandler_t)(int);

int istat;				/* look in "UNIX programming" */
int iflag = 0;				/* interrupt flag used by mdisplay */

#define BUFFER	0
#define CHARHT	22*21
#define PLTCMD	0200
#define CMDMSK	0177
#define INCHDOT	(8.0/1024.)
#define STDSIZ 16384	/* standard coordinate range */
#define STDSCL 	1.	/* STDSIZ/780 to scale standard plot */
#define TWMUL 	64	/* # of pixels in standard text width size spec. */
#define THMUL 	96	/* 3/2*TWMUL */
#define STDCW 	5	/* standard text size (about 15 x 23) */
#define XMAX 	32000	/* maximum x value */
#define XMIN 	(-16384)	/* minimum x value */
#define YMAX 	32000	/* maximum y value (manual says 779, but 780 fits) */
#define YMIN 	(-16384)	/* minimum y value */
#define XOFFS   (STDSIZ/2)  /* offset of standard picture */
#define YOFFS   (STDSIZ/2)  /* offset of standard picture */
#define XCENT   XOFFS	/* center of plotter platen */
#define YCENT	YOFFS
#define SCAL 	(STDSIZ /10000.)  /* scale for xcent, ycent */

#define ENDRAS 	(-1)	/* end of raster data flag */
#define ERASE	014
#define ESC	033
#define CSZ	032
#define GS	035
#define CUS	037
#define GRKFONT 'g'	/* greek character set */
#define NPENS 8		/* max number of plotter pens */

int xmax = XMAX;	/* plot window parameters defaulted */
int xmin = XMIN;	/* to maximum size 		*/
int ymax = YMAX;
int ymin = YMIN;
int xnew,ynew;		/* new pen location */
int xold,yold;		/* old pen location */
double scale;
int separate;		/* output redirected */
int intract, quitact;	/* old signal actions */
int pen,npen,speed;	/* pen, speed for plotter */
double xcent,ycent;	/* position of center of picture */
int makfil;		/* make output files from input file names */

struct inode
  {	int dev;
	int inum;
	int flags;
  };
struct
   {
	char ispeed, ospeed;
	char erase, kill;
	int mode;
   } ttystat;

int fixtty();
FILE *pltout, *pltin, *prinout, *temp;
char inbuf[BUFSIZ], outbuf[BUFSIZ];
char *strrchr(const char *, int);
int atoi(char*);

int getwsx();
void putwsx();
void doplot(void);

/*-------------------------------------*/

void onintr(void)
{
  extern int iflag;

  iflag = 1; 
  signal (2,(sighandler_t)onintr);		/* go to "onintr" at ^C */
}

/*-------------------------------------*/

int main(int argc, char **argv)
{
	FILE *freopen(const char *, const char *, FILE *);
	double atof(const char *);
	float fscale;
	register char *cptr;
	char outnam[20], *namend;

	iflag = 0;
   	istat = (int) signal (2, (sighandler_t) 1);	/* save original status */
   	if (istat != 1)
  	   signal (2,(sighandler_t) onintr);		/* go to "onintr" at ^C */

/*	fstat(1, inbuf);
	fstat(2, outbuf);
	separate = (inbuf->st_ino == outbuf->st_ino)? 0: 1; */
	pltin= stdin;
	pltout= stdout;
	prinout= stderr;
	npen = 0;
	makfil = 0;
/*	setbuf(pltin,inbuf);
	if(BUFFER) setbuf(pltout,outbuf); */
	scale= STDSCL;
	xcent = XCENT;
	ycent = YCENT;

	do	/* loop over arg list */
	   {
		argc--; argv++;
		cptr= *argv;
		if(argc)
		 {
			if (*cptr == '-')
		         {
				cptr++;
				switch(*cptr)
				 {
					case 'i':
						signal (2,(sighandler_t)1);
						break;
					case 'f':		/* output files */
						makfil = 1;
						break;
					case 'm':	/* mag */
						argv++; argc--;
						fscale= atof(*argv);
						if(fscale<=0.0)
						  scale= STDSCL;
					        else
						  scale= STDSCL/fscale;
						break;
					 case 'p':	/* pen override */
						npen = atoi(*(++argv));
						argc--;
						break;
					 case 'x':	/* x center */
						xcent = atof(*(++argv));
						xcent *= SCAL;
						argc--;
						break;
					 case 'y':	/* y center */
						ycent = atof(*(++argv));
						ycent *= SCAL;
						argc--;
						break;
					default:
					break;  /* default */

				  }		/* switch */
			  }			/* if *cptr */
			else
			 {			/* *cptr not= '-' */

			   if((temp=freopen(cptr,"r",pltin)) == NULL)
			      {
				fprintf(stderr,"posit:cannot open %s\n",cptr);
				fflush(stderr);
				continue;
			      }
			   else	pltin= temp;
       			if (makfil)			/* make file mode */
			 {
	  		  namend = strrchr (cptr,'.');	/* find last dot in name */
	  		  if (namend) 
	    		   {
	     		     strncpy (outnam,cptr,namend-cptr);
	     		     *(outnam+(int)(namend-cptr)) = 0;
	    		   }
	  		  else strcpy (outnam,cptr);	/* '.' not found in name */
	  		 strcat (outnam,".p");		/* append '.p' to name */
          		if((freopen(outnam,"w",stdout))==NULL)
           		  {
             		     fprintf(stderr,"posit: cannot open %s\n",outnam);
	     		     fflush (stderr);
             		     continue;
           		  }
		       }
		     doplot();
		     if (argc <= 1) break;
		   }
		 }				/* if (argc) */
		else doplot();
		if (iflag) break;
  	} while(argc>0);
}

void doplot(void)
{
	register int i;
	register int c;
	int nlines,nesc,offset,count,header;	

	pltout = stdout;
	while((c=getc(pltin))!= EOF)
	   {
		if (iflag) break;
		if ((c&PLTCMD) == 0 || ((c&CMDMSK) < 'a' || (c&CMDMSK) > 'z'))
		   {
		/*	if(!separate)
			   {
				putc(CUS,pltout);
				fflush(pltout); 
			   }
		*/
			do
			   {
				putc(c,pltout);
			   }  while(((c=getc(pltin))&PLTCMD) == 0);
			if( c == EOF ) break;
			ungetc(c,pltin);
			continue;
		   }

		putc(c,pltout);		/* send the command char */
		switch (c& CMDMSK)	/* command list */
		   {
			case 'm':		/* move */
			   	putwsx ((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
			   	putwsx ((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				break;
			case 'd':		/* draw */
			   	putwsx ((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
			   	putwsx ((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				break;
			case 'k':               /* circle */
				putwsx((int)(getwsx(pltin)/scale),pltout); /* radius */
				putc(getc(pltin),pltout);
				break;
			case 'l':               /* rectangle */
				putc(getc(pltin),pltout);		/* fill */
				putc(getc(pltin),pltout);		/* "m" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				putc(getc(pltin),pltout);		/* "d" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				putc(getc(pltin),pltout);		/* "d" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				putc(getc(pltin),pltout);		/* "d" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				break;
   
			case 'o':               /* fill */
				putc(getc(pltin),pltout);	/* fill */
				break;
			case 'v':               /* triangle */
				putc(getc(pltin),pltout);	/* fill */
				putc(getc(pltin),pltout);	/* "m" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				putc(getc(pltin),pltout);	/* "d" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				putc(getc(pltin),pltout);	/* "d" */
				putwsx((int)((getwsx(pltin)-XOFFS)/scale+xcent),pltout);
				putwsx((int)((getwsx(pltin)-YOFFS)/scale+ycent),pltout);
				break;
   
			case 's':		/* set up Text modes */
				putwsx(getwsx(pltin),pltout);
				putwsx(getwsx(pltin),pltout);
				break;
			case 't':		/* text */
				while ((c=getc(pltin)))	/* get string */
					putc(c,pltout);
				putc(c,pltout);		/* send the zero at end */
				break;
			case 'z':		    /* change dash pattern */
				putc(getc(pltin),pltout);
				break;
			case 'c':			/* change pen */
				if (npen) {	 	/* ignore if -p nonzero */
				   putc(npen,pltout);
				   getc(pltin);
				 }
				else
				   putc(getc(pltin),pltout);
				fflush (pltout);
				break;
			case 'b':		/* break */
			case 'e':		/* erase */
				fflush(pltout);
			case 'f':		/* fat */
				putwsx(getwsx(pltin),pltout);
				break;
			case 'r':		/* raster data */
				putc(getc(pltin),pltout);
				while (header=getwsx(pltin) != ENDRAS)
				  {
				    putwsx(header,pltout);
				    count = header & 0377;
				    offset = (header>>8)&0377;
				    if (count == 0377)		/* escape */
				      {  switch (offset)
					 {
					  case 't':
					  case 'T':
					   for(count=0;putc(getc(pltin),pltout)>0;)
						count++;
					   if (count&1)  putc(getc(pltin),pltout);
					   break;
					 }
				        nesc++;
				      }
				     else		/* raster line */
					{
					  while(count--) putwsx(getwsx(pltin),pltout);
					}
				  }
				break;
			case 'w':		/* window */
				putwsx(getwsx(pltin),pltout);
				putwsx(getwsx(pltin),pltout);
				putwsx(getwsx(pltin),pltout);
				putwsx(getwsx(pltin),pltout);
				break;
			case 'x':		/* display x-hairs */
				break;
			case 'p':		/* purge pltout buffers */
				fflush(pltout);
				break;
			case 'a':		/* text mode */
			case 'g':		/* graphics mode */
			case 'h':		/* hidden line commands */
			case 'i':
			case 'j':
			case 'n':		/* no op */
				break;
			default: 		/* error */
			   fprintf(stderr,"\nposit: invalid plot command %c\n",c);
			   fflush(stderr);
			   break;

			}	/* switch */
		}		/* while */
	fflush(pltout);
}


