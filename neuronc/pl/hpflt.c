#
/* Filter for HP 7221A plotter.  R.G.Smith 
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
#include	<sys/types.h>
#include	<sys/stat.h>

int istat;				/* look in "UNIX programming" */
int iflag = 0;				/* interrupt flag used by mdisplay */

#define BUFFER	0
#define CHARHT	22*21
#define PLTCMD	0200
#define CMDMSK	0177
#define INCHDOT	(8.0/1024.)
#define STDSIZ 16384	/* standard coordinate range */
#define STDSCL 	1.6384	/* STDSIZ/780 to scale standard plot */
#define TWMUL 	64	/* # of pixels in standard text width size spec. */
#define THMUL 	96	/* 3/2*TWMUL */
#define STDCW 	5	/* standard text size (about 15 x 23) */
#define XMAX 	15200	/* maximum x value */
#define XMIN 	0	/* minimum x value */
#define YMAX 	10000	/* maximum y value (manual says 779, but 780 fits) */
#define YMIN 	0	/* minimum y value */
#define XOFFS   (STDSIZ/2)  /* offset of standard picture */
#define YOFFS   (STDSIZ/2)  /* offset of standard picture */
#define XCENT   7500	/* center of plotter platen */
#define YCENT	5000
#define ENDRAS 	-1	/* end of raster data flag */
#define ERASE	014
#define ESC	033
#define CSZ	032
#define GS	035
#define CUS	037
#define GRKFONT 'g'	/* greek character set */
#define NPENS 16	/* max number of plotter pens */

int grkflag;		/* flag to return ascii chars */
int fatbase,fat;	/* line thickness parameters */
int xmax=XMAX;		/* plot window parameters defaulted */
int xmin=XMIN;		/* to maximum size 		*/
int ymax=YMAX;
int ymin=YMIN;
int xnew,ynew;		/* new pen location */
int xold,yold;		/* old pen location */
int xorigin,yorigin;	/* origin for text mode */
unsigned tangle, tsize;	/* new text angle, size */
char tfont;		/* future text font */
double scale;
int stop=0; /* 3;	/* if 1, stop before 'e', 'b' */
int ty;			/* counter for y pos. of straight text */
int separate;		/* output redirected */
int intract, quitact;	/* old signal actions */
int pen,npen,speed;	/* pen, speed for plotter */
int xcent,ycent;	/* position of center of picture */
int wait;		/* wait = 1 means wait for keypress before each pen */
int savef;		/* savef = 1 means save pen files on exit */
int rotflag;		/* rotflag = 1 means rotate picture sideways */

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

static char swnam[NPENS][20];		/* file names for pen sort bins */
static FILE *swfil[NPENS];		/* file pointers for pen sort bins */
static int pensort=0;			/* pen sort flag */

/*-------------------------------------*/

onintr()

{
  extern int iflag;

  iflag = 1; 
  signal (2,onintr);		/* go to "onintr" at ^C */
}

/*-------------------------------------*/

main(argc,argv)
int argc; char **argv;
   {
	FILE *freopen();
	double atof();
	float fscale;
	register char *cptr;

	iflag = 0;
   	istat = signal (2, 1);	/* save original status */
   	if (istat != 1)
  	   signal (2,onintr);		/* go to "onintr" at ^C */

/*	fstat(1, inbuf);
	fstat(2, outbuf);
	separate = (inbuf->st_ino == outbuf->st_ino)? 0: 1; */
	pltin= stdin;
	pltout= stdout;
	prinout= stderr;
	speed = 4;
	pen = 1;
	npen = 0;
	wait = 0;
	savef = 0;
	rotflag = 0;
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
					case 'c':
						pensort = !pensort;
						break;
					case 's':		/* speed */
						argv++; argc--;
						speed= atoi(*argv);
						if (speed <= 0) speed = 4;
						else if (speed > 30) speed = 30;
						break;
					case 'f':
						fatbase = atoi(++cptr);
						if (fatbase < 0)
						fatbase = 0;
						fat = fatbase * 10;
						break;
					case 'i':
						signal (2,1);
						break;
					case 'm':	/* mag */
						argv++; argc--;
						fscale= atof(*argv);
						if(fscale<=0.0)
						  scale= STDSCL;
					        else
						  scale= STDSCL/fscale;
						break;
					case 'x':	/* x center */
						xcent = atoi(*(++argv));
						argc--;
						break;
					case 'y':	/* y center */
						ycent = atoi(*(++argv));
						argc--;
						break;
					case 'p':	/* pen */
						argv++; argc--;
						npen = atoi (*argv);
						break;
					case 'r':  /* rotate pic sideways */ 
						rotflag = !rotflag;
						break;
					case 'w':  /* wait before each pen */
						wait = !wait;
						break;
					case 'z':  /* save pen files on exit */
						savef = !savef;
						break;
					default:
					break;  /* default */

				  }		/* switch */
			  }			/* if *cptr */
			else
			 {			/* *cptr not= '-' */

			   if((temp=freopen(cptr,"r",pltin)) == NULL)
			      {
				fprintf(stderr,"hpflt:cannot open %s\n",cptr);
				fflush(stderr);
				continue;
			      }
			   else	pltin= temp;
			   if (wait) filwait(cptr);	/* wait for keypress */
			   doplot();
			   if (argc <= 1) break;
			 }
		 }				/* if (argc) */
		else doplot();
		if (iflag) break;
  	} while(argc>0);
}

doplot()
   {
	register int i;
	register int c;
	static int ch;
	static int penused[NPENS+1]; 
	

	pltout = stdout;
	if (pensort)
	  {
	   for (i=0; i<NPENS+1; i++)
	     penused[i] = 0;
	   switinit();
	  }

	pset2();		/* set limits and grid size */
	pveloc(speed);		/* source in "plotsub.c" */
	if (npen) ppen(npen);
	else      ppen(pen);
	ty= YMAX-CHARHT;
	while((c=getc(pltin))!= EOF)
	   {
		if (iflag) break;
		if ((c&PLTCMD) == 0 || ((c&CMDMSK) < 'a' || (c&CMDMSK) > 'z'))
		   {
			if(!separate)
			   {
				plot(0, ty, 1);
		/*		putc(CUS,pltout); */
				fflush(pltout);
				stop &= 1;
			   }
			do
			   {
				if( c == '\n') ty -= CHARHT;
				putc(c,prinout);
			   }  while(((c=getc(pltin))&PLTCMD) == 0);
			if( c == EOF ) break;
			ungetc(c,pltin);
			continue;
		   }
		switch (c& CMDMSK)	/* command list */
		   {
			case 'm':		/* move */
			   xold= xorigin= (getwsx(pltin)-XOFFS)/scale;
			   yold= yorigin= (getwsx(pltin)-YOFFS)/scale;
				break;
			case 'd':		/* draw */
				xnew=(getwsx(pltin)-XOFFS)/scale;
				ynew=(getwsx(pltin)-YOFFS)/scale;
				if(fat)
				  {
					if(abs(xnew-xold) >= abs(ynew-yold))
						for(i= -(fat/2);i<=(fat+1)/2;i+=10)
						  dump(xold,yold+i,xnew,ynew+i);
					  else
						for(i= -(fat/2);i<=(fat+1)/2;i+=10)
						  dump(xold+i,yold,xnew+i,ynew);
				  }
				  else  dump(xold,yold,xnew,ynew);
				xold=xnew;
				yold=ynew;
				break;
			case 's':		/* set up Text modes */
				tsize = getc(pltin)&0377;
				tfont = getc(pltin)&0377;
				tangle = getwsx(pltin);
				if (tfont == GRKFONT) grkflag = 1;
				else grkflag = 0;
				break;
			case 't':		      /* text */
				putsym();
				break;
			case 'c':		      /* change pen */
				if (npen) break;      /* ignore if -p nonzero */
				pen = getc (pltin);
				pen = ((pen-1) % NPENS) + 1;
				if (pensort)		/* sort on pen number */
			         {
				  if (pen != 0)
				   {
				    penswitch (pen);
				    if (! penused[pen])
				     {
				/*	 penused[pen] = 1;  */
				         ppen (pen);
				     }
				   }
				  else		/* pen = 0 */
				   ppen (pen);
				 }
				else		/* not pensort */
				  ppen (pen);
				fflush (pltout);
				break;
			case 'b':		/* break */
			case 'e':		/* erase */
				fflush(pltout);
				if (stop == 1)
					read(2, outbuf, BUFSIZ);
				/* putc(ESC,pltout);
				putc(ERASE,pltout); */
				fflush(pltout);
			/*	sleep(2); */
				ty= YMAX - CHARHT;
				break;
			case 'f':		/* fat */
				fat= fatbase+ getwsx(pltin);
				fat *= 10;
				break;
			case 'r':		/* raster data */
				putras();
				break;
			case 'w':		/* window */
				xmin=((unsigned) getwsx(pltin))/scale;
				if(xmin<XMIN) xmin=XMIN;
				xmax=((unsigned) getwsx(pltin))/scale;
				if(xmax>XMAX) xmax=XMAX;
				ymin=((unsigned) getwsx(pltin))/scale;
				if(ymin<YMIN) ymin=YMIN;
				ymax=((unsigned) getwsx(pltin))/scale;
				if(ymax>YMAX) ymax=YMAX;
				break;
			case 'x':		/* display x-hairs */
				if (1) break;	/* ignore for 7221A (RGS) */
			
				fflush(pltout);
				gtty(fileno(pltout),&ttystat);
/* used to be:	ttystat.mode =| 040;	/* raw mode */
				if ((intract = (signal(2, 1))&1) == 0)
					signal(2, fixtty);	/* interrupts */
				if ((quitact = (signal(3, 1))&1) == 0)
					signal(3, fixtty);	/* quits */
				ttystat.mode &= ~010;	/* no echo */
				stty(fileno(pltout),&ttystat);
			    	do
				  {	putc(ESC,pltout);
					putc(CSZ,pltout);
					fflush(pltout);
					c= getc(pltin);	/* get next cmd */
				  } while (c == ('x'|PLTCMD));
/* used to be:	putc('\r', pltout);	/* get out of funny mode */
				fflush(pltout);
				ungetc(c,pltin);
		 		ttystat.mode &= ~040;	/* turn off raw */
				ttystat.mode |= 010;
				stty(fileno(pltout),&ttystat);
				signal(2, intract);
				signal(3, quitact);
				break;
			case 'p':		/* purge pltout buffers */
				fflush(pltout);
				break;
			case 'z':		/* dash command */
				c= getc(pltin);	/* get next cmd */
				break;
			case 'h':		/* hidden line commands */
			case 'i':
			case 'j':
			case 'n':		/* no op */
			case 'g':
			case 'a':		/* no op */
				break;
			default: 		/* error */
				fprintf(stderr,"\nhpflt: invalid plot command %c\n",c);
				fflush(stderr);

			}	/* switch */
		}		/* while */
/*
 * finish up
 */
	if (pensort)
	 {
	  switclos();
	  for (i=0; i<NPENS; i++)
	   {				/* copy all four pen files to stdout */
	    if (wait) penwait(i);	/* wait for user to type key  */
	    pltin = swfil[i];
	    while ((ch = getc(pltin)) != EOF) putc(ch,stdout);
	    if (wait)
	     {				/* after each pen is used: */
	       ppen (0);		/* put pen back in stall */
	       papnt (15200,10000);	/* go to corner */
	       ppen (0);		/* dummy pen command to flush output */
	     }
	    ppflush();
	    fflush(pltout);
	   }
	  switend();
	 }

        ppen (0);			/* put the pen back */
	papnt (15200,10000);		/* go to corner */
	ppen (0);			/* dummy pen command to flush output */
        ppflush();
	fflush(pltout);
}


filwait(fil)
     char *fil;
{
   FILE *ttyfil;

  fprintf (stderr,"Please type <RET> for file '%s'\n",fil);
  ttyfil = fopen ("/dev/tty","r");
  getc (ttyfil);
  fclose (ttyfil);
}

penwait(pen)
     int pen;
{
   FILE *ttyfil;

  fprintf (stderr,"Please type <RET> for pen %d\n",pen+1);
  ttyfil = fopen ("/dev/tty","r");
  getc (ttyfil);
  fclose (ttyfil);
}


dump(x1,y1,x2,y2)	/* send a vector, if possible */
int x1,y1,x2,y2;
   {
	static int xlst, ylst, temp;
	extern int xcent,ycent;

	if (rotflag)
	  {
		temp = x1;
		x1 = -y1;
		y1 = temp;
		temp = x2;
		x2 = -y2;
		y2 = temp;
	  } 
	x1 += xcent;
	y1 += ycent;
	x2 += xcent;
	y2 += ycent;
	if( clip(&x1,&y1,&x2,&y2) ) return;
	if( (x1 != xlst) || (y1 != ylst) )
		plot(x1,y1,1);
	plot(x2,y2,0);
	xlst= x2;
	ylst= y2;
   }


plot(x,y,icode)
int x,y,icode;
   {
	if (icode)	/* move pen to abs loc. */
	   {
		papnt (x,y);
	   }
	 else		/* draw line to abs loc. */
	   {
		pdraw (x,y);
	   }
   }


/* code for clip */
#define code(x,y) (x<xmin?1:(x>xmax?2:0))|(y<ymin?4:(y>ymax?8:0))

clip(x1,y1,x2,y2)	/* window the plot */
int *x1,*y1,*x2,*y2;
   {
	register int c1,c2,temp;
	int swap;
	c1=code(*x1,*y1);
	c2=code(*x2,*y2);
	swap=0;
	if(!(c1||c2)) return(0); /* line completely in bounds */
	while(c1|c2)
	   {
		if( c1&c2 ) return(1);  /* line completely out of bounds */
		if(!c1)	/* interchange endpoints */
		   {
			temp= *x1;*x1= *x2;*x2=temp;
			temp= *y1;*y1= *y2;*y2=temp;
			temp=c1;c1=c2;c2=temp;
			swap= ~swap;
		   }
		if(c1<4)	/* move endpoint in x */
		   {
			temp=(c1&2?xmax:xmin);
			*y1= solve(temp,*x1,*y1,*x2,*y2);
			*x1=temp;
		   }
		  else		/* move endpoint in y */
		   {
			temp=(c1&8?ymax:ymin);
			*x1= solve(temp,*y1,*x1,*y2,*x2);
			*y1=temp;
		   }
		c1=code(*x1,*y1);
	   }
	if( swap )	/* put endpoints in order */
	   {
		temp= *x1; *x1= *x2; *x2=temp;
		temp= *y1; *y1= *y2; *y2=temp;
	   }
	return(0);
   }




solve(pnot,p1,q1,p2,q2) /* solve linear eqn in integers */
int pnot,p1,q1,p2,q2;
   {
	register int pmid,qmid;
	if(p1>p2)
	   {
		pmid=p1;p1=p2;p2=pmid;
		qmid=q1;q1=q2;q2=qmid;
	   }
	if(pnot<=p1) return(q1);
	if(pnot>=p2) return(q2);
	while(1)  /* iterate until convergence */
	   {
		pmid= (p1+p2)>>1;
		qmid= (q1+q2)>>1;
		if(pmid<pnot) {p1=pmid;q1=qmid;}
		 else if(pmid>pnot) {p2=pmid;q2=qmid;}
			else return(qmid);
	   }
   }


/*
 * plot maximum enclosing box for raster data
 */
putras()
   {
	register int header,count,y;
	int rxmin, rymin, rxmax, rymax, offset, dotsize, raster;

	dotsize = getc(pltin);
	if (dotsize == 0)
		dotsize = STDSIZ/2048;	/* mp resolution */
	raster = 0;
	rxmin= rxmax= xold;
	rymin= 32767;
	rymax= yold;
	while( (header=getwsx(pltin)) != ENDRAS )
	  {
		count= header&0377;
		offset= (header>>8)&0377;
		if (count == 0377)	/* escape */
		  {	switch (offset)
			  {
			  case 't':
			  case 'T':
				for (count=0; getc(pltin)>0; count++)
					;
				if (count&1)
					getc(pltin);
				break;
			  }
			continue;
		  }
		y = yold + (offset*16*dotsize)/STDSCL;
		if(rymin > y) rymin=y;
		y += (count*16*dotsize)/STDSCL;
		if(rymax < y) rymax=y;
		while(count--) getwsx(pltin);
		raster++;
	  }
	rxmax += (raster*dotsize)/STDSCL + 1;
	dump(rxmin,rymin,rxmin,rymax);
	dump(rxmin,rymax,rxmax,rymax);
	dump(rxmax,rymax,rxmax,rymin);
	dump(rxmax,rymin,rxmin,rymin);
   }

/*
 * Restore tty mode and exit.
 */
fixtty()
  {
	ttystat.mode &= ~040;	/* turn off raw */
	ttystat.mode |= 010;	/* turn on echo */
	stty(fileno(pltout),&ttystat);
	exit(1);
  }

 psend(ch)
   char ch;

/* Send char to printer routines.
  Needed by "PLOTSUB" library. (RGS) */

  {
    putc (ch, pltout);
  } 

/*-------------------------------------*/

switinit()

/* open four files, one for each pen */

{
   static char template[] = {"/tmp/penXXXXXX"};
   static int i;

 for (i=0; i<NPENS; i++)
  {
   strcpy (swnam[i], template);
   swfil[i] = fopen (mktemp (swnam[i]),"w");
  }
}

/*-------------------------------------*/

switclos()

/* close files then reopen for reads */

{
  int i;

 for (i=0; i<NPENS; i++)
  {
   fclose (swfil[i]);
   swfil[i] = fopen (swnam[i],"r");
  }
 pltout = stdout;
}

/*-------------------------------------*/

switend()

/* erase files when done */

{
  static int i;

 for (i=0; i<NPENS; i++)
  {
   if (!savef) {		/* erase files except when -z flag set */
      fclose (swfil[i]);
      unlink (swnam[i]);
    }
  }
 pltout = stdout;
}

/*-------------------------------------*/

penswitch (p)
    int p;

/* switch "pltout" to a different
output file */

{
  extern FILE *pltout;
  extern FILE *swfil[NPENS];

 if (p < 1) p = 1;
 if (p > NPENS) p = ((p-1) % NPENS)+1;
 p--;				/* pen range is 1-8; files are 0-7 */
 fflush (pltout);
 pltout = swfil[p];
}

/*-------------------------------------*/

#include "putsym.h"
