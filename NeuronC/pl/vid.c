#
#
/*
 * vector plot filter for Tektronix 4012
 *
 * author- Rob Clayton, Dept. of Geophysics, Stanford Univ.
 * Modified for Version 7, 12/3/79
 * Mitchell Gart, Ampex Corp.
 *
 * Modified for Venix on IBM-PC	with either:
 *   Monochrome, Enhanced graphics, or Hercules boards.
 *
 *	Oct-85   R.G.Smith
 *
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<unistd.h>
#include	<signal.h>

#include	"colors.h"
#include 	"vdef.h"

#define BUFFER	1
#define CHARHT	22
#define PLTCMD	0200
#define CMDMSK	0177
#define INCHDOT	(8.0/1024.)
#define HPSCAL  (10000./16384.)	/* scale between HP plotter and square scrn */
#define STDSIZ  16384	/* standard coordinate range */
#define STDSCL  1	/* STDSIZ/350 to scale standard plot */
#define XOFFS  (STDSIZ/2) /* offset of standard picture */
#define YOFFS  (STDSIZ/2) /* offset of standard picture */
#define TWMUL 	64	/* # of pixels in standard text width size spec. */
#define THMUL 	96	/* 3/2*TWMUL */
#define STDCW 	5	/* standard text size (about 15 x 23) */
#define XMAX    22599	/* maximum x value */
#define XMIN 	0	/* minimum x value */
#define YMAX 	16383	/* maximum y value (man says 779,but 780 fits) */
#define YMIN 	0	/* minimum y value */
#define XCENT	8192	/* X center of picture  (480/348) * 8192 */
#define YCENT 	8192	/* Y center of picture */
#define ENDRAS 	-1	/* end of raster data flag */
#define ERASE	014
#define ESC	033
#define CSZ	032
#define GS	035
#define CUS	037
#define GRKFONT 'g'	/* greek font for 's' command */

int grkflag=0;		/* flag to return ascii chars */
int fatbase,fat;	/* line thickness parameters */
int xmax=XMAX;		/* plot window parameters defaulted */
int xmin=XMIN;		/* to maximum size 		*/
int ymax=YMAX;
int ymin=YMIN;
int xnew,ynew;		/* new pen location */
int xold,yold;		/* old pen location */
static int xlst, ylst;  /* last pos for dump() */

int xorigin,yorigin;	/* origin for text mode */
unsigned tangle, tsize;	/* new text angle, size */
char tfont;		/* future text font */
double scale=1.0;
int ty;			/* counter for y pos. of straight text */
int separate;		/* output redirected */
int intract, quitact;	/* old signal actions */
char *display_name=NULL;/* the display name for X windows */

static int xcent=0;
static int ycent=0;	/* center of picture */
static int erasflg=1;	/* erasflg = 1 means erase before displaying pict */
static int clear=1;	/* clear = 1 means allow picture clear by file */
static int clst=0;	/* last color for dump() */

char *pagefilenam = 0;		/* file name for mprintc */
char *comment = 0;		/* comment for mprinta */
int multfl = 0;			/* =1 -> multiple displays allowed (mprintx) */
int scrrev = 0;			/* =1 -> screen display upside down */
int noexit = 0;			/* =1 -> don't exit at end of picture */ 
double rastres = 150.0;		/* raster resolution in dpi */
double scrsiz = 1.0;		/* size of X-win windows display */
double fscale=1.0;		/* scale from "-m" command line switch */
double pageinterval;		/* interval for showing movie frames */
int radius=0;			/* radius for circle */
int scrntyp = SCRNTYP;		/* Monochrome, Enhanced or Hercules */
int scrcolor = WHITE;  		/* color is initially set to white */
int rotflag=0;			/* rotflag = 1 means rotate picture sideways */
static int color = BLACK;  	/* color is initially set to black */
static int pen = -1;		/* color override */
int backgr = -1;		/* background color override */
int bitres = 10;		/* size of pixel array for "mprinttf.c" */

static char tbuf[256] = {0};	/* temp buffer for printing string */
static char *ptbuf;
static char *progname = {0};

struct
   {
	char ispeed, ospeed;
	char erase, kill;
	int mode;
   } ttystat;

int fixtty(void);
FILE *pltout, *pltin, *prinout, *temp;
FILE *fopen(const char *, const char *);
int pargc; char **pargv;
void printhelp();
void dump(int x1, int y1, int x2, int y2);
void plot(int x, int y, int icode);
int xclip(int *x1, int *y1, int *x2, int *y2);
void putras(void);;
void doplot(void);
int solve(int pnot, int p1, int q1, int p2, int q2);

void txtevent(int arg);
void txtmode (int arg);
int txsize (void), tysize(void);
void txtclr(void);
void txtcolr(int color);
int getwsx(FILE *stream);
void tdrcirc(int x, int y, int rad, int fill);
void tdrrect(int x1,int y1,int x2,int y2,int x3,int y3,int x4,int y4,int fill);
void tfill(int fill);
void tdrtri(int x1,int y1,int x2,int y2,int x3,int y3,int fill);
void tmove(int x, int y);
void tdraw(int x, int y);
void txtstr(char *str);
void txtpage(void);
void tgrdump(void);

/*-----------------------------------------*/

int main(int argc, char **argv)
{
	FILE *freopen(const char *, const char *, FILE *);
	double atof(const char *);
	register char *cptr;

	pargc = argc;
	pargv = argv;
	pltin= stdin;
	pltout= stdout;
	prinout= stdout;
	progname = argv[0];
	/*
	 * the following line is needed for an adm-3a with a RG-512
	 * card.  We have to take the adm out of "adm alpha mode."
	 */
/*	putc(GS, pltout); */
	scale= STDSCL;

	do	/* loop over arg list */
	   {
		argc--; argv++;
		cptr= *argv;
		if(argc) {

		 if (*cptr == '-') {
			cptr++;
			switch(*cptr)
			  {
			  case 'a':
				scrntyp = POSTSCRIPT;
				break;
			  case 'b': 
				scrntyp = DASHPS;
				break;
			  case 'c': 
				scrntyp = COLORPS;
				break;
			  case 'd': 
				display_name = *(++argv); argc--;
				break;
			  case 'e': 
				scrntyp = ENHANCED;
				break;
			  case 'F': 
				scrntyp = XFIG;
				break;
			  case 'g': 
				scrntyp = VGA;
				break;
			  case 'h': printhelp(); /* print help table */
				break;
			  case 'H': 
				scrntyp = HERCULES;
				break;
			  case 'L':	
				comment = *(++argv); argc--;
				break;
			  case 'X':	
				scrntyp = XWIN;
				break;
			  case 'M':	
				multfl = !multfl;
				break;
			  case 'j':		/* screen upside down */
				noexit = ! noexit;
				break;
			  case 'k':		/* screen upside down */
				scrrev = ! scrrev;
				break;
			  case 'l':
				scrntyp = LASERJET;
				break;
			  case 't':
				scrntyp = TX4014;
				break;
			  case 's':		
				++argv; argc--;
				bitres = atoi(*argv);
				break;
			  case 'f':
				fatbase = atoi(++cptr);
				if (fatbase < 0)
					fatbase = 0;
				fat = fatbase;
				break;
			  case 'i':
				/* signal(SIGINT, ((__sighandler_t)1));  */
				break;
			  case 'm':
				cptr = *(++argv); argc--;
				fscale= atof(cptr);
				if(fscale<=0.0) {
					fscale = 1.0;
					scale= STDSCL;
				}
				else
					scale= STDSCL/fscale;
				break;
			   case 'x':	/* x center */
				xcent = atoi(*(++argv)) / HPSCAL;
				argc--;
				break;
			   case 'y':	/* y center */
				ycent = atoi(*(++argv)) / HPSCAL;
				argc--;
				break;
			   case 'p':	/* pen */
				pen = atoi (*(++argv));
				argc--;
				break;
			   case 'B':	/* background color */
				backgr = atoi (*(++argv));
				argc--;
				break;
			   case 'r':	/* rotate */
				rotflag = !rotflag;
				break;
			   case 'P':
				pagefilenam = *(++argv); argc--;
				break;
			   case 'R':	/* raster resolution, i.e 150 dpi */
				rastres = atof(*(++argv));
				argc--;
				break;
			   case 'S':	/* interval for showing movie frames */
				pageinterval = atof(*(++argv));
				argc--;
				break;
			   case 'v':	/* erase */
				erasflg = !erasflg;
				break;
			   case 'w':	/* window size  */
				scrsiz *= atof(*(++argv));
				argc--;
				break;
			   case 'z':	/* no clear */
				clear = !clear;
				break;
			  default:
				break;
			  }	/* switch */

			continue;
		   }	/* if (*cptr) */

		 else {
		    if((temp=fopen(cptr,"r")) == NULL) {
			fprintf(stderr,"vid:cannot open %s\n",cptr);
			continue;
		    }
		    else pltin= temp;

		    doplot();
		    if (argc <= 1) break;
		 }
	     }			/* if (argc) */

	     else {		/* if (! argc) */
		doplot();			/* no files, only stdin */

	    }		/* else if (!argc) */
	  } while(argc);

	switch (scrntyp) {

	  case XWIN:    txtevent(1); break;

	  case XFIG:    
	  case DASHPS:    
	  case COLORPS:    
	  case POSTSCRIPT: txtmode (STDTXT); break;

		  default: break;
	}
        if (noexit) while (1);
	exit(0);
   }

void printhelp()

{
  printf ("%s version #%-4.3g\n",progname,4.05);
  printf ("vid -a     make 'PostScript' file on output: 'vid -a >file'\n");
  printf ("    -b     make dashed PS file suitable for xeroxing\n");
  printf ("    -c     make color PS plot: 'vid -c >file'\n");
  printf ("    -h     display this help message.\n");
  printf ("    -l     make 'PostScript raster' file in stdout.\n");
  printf ("    -e     display to Enhanced Graphics Adapter (EGA).\n");
  printf ("    -H     display to monochrome Hercules board.\n");
  printf ("    -t     display to Tektronix 4014 type display terminal.\n");
  printf ("    -X     use X-windows for graphics display.\n");
  printf ("    -d name define name of X-windows screen to display on.\n");
  printf ("    -k     invert screen (upside down).\n");
  printf ("    -m n   use 'n' as mag, default 1.0.\n");
  printf ("    -p n   use 'n' as color (16 colors)\n");
  printf ("    -B n   use 'n' as background color\n");
  printf ("    -P name  make individual .ps frames, use with -c \n");
  printf ("    -r     rotate picture sideways 90 deg.\n");
  printf ("    -S n   set page display inteval (sec)\n");
  printf ("    -w n   use 'n' as size of X-window (1.0 = half of screen).\n");
  printf ("    -x n   use 'n' as X offset (default 8192). Defines center.\n");
  printf ("    -y n   use 'n' as Y offset (default 8192). Defines center.\n");
  printf ("    -v     don't erase screen when first starting graphics.\n");
  exit(0);
}

void doplot(void)
{
	register int i;
	register int c;
	int fill;

	txtmode (GRAPHICS); 
	if (! xcent) xcent = txsize() / 2;
	if (! ycent) ycent = tysize() / 2;
	if (erasflg) txtclr(); 
        erasflg = 0;
	if (pen>=0) color = pen;	/* default color override */
	txtcolr(color);
	ty= YMAX-CHARHT;
	while((c=getc(pltin))!= EOF)
	   {
		if ((c&PLTCMD) == 0 || ((c&CMDMSK) < 'a' || (c&CMDMSK) > 'z'))
		   {
			do
			   {
				if( c == '\n') ty -= CHARHT;
/*				putc(c,prinout);  */
			   }  while(((c=getc(pltin)) & PLTCMD) == 0);
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
						for(i= -(fat/2);i<=(fat+1)/2;i++)
						  dump(xold,yold+i,xnew,ynew+i);
					  else
						for(i= -(fat/2);i<=(fat+1)/2;i++)
						  dump(xold+i,yold,xnew+i,ynew);
				  }
				  else  dump(xold,yold,xnew,ynew);
				xold=xnew;
				yold=ynew;
				break;
			case 'k':		/* circle */
			        radius = (getwsx(pltin))/scale;
				fill = getc(pltin)&0377;
				tdrcirc (xold+xcent,yold+ycent,radius,fill);
				xlst = xold;	
				ylst = yold;	
				break;
			case 'l':		/* rectangle */
				{
				int x1,y1,x2,y2,x3,y3,x4,y4,temp;
 
				fill = getc(pltin)&0377;
				temp = getc(pltin)&0377;     /* "m" */
				x1 = (getwsx(pltin)-XOFFS)/scale;
				y1 = (getwsx(pltin)-YOFFS)/scale;
				temp = getc(pltin)&0377;     /* "d" */
				x2 = (getwsx(pltin)-XOFFS)/scale;
				y2 = (getwsx(pltin)-YOFFS)/scale;
				temp = getc(pltin)&0377;     /* "d" */
				x3 = (getwsx(pltin)-XOFFS)/scale;
				y3 = (getwsx(pltin)-YOFFS)/scale;
				temp = getc(pltin)&0377;     /* "d" */
				x4 = (getwsx(pltin)-XOFFS)/scale;
				y4 = (getwsx(pltin)-YOFFS)/scale;
				x1 += xcent; y1 += ycent;
				x2 += xcent; y2 += ycent;
				x3 += xcent; y3 += ycent;
				x4 += xcent; y4 += ycent;
				tdrrect (x1,y1,x2,y2,x3,y3,x4,y4,fill);
				xlst = x1;
				ylst = y1;
				}

				break;
			case 'o':		/* fill */
				fill = getc(pltin)&0377;
				tfill (fill);
				xlst = ylst = 0;
				break;
			case 'v':		/* triangle */
				{
				int x1,y1,x2,y2,x3,y3,temp;
 
				fill = getc(pltin)&0377;
				temp = getc(pltin)&0377;     /* "m" */
				x1 = (getwsx(pltin)-XOFFS)/scale;
				y1 = (getwsx(pltin)-YOFFS)/scale;
				temp = getc(pltin)&0377;     /* "d" */
				x2 = (getwsx(pltin)-XOFFS)/scale;
				y2 = (getwsx(pltin)-YOFFS)/scale;
				temp = getc(pltin)&0377;     /* "d" */
				x3 = (getwsx(pltin)-XOFFS)/scale;
				y3 = (getwsx(pltin)-YOFFS)/scale;
				x1 += xcent; y1 += ycent;
				x2 += xcent; y2 += ycent;
				x3 += xcent; y3 += ycent;
				tdrtri (x1,y1,x2,y2,x3,y3,fill);
				xlst = x1;
				ylst = y1;
				}

				break;
			case 's':		/* set up Text modes */
				tsize = (getc(pltin)&0377) * fscale;
				tfont = getc(pltin)&0377;
				tangle = getwsx(pltin);
				if (tfont == GRKFONT) grkflag=1;
				else grkflag = 0;
				break;
			case 't':		/* text */
                        	for(ptbuf=tbuf; c=getc(pltin); ) {
				  *ptbuf++ = c;
                          	}
				*ptbuf++ = 0;
				{
				  int temp,x1,y1;

			 	  x1 = xold;
				  y1 = yold;
				  if (rotflag && (scrntyp !=POSTSCRIPT)
					      && (scrntyp !=DASHPS) 
					      && (scrntyp !=COLORPS) )
	  				{
						temp = y1;
						x1 = - y1;
						y1 = temp;
	  				} 
				  x1 += xcent;
				  y1 += ycent;
			 	  tmove(x1,y1);
				}
				txtstr(tbuf);
				break;
			case 'c':		/* change plotter pen */
				color = getwsx(pltin); /* 2 char color */
				/* else {		/* color is 1 char */
				/*   char *lowcol,col;
				   lowcol = (char *)&color;
				   col = *lowcol;
				   color = col;
				} */
				txtcolr (color);
				break;
			case 'g':if (!(scrntyp&(XWIN|POSTSCRIPT|DASHPS
					   |COLORPS|LASERJET|XFIG))) {
				   txtmode(GRAPHICS);  /* graphics mode */
			        }
				break;
			case 'a':if (!(scrntyp&(XWIN|POSTSCRIPT|DASHPS
					   |COLORPS|LASERJET|XFIG))) {
				   txtmode(STDTXT);    /* change to text mode*/
				}
				break;
			case 'h':		/* hidden line commands */
			case 'i':
			case 'j':
			case 'b':		/* break */
				break;
			case 'e':		/* erase */
				if (!clear) break;
				txtclr();
		/*		fflush(pltout);  */
		/*		putc(ESC,pltout);
				putc(ERASE,pltout);
				fflush(pltout);  */
			/*	sleep(2); */
				ty= YMAX - CHARHT;
				break;
			case 'f':		/* fat */
				fat= fatbase+ getwsx(pltin);
				break;
			case 'r':		/* raster data */
				putras();
				break;		/* */

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
		/*		fflush(pltout);  */
		/*		gtty(fileno(pltout),&ttystat);  */
/* used to be:	ttystat.mode =| 040;	/* raw mode */
		/*		if ((intract = (signal(2, 1))&1) == 0)
					signal(2, fixtty);	/* interrupts */
		/*		if ((quitact = (signal(3, 1))&1) == 0)
					signal(3, fixtty);	/* quits */
		/*		ttystat.mode &= ~010;	/* no echo */
		/*		stty(fileno(pltout),&ttystat);  */
		/*	    	do
				  {	putc(ESC,pltout);
					putc(CSZ,pltout);
					fflush(pltout);
					c= getc(pltin);	/* get next cmd */
		/*		  } while (c == ('x'|PLTCMD)); */
/* used to be:	putc('\r', pltout);	/* get out of funny mode */
		/*		fflush(pltout);  */
				ungetc(c,pltin);
		 		ttystat.mode &= ~040;	/* turn off raw */
				ttystat.mode |= 010;
		/*		stty(fileno(pltout),&ttystat); */
		/*		signal(2, intract);
				signal(3, quitact); */
				break;
			case 'p':		/* purge pltout buffers */
				if (scrntyp & TX4014) {
 					putc (CUS,pltout);
					fflush(pltout);
				}
				clst = 0;	/* reset polyline in dump */
				break;
			case 'z':		/* no op */
				 getc(pltin);
			case 'n':
				txtpage();	/* new page */
				if (pageinterval>0) usleep ((long int)(pageinterval*1e6));
				break;
			default: 		/* error */
			   fprintf(stderr,"\nvid: invalid plot command %c\n",c);
				fflush(stderr);
			/*	exit(1); */
			}
 		if (scrntyp==XWIN) txtevent(0); 
		}
	tgrdump();   		/* output to printer */
}


void dump(int x1, int y1, int x2, int y2)
{
    static int temp;

	if (rotflag && (! (scrntyp&(POSTSCRIPT | COLORPS | DASHPS | XFIG))))
	  {
		temp = x1;
		x1 = - y1;
		y1 = temp;
		temp = x2;
		x2 = - y2;
		y2 = temp;
	  } 
	x1 += xcent;
	y1 += ycent;
	x2 += xcent;
	y2 += ycent;


	if( xclip(&x1,&y1,&x2,&y2) ) return;  
   	 if( (x1 != xlst) || (y1 != ylst) || (color != clst) ) /* */
	  plot(x1,y1,1);
	plot(x2,y2,0);
	xlst= x2;
	ylst= y2;
        clst = color;
   }


void plot(int x, int y, int icode)
{
	/*
	 * this code takes advantage of
	 * vector addresses.
	 */
	register char tempx, tempy;
	static char xhi,yhi,ylo;

	   if (icode)			/* move pen */
	        tmove (x,y);
   	   else tdraw (x,y);
}


/* code for clip */
#define code(x,y) (x<xmin?1:(x>xmax?2:0))|(y<ymin?4:(y>ymax?8:0))

int xclip(int *x1, int *y1, int *x2, int *y2)	/* window the plot */
                    
   {
	int c1,c2,temp;
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




int solve(int pnot, int p1, int q1, int p2, int q2) /* solve linear eqn in integers */
                     
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
		pmid= ((long)p1+p2)>>1;
		qmid= ((long)q1+q2)>>1;
		if(pmid<pnot) {p1=pmid;q1=qmid;}
		 else if(pmid>pnot) {p2=pmid;q2=qmid;}
			else return(qmid);
	   }
   }


/*
 * plot maximum enclosing box for raster data
 */
void putras(void)
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
int fixtty(void)
{
/*	ttystat.mode =& ~040;	/* turn off raw */
/*	ttystat.mode =| 010;	/* turn on echo */
/*	stty(fileno(pltout),&ttystat);
	exit(1); */
  }

/* #include "putsym.h" */

