/* Module mprinta in program vid */

/* Black-and-White Post-Script device driver.  

This module with  device-independent subroutines in "mvid.c".
Nothing real fancy here.  All this program will do
for now is dump out on standard output.  We can redirect
at a later stage.  Let's just see it work for now.

To reduce output filesize, I've added some macros for move and draw.
Feb, 12, 1992 --> Andrew

Modified for black-and-white Postscript, R.G.Smith   July 1992
 Picture is identical to that produced by "mprintc.c"
 except all lines are black.

*/

#include	<stdio.h>
#include	<string.h>
#include	<time.h>
#include	<sys/types.h>
/* #include	<sys/timeb.h> */

/* Define image to be 7.5 x 7.5; User can imbed postscript commands to shrink
   image to any size if needed.
*/

#define XSCALE 30.339    /* STDSIZ/(7.5*72); LW is 72 dpi */
#define YSCALE 30.339
#define XHHI 750
#define XHI 750
#define YHI 675

#define STDSIZ 16383.
#define GRKFONT 'g'	/* greek font for 's' command */
#define NCOLORS 256

static double xpos, ypos;
static double xscale = XSCALE;
static double yscale = YSCALE;
static double yhi = YHI;

extern int rotflag;
extern char tfont;
extern unsigned int tsize,tangle;
extern char *comment;
extern int scrrev;
extern int backgr;
extern int pargc;
extern char **pargv;
int gmodea(int argc, char **argv);
int setccol(int col);

/* --------------------------------- */

int tainit(void)
{
gmodea(pargc,pargv);
}

void taexit(void)
{
   printf("stroke showpage\n");
   printf("restore\n");
/*   printf("clear\n"); */
}
void clrca(void)
{}
void xycona(void)
{}

void xyvida(void)
{}

void adot(void)
{}

int setacol(int pen)
{
    double grayval;
    static int oldcolor,scrcolor;

   oldcolor = scrcolor;
   scrcolor = pen;

    printf("stroke\n");
    printf("0 setgray\n");
    switch (pen) {
	case 0:
	    printf("1 setlinewidth\n");			/* black */
	    break;
	case 1:
	    printf("1 setlinewidth\n");			/* blue */
	    break;
	case 2:
	    printf("1 setlinewidth\n");			/* green */
	    break;
	case 3:
	    printf("1 setlinewidth\n");			/* cyan */
	    break;
	case 4:
	    printf("1 setlinewidth\n");			/* red */
	    break;
	case 5:
	    printf("1 setlinewidth\n");			/* magenta */
	    break;
	case 6:
	    printf("1 setlinewidth\n");			/* brown */
	    break;
	case 7:
	    printf("1 setlinewidth\n");			/* white */
	    break;
	case 8:
	    printf("1.5 setlinewidth\n");		/* gray */
	    break;
	case 9:		
	    printf("1.5 setlinewidth\n");		/* light blue */
	    break;
	case 10:
	    printf("1.5 setlinewidth\n");		/* light green */
	    break;
	case 11:
	    printf("1.5 setlinewidth\n");		/* light cyan */
	    break;
	case 12:
	    printf("1.5 setlinewidth\n");		/* light red */
	    break;
	case 13:
	    printf("1.5 setlinewidth\n");		/* light magenta */
	    break;
	case 14:
	    printf("1.5 setlinewidth\n");		/* yellow */
	    break;
	case 15:
	    printf("1.5 setlinewidth\n");		/* bright white */
	    break;
	default :
	    grayval = ((double)pen)/NCOLORS;
            printf("%g setgray\n",1.0-grayval);
	    break;	
	}
  return (oldcolor);
}

void movea(int ix, int iy)
{
 extern int scrrev;

/* if (rotflag) {
   xscale = 22.754;  /* STDSIZ/(10 * 72) */
/*   yscale = 22.754; 
 }
*/

 xpos = ix/xscale;
 ypos = iy/yscale;
 if (xpos >XHHI) xpos=XHHI;
 if (xpos <0) xpos=0;
 if (ypos >yhi) ypos=yhi;
 if (ypos <0) ypos=0;
 if (scrrev) ypos = yhi-ypos;
 printf("s\n%-5.4g %-5.4g m\n",xpos,ypos);
}

void drawa(int ix, int iy)
{
 extern int scrrrev;
 double x,y;

 x = ix/xscale;
 y = iy/yscale;
 if  (scrrev) y=yhi-y;
 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;
 printf("%-5.4g %-5.4g l\n",x,y);
 xpos=x;
 ypos=y;
}

void tadrcirc(int ix, int iy, int irad, int fill)
{
 double rad,x,y;

 x = ix/yscale;
 y = iy/yscale;
 rad = irad/yscale;
 if  (scrrev) y=yhi-y;
 if (x>XHHI) x=XHHI;
 if (x<0) x=0;
 if (y>yhi) y=yhi;
 if (y<0) y=0;
 printf("s %-5.4g %-5.4g %-5.4g 0 360 arc",x,y,rad);
 if (!fill) printf(" s\n");
 else       printf(" fill\n");
 xpos = x;
 ypos = y;
}

void tadrrect(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int ix4, int iy4, int fill)
{
 extern int scrrev;
 double x1,y1,x2,y2,x3,y3,x4,y4;
 int color;

 x1 = ix1/xscale;
 y1 = iy1/yscale;
 x2 = ix2/xscale;
 y2 = iy2/yscale;
 x3 = ix3/xscale;
 y3 = iy3/yscale;
 x4 = ix4/xscale;
 y4 = iy4/yscale;

 if  (scrrev) {
   y1=yhi-y1;
   y2=yhi-y2;
   y3=yhi-y3;
   y4=yhi-y4;
 }
 if (!fill && backgr>=0) {
   color = setccol(backgr);
   printf("%-5.4g %-5.4g m\n",x1,y1);
   printf("%-5.4g %-5.4g l\n",x2,y2);
   printf("%-5.4g %-5.4g l\n",x3,y3);
   printf("%-5.4g %-5.4g l\n",x4,y4);
   printf("%-5.4g %-5.4g l\n",x1,y1);
   printf(" fill\n");
   setccol(color);
 }
 printf("%-5.4g %-5.4g m\n",x1,y1);
 printf("%-5.4g %-5.4g l\n",x2,y2);
 printf("%-5.4g %-5.4g l\n",x3,y3);
 printf("%-5.4g %-5.4g l\n",x4,y4);
 printf("%-5.4g %-5.4g l\n",x1,y1);
 if (!fill) printf("s\n");
 else       printf("fill\n");
 xpos = x1;
 ypos = y1;
}

void tadrtri(int ix1, int iy1, int ix2, int iy2, int ix3, int iy3, int fill)
{
 extern int scrrev;
 double x1,y1,x2,y2,x3,y3;
 int color;

 x1 = ix1/xscale;
 y1 = iy1/yscale;
 x2 = ix2/xscale;
 y2 = iy2/yscale;
 x3 = ix3/xscale;
 y3 = iy3/yscale;

 if  (scrrev) {
   y1=yhi-y1;
   y2=yhi-y2;
   y3=yhi-y3;
 }

 if (!fill && backgr>=0) {
   color = setccol(backgr);
   printf("%-5.4g %-5.4g m\n",x1,y1);
   printf("%-5.4g %-5.4g l\n",x2,y2);
   printf("%-5.4g %-5.4g l\n",x3,y3);
   printf("%-5.4g %-5.4g l\n",x1,y1);
   printf(" fill\n");
   setccol (color);
 }
 printf("%-5.4g %-5.4g m\n",x1,y1);
 printf("%-5.4g %-5.4g l\n",x2,y2);
 printf("%-5.4g %-5.4g l\n",x3,y3);
 printf("%-5.4g %-5.4g lt\n",x1,y1);
 if (fill) printf(" s\n");
 else      printf(" fill\n");
 xpos = x1;
 ypos = y1;
}

void tafill (int fill)

{
 if (fill) printf("fill\n");
}

double xscalea (double size)
{
  double temp;

  temp = xscale;
/*  xscale = STDSIZ / size; */
  return (temp);
}

double yscalea (double size)
{ 
  double temp;

  temp = yscale;
/*  yscale = STDSIZ / size; */
  return (temp);
}

int taxsize(void)
{ return STDSIZ; }

int taysize(void)
{ return STDSIZ; }

void putstra(char *str)
{
   static int oldtsize=0,oldtfont=0;

  if (tangle) printf ("%-5.3g rotate\n",(tangle*360./32768.));
  if (tsize!=oldtsize || oldtfont!=tfont) {
    if (tsize==0) tsize=10;
    if (tfont == GRKFONT) printf ("/Symbol findfont\n");
    else                  printf ("/Times-Roman findfont\n");
    printf ("%d scalefont setfont\n",(int)(tsize*13500./STDSIZ));
  }
  printf ("(%s) show\n",str);
  if (tangle) printf ("-%-5.3g rotate\n",(tangle*360./32768.));
  oldtsize=tsize;
  oldtfont=tfont;
}

int gmodea(int argc, char **argv)
{
/* Variables declared here.... */
int	x,y;
long	clock;
char	date[30];
extern int rotflag;

    time(&clock);
    strcpy(date, ctime(&clock));

    printf("%%!PS-Adobe1.1\n");
    printf("%%%%Creator: Vid - Display routine by R.G. Smith and");
    printf(" A. Hsu.\n");
    printf("%%%%Title: %s\n",argv[0]);
    printf("%%%%CreationDate: %s",&date[0]);
    printf("%%%%Pages: 1\n");
    printf("%%%%BoundingBox: 50 108 600 680\n");
    printf("%%%%EndComments\n");
    printf("%%%%BeginProlog\n");
    printf("/m {moveto} def\n");
    printf("/l {lineto} def\n");
    printf("/s {stroke} def\n");
    printf("%%%%EndProlog\n");
    printf("%%%%Page: 0 1\n");
    printf("\nnewpath\n");
    printf("save\n");
    printf("%% This is for text at the top of page..\n");
    printf("/CommentFont\n  /Helvetica findfont 10 scalefont def\n");
    printf("/Date (%s) def\n", &date[0]);
/*    if (comment)  printf("/Comment (%s) def\n",comment);
    else          printf("/Comment ([         ]) def\n");
*/
    printf("%% This is where you can scale, rotate and translate the image\n");
    printf("%%     if you feel the need to manipulate the image further...\n");
    if (rotflag) {
      printf("%6.3f rotate\n",rotflag*90.0);
      printf("30 -580  translate\n");
    } 
    else {
      printf("0 rotate\n");
      printf("50 108 translate\n");
    }
    printf("1 1 scale\n");
    /* Let's try drawing an 8 x 8 box for a frame...*/
    printf("%% Border Box Begins here\n");
    printf("%% 0 setgray\n");
    printf("%% 0 0 moveto\n");
    printf("%% 650 0 rlineto\n");
    printf("%% 0 800 rlineto\n");
    printf("%% -650 0 rlineto\n");
    printf("%% closepath\n");
    printf("%% fill\n");
    printf("%% 1 setlinewidth\n");
    printf("%% stroke\n");
    printf("0 650 moveto\n");

    printf("CommentFont setfont\n");
    printf("Date show\n");
/*    printf("0 -72 moveto\n");
    printf("Comment show\n");
*/

if (backgr>=0) {
    setacol(backgr);
    printf("\nnewpath\n");
    printf("-100 -100 moveto\n");
    printf("700 0 rlineto\n");
    printf("0 800 rlineto\n");
    printf("-700 0 rlineto\n");
    printf("closepath\n");
    printf("fill\n\n");
}

    printf("%% Vid data begins here...\n");
    printf(".5 setlinewidth\n");
}
