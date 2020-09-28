/* Segment PUTSYMX in program montage */

/* Drivers for software-generated text
   to be used with any graphics card
*/

#include	<stdio.h>
#include	"../h/grdef.h"

#include        "chars.h"
#define DROPBIT 02000
#define BS      010
#define TAB     011
#define CR      015
#define NL      012
#define ESC     033
#define US      'U'
#define DS      'D'
#define USHALF  'u'
#define DSHALF  'd'
#define GREEK   'g'
#define SIZEP   'S'
#define SIZEM   's'
#define FATP    'F'
#define FATM    'f'
#define BSPACE  'b'
#define QROT    8192    /* text angle for a quarter rotation (90 degrees) */

#define sq(x) ((double)((x)*(x)))

#define LINHT   1.5	/* dist to move for next line */
#define STDSCL  1
#define TWMUL 	64	/* # of pixels in standard text width size spec. */
#define THMUL 	120	/* 96	/* 3/2*TWMUL */
#define STDCW 	5	/* standard text size (about 15 x 23) */
#define XMAX 	(STDSIZ * 3/2)	/* maximum x value */
#define XMIN 	0		/* minimum x value */
#define YMAX 	STDSIZ	/* maximum y value (man says 779,but 780 fits) */
#define YMIN 	0		/* minimum y value */
#define GRKFONT 'g'	/* greek font for 's' command */

int grkflag=0;		/* flag to return ascii chars */
int fatbase,fat;	/* line thickness parameters */
int xmax=XMAX;		/* plot window parameters defaulted */
int xmin=XMIN;		/* to maximum size 		*/
int ymax=YMAX;
int ymin=YMIN;

int xnew,ynew;		/* new pen location */
int xold,yold;		/* coords of char */
int xorigin,yorigin;	/* coords for origin of a line of text */
unsigned tangle, tsize;	/* new text angle, size */
static double scale=1;
static int xcent,ycent;	/* center of picture */

/*------------------------------------------*/

setorg(x,y)
   int x,y;

/* Set the origin coords and char coords
   for a line of text. */

{
  xorigin = xold = x;
  yorigin = yold = y;
}

/*------------------------------------------*/

dump(x1,y1,x2,y2)
int x1,y1,x2,y2;
   {
	static int xlst, ylst, temp;
	extern int xcent,ycent;

	x1 += xcent;
	y1 += ycent;
	x2 += xcent;
	y2 += ycent;
	if( (x1 != xlst) || (y1 != ylst) )
		plot(x1,y1,1);
	plot(x2,y2,0);
	xlst= x2;
	ylst= y2;
   }

/*------------------------------------------*/

plot(x,y,icode)
int x,y,icode;

{
   if (icode) tmove(x,y);
   else	      tdraw(x,y);
}

/*------------------------------------------*/

putsymc(ch)

{
  static char str[2] = 0;

  str[0] = ch;
  putsymx(str);
}

/*------------------------------------------*/

int *symadd;
char *symbase;
double c, s;

putsymx(str)
    char *str;

/*
 * interpret characters into vectors
 */
   {
        int tfat;
        int size,move,vmove,add,drop,xp,yp,i;
        char sym,xyw,*symptr;
        register int orient,a,b;
 /*       double dist;  					*/
 /*       double arg, cos(), sin(), fabs(), floor(), sqrt(); 	*/

        symadd= ascii.saddr;
        symbase= ascii.svec;
	if (grkflag)
	   {
		symadd = greek.saddr;
		symbase = greek.svec;
	   }
        size = (tsize)? tsize: STDCW;
        orient = 0;
        tfat = 0;
/*        if( orient%QROT )       /* compute transforms for non-90 angles */
/*           {
                arg = orient*3.14159265/(QROT*2.);
                c= cos( arg );
                s= sin( arg );
           }
*/
        move = size*TWMUL/scale;
        vmove = size*THMUL/scale;
/*	while((sym=getc(pltin))) */
	while ((sym = *str++))
           {
                if(sym < 040)
                  {
                  switch (sym)  /* standard carriage controls */
                    {
                    case -1:
                        return;
                    case BS:
                        movx(-move,0,orient);
                        break;
                    case TAB:
                        switch (orient)
                          {
                          case 0:
                          case 2*QROT:
                                do movx(move, 0, orient);
                                  while ((abs(xold - xorigin)/move)%8);
                                break;
                          case QROT:
                          case 3*QROT:
                                do movx(move, 0, orient);
                                  while ((abs(yold - yorigin)/move)%8);
                                break;
                          default:
                     /*           do
                                  {     movx(move, 0, orient);
                                        dist = sqrt(sq(xold - xorigin) +
                                                    sq(yold - yorigin));
                                  } while ((int)floor(dist/move + .5) % 8);
                      */          break;
                          }
                        break;
                    case NL:
                        switch(orient)
                          {
                          case 0:
                                yorigin -= vmove; break;
                          case QROT:
                                xorigin += vmove; break;
                          case 2*QROT:
                                yorigin += vmove; break;
                          case 3*QROT:
                                xorigin -= vmove; break;
                          default: 
                   /*      yorigin = floor(yorigin - c*vmove + 0.5);
                           xorigin = floor(xorigin + s*vmove + 0.5);
                    */            break;
                          }
                    case CR:
                        xold = xorigin;
                        yold = yorigin;
                        break;
                    case ESC:   /* special controls */
                        switch(sym = *str++)
                          {
                          case -1:
                          case 0:
                                return;
                          case DS:
                                movx(0,-vmove,orient);
                                break;
                          case DSHALF:
                                movx(0,-(vmove/2),orient);
                                break;
                          case US:
                                movx(0,vmove,orient);
                                break;
                          case USHALF:
                                movx(0,vmove/2,orient);
                                break;
                          case BSPACE:
                                movx(-move,0,orient);
                                break;
                          case GREEK:
                                symadd= greek.saddr;
                                symbase= greek.svec;
                                grkflag=1;
                                break;
                          case SIZEP:
                                size++;
                                move  = size*TWMUL/STDSCL;
                                vmove = size*THMUL/STDSCL;
                                break;
                          case SIZEM:
                                size--;
                                if(size<1)      size=1;
                                move    = size*TWMUL/STDSCL;
                                vmove   = size*THMUL/STDSCL;
                                break;
                          case FATP:
                                tfat++;
                                break;
                          case FATM:
                                tfat--;
                                if(tfat<0)      tfat=0;
                                break;
                          }
                        break;
                    default:
                        goto    moveahx;
                    }
                  continue;
                  }
                if(sym == 0x20) goto moveahx;
                add = symadd[sym-0x20];
                symptr = symbase + (add & 0x3ff);
                drop = (add & DROPBIT ? 2 : 0);
                xnew = xold;
                ynew = yold;
                do
                   {
                        xyw = *symptr++;
                        a = ((xyw & 0160)>>4) * move/6;
                        b = ((xyw & 07) - drop) * vmove/9;
                        switch(orient)
                          {
                          case 0:
                                xp=xold+a; yp=yold+b; break;
                          case QROT:
                                xp=xold-b; yp=yold+a; break;
                          case 2*QROT:
                                xp=xold-a; yp=yold-b; break;
                          case 3*QROT:
                                xp=xold+b; yp=yold-a; break;
                          default:
                   /*             xp= floor(xold+ c*a - s*b +0.5);
                                yp= floor(yold+ s*a + c*b +0.5);
                    */            break;
                          }
                        if( !(xyw&0200) )
                                if(tfat)
                                        if(abs(xp-xnew) >= abs(yp-ynew) )
                                           for(i=  -(tfat/2);i<=(tfat+1)/2;i++)
                                              dump(xnew,ynew+i,xp,yp+i);
                                          else
                                            for(i=  -(tfat/2);i<=(tfat+1)/2;i++)
                                              dump(xnew+i,ynew,xp+i,yp);
                                 else  dump(xnew,ynew,xp,yp);
                        xnew=xp;
                        ynew=yp;
                   }     while( !(xyw&010) );
        moveahx:
                switch(orient)
                  {
                  case 0:
                        xold += move; break;
                  case QROT:
                        yold += move; break;
                  case 2*QROT:
                        xold -= move; break;
                  case 3*QROT:
                        yold -= move; break;
                  default:
             /*            xold = floor(c*move + xold + 0.5);
                        yold = floor(s*move + yold + 0.5);
              */          break;
                  }
  /*         if(grkflag)
              {
                    symadd= ascii.saddr;
                    symbase= ascii.svec;
                    grkflag=0;
              }  */
           }
   }

/***************************************************/

movx(hadd,vadd,orient)
int hadd,vadd,orient;
   {
        switch(orient)
          {
          case 0:
                xold += hadd; yold += vadd; break;
          case QROT:
                xold -= vadd; yold += hadd; break;
          case 2*QROT:
                xold -= hadd; yold -= vadd; break;
          case 3*QROT:
                xold += vadd; yold -= hadd; break;
          default:
        /*        xold = floor(c*hadd - s*vadd + xold + 0.5);
                yold = floor(s*hadd + c*vadd + yold + 0.5);
         */       break;
          }
   }

/***************************************************/


