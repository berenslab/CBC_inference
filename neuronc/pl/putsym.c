/* module putsym.c in program vid */
/* Draws vector font at any orientation on video card */

#include        "chars.h"
#include        "mdef.h"

#define DROPBIT 02000
#define BS      010
#define TAB     011
#define CR      015
#define NL      012
#ifndef ESC
#       define ESC      033
#endif
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

#define STDSIZ  16384	/* standard coordinate range */
#define STDSCL  1	/* STDSIZ/350 to scale standard plot */
#define TWMUL 	50	/* 64/* # of pixels in standard text width size spec. */
#define THMUL 	96	/* 3/2*TWMUL */
#define STDCW 	5	/* standard text size (about 15 x 23) */

#define XMAX    22599		/* maximum x value */
#define XMIN 	0		/* minimum x value */
#define YMAX 	16383	/* maximum y value (man says 779,but 780 fits) */
#define YMIN 	0		/* minimum y value */

extern int grkflag;		/* flag to return ascii chars */
extern int fatbase,fat;		/* line thickness parameters */
extern int xmax;		/* plot window parameters defaulted */
extern int xmin;		/* to maximum size 		*/
extern int ymax;
extern int ymin;
extern int xnew,ynew;		/* new pen location */
extern int xold,yold;		/* old pen location */
extern int xorigin,yorigin;	/* origin for text mode */
extern unsigned tangle, tsize;	/* new text angle, size */
extern char tfont;		/* future text font */
extern double scale;

int *symadd;
char *symbase;
double c, s;
double cos(), sin(), fabs(), floor(), sqrt();
void mov(int hadd,int vadd,int orient);
void dump(int x1,int y1,int x2,int y2);

/*------------------------------------------*/

void setorg(x,y)
   int x,y;

/* Set the origin coords and char coords
   for a line of text. */

{
  xorigin = xold = x;
  yorigin = yold = y;
}

/*------------------------------------------*/

void putsym(str)
   char *str;
/*
 * interpret characters into vectors
 */
   {
        int tfat,xyw;
        int size,move,vmove,add,drop,xp,yp,i;
        char sym,*symptr;
        register int orient,a,b;
        double dist;
        double arg;

        symadd= ascii.saddr;
        symbase= ascii.svec;
	if (grkflag)
	   {
		symadd = greek.saddr;
		symbase = greek.svec;
	   }
      /*  size = (tsize)? tsize: STDCW; */
        size = (tsize)?tsize/4: STDCW;	 /* changed because CWMIN in textf.c */
  	orient= tangle;		/*  was divided by factor of 4 */
        tfat= fat;
        if( orient%QROT )       /* compute transforms for non-90 angles */
           {
                arg = orient*3.14159265/(QROT*2.);
                c= cos( arg );
                s= sin( arg );
           }
        move = size*TWMUL/scale;		/* 0.8 narrows font */
        vmove = size*THMUL/scale;
/*      while((sym=getc(pltin))) */
	while ((sym = *str++))
           {
                if(sym < 040)
                  {
                  switch (sym)  /* standard carriage controls */
                    {
                    case -1:
                        return;
                    case BS:
                        mov(-move,0,orient);
                        break;
                    case TAB:
                        switch (orient)
                          {
                          case 0:
                          case 2*QROT:
                                do mov(move, 0, orient);
                                  while ((abs(xold - xorigin)/move)%8);
                                break;
                          case QROT:
                          case 3*QROT:
                                do mov(move, 0, orient);
                                  while ((abs(yold - yorigin)/move)%8);
                                break;
                          default:
#define sq(x) ((double)((x)*(x)))
                                do
                                  {     mov(move, 0, orient);
                                        dist = sqrt(sq(xold - xorigin) +
                                                    sq(yold - yorigin));
                                  } while ((int)floor(dist/move + .5) % 8);
                                break;
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
                                yorigin = floor(yorigin - c*vmove + 0.5);
                                xorigin = floor(xorigin + s*vmove + 0.5);
                                break;
                          }
                    case CR:
                        xold=   xorigin;
                        yold=   yorigin;
                        break;
                    case ESC:   /* special controls */
          /*            switch(sym=getc(pltin))   */
                        switch(sym = *str++)
                          {
                          case -1:
                          case 0:
                                return;
                          case DS:
                                mov(0,-vmove,orient);
                                break;
                          case DSHALF:
                                mov(0,-(vmove/2),orient);
                                break;
                          case US:
                                mov(0,vmove,orient);
                                break;
                          case USHALF:
                                mov(0,vmove/2,orient);
                                break;
                          case BSPACE:
                                mov(-move,0,orient);
                                break;
                          case GREEK:
                                symadd= greek.saddr;
                                symbase= greek.svec;
                                grkflag=1;
                                break;
                          case SIZEP:
                                size++;
                                move    = size*TWMUL/STDSCL;
                                vmove   = size*THMUL/STDSCL;
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
                        goto    moveah;
                    }
                  continue;
                  }
                if(sym == 040) goto moveah;
                add= symadd[sym-040];
                symptr= symbase +(add&01777);
                drop= (add&DROPBIT ? 2 : 0);
                xnew=xold;
                ynew=yold;
                do
                   {
                        xyw= *symptr++;
                        a= ((xyw&0160)>>4)*move/6;
                        b= ((xyw&07) - drop)*vmove/9;
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
                                xp= floor(xold+ c*a - s*b +0.5);
                                yp= floor(yold+ s*a + c*b +0.5);
                                break;
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
        moveah:
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
                        xold = floor(c*move + xold + 0.5);
                        yold = floor(s*move + yold + 0.5);
                        break;
                  }
  /*         if(grkflag)
              {
                    symadd= ascii.saddr;
                    symbase= ascii.svec;
                    grkflag=0;
              }  */
           }
   }
void mov(int hadd,int vadd,int orient)
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
                xold = floor(c*hadd - s*vadd + xold + 0.5);
                yold = floor(s*hadd + c*vadd + yold + 0.5);
                break;
          }
   }

/*************************/

void getln(void)
{}

