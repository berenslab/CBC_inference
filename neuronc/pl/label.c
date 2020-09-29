/* Segment LABEL in program LABELS */

/* Routine to write labels on pictures. */

/* Latest mod   11-Mar-84		R.G.Smith	*/


/* 	Input statement formats:

	m x y			move to x,y
	d x y			draw to x,y
	l x1 y1 x2 y2		draw line segment from x1,y1 to x2,y2
	rm x y			rel move to x,y
	rd x y			rel draw to x,y
	rl x1 y1 x2 y2		rel draw line segment from x1,y1 to x2,y2
	ro x			rotate whole frame x degrees.
	rt x			rotate text x degrees.
	t x y "text"		print text at loc x,y.
	u "text"		print text at current loc.
	c x 			change char size to x (fract of 1)
	ci x			make cirle with radius x
	p x			change pen to color x (0 = put back)
	s x			make picture size = x
	e			erase screen
	o x			make dashed line number x
	or x y			sets origin of frame to x,y.
	w x1 x2 y1 y2		window from (xmin,xmax,ymin,ymax)
	fo font			change text to font
	fr name			move to frame name
	fl shape size		draw flag like filter in "mpict"
	i name			include file name
	q,x			exit


The picture is drawn on a grid
of 0 to 1 for both X and Y,
except when scaled.

*/ 	
	

#include <stdio.h>
#include <signal.h>
#include <string.h>
#include "stdplt.h"
#include "mdef.h"

typedef void (*sighandler_t)(int);
int system(const char *);

int istat;
int iflag;

void ctext(void);
void cgraphics(void);
void drtri();
void drsqr();
void drrect();
void drast();
void drchar();
void drcirc();
void drcros();
void erase();
char cfont();
double crotate();

/*-------------------------------*/

void onintr(void)

{
  extern int iflag;

  iflag = 1; 
  signal (2,(sighandler_t)onintr);		/* go to "onintr" at ^C */
}

/*-------------------------------*/

void label(FILE *instream)

{
   static double x,y,x2,y2,size=0;
   static TEXT labl[1000], str[80], comnd[80], comch[10];
   int  ch,dashnum,flsize,fill,aspect;
   char delim,filnam[20];
   FILE *tempin;

 iflag = 0;
 istat = (int) signal (2, (sighandler_t)1);	/* save original status */
 if (istat != 1)
   signal (2,(sighandler_t)onintr);		/* go to "onintr" at ^C */

 size = 1.0;
 frame ("zzyz");
 while ((ch = fscanf (instream,"%s",comnd)) > 0)
   {
    if (iflag) break;
    x = y = 0.0;
    *labl = 0; *str = 0;
    switch (*comnd & 0X7F)
     {

    case '!':	ctext();		/* change back to text screen */
	 	fflush(stdout);
	  	fgets(labl,1000,stdin);		/* run system command */
		sprintf (str,"csh -c '%s'>/dev/tty",labl);
		system(str);		/* do system command, stop */
	   	fprintf(stderr,"\nType <RET> to continue...");
	   	fscanf(instream,"\n");
		cgraphics();		/* change back to graphics screen */
	 	fflush(stdout);
	   	break;

    case 'T':					/* Text at pos. x,y */
    case 't': fscanf (instream," %lf%*[ ,]%lf", &x,&y);
    	      fscanf (instream,"%*[ ,]%c",&delim);
	      sprintf (str,"%%[^%c]%c",delim,delim);
	      fscanf (instream,str,labl);
	      move (x,y);
	      textf (labl);
	      purge();
	     break;

    case 'U':					/* Text at current loc. */
    case 'u': fscanf (instream," %c",&delim);
	      sprintf (str,"%%[^%c]%c",delim,delim);
	      fscanf (instream,str,labl);
	      textf (labl);
	      purge();
	     break;

    case 'M':					/* move */
    case 'm': fscanf (instream," %lf%*[ ,]%lf", &x,&y);
 	      move (x,y);
	     break;

    case 'D':					/* draw */
    case 'd': fscanf (instream," %lf%*[ ,]%lf", &x,&y);
	      draw (x,y);
	      purge();
	     break;

    case 'F':	
    case 'f': switch (*(comnd+1) & 0x7f)
		{
        case 'L':				/* flag */
        case 'l': fscanf (instream," %s",comnd);	/* shape */
        	  if (strncmp (comnd,"ch",2)==0)
		    fscanf (instream," %c",comch);	/* char to print */
        	  if (strncmp (comnd,"rect",2)==0)
		    fscanf (instream," %d %d",aspect,fill);/* char to print */
        	  if (strncmp (comnd,"sqr",2)==0)
		    fscanf (instream," %d",fill);	/* char to print */
        	  if (strncmp (comnd,"circ",2)==0)
		    fscanf (instream," %d",fill);	/* char to print */
        	  fscanf (instream," %d",&flsize);	/* size */
	          switch (*comnd)
		    {
			case 't': drtri(flsize); break;   /* triangle */
			case 's': drsqr(flsize,fill); break;   /* square */
			case 'r': drrect(flsize,aspect,fill); break; /* rect */
			case 'a': drast(flsize); break;   /* asterisk */
			case 'c':
			 switch (*(comnd+1))
			  {					    /* char */
				case 'h':  drchar(*comch,flsize); break;
				case 'i':  drcirc(flsize,fill); break;   /* circle */
				case 'r':  drcros(flsize); break;   /* cross */
			  }
			 break;
		    }
	     break;

	case 'R':				/* frame */
	case 'r': fscanf (instream," %s", labl);
	          frame (labl);
	          purge();
	         break;

        case 'O':				/* font */
        case 'o': fscanf (instream," %s",comnd);	/* font (A=default) */
	      ch = cfont (*comnd);
	     break;

		}		/* 'f' */
 
    case 'R':
    case 'r': switch (*(comnd+1) & 0x7f)
		  {

 	case 'M':					/* rel move */
	case 'm': fscanf (instream," %lf%*[ ,]%lf", &x,&y);
 	      rmove (x,y);
	      purge();
	     break;

	case 'D':					/* rel draw */
	case 'd': fscanf (instream," %lf%*[ ,]%lf", &x,&y);
	      rdraw (x,y);
	      purge();
	     break;

	case 'L':					/* rel line */
	case 'l': fscanf(instream," %lf%*[ ,]%lf%*[ ,]%lf%*[ ,]%lf", &x,&y,&x2,&y2);
	      rmove (x,y);
	      rdraw (x2,y2);
	      purge();
	     break;

	case 'O':
	case 'o': fscanf (instream," %lf",&x);
		  rotate (x * DEG);
	  	  purge();
		break;

	case 'T':
	case 't': fscanf (instream," %lf",&x);
		  crotate (x * DEG,".");
	  	  purge();
		break;

	     }	/* switch ()  relative */
	    break;

    case 'L':					/* line */
    case 'l': fscanf (instream," %lf%*[ ,]%lf%*[ ,]%lf%*[ ,]%lf", &x,&y,&x2,&y2);
	      move (x,y);
	      draw (x2,y2);
	      purge();
	     break;

    case 'C':
    case 'c': switch (*(comnd+1) & 0x7f)
           {    
     
     case '\n':
     case '\0':
     case 'h': fscanf (instream," %lf", &x);		/* Char size */
	      cwidth (x,".");
	     break;

     case 'i': fscanf (instream," %lf", &x2);
	      _drcirc(x2,0);
	      purge();
	     break;
           }
         break;

    case 'P':					/* Pen number */
    case 'p': fscanf (instream," %lf", &x);
	      cpen ((int)x);
	     break;

    case 'S':
    case 's': fscanf (instream," %lf",&size);		/* Size of picture */
	      scale (size,size);
	     break;

    case 'O':
    case 'o': switch (*(comnd+1) & 0x7f)
		  {

	case 'R':
	case 'r': fscanf (instream," %lf%*[ ,]%lf",&x,&y);
		  origin (x,y);
		  purge ();
		 break;

	default: fscanf (instream," %d",&dashnum);	   /* dashed line */
	        dash (dashnum);
	        break;

	        }		/* O */
	       break;
     case 'W':					/* window */
     case 'w': fscanf (instream," %lf%*[ ,]%lf%*[ ,]%lf%*[ ,]%lf", &x,&x2,&y,&y2);
	     window (x,x2,y,y2);
	     purge();
	    break;
 

	      
    case 'I':
    case 'i': tempin = instream;			/* include file */
	      fscanf (instream," %s",filnam);
	      if ((tempin = fopen (filnam,"r")) == NULL)
		{
		  fprintf (stderr,"Label: cannot open '%s'\n",filnam);
		  fflush  (stderr);
	          break;
	        }
	      label (tempin);
	      fclose (tempin);
	     break; 

    case 'E':
    case 'e':   erase();
		purge();
	       break;
    case 'q':
    case 'Q':
    case 'x':
    case 'X': 
	     purge();
	     return;				/* exit */

    case '#':					/* comment */
     default: break;

     }

   fflush (stdout);

/*    fprintf (stderr,"%c %d %d '%s'\r\n",ch,x,y,labl);
*/
   }

  purge ();
}

/*******************************/

/* not currently used */

/*

drcircl(rad)
    double rad;
{
    static int i,incr;
    static double theta,dtheta,arclen;
    double cos(double), sin(double);

  rmove (rad,0.0);
  incr = 200;
  if (rad > .2) incr = 400;
  if (rad > .3) incr = 1000;
  dtheta = 2. * PI / incr;
  arclen = dtheta * rad;
  for (theta=0.0,i=0; i<=incr; i++)
   {
     theta += dtheta;
     rdraw (-arclen * sin(theta), arclen * cos(theta));
  }
  rmove(-rad,0.0);
}

*/

