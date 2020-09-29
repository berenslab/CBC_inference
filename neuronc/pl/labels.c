/* Program LABELS */

/* Labels pictures on standard output */

/* Latest Mod	Mar 84		R.G.Smith */

#include <stdio.h>
#include "mdef.h"


FILE *pictin;
FILE *stdplt;
char inbuf[BUFSIZE];

void label(FILE *);
void signal(int,int);

/****************************************/

void run(void)

/* get user input from stdin */


{
 label(pictin);
}


/****************************************/
void main(argc,argv)
   int argc;
   char **argv;

{
   FAST TEXT *cptr;
   int i;
   FILE *temp,*freopen();
	 
 stdplt = stdout;
 pictin = stdin;
 if (argc==1)			/* if user needs help */
    run();
 else
 do					/* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {

          case 'i':
		signal (2,1);
		break;

        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"Labels: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         }
       run();
       if (argc <= 1) break;
     }
    }
   else run();
  }
 while (argc > 0);
}

/****************************************/
