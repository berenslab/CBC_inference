/* Bins */

/* Prints a histogram list on standard output */

/* Standard input is a list of numbers;
   output is a list of x,y points which describe the
   distribution of the original list.
*/

#include <stdio.h>


#define BINSIZ 5000
#define BINWID 1.0

#define STRSIZ 80
int bin[BINSIZ];

double binwid,offset;
double data;
static int histofl = 0;
static int nocentfl = 0;


double atof();

FILE *pictin;

/****************************************/

main(argc,argv)
   int argc;
   char **argv;

{
   char *cptr;
   int i;
   FILE *temp,*freopen();
   double atof();
	 
 binwid = BINWID;
 offset = 0.0;
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

      case 'h':   histofl = 1;		/* generate histogram bars */
		 break;
      
      case 'n':  nocentfl = 1; 		/* bins not centered on even binwidths*/
		 break;
      
      case 'o':  argc--; argv++;	/* offset values for array */
		 offset = atof (*argv);
		 break;
      
      case 'w':  argc--; argv++;	/* binwidth */ 
		 binwid = atof (*argv);
		 if (binwid <= 0.0) binwid = BINWID;
		 break;

      default:	fprintf (stderr,"bin -h -n -w width -o offset\n");
 		 return;


        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"Bins: cannot open %s\n",cptr);
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

/*------------------------------------*/

run()

{
   double data;
   int bini;
   int lolimit, hilimit;
   int i,j,n,found;
   static char str[STRSIZ]={0};

 lolimit = BINSIZ;
 hilimit = 0;
 for(n=found=0; !found; ) {
   if (! fgets (str,STRSIZ,pictin)) break;
   if (*str == '#') continue;
   if ((j=sscanf(str,"%lf",&data)) > 0) {	 /* read in data */
      if (! nocentfl) data += binwid / 2;
      data /= binwid;
      bini = data - offset;			/* bini = index into array */
      if (bini > BINSIZ) bini = BINSIZ;
      if (bini < 0) bini = 0; 
      if (bini < lolimit) lolimit = bini;
      if (bini > hilimit) hilimit = bini;
      bin[bini]++;
      n++;
     }
     else found=1;
   }

  for (i=lolimit; i<=hilimit; i++)		/* output list of values */
   {
    if (histofl) {
      if (nocentfl)
       {
	if (i==lolimit)
        printf ("%g  %d\n",(i * binwid          + offset), 0);
        printf ("%g  %d\n",(i * binwid          + offset), bin[i]);
        printf ("%g  %d\n",(i * binwid + binwid + offset), bin[i]);
	if (i==hilimit)
        printf ("%g  %d\n",(i * binwid + binwid + offset + binwid*.01), 0);
       }
      else	/* normal bin centering */
       {
	if (i==lolimit)
        printf ("%g  %d\n",(i * binwid - binwid / 2 + offset), 0);
        printf ("%g  %d\n",(i * binwid - binwid / 2 + offset), bin[i]);
        printf ("%g  %d\n",(i * binwid + binwid / 2 + offset), bin[i]);
	if (i==hilimit)
        printf ("%g  %d\n",(i * binwid + binwid / 2 + offset), 0);
       }
  }
    else
      printf ("%g  %d\n",(i * binwid + offset), bin[i]);
   }

}
