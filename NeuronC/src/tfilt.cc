/* 	Program TFILT for filter of list of numbers */

/*   	R.G.Smith Nov, 1999	*/


#ifdef __cplusplus
extern "C" {
#endif

#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif
#include <stdio.h>
double atof(const char*str);

#ifdef __cplusplus
}
#endif

FILE *pictin;

int fflg=1;

static char *filnam = 0;
static double tau=1.0;
static double timestep=1.0;

void run(void);

/****************************************/

main(int argc, char **argv)

{
   int i;
   FILE *temp,*freopen();
   char *cptr=0;
	 
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

          case 'f':
		fflg = !fflg;	
		break;
     
          case 's':
		argc--; argv++;		/* time step */
		timestep = atof(*argv);
		break;
     
          case 't':
		argc--; argv++;		/* filter tau */
		tau = atof(*argv);
		break;
     
	  default:
		ncfprintf (stderr,"tfilt: unknown switch '%s'\n",*argv);
		exit(2);

        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           ncfprintf(stderr,"var: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         }
       else {
         filnam = cptr;
         run();
       }
       if (argc <= 1) break;
     }
    }
    else run();
  }
 while (argc > 0);
}

/****************************************/

void run(void)

/* skip blank lines and read list of
   x,y values until either a blank line or the
   end of text */

#define STRSIZ 80
 
{
    static int  i;
    static int n,nfields;
    static char str[STRSIZ]={0};
    static int found;
    static double x,y,ndata;
    static double olddata=0;
 
    double k = 1.0;
    double data = 0.0;		/* data entry */

 k = 1.0 - exp (-timestep / tau);	/* fix error for small tau */ 
 if (k > 1.0) k = 1.0;			/* limit freq to reasonable values */
 if (k <= 0) k = 1e-6;

 n = 0;
 for(found=0; !found; )			/* ignore leading blank lines */
  {					/* ignore lines starting with '#' */
   if (! fgets (str,STRSIZ,pictin)) return;
   if (*str == '#') continue;
   if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {found=1; break;}
  }

 for(found=0; !found; ) {
   if (*str == '#') continue;
   if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {
      switch (nfields) {
        case 1: data = x; break;
        case 2: data = y; break;
      }
      ndata = (data - olddata) * k + olddata;
      olddata = ndata;
      switch (nfields) {
        case 1: printf ("%g\n",ndata);
		break;
        case 2: printf ("%g %g\n",x,ndata);
		break;
      }
   }
   else found = 1;
   if (! fgets (str,STRSIZ,pictin)) break;
  }

}

