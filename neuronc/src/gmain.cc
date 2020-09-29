/* Module gmain in program gausnn */

/* Makes a gaussian nearest-neighbor distance */

/* R.G.Smith */

#include <stdio.h>
#include "ncio.h"

double sdensity = 0;
double smean = 0;
double sstdev = 0;
double sms = 0;
double sframex = 0;
double sframey = 0;
double scenterx = 1e38;
double scentery = 1e38;
int snumcells = 0;
int rseed = 1234;
int printfl = 0;
int first_center = 0;
double *x=0,*y=0;

int cumrand=0;

#ifdef __cplusplus
extern "C" {
#endif

double atof(const char *);

#ifdef __cplusplus
}
#endif

void setrand(int n);
int gausnn(double mean, double stdev, double density, double ms, int grseed,
	double framex, double framey, double xcent, double ycent, 
	int numcells, double **xarr, double **yarr, int first_center, 
	int filout, int textfl, int printfl);
                                                            
extern FILE *outfil;

/* -------------------------------------------------------------- */

int main(int argc, char **argv)
{
   char *cptr;
   FILE *freopen(const char *, const char *, FILE *);
         
 outfil = (FILE *)NULL;
 if (argc==1) {                   /* if user needs help */
   ncfprintf (stderr,"    Usage: gausnn -d density -t regularity file\n");
   ncfprintf (stderr,"    Usage: gausnn -m mean -s stdev file\n");
   ncfprintf (stderr,"     Other options: -n ncells -r seed -p (debug)\n");
   ncfprintf (stderr,"                    -x xsize  -y ysize\n");
   return (0);
 }
 else
 do                                     /* if there are any switches */
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
     
          case 'd': 
                argv++; argc--;
                sdensity = atof(*argv);
                break;

          case 'f': 
                argv++; argc--;
                first_center = atof(*argv);
                break;

          case 'm': 
                argv++; argc--;
                smean = atof(*argv);
                break;

          case 'n': 
                argv++; argc--;
                snumcells = (int)atof(*argv);
                break;

          case 'p': 
                printfl = !printfl;
                break;

          case 'r': 
                argv++; argc--;
                rseed = (int)atof(*argv);
                break;

          case 's': 
                argv++; argc--;
                sstdev = atof(*argv);
                break;

          case 't': 
                argv++; argc--;
                sms = atof(*argv);
                break;

          case 'x': 
                argv++; argc--;
                sframex = atof(*argv);
                break;

          case 'X': 
                argv++; argc--;
                scenterx = atof(*argv);
                break;

          case 'y': 
                argv++; argc--;
                sframey = atof(*argv);
                break;

          case 'Y': 
                argv++; argc--;
                scentery = atof(*argv);
                break;

          default:
                ncfprintf (stderr,"gausnn: unknown switch '%s'\n",*argv);
                return(1);

        }  /* switch */
      }    /* if */
     else
      {
       if((outfil=fopen(cptr,"w"))==NULL)
         {
           ncfprintf(stderr,"gausnn: cannot open %s\n",cptr);
           fflush (stderr);
           continue;
         }
       if (argc <= 1) break;
     }
    }
  }
 while (argc > 0);

 if (rseed) setrand(rseed);	/* initialize random number generator */

 gausnn(smean,sstdev,sdensity,sms,0,sframex,sframey,scenterx,scentery,
			snumcells,&x,&y,first_center,1,1,printfl);
}

/* ---------------------------------------------------------- */

