/* 	Program VAR for variance of list of numbers */

/*   	Latest mod	28-Apr-88	R.G.Smith

*/

#include <stdio.h>
#include <stdlib.h>
FILE *pictin;

int sflg=0;
int fflg=1;
static char *filnam = 0;

double sqrt();
double atof();
void run();

/****************************************/

int main(int argc, char **argv)

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

          case 'i':
		// signal (2,1);
		break;
     
          case 'f':
		fflg = !fflg;	
		break;
     
          case 's':
		sflg = !sflg;	
		break;
     
	  default:
		fprintf (stderr,"var: unknown switch '%s'\n",*argv);
		exit(2);

        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"var: cannot open %s\n",cptr);
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

void run()

/* skip blank lines and read list of
   x,y values until either a blank line or the
   end of text */

#define STRSIZ 80
 
{
    static int  i;
    static int n,nfields;
    static char str[STRSIZ]={0};
    static int found;
    static double x,y;
 
    double k = 1.0;
    double data = 0.0;		/* data entry */
    double adj=0.0, adjsq=0.0;
    double sum=0.0, sumsq=0.0;	/* sum of adjusted data */
    double var=0.0, stdev=0.0,mean=0.0;
    double var2=0.0, stdev2=0.0;
    double stdevx=0.0, meanstd=0.0;
 
 n = 0;
 for(found=0; !found; )			/* ignore leading blank lines */
  {					/* ignore lines starting with '#' */
   if (! fgets (str,STRSIZ,pictin)) return;
   if (*str == '#') continue;
   if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {found=1; break;}
  }

 if (found) {
      switch (nfields) {
        case 1: data = x; break;
        case 2: data = y; break;
      }
      adj = data * k;
      adjsq = adj * adj;
      sum += adj;
      sumsq += adjsq;
      n++;
 }
 for(found=0; !found; ) {
   if (! fgets (str,STRSIZ,pictin)) break;
   if (*str == '#') continue;
   if ((nfields=sscanf(str,"%lf%*[ ,\t]%lf",&x,&y)) > 0) {
      switch (nfields) {
        case 1: data = x; break;
        case 2: data = y; break;
      }
      adj = data * k;
      adjsq = adj * adj;
      sum += adj;
      sumsq += adjsq;
      if (++n >= 1000000) break;
   }
   else found = 1;
  }

  if (n == 0) n = 1;
  mean = sum / n;
  var = (sumsq - (sum * mean)) / n;
  if (n==1) n = 2;
  var2 = (sumsq - (sum * mean)) / (n-1);
  stdev = sqrt (var);
  stdev2 = sqrt (var2);
  stdevx = stdev;
  if (stdevx == 0) stdevx = 1;
  meanstd = mean / stdevx;

  if (sflg) {
    if (fflg && filnam) printf ("%s: ",filnam);
    printf ("%9.3g %9.3g %9.3g\n",mean,stdev,meanstd);
  }
  else {
    if (fflg && filnam) {
       printf ("%s: n %d sum %g mean %g var %g var(n-1) %g \n",
		filnam,n,sum,mean,var,var2);
       printf ("    stdev %g stdev(n-1) %g m/s %g\n",
		stdev,stdev2,meanstd);
    }
    else {
      printf ("n %d sum %g mean %g var %g var(n-1) %g ",
		n,sum,mean,var,var2);
      printf ("stdev %g stdev(n-1) %g m/s %g\n",
		stdev,stdev2,meanstd);
    }
  }
}

