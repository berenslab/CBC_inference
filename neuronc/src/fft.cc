/* Program fft */

/* computes fast fourier transform */

/* prints transform on standard output */

/* Latest Mod	Nov 87		R.G.Smith */

#include <stdio.h>
// #include "mdef.h"

FILE *pictin;
FILE *ncstdout = stdout;
FILE *ncstderr = stderr;

#define PI 3.14159265358979323846264

static int k;
static int w;

#ifdef __cplusplus
extern "C" {
#endif

double atof(char *str);
double cos(double x);
double sin(double x);

#ifdef __cplusplus
}
#endif

double mulreal(double a,double b,double c,double d);
double mulimag(double a,double b,double c,double d);
void run(void);
int readin(double *x,double *y);
void fftstp (double *zinr,double *zini,int after,int now,int before,double *zoutr,double *zouti);

/****************************************/

int main(int argc, char **argv)

{
   char *cptr;
   int i;
   FILE *temp,*freopen();
	 
 pictin = stdin;
 if (argc==1)			/* if user needs help */
 {
 /*   fprintf (ncstdout,"\nfft program Nov 87\r\n"); */
    run();
 }
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
		//signal (2,1);
		break;
     
      	  case 'w': 
		argv++; argc--;
		w = atof(*argv);
		break;

      	  case 'k': 
		argv++; argc--;
		k = atof(*argv);
		break;

       }  /* switch */
     }	   /* if */
     else
     {
       if((pictin=fopen(cptr,"r"))==NULL)
       {
           fprintf(ncstderr,"fft: cannot open %s\n",cptr);
	   fflush (ncstderr);
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

void run(void)

/* On input, z1 and z2 are complex n-vectors,
   n is length of z1 and z2,
   inzee indicates whether z1 or z2 is to be transformed.

   Both z1 and z2 are used as work arrays.
   
   On output, z1 or z2 contains the desired transform
   in the correct order. 
   Inzee indicates whether z1 or z2 contains the transform.
   inzee = 1  => z1 contains transform
   inzee = 2  => z2 contains transform

   Method:  

     The integer n is divided into its prime factors (up to
   a point).  For each such factor p, the p-transform of 
   appropriate p-subvectors of z1 (or z2) is calculated
   in fftstp and stored in a suitable way in z2 (or z1).

   Adapted from fortran routine in:

	 Conte and de Boor,
	 'Elementary Numerical Analysis' 1980
		Chapter 6, "Fast Fourier Transforms", p 283
*/

#define FFTSIZ 2048
 
{
    static int i,j,n;
    static double z1real[FFTSIZ] = {0};
    static double z1imag[FFTSIZ] = {0};
    static double z2real[FFTSIZ] = {0};
    static double z2imag[FFTSIZ] = {0};
    static int prime[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    static int inzee,after,before,next,nextmx,now;

 inzee = 1;
 nextmx = 12;
 n = readin(z1real,z1imag); 		/* read first list */
 after = 1;
 before = n;
 next = 1;
 while (before != 1) {
   if ((before / prime[next-1]) * prime[next-1] < before) {
    next++;
     if (next <= nextmx) continue;
     else {
	now = before;
	before = 1;
     }
   }
   else {
     now = prime[next-1];
     before = before / prime[next-1];
   }

   if (inzee == 1)
      fftstp(z1real, z1imag, after, now, before, z2real, z2imag);
   else
      fftstp(z2real, z2imag, after, now, before, z1real, z1imag);
   inzee = 3 - inzee;
   if (before != 1) after *= now;

 }    /* while (before ) */

if (inzee==1) {
   for (i=0; i<n; i++) 
     printf("%10g %10g\n",z1real[i],z1imag[i]);
}
else {
   for (i=0; i<n; i++)
     printf("%10g %10g\n",z2real[i],z2imag[i]);
}
}

/****************************************/

void fftstp (double *zinr,double *zini,int after,int now,int before,double *zoutr,double *zouti)
{
    int ia,ib,in,j;
    static double angle;
    static double argreal,argimag;
    static double omegreal,omegimag;
    static double valreal,valimag;
    static double tvalreal,targreal;
    int pnt;

 angle = 2 * PI / (now * after);
 omegreal = cos(angle);
 omegimag = -sin(angle); 
 argreal  = 1.0;
 argimag  = 0.0;
 for (j=0; j<now; j++) {
   for (ia=0; ia<after; ia++) {
     for (ib=0; ib<before; ib++) {
       pnt = ((now-1)*before+ib)*after+ia;
       valreal = *(zinr+pnt);
       valimag = *(zini+pnt);
       for (in=now-2; in>=0; in--)  {
	 pnt = (in*before+ib)*after+ia;
	 tvalreal = mulreal(valreal,valimag,argreal,argimag);
	 tvalreal = tvalreal + *(zinr+pnt);
	 valimag = mulimag(valreal,valimag,argreal,argimag);
	 valimag = valimag + *(zini+pnt);
	 valreal = tvalreal;
      }
      pnt = (ib*now+j)*after+ia;
      *(zoutr+pnt) = valreal;
      *(zouti+pnt) = valimag;
    }
    targreal = mulreal(argreal,argimag,omegreal,omegimag);
    argimag = mulimag(argreal,argimag,omegreal,omegimag);
    argreal = targreal;
  }
 }
}

/****************************************/

double mulreal(double a,double b,double c,double d)

{
 return (a*c - b*d);
}

/****************************************/

double mulimag(double a,double b,double c,double d)

{
 return (b*c + a*d);
}

/****************************************/

int readin(double *x,double *y)

/* skip blank lines and read list of
   x,y values until either a blank line or the
   end of text */

#define STRSIZ 80
 
{
    static int  i,xx;
    static int np;
    static char str[STRSIZ]={0};
    static int found;
 
 np = 0;
 for(found=0; !found; )			/* look for leading blank lines */
 {
   if (not fgets (str,STRSIZ,pictin)) return (0);
   if (sscanf(str,"%f%*[ ,]%f",x,y) > 0) break;
 }

 x++; y++; np++;
 for(found=0; !found; )
 {
   if (not fgets (str,STRSIZ,pictin)) break;
   if ((xx=sscanf(str,"%f%*[ ,]%f",x,y)) > 0)
   {
     x++,y++,np++;
     if (np >= FFTSIZ) break;
   }
   else found = 1;
 }

  return (np);
}

/*****************************************/
   
