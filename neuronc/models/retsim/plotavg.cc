/* plotavg  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ncfuncs.h"
#include "ncinit.h"

/* Split .r files with repeating zero time into separate traces */

int ntraces;
int sttrace;
int tracelen;
int plotn;
int plotd;
int t, t1, t2, t3, t4;
int ntimes;
int navg;

extern int yprecision;

int ninfo;

double *data_in=NULL;
double *data_out=NULL;


int main (int argc, char **argv) 
{
    int k,n,x,y,nrows, ncols;
    int xstart, xend, len;
    double zeronum,val;
    const char *progvar;
    double *data_avg;
    int *times;

   progvar = argv[0];

   setptr("ntraces",    &ntraces);
   setptr("sttrace",    &sttrace);
   setptr("tracelen",   &tracelen);
   setptr("t1",         &t1);
   setptr("t2",         &t2);
   setptr("t3",         &t3);
   setptr("t4",         &t4);
   setptr("navg",       &navg);
   setptr("ninfo",      &ninfo);
   setptr("yprecision", &yprecision);

   setvars(argc, argv);

   info = 0;
   if (notinit(ninfo)) ninfo = 0;
   if (notinit(yprecision)) yprecision = 18; /* precision of y values in printout */

   if (notinit(ntraces))   ntraces = 0;		// number of traces (voltage clamp steps), 0 -> auto)
   if (notinit(sttrace))   sttrace = 0;		// starting trace (voltage clamp steps), 0=first)
   if (notinit(tracelen)) tracelen = 0;		// number of rows (length) in each trace, 0 -> auto
   if (notinit(plotn))       plotn = 5;		// the plot (column) number to split (0=time column)
   if (notinit(plotd))       plotd = -1;	// the plot (column) number to subtract before split
   if (notinit(navg))         navg = 10;       	// how many points to average for each output point

   data_in = fread ("stdin",&nrows,&ncols);

   if (ninfo>=1) fprintf(stderr,"ncols %d nrows %d\n",ncols, nrows);

   if (ntraces==0) ntraces = ncols;
   if (tracelen==0)  tracelen = nrows;


   ntimes = 3;
   // if (notinit(t4))  ntimes--;	
   if (notinit(t3))  ntimes--;	
   if (notinit(t2))  ntimes--;	
   if (notinit(t1))  ntimes--;	
   if (ntimes==0) fprintf(stderr,"plotavg: no times selected.\n");
    
   if (ntimes >= 4) 
      if (t4 < 0 || t4 >= tracelen) {
	   fprintf (stderr, "plotavg: t4 is out of range.\n");
      }
   if (ntimes >= 3) 
      if (t3 < 0 || t3 >= tracelen) {
	   fprintf (stderr, "plotavg: t3 is out of range.\n");
      }
   if (ntimes >= 2)
      if (t2 < 0 || t2 >= tracelen) {
	   fprintf (stderr, "plotavg: t2 is out of range.\n");
      }
   if (ntimes >= 1)
      if (t1 < 0 || t1 >= tracelen) {
	   fprintf (stderr, "plotavg: t1 is out of range.\n");
      }


   data_avg = (double *)emalloc(ntimes * ntraces * sizeof(double));

   // get first bunch of points

   times = (int *)emalloc(ntimes*sizeof(int));
   times[0] = t1;
   times[1] = t2;
   times[2] = t3;
   times[3] = t4;

   if (ninfo>=1) {
      fprintf(stderr,"times ");
      for (t=0; t<ntimes; t++) {
	fprintf(stderr,"%d ",times[t]);
      }
      fprintf(stderr,"\n");
   }

   for (t=0; t<ntimes; t++) {
        for (y=0; y<ncols; y++) {
           data_avg[t*ncols+y] = 0;
        }
    }

   if (ninfo>=1) fprintf(stderr,"ntraces %d tracelen %d navg %d\n", ntraces, tracelen, navg);

   for (t=0; t<ntimes; t++) {
      xend = times[t] + navg;
      for (x=times[t]; x<xend; x++) {
        for (y=0; y<ncols; y++) {
          data_avg[t*ncols+y] += data_in[x*ncols+y];
        };
      };
      for (y=0; y<ncols; y++) {
        printf ("%-11.*2$g\t",data_avg[t*ncols+y] / navg,yprecision);
      }
      printf ("\n");
   }

}

