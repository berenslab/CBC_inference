/* plotsplitf  */

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
int timetrace;
int stimtrace;

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

   progvar = argv[0];

   setptr("ntraces",    &ntraces);
   setptr("sttrace",    &sttrace);
   setptr("tracelen",   &tracelen);
   setptr("plotn",      &plotn);
   setptr("plotd",      &plotd);
   setptr("ninfo",      &ninfo);
   setptr("yprecision", &yprecision);
   setptr("timetrace",  &timetrace);
   setptr("stimtrace",  &stimtrace);

   setvars(argc, argv);

   info = 0;
   if (notinit(ninfo)) ninfo = 0;
   if (notinit(yprecision)) yprecision = 18; /* precision of y values in printout */

   if (notinit(ntraces))   ntraces = 0;		// number of traces (voltage clamp steps), 0 -> auto)
   if (notinit(sttrace))   sttrace = 0;		// starting trace (voltage clamp steps), 0=first)
   if (notinit(tracelen)) tracelen = 0;		// number of rows (length) in each trace, 0 -> auto
   if (notinit(plotn))       plotn = 5;		// the plot (column) number to split (0=time column)
   if (notinit(plotd))       plotd = -1;	// the plot (column) number to subtract before split
   if (notinit(timetrace))   timetrace = -1;	// the time column number, -1 = none
   if (notinit(stimtrace))   stimtrace = -1;	// the stimulus plot number for 1st column, -1 = none

   data_in = fread ("stdin",&nrows,&ncols);

   if (ninfo>=1) fprintf(stderr,"ncols %d nrows %d\n",ncols, nrows);

   zeronum = 1e-8;

   // find start of data

   for (x=0; x<nrows; x++) {
      val = data_in[x*ncols+0];
      if (fabs(val) <= zeronum) break;
   }

   xstart = x;

   // find length of first trace

   for (x++; x<nrows; x++) {
      val = data_in[x*ncols+0];
      if (fabs(val) <= zeronum) break;
   };
   xend = x;

   // count traces

   for (x=n=0; x<nrows; x++) {
      if (fabs(data_in[x*ncols+0]) <= zeronum) n++;
   };

   if (ntraces==0)  ntraces = n;
   if (tracelen==0) tracelen = xend - xstart - 1;

   if (ninfo>=1) fprintf(stderr,"ntraces %d tracelen %d xstart %d xend %d\n",ntraces, tracelen, xstart, xend);

   if ((data_out=(double *)emalloc((tracelen+1)*ntraces*sizeof(double)))==NULL) {
	   fprintf (stderr, "plotsplit: can't allocate output array.\n");
   }

   // remove predur

   len = nrows - xstart;
   for (x=0; x<len; x++) {
     for (y=0; y<ncols; y++) {
       data_in[x*ncols+y] = data_in[(x+xstart)*ncols+y];
     };
   };

   // split the traces

   for (x=k=y=0; x<len; x++) {
        if (plotd>0) 
         data_out[k*ntraces+y] = data_in[x*ncols+plotn] - data_in[x*ncols+plotd];
        else  
         data_out[k*ntraces+y] = data_in[x*ncols+plotn];
        if (k++>=tracelen) {
          k = 0;
          y++;
        };
        if (y>=ntraces) break;
   };

   for (x=0; x<=tracelen; x++) {
     if (timetrace >= 0) {
          val = data_in[x*ncols+timetrace];
	  printf ("%-11.*2$g",val,yprecision);
     }
     if (stimtrace >= 0) {
	  if (timetrace>=0) printf("\t");
          val = data_in[x*ncols+stimtrace];
	  printf ("%-11.*2$g",val,yprecision);
     }
     for (y=sttrace; y<sttrace+1; y++) {
        if (timetrace>=0 || stimtrace>=0) printf("\t");
       val = data_out[x*ntraces+y];
       printf ("%-11.*2$g",val,yprecision);
     };
     for ( ; y<ntraces; y++) {
       val = data_out[x*ntraces+y];
       printf ("\t%-11.*2$g",val,yprecision);
     };
     printf ("\n");
   };
}

