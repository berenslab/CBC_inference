/* modelfit */
/* curve fitting to a model */

/* R.G. Smith   Mar, 2014 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "ncfuncs.h"
#include "ncinit.h"

#include "lmmin.h"
#include "lm_eval.h"

#include "ncio.h"
#include "lm_funcs.h"

double c1, c1_max, c1_min;
double c2, c2_max, c2_min;
double c3, c3_max, c3_min;
double c4, c4_max, c4_min;
double c5, c5_max, c5_min;
double c6, c6_max, c6_min;
double c7, c7_max, c7_min;
double c8, c8_max, c8_min;
double c9, c9_max, c9_min;
double c10, c10_max, c10_min;
double c11, c11_max, c11_min;
double c12, c12_max, c12_min;

const char *p1;
const char *p2;
const char *p3;
const char *p4;
const char *p5;
const char *p6;
const char *p7;
const char *p8;
const char *p9;
const char *p10;
const char *p11;
const char *p12;

int nparms;		// number of parameters to be passed to the simulator
int filesize;
int filesizet;
int filesize2d;		// data file
int filesize2t;		// template file

const char *expt_string;
const char *data_split;
char data_split_string[200];
char data_split_cmd[200];
char exptcmd[5000];
char runcmd[5000];
int tracelen;		// length of each trace in data file, 0 -> auto
int plotn;		// trace number from .r file to fit
double waitsecs;	// number of seconds to wait after spawning model in background 

#define FILENAMESIZ 100
const char *model_output_file;
char model_output_data_file[FILENAMESIZ];
const char *data_file;
const char *template_file;
const char *template_file2;
char *xsystem (const char *str);

/*-----------------------------------------------------------------------------*/

double *model_data_arr=NULL;
int    *addr_arr=NULL;
double *template_arr=NULL;
double *template_arr2=NULL;

double run_model (double *coeff, int n) 

/* Run a retsim model, fit to original data.
/* If n > 0, model is run in background 

  Coefficients are: 
	coeff[0] = C1 (amplitude)
	coeff[1] = C2 (radius)
	coeff[2] = C3 (amplitude)
	coeff[3] = C4 (radius)
	coeff[4] = C5 offset
	coeff[5] = C6 offset
	coeff[6] = C7 offset
	coeff[7] = C8 offset
	coeff[8] = C9 offset
	coeff[9] = C10 offset
*/

#define _C1 0
#define _C2 1
#define _C3 2
#define _C4 3
#define _C5 4
#define _C6 5
#define _C7 6
#define _C8 7
#define _C9 8
#define _C10 9
#define _C11 10
#define _C12 11
#define NCOEFF 12

{
    int i, xval;
    double lc1, lc2, lc3, lc4, lc5, lc6, lc7, lc8, lc9, lc10, lc11, lc12;
    double retval;
    char *runmodel_cmd;

  //  fprintf (stderr,"x %18.18g c1 %18.18g\n",x, coeff[0]);
      lc1 = coeff[_C1];
      //if (!notinit(c1_max) && lc1 > c1_max) lc1 = c1_max;
      //if (!notinit(c1_min) && lc1 < c1_min) lc1 = c1_min;
      lc2 = coeff[_C2];
      //if (!notinit(c2_max) && lc2 > c2_max) lc2 = c2_max;
      //if (!notinit(c2_min) && lc2 < c2_min) lc2 = c2_min;
      lc3 = coeff[_C3];
      //if (!notinit(c3_max) && lc3 > c3_max) lc3 = c3_max;
      //if (!notinit(c3_min) && lc3 < c3_min) lc3 = c3_min;
      lc4 = coeff[_C4];
      //if (!notinit(c4_max) && lc4 > c4_max) lc4 = c4_max;
      //if (!notinit(c4_min) && lc4 < c4_min) lc4 = c4_min;
      lc5 = coeff[_C5];
      //if (!notinit(c5_max) && lc5 > c5_max) lc5 = c5_max;
      //if (!notinit(c5_min) && lc5 < c5_min) lc5 = c5_min;
      lc6 = coeff[_C6];
      //if (!notinit(c6_max) && lc6 > c6_max) lc6 = c6_max;
      //if (!notinit(c6_min) && lc6 < c6_min) lc6 = c6_min;
      lc7 = coeff[_C7];
      //if (!notinit(c7_max) && lc7 > c7_max) lc7 = c7_max;
      //if (!notinit(c7_min) && lc7 < c7_min) lc7 = c7_min;
      lc8 = coeff[_C8];
      //if (!notinit(c8_max) && lc8 > c8_max) lc8 = c8_max;
      //if (!notinit(c8_min) && lc8 < c8_min) lc8 = c8_min;
      lc9 = coeff[_C9];
      //if (!notinit(c9_max) && lc9 > c9_max) lc9 = c9_max;
      //if (!notinit(c9_min) && lc9 < c9_min) lc9 = c9_min;
      lc10 = coeff[_C10];
      //if (!notinit(c10_max) && lc10 > c10_max) lc10 = c10_max;
      //if (!notinit(c10_min) && lc10 < c10_min) lc10 = c10_min;
      lc11 = coeff[_C11];
      //if (!notinit(c11_max) && lc11 > c11_max) lc11 = c11_max;
      //if (!notinit(c11_min) && lc11 < c11_min) lc11 = c11_min;
      lc12 = coeff[_C12];
      //if (!notinit(c12_max) && lc12 > c12_max) lc12 = c12_max;
      //if (!notinit(c12_min) && lc12 < c12_min) lc12 = c12_min;

      //fprintf (stderr,"%s %s %g >& %s; %s\n", expt_string, p1, c1, model_output_file,data_split_string);
       
      switch (nparms) {
	   case 1: sprintf (exptcmd,"%s %s %18.18g",
		       		expt_string, p1, lc1); 
		   		break;
	   case 2: sprintf (exptcmd,"%s %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2);
		  		break;

	   case 3: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3);
		  		break;

	   case 4: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4);
		  		break;

	   case 5: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5);
		  		break;

	   case 6: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6);
		  		break;

	   case 7: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6, p7, lc7);
		  		break;

	   case 8: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6, p7, lc7, p8, lc8);
		  		break;

	   case 9: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6, p7, lc7, p8, lc8, p9, lc9);
		  		break;
	   case 10: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6, p7, lc7, p8, lc8, p9, lc9, p10, lc10);
		  		break;
	   case 11: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6, p7, lc7, p8, lc8, p9, lc9, p10, lc10, p11, lc11);
		  		break;
	   case 12: sprintf (exptcmd,"%s %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g %s %18.18g",
		       		expt_string, p1, lc1, p2, lc2, p3, lc3, p4, lc4, p5, lc5, p6, lc6, p7, lc7, p8, lc8, p9, lc9, p10, lc10, p11, lc11, p12, lc12);
		  		break;
      }
      if (n==0)
        sprintf (runcmd,"%s >& %s_%d\n", exptcmd, model_output_file, n);
      else 
        sprintf (runcmd,"%s >& %s_%d &\n", exptcmd, model_output_file, n);	// run multi proc in background

      //fprintf (stderr,"%s\n",expt_string);
      //fprintf (stderr,"%s\n",exptcmd);
      printf ("%s",runcmd);
      fflush (stdout);
      // fprintf (stderr,data_split_string);
      system (runcmd);
      if (n>0 && waitsecs>0) usleep(int(waitsecs*1e6));				// if run multi, wait
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int wait_for_model_data(int n)

/* Poll the running model and return 1 when it is done. */
/* If n==0, return immediately */
 
{
	int retval;
        char tailcmd[500] = {0};

   sprintf (tailcmd,"tail -n 1 %s_%d\n", model_output_file, n);
   for (retval=0; retval==0;) {		// poll for all models done
      if (strstr(xsystem(tailcmd),"done")!=NULL) retval = 1;
      if (retval==0) simwait (5);	// if not ready, sleep for 5 seconds before polling again
   }
   return retval;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double get_model_data(double x, double *coeff, int n)

/* Read model output file, process it into columns, then return yval associated with xval */
/* For n=0, models are not run in background */

{
      int i, xval, nrows, ncols, addr_val;
      static char model_output_data_file_n[FILENAMESIZ];
      double retval;

  xval = x;
  if (n<0) {
      run_model(coeff,-n);		// call with negative n starts the models running
      return 0;
  }
  if (xval==0) {
     if (n==0) {
       run_model(coeff,n);		// call with zero n starts model running, don't wait
     }
     else  wait_for_model_data(n);  // n > 0
     sprintf (data_split_cmd, data_split_string, n,n);
     xsystem (data_split_cmd);
     if (model_data_arr !=NULL) efree (model_data_arr);

     sprintf (model_output_data_file_n, model_output_data_file, n);
     if ((model_data_arr=fread (model_output_data_file_n, &nrows, &ncols))==NULL) {
	fprintf (stderr,"modelfit: can't read model output data file '%s'\n",model_output_data_file);
	exit (2);
     }
     filesize2d = nrows * ncols;
     if (template_arr2 != NULL && filesize2t >= filesize2d) { 
        for (i=0; i<filesize2d; i++) {
            model_data_arr[i] *= template_arr2[i]; // mask the data to limit the number of points.
        } 
     } 
  }
  if (model_data_arr!=NULL) {
       if (addr_arr!=NULL) addr_val = addr_arr[xval]; 	// locate xval from template address array
       else addr_val = xval;
       retval = model_data_arr[addr_val];
  }
  else retval = 0;
  // fprintf (stderr,"x %d retval %g addr_val %d\n",xval,retval, addr_val);
  return retval;
}

/*-----------------------------------------------------------------------------*/

double coeff[NCOEFF]     = {0};
double *coeffc;

main (int argc, char **argv) 
{
    int i, j, datasize, maxdata, fsize;
    int nlong, nwid, nlong2, nwid2;
    double *data_arr;
    double *xydata;
    const char *progname;
    FILE *fp;
    
  progname = argv[0];
  if (argc<=1) {
     fprintf (stderr,"Usage: modelfit <parameters>\n");
     fprintf (stderr,"	   --p1 xxx  --c1 n       coeff 1 (1.0)\n");
     fprintf (stderr,"	             --c1_max n   coeff 1 max \n");
     fprintf (stderr,"	             --c1_min n   coeff 1 min \n");
     fprintf (stderr,"	   --p2 xxx  --c2 n       coeff 2\n");
     fprintf (stderr,"	   --p3 xxx  --c3 n       coeff 3\n");
     fprintf (stderr,"	   --p4 xxx  --c4 n       coeff 4\n");
     fprintf (stderr,"	   --p5 xxx  --c5 n       coeff 5\n");
     fprintf (stderr,"	     --data_file xxx      file name (file.dat)\n");
     fprintf (stderr,"	     --template_file xxx  file name (file.dat)\n");
     fprintf (stderr,"	     --waitsecs n         wait time (1)\n");
     exit(0);
     }

  setptr("c1",      &c1);
  setptr("c1_max",  &c1_max);
  setptr("c1_min",  &c1_min);
  setptr("c2",      &c2);
  setptr("c2_max",  &c2_max);
  setptr("c2_min",  &c2_min);
  setptr("c3",      &c3);
  setptr("c3_max",  &c3_max);
  setptr("c3_min",  &c3_min);
  setptr("c4",      &c4);
  setptr("c4_max",  &c4_max);
  setptr("c4_min",  &c4_min);
  setptr("c5",      &c5);
  setptr("c5_max",  &c5_max);
  setptr("c5_min",  &c5_min);
  setptr("c6",      &c6);
  setptr("c6_max",  &c6_max);
  setptr("c6_min",  &c6_min);
  setptr("c7",      &c7);
  setptr("c7_max",  &c7_max);
  setptr("c7_min",  &c7_min);
  setptr("c8",      &c8);
  setptr("c8_max",  &c8_max);
  setptr("c8_min",  &c8_min);
  setptr("c9",      &c9);
  setptr("c9_max",  &c9_max);
  setptr("c9_min",  &c9_min);
  setptr("c10",      &c10);
  setptr("c10_max",  &c10_max);
  setptr("c10_min",  &c10_min);
  setptr("c11",      &c11);
  setptr("c11_max",  &c11_max);
  setptr("c11_min",  &c11_min);
  setptr("c12",      &c12);
  setptr("c12_max",  &c12_max);
  setptr("c12_min",  &c12_min);

  setptr("p1",      &p1);
  setptr("p2",      &p2);
  setptr("p3",      &p3);
  setptr("p4",      &p4);
  setptr("p5",      &p5);
  setptr("p6",      &p6);
  setptr("p7",      &p7);
  setptr("p8",      &p8);
  setptr("p9",      &p9);
  setptr("p10",     &p10);
  setptr("p11",     &p11);
  setptr("p12",     &p12);

  setptr("lm_multi_proc",  &lm_multi_proc);

  setptr("expt_string",   &expt_string);
  setptr("data_split",    &data_split);
  setptr("data_file",     &data_file);
  setptr("template_file", &template_file);
  setptr("template_file2", &template_file2);

  setptr("info", &info);
  setptr("filesize", &filesize);
  setptr("filesizet", &filesizet);
  setptr("tracelen", &tracelen);
  setptr("plotn",    &plotn);
  setptr("waitsecs", &waitsecs);

  setptr("model_output_file",  &model_output_file);

  /*
  c1 = 1e-9; c1_max = 1e9; c1_min = 0;
  c2 = 1e-9; c2_max = 1e9; c2_min = 0;
  c3 = 1e-9; c3_max = 1e9; c3_min = 0;
  c4 = 1;    c4_max = 1e9; c4_min = 0;
  c5 = 1;    c5_max = 1e9; c5_min = 0;
  c6 = 1;    c6_max = 1e9; c6_min = 0;
  c7 = 1;    c7_max = 1e9; c7_min = 0;
  c8 = 1;    c8_max = 1e9; c8_min = 0;
  c9 = 1;    c9_max = 1e9; c9_min = 0;
  c10 = 1;   c10_max = 1e9;c10_min = 0;
  c11 = 1;   c11_max = 1e9;c11_min = 0;
  c12 = 1;   c12_max = 1e9;c12_min = 0;

  c1_min = 0;
  c2_min = 0;
  c3_min = 0;
  c4_min = 0;
  c5_min = 0;
  c6_min = 0;
  c7_min = 0;
  c8_min = 0;
  c9_min = 0;
  c10_min = 0;
  c11_min = 0;
  c12_min = 0;
  */

  setvars(argc,argv);		// set variables from command line

  if (notinit(model_output_file)) {
       model_output_file = emalloc(FILENAMESIZ*sizeof(char));
       sprintf ((char *)model_output_file,"modelfit%06d",getpid());
  }

  if (notinit (lm_multi_proc)) lm_multi_proc = 1;

  nparms = NCOEFF;
  if (notinit(p12)) nparms--;
  if (notinit(p11)) nparms--;
  if (notinit(p10)) nparms--;
  if (notinit(p9)) nparms--;
  if (notinit(p8)) nparms--;
  if (notinit(p7)) nparms--;
  if (notinit(p6)) nparms--;
  if (notinit(p5)) nparms--;
  if (notinit(p4)) nparms--;
  if (notinit(p3)) nparms--;
  if (notinit(p2)) nparms--;

  if (notinit(expt_string)) expt_string = "retsim --expt dsgc_chans";
  if (notinit(p1)) p1 = "--ampa_cond";
  //if (notinit(p2)) p2 = "--nmda_cond";
  //if (notinit(p3)) p3 = "--gaba_cond";
  //if (notinit(p4)) p4 = "--dri";
  if (notinit(tracelen)) tracelen = 0;
  if (notinit(plotn)) plotn = 6;
  if (notinit(waitsecs)) waitsecs = 0.5;
//  if (notinit(data_split)) data_split = "plotmod -t -s 18 %s_%%d | plotsplit --plotn 18 --tracelen %d --info 0 > %s_%%d.data\n";
// if (notinit(data_split)) data_split = "plotmod -t -s 18 %s_%%d | plotsplit --plotn 6 --tracelen %d --info 0 > %s_%%d.data\n";
  if (notinit(data_split)) data_split = "plotmod -t -s 18 %s_%%d | plotsplit --plotn %d --tracelen %d --info 0 > %s_%%d.data\n";
  sprintf (model_output_data_file,"%s_%%d.data",model_output_file);
  sprintf (data_split_string,data_split,model_output_file, plotn, tracelen, model_output_file);

  if (notinit(data_file)) data_file = "file.data";	// original data file to be fitted
  if ((data_arr=fread (data_file, &nlong, &nwid))==NULL) {
	fprintf (stderr,"modelfit: can't read data file '%s'\n",data_file);
	exit (2);
  }	
  if (notinit(filesize)) filesize = nlong * nwid;

  if (!notinit(template_file)) {	// if template file specified for original data file 
    if ((template_arr=fread (template_file, &nlong, &nwid)) == NULL) {
	fprintf (stderr,"modelfit: can't read template file '%s'\n",template_file);
	exit (2);
    }	
    if (notinit(filesizet)) filesizet = nlong * nwid;
    for (i=j=0; i<filesizet; i++) {
       if (template_arr[i] > 0) j++; // count nonzero template weighting values.
    }
    datasize = j;	// limit size of data array to non-zero template weighting values
    fprintf (stderr,"datasize %d\n",datasize);
  }
  else  datasize = filesize;	     // no template file
  if ((xydata = (double *)emalloc(2*datasize*sizeof(double)))==NULL) 
	fprintf (stderr,"modelfit: can't allocate memory for array\n");

  if (!notinit(template_file)) {	// if data template file specified

    fsize = min(filesize,filesizet);
    for (i=j=0; i<fsize; i++) {
       if (template_arr[i] > 0) {
         xydata[j] = j;			// xydata[0][0..n] -> X values (2d matrix condensed into 1d)
         xydata[datasize + j] = data_arr[i] * template_arr[i]; // xydata[1][0..n] -> Z values to fit with ltsq
	 j++;
       }
    } 
  } else {				// no template_file
    for (i=0; i<datasize; i++) {
       xydata[i] = i;			// xydata[0][0..n] -> X values (2d matrix condensed into 1d)
       xydata[datasize+i] = data_arr[i];  // xydata[1][0..n] -> Z values to fit with least-squares 
    }
  }

  if (!notinit(template_file2)) {	// if data template file specified

    if ((addr_arr = (int *)emalloc(datasize*sizeof(int)))==NULL) {  // make address array
	fprintf (stderr,"modelfit: can't allocate memory for addr array\n");
        exit(2);
    }
    if ((template_arr2=fread (template_file2, &nlong2, &nwid2))==NULL) { // read fitting template 
	fprintf (stderr,"modelfit: can't read template2 file '%s'\n",template_file2);
	exit (2);
    }
    filesize2t = nlong2 * nwid2;

    for (i=j=0; i<filesize2t; i++) {
       if (template_arr2[i] > 0) {
	 addr_arr[j] = i;		// save the location of the non-zero fitting template values
	 j++;
       }
    }
    if (datasize != j) {
        fprintf (stderr,"modelfit: number of points in template files doesn't match %d %d\n",datasize,j);
        exit(2);
    }
  }

  // Starting coefficients:
  
  coeff[_C1] = c1;
  coeff[_C2] = c2;
  coeff[_C3] = c3;
  coeff[_C4] = c4;
  coeff[_C5] = c5;
  coeff[_C6] = c6;
  coeff[_C7] = c7;
  coeff[_C8] = c8;
  coeff[_C9] = c9;
  coeff[_C10] = c10;
  coeff[_C11] = c11;
  coeff[_C12] = c12;
 

  if (!notinit(c1_min)  || !notinit(c1_max) ||		// if any constraints have been set
      !notinit(c2_min)  || !notinit(c2_max) ||
      !notinit(c3_min)  || !notinit(c3_max) ||
      !notinit(c4_min)  || !notinit(c4_max) ||
      !notinit(c5_min)  || !notinit(c5_max) ||
      !notinit(c6_min)  || !notinit(c6_max) ||
      !notinit(c7_min)  || !notinit(c7_max) ||
      !notinit(c8_min)  || !notinit(c8_max) ||
      !notinit(c9_min)  || !notinit(c9_max) ||
      !notinit(c10_min) || !notinit(c10_max) ||
      !notinit(c11_min) || !notinit(c11_max) ||
      !notinit(c12_min) || !notinit(c12_max)   ) { 

  if ((coeffc=(double*)emalloc (NCOEFF*2*sizeof(double)))==NULL) {
	fprintf (stderr,"modelfit: can't allocate memory for constraints array\n");
	exit(2);
  }

#define LC 0
#define UC 1

  coeffc[_C1*2+LC] = c1_min;
  coeffc[_C1*2+UC] = c1_max;
  coeffc[_C2*2+LC] = c2_min;
  coeffc[_C2*2+UC] = c2_max;
  coeffc[_C3*2+LC] = c3_min;
  coeffc[_C3*2+UC] = c3_max;
  coeffc[_C4*2+LC] = c4_min;
  coeffc[_C4*2+UC] = c4_max;
  coeffc[_C5*2+LC] = c5_min;
  coeffc[_C5*2+UC] = c5_max;
  coeffc[_C6*2+LC] = c6_min;
  coeffc[_C6*2+UC] = c6_max;
  coeffc[_C7*2+LC] = c7_min;
  coeffc[_C7*2+UC] = c7_max;
  coeffc[_C8*2+LC] = c8_min;
  coeffc[_C8*2+UC] = c8_max;
  coeffc[_C9*2+LC] = c9_min;
  coeffc[_C9*2+UC] = c9_max;
  coeffc[_C10*2+LC] = c10_min;
  coeffc[_C10*2+UC] = c10_max;
  coeffc[_C11*2+LC] = c11_min;
  coeffc[_C11*2+UC] = c11_max;
  coeffc[_C12*2+LC] = c12_min;
  coeffc[_C12*2+UC] = c12_max;
#undef LC
#undef UC
  }
  else coeffc = NULL;
  
 lmfit (get_model_data, datasize, xydata, nparms, coeff, coeffc);
}
