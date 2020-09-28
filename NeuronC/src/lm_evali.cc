
/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lm_eval.c
 *
 * Contents: Default implementation of 
 *           - user-defined data evalution function,
 *           - user-defined printout routine.
 *           Declarations are in lm_eval.h.
 *
 * Usage:    Least-squares fitting of a simple one-dimensional data set
 *           is shown in lm_test.c. For other applications, the data 
 *           must be modified.
 *
 * Author:   Joachim Wuttke 2004-8 
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Public domain.
 * Version:  lmfit-2.4
 */
 
#include "nc.h"
#include "lmmin.h"
#include "lm_eval.h"
#include "control.h"
#include <stdio.h>
#include <math.h>

double callfuncxp  (Symbol *funcp, double xval, ncarray *par, int model_num);
double callfuncxp2d(Symbol *funcp, double xval, double yval, ncarray *par, int model_num);

/*-------------------------------------------------------------------------*/

void lm_evaluate_interp(double *par, int m_dat, double *fvec,
		                         void *data, int *info, int model_num)

{
    int i;
    lm_data_type *mydata;
    mydata = (lm_data_type *) data;
    ncarray coeff;

   // fprintf (stderr,"B, %d  %x\n",mydata->user_n_p, mydata);
			    
   coeff.arr = par;
   coeff.ndim = 1;
   coeff.dim[0] = mydata->user_n_p;
   for (i = 0; i < m_dat; i++) {
       fvec[i] = mydata->user_y[i] -
       callfuncxp(mydata->user_funcp, mydata->user_t[i], &coeff, model_num);
    }
    *info = *info; /* to prevent a 'unused variable' warning */
}
			    
/*-------------------------------------------------------------------------*/

void lm_evaluate_2d_interp(double *par, int m_dat, double *fvec,
		                         void *data, int *info, int model_num)

{
    int i,j,xysiz;
    lm_data_type *mydata;
    mydata = (lm_data_type *) data;
    ncarray coeff;

   // fprintf (stderr,"B, %d  %x\n",mydata->user_n_p, mydata);
			    
   coeff.arr = par;
   coeff.ndim = 1;
   coeff.dim[0] = mydata->user_n_p;

   xysiz = sqrt(m_dat);
   for (i = 0; i < xysiz; i++) {
     for (j = 0; j < xysiz; j++) {
       fvec[i] = mydata->user_y[i*xysiz+j] -
       callfuncxp2d(mydata->user_funcp2d, mydata->user_t[i], mydata->user_t[xysiz+j], &coeff, model_num);
     }
   }
   *info = *info; /* to prevent a 'unused variable' warning */
}
			    
/*-------------------------------------------------------------------------*/

void lm_print_interp(int n_par, double *par, int m_dat, double *fvec,
		      void *data, int iflag, int iter, int nfev)
/*
 *       data  : for soft control of printout behaviour, add control
 *                 variables to the data struct
 *       iflag : 0 (init) 1 (outer loop) 2(inner loop) -1(terminated)
 *       iter  : outer loop counter
 *       nfev  : number of calls to *evaluate
 */
{
    double f, y, t;
    int i;
    lm_data_type *mydata;
    mydata = (lm_data_type *) data;
    ncarray coeff;

  if (info >=2) {
    if (iflag == 2) {
	printf("trying step in gradient direction\n");
    } else if (iflag == 1) {
	printf("determining gradient (iteration %d)\n", iter);
    } else if (iflag == 0) {
	printf("starting minimization\n");
    } else if (iflag == -1) {
	printf("terminated after %d evaluations\n", nfev);
    }

    printf("  par: ");
    for (i = 0; i < n_par; ++i)
	printf(" %12g", par[i]);
    printf(" => norm: %12g\n", lm_enorm(m_dat, fvec));

    if (iflag == -1) {
	printf("  fitting data as follows:\n");
        coeff.arr = par;
        coeff.ndim = 1;
        coeff.dim[0] = mydata->user_n_p; 
	for (i = 0; i < m_dat; ++i) {
	    t = (mydata->user_t)[i];
	    y = (mydata->user_y)[i];
	    f = callfuncxp(mydata->user_funcp, t, &coeff, 0);
	    printf("    t[%2d]=%12g y=%12g fit=%12g residue=%12g\n",
		   i, t, y, f, y - f);
	}
    }
  }
}

/*-------------------------------------------------------------------------*/

void lm_print_2d_interp(int n_par, double *par, int m_dat, double *fvec,
		      void *data, int iflag, int iter, int nfev)
/*
 *       data  : for soft control of printout behaviour, add control
 *                 variables to the data struct
 *       iflag : 0 (init) 1 (outer loop) 2(inner loop) -1(terminated)
 *       iter  : outer loop counter
 *       nfev  : number of calls to *evaluate
 */
{
    double f, x, y, z;
    int i, j, xysiz;
    lm_data_type *mydata;
    mydata = (lm_data_type *) data;
    ncarray coeff;

  if (info >=2) {
    if (iflag == 2) {
	printf("trying step in gradient direction\n");
    } else if (iflag == 1) {
	printf("determining gradient (iteration %d)\n", iter);
    } else if (iflag == 0) {
	printf("starting minimization\n");
    } else if (iflag == -1) {
	printf("terminated after %d evaluations\n", nfev);
    }

    printf("  par: ");
    for (i = 0; i < n_par; ++i)
	printf(" %12g", par[i]);
    printf(" => norm: %12g\n", lm_enorm(m_dat, fvec));

    if (iflag == -1) {
	printf("  fitting data as follows:\n");
	printf("------------   x     y      z                f                  (z-f)\n");
        xysiz = sqrt(m_dat);
        coeff.arr = par;
        coeff.ndim = 1;
        coeff.dim[0] = mydata->user_n_p; 
	for (i = 0; i < m_dat; ++i) {
	    x = (mydata->user_t)[i];
	    y = (mydata->user_t)[xysiz+j];
	    z = (mydata->user_y)[i*xysiz+j];
	    f = callfuncxp2d(mydata->user_funcp, x, y, &coeff, 0);
	    printf("    t[%2d,%2d]=%5g,%5g z=%12g fit=%12g residue=%12g\n",
		   i, j, x, y ,z, f, y - f);
	}
    }
  }
}

/*-------------------------------------------------------------------------*/
