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
#include <stdio.h>
#include <math.h>

/*-------------------------------------------------------------------------*/

void lm_evaluate_default(double *par, int m_dat, double *fvec,
			 void *data, int *info, int model_num)
/* 
 *	par is an input array. At the end of the minimization, it contains
 *        the approximate solution vector.
 *
 *	m_dat is a positive integer input variable set to the number
 *	  of functions.
 *
 *	fvec is an output array of length m_dat which contains the function
 *        values the square sum of which ought to be minimized.
 *
 *	data is a read-only pointer to lm_data_type, as specified by lm_eval.h.
 *
 *      info is an integer output variable. If set to a negative value, the
 *        minimization procedure will stop.
 *
 *      model_num is an integer that represents which param is being tested
 *        for the purpose of tracking the models running in parallel in the bknd.
 */
{
    int i;
    lm_data_type *mydata;
    mydata = (lm_data_type *) data;

    printf ("evaluate: lm_multi %d model_num %d\n",lm_multi_proc,model_num);

    if (lm_multi_proc && model_num < 0) mydata->user_func(mydata->user_t[0], par, model_num);
    else {
      for (i = 0; i < m_dat; i++)
	fvec[i] = mydata->user_y[i]
	    - mydata->user_func(mydata->user_t[i], par, model_num);

    *info = *info;		/* to prevent a 'unused variable' warning */
    /* if <parameters drifted away> { *info = -1; } */
    }
}

void lm_evaluate_default(double *par, int m_dat, double *fvec,
			 void *data, int *info)
{
     int model_num = 0;

    lm_evaluate_default(par, m_dat, fvec, data, info,  model_num);
}

/*-------------------------------------------------------------------------*/

void lm_evaluate_2d(double *par, int m_dat, double *fvec,
			 void *data, int *info, int model_num)
/* 
 *      2D (x,y) evaluate function.
 *
 *	par is an input array. At the end of the minimization, it contains
 *        the approximate solution vector.
 *
 *	m_dat is a positive integer input variable set to the number
 *	  of function values.
 *
 *	fvec is an output array of length m_dat which contains the function
 *        values the square sum of which ought to be minimized.
 *
 *	data is a read-only pointer to lm_data_type, as specified by lm_eval.h.
 *
 *      info is an integer output variable. If set to a negative value, the
 *        minimization procedure will stop.
 */
{
    int i,j,xysiz;
    double val;
    lm_data_type *mydata;
    mydata = (lm_data_type *) data;

  if (lm_multi_proc && model_num < 0) mydata->user_func(mydata->user_t[0], par, model_num);
  else {
    xysiz = sqrt(m_dat);
    for (i = 0; i < xysiz; i++) {
      for (j = 0; j < xysiz; j++) {
	val = mydata->user_y[i*xysiz+j]
	    - mydata->user_func2d(mydata->user_t[i],mydata->user_t[xysiz+j], par, model_num);
	fvec[i*xysiz+j] = val;
      }
    }
    *info = *info;		/* to prevent a 'unused variable' warning */
    /* if <parameters drifted away> { *info = -1; } */
  }
}

void lm_evaluate_2d(double *par, int m_dat, double *fvec,
			 void *data, int *info)
{
     int model_num = 0;

    lm_evaluate_2d(par, m_dat, fvec, data, info,  model_num);
}


/*-------------------------------------------------------------------------*/

void lm_print_default(int n_par, double *par, int m_dat, double *fvec,
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
	for (i = 0; i < m_dat; ++i) {
	    t = (mydata->user_t)[i];
	    y = (mydata->user_y)[i];
	    f = mydata->user_func(t, par, 0);
	    printf("    t[%2d]=%12g y=%12g fit=%12g residue=%12g\n",
		   i, t, y, f, y - f);
	}
    }
}

/*-------------------------------------------------------------------------*/

void lm_print_2d(int n_par, double *par, int m_dat, double *fvec,
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
	for (i = 0; i < xysiz; i++) 
	  for (j = 0; j < xysiz; j++) {
	    x = (mydata->user_t)[i];
	    y = (mydata->user_t)[xysiz+j];
	    z = (mydata->user_y)[i*xysiz+j];
	    f = mydata->user_func2d(x, y, par, 0);
	    printf("  t[%2d,%2d]=%5g,%5g z=%12g fit=%12g residue=%12g\n",
		   i, j, x, y, z, f, z - f);
	}
    }
}

/*-------------------------------------------------------------------------*/
