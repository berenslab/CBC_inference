
/* lm_funcs.cc */

#include <stdio.h>

#include "ncio.h"
#include "lmmin.h"
#include "lm_eval.h"

extern int info;

/*------------------------------------------------------*/

void lmfit(double (*user_func) (double user_x_point, double *coeff, int model_num), 
		int m_dat, double *xydata, int n_p, double *coeff, double *coeffc)

/* Levenberg-Marquardt curve fitting on an array. */
/* Adapted from lmfit-2.4 from sourceforge.net */

/* Fit 1D function y = user_func2d(x,coeff).
   xydata array must be: arr[2][m_dat], arr[0] is x, and arr[1] is y.
      coeff array must be: coeff[n_p], which is starting value of points.  
*/
{
	lm_control_type control;
	lm_data_type data;

   lm_initialize_control(&control);

   data.user_func = user_func;		/* the user-defined fit function */
   data.user_t = xydata;
   data.user_y = xydata+m_dat;
   data.user_n_p = n_p;

   info=2;
   lm_minimize (m_dat,n_p,coeff,coeffc,lm_evaluate_default,lm_print_default,&data,&control);
   if (info>=1) ncfprintf(stdout,"status: %s after %d evaluations\n",
		                 lm_shortmsg[control.info], control.nfev);

}

/*------------------------------------------------------*/

void lmfit2d(double (*user_func2d) (double user_x_point, double user_y_point, double *coeff, int model_num), 
			int m_dat, double *xydata, int n_coeff, double *coeff, double *coeffc)

/* Levenberg-Marquardt curve fitting on an array. */
/* Adapted from lmfit-2.4 from sourceforge.net */

/* Fit 2D function z = user_func2d(x,y,coeff).
   xydata array must be: arr[2][m_dat], arr[0] is x [0...xysiz], 
   followed by y [xysiz+0...xysiz*2], where xysiz=sqrt(m_dat), and arr[1] is z[0...m_dat].
      coeff array must be: coeff[n_p], which is starting value of points.  
*/
{
	lm_control_type control;
	lm_data_type data;

   lm_initialize_control(&control);

   data.user_func2d = user_func2d;		/* the user-defined fit function */
   data.user_t = xydata;
   data.user_y = xydata+m_dat;
   data.user_n_p = n_coeff;

   lm_minimize (m_dat,n_coeff,coeff,coeffc,lm_evaluate_2d,lm_print_2d,&data,&control);
   if (info>=1) ncfprintf(stdout,"status: %s after %d evaluations\n",
		                 lm_shortmsg[control.info], control.nfev);
}

