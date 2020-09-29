// cone_ode
//
//     Program accompanying 'Simulating human cones from mid-mesopic
//     up to high-photopic luminances' (J.H. van Hateren, 2006)
//
//     Calculates normalized output of human cones (>response) to input time   
//     series (with intensities in trolands at 1 ms resolution)
//
//     Implemented in matlab using standard ODE integration routines
//     Uses: steady_ode.m, cone_func.m
//
//     Usage: cone_ode < stimulus; (plots results of ODE)
//            see further README file and article
//
//     Adapted from:
//
//     J.H van Hateren, 28/8/06
//
//     van Hateren JH, Snippe HP (2006) Simulating human cones from 
//       mid-mesopic up to high-photopic luminances. Journal of Vision 7(4):1-11.
//
//     R.G. Smith  2011-09-29

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double tau_r=0;
static double tau_b=0;
static double tau_c=0;
static double tau_vc=0;
static double tau_e=0;
static double rk_b=0;
static double cn=0;
static double beta_e=0;
static double beta_e_max=0;
static double a_c=0;
static double rnc=0;
static double c_beta=0;
static double rk_beta=0;
static double rnx=0;

static double a_is=0;
static double gamma_is=0;
static double tau_is=0;

static double resp_r=0;
static double resp_b=0;
static double resp_e=0;
static double resp_x=0;
static double resp_c=0;
static double resp_os=0;
static double resp_is=0;
static double atten_i=0;

static double stim=0;

/*----------------------------------------------------------------------*/

void cone_func(double stim)

{
   double alpha, beta;
   double resp_im;

  resp_r  += (stim* (1-resp_b-cn*resp_r) - resp_r)/tau_r;
  resp_b  += cn*resp_r/tau_r - rk_b/tau_b * resp_b/(resp_b+rk_b);
  resp_e  += (resp_r-resp_e)/tau_e;

  beta     = c_beta + rk_beta*resp_e;
  beta_e   = beta / (1 + beta/beta_e_max);
  alpha    = 1 / (1 + pow(a_c*resp_c,rnc));

  resp_x  += alpha - beta_e*resp_x;
  resp_os  = pow(resp_x,rnx);

  resp_c  += (resp_os - resp_c)/tau_c;

  resp_is += (resp_os/atten_i - resp_is)/tau_vc;
  resp_im  = pow(resp_is,gamma_is);

  atten_i += (a_is*resp_im - atten_i)/tau_is;

 //printf ("%g %g %g %g %g %g %g %g\n",stim,resp_r, resp_b, resp_e, beta_e, resp_os, resp_c, resp_is);
}

/*----------------------------------------------------------------------*/

// steady_ode
//
// finds steady-state: steady=0 if x=steady-state of resp_c
//
//     Function accompanying 'Simulating human cones from mid-mesopic
//     up to high-photopic luminances' (J.H. van Hateren, 2006)
//
//     Used by: cone_ode.m
//
//     Author: J.H. van Hateren, 28/8/06
//

double steady_ode (double x)

{
    double rval, beta, stimn;

  if (stim==0) resp_e=0;
  else {
   stimn=stim*cn;
   resp_b=0.5 * (1-rk_b-tau_r/tau_b*rk_b*(stimn+1)/stimn +
        sqrt(pow(1-rk_b-tau_r/tau_b*rk_b*(stimn+1)/stimn,2)+4*rk_b));
   resp_r=(1-resp_b)*stim/(1+stimn);
   resp_e=resp_r;
  }
  beta=c_beta+rk_beta*resp_e;
  beta_e=beta/(1+beta/beta_e_max);
  rval = x - pow((1/(1+pow(a_c*x,rnc)))/beta_e,rnx);
  return rval;
}

/*----------------------------------------------------------------------*/

double fzero (double (*func)(double), double near_val)
{
    double x0, x1, x2;
    double crit = 1e-6;
  
    for (x0=near_val,x1=x0+1; fabs(x1-x0) > crit; ) {
	 // printf ("      x0 %g %g\n",x0,x1);
         x2 = x1 - func(x1) * (x1 - x0) / (func(x1) - func(x0));
	 x0 = x1; x1 = x2;
    }
  // printf ("final x0 %g %g\n",x0,x1);
  return x0;
}	


/*----------------------------------------------------------------------*/

int main (int argc, char **argv)
{
   int i;
   double stimn, tau=1.3;
   double beta, resp_q, gain_x, resp_im;
   double dark_resp_os,dark_resp_is;

  tau_r=3.4*tau;          // R* lifetime (ms)
  cn=4.1e-9*tau_r/3.4;    // normalization constant for cone-
                          // bleaching (4.1e-9 if tau_r=3.4 ms)
  tau_b=25000*tau;        // time constant of bleaching recovery
  tau_e=8.7*tau;          // E* lifetime
  tau_c=3*tau;            // time constant of Ca2+ extrusion
  tau_vc=4*tau;           // membrane time constant
  tau_is=90*tau;          // time constant of membrane nonlinearity
  rk_b=0.2*tau;           // parameter of bleaching recovery
  c_beta=2.8e-3/tau;      // dark PDE activity (=1/tau_D)
  rk_beta=1.4e-4/tau;     // E* dependence of PDE activity
  beta_e_max=4/tau;       // parameter diffusion-limited cGMP hydrolysis
  rnx=1;                  // apparent Hill coefficient CNG channels
  rnc=4;                  // Hill coefficient GC activation
  a_c=0.23;               // scaling constant of GC activation
  gamma_is=0.7;           // parameter of membrane nonlinearity
  a_is=2.9e-2;            // parameter of membrane nonlinearity

  // get dark value of resp_is
 
  stim = 0;
  resp_c = fzero(steady_ode,1);
  // printf ("# dark resp_c %g\n",resp_c);

  beta    = c_beta;
  beta_e  = beta / (1+beta/beta_e_max);
  resp_q  = 1/beta_e;
  gain_x  = 1/(1 + pow(a_c*resp_c,rnc));
  resp_x  = gain_x*resp_q;
  dark_resp_os = pow(resp_x,rnx);
  dark_resp_is= pow(dark_resp_os/a_is,(1/(1+gamma_is)));  
 
  // printf ("# dark resp_is %g\n",dark_resp_is);

  // determine adaptive state

  stim=10000;

  resp_c = fzero(steady_ode,1);
  printf ("# adapted resp_c %g\n",resp_c);

  if (stim==0) {
    resp_b  = 0;
    resp_r  = 0;
    resp_e  = 0;
  }
  else {
    stimn  = stim * cn;
    resp_b = 0.5 * (1-rk_b-tau_r/tau_b*rk_b*(stimn+1)/stimn +
         sqrt(pow(1-rk_b-tau_r/tau_b*rk_b*(stimn+1)/stimn,2)+4*rk_b));
   resp_r = (1-resp_b)*stim/(1+stimn);
   resp_e = resp_r;
  }
  beta    = c_beta + rk_beta*resp_e;
  beta_e  = beta/(1 + beta/beta_e_max);
  resp_q  = 1/beta_e;
  gain_x  = 1/(1+pow(a_c*resp_c,rnc));
  resp_x  = gain_x*resp_q;
  resp_os = pow(resp_x,rnx);
  resp_is = pow(resp_os/a_is,(1/(1+gamma_is)));
  resp_im = pow(resp_is,gamma_is);
  atten_i = a_is*resp_im;

  while (scanf("%lg",&stim) > 0) {
     cone_func(stim);
     printf ("%g\n",(resp_is-dark_resp_is)/dark_resp_is);
     // printf ("%g\n",resp_os);
  }
    
}

