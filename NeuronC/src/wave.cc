/* wave.cc: a program to produce relative spectral sensitivity
   arrays for the three cone pigment types.
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif
void exit (int n);

#ifdef __cplusplus
}
#endif

#include "ncsub.h"
#include "ncomp.h"

/* double log(double), exp(double); */
double lamm;

#define LN10	2.30258509299404568402
#define LAMM 561.0

#define TURTLECONES 1	/* =1 -> add pigment types 4,5,6 for turtle cones */
#define GOLDFISHCONES 1	/* =1 -> add pigment types 4,5,6 for goldfish cones */

void calcwav(int pigm);
void copywavt(void);
void copywavg(void);
void wavout(double val);
void arrfill (double val, int pigm);

static int p=0;

double pigmarr[TOTREC][PIGMSIZ] = {0};

double lumlight [LUMINS][LIGHTS] = {0};

double pigmavg[TOTREC][LIGHTS] = {0};

			/* OD values taken from Bowmaker et al 1980 *:

			/* Bowmaker JK, Dartnall HJ, Mollon JD. (1980) Microspectrophotometric demonstration of  */
			/* four classes of photoreceptor in an old world primate, Macaca fascicularis.  */
			/* J Physiol. 298:131-143. */

//  			/* The old method: action spectrum is for absolute OD */ 
//
//			/* This table gives log10 of abs sensitivity at lmax: */
//
//double pigmconst[TOTREC] = {	
//				-0.19031465,	/* rod log10(1-exp10(-.450)) */
//				-0.153996,	/* red                .525   */
//				-0.153996,	/* green	      .525   */
//				//-0.23784418,	/* blue		      .375   */
//				-0.153996,	/* blue		      .525   */
//				-0.15490196,	/* turtle red         .5228  */
//				-0.15490196,	/* turtle green       .5228  */
//				-0.15490196,	/* turtle blue        .5228  */
//				-0.153996,	/* red   without cascade */	
//				-0.153996,	/* green without cascade */	
//				-0.23784418,	/* blue  without cascade */	
//				-0.37478641,	/* rabbit rod         0.238   */
//				-0.88767125,	/* rabbit M cone      0.06024 */
//				-0.88767125,	/* rabbit S cone      0.06024 */
//				-0.37478641,	/* guinea pig rod     0.238   */
//				-0.88767125,	/* guinea pig M cone  0.06024 */
//				-0.88767125	/* guinea pig S cone  0.06024 */
//				-0.88767125	/* goldfish L cone    0.06024 */
//				-0.88767125	/* goldfish M cone    0.06024 */
//				-0.88767125	/* goldfish S cone    0.06024 */
//				-0.88767125	/* van hateren L cone 0.06024 */
//				-0.88767125	/* van hateren M cone 0.06024 */
//				-0.88767125	/* van hateren S cone 0.06024 */
//				-0.19031465,	/* invergo rod log10(1-exp10(-.450)) */
//				-0.19031465,	/* ChR         log10(1-exp10(-.450)) */
//			    } ;	

			/* New method: action spectrum is for spec OD */ 
			/* This table gives peak spec OD (per um) for each pigment at lmax: */
double peakspecod[TOTREC] = {	
				0.018,		/* rod log10(1-exp10(-0.018)) */
				0.016,		/* 1 L red                0.016   */
				0.016,		/* 2 M green	      0.016   */
				0.016,		/* 3 S blue	      0.016   */
				0.015,		/* 4 turtle red         0.015  */
				0.015,		/* 5 turtle green       0.015  */
				0.015,		/* 6 turtle blue        0.015  */
				0.016,		/* 7 red   without cascade */	
				0.016,		/* 8 green without cascade */	
				0.016,		/* 9 blue  without cascade */	
				0.017,		/* 10 rabbit rod         0.017    */
				0.008032,	/* 11 rabbit M cone      0.008032 */
				0.008032,	/* 12 rabbit S cone      0.008032 */
				0.017,		/* 13 guinea pig rod     0.017    */
				0.008032,	/* 14 guinea pig M cone  0.008032 */
				0.008032,	/* 15 guinea pig S cone  0.008032 */
				0.016,		/* 16 goldfish L         0.016    */
				0.016,		/* 17 goldfish M         0.016    */
				0.016,		/* 18 goldfish S         0.016    */
				0.016,		/* 19 van hateren L      0.016    */
				0.016,		/* 20 van hateren M      0.016    */
				0.016,		/* 21 van hateren S      0.016    */
				0.018,		/* 22 invergo mouse rod  0.018    */
				0.018		/* 23 ChR                0.018    */
			    } ;	
double pigmlen[TOTREC] = {	
				25.0,		/* rod path length */
				35.0,		/* red cone path length */
				35.0,		/* green cone path length */
				35.0,		/* blue cone path length */
				35.0,		/* turtle red cone path len */
				35.0,		/* turtle green cone path len */
				35.0,		/* turtle blue cone path len */
				35.0,		/* red   without cascade */	
				35.0,		/* green without cascade */	
				35.0,		/* blue  without cascade */	
				14.0,		/* rabbit rod path len */
				 7.5,		/* rabbit M cone path len */
				 7.5,		/* rabbit S cone path len */
				14.0,		/* guinea pig rod path len */
				 7.5,		/* guinea pig M cone path len */
				 7.5,		/* guinea pig S cone path len */
				 7.5,		/* goldfish L cone path len */
				 7.5,		/* goldfish M cone path len */
				 7.5,		/* goldfish S cone path len */
				35.0,		/* van hateren L cone path len */
				35.0,		/* van hateren M cone path len */
				35.0,		/* van hateren S cone path len */
				25.0,		/* invergo mouse rod path length */
				25.0		/* ChR path length */
			    } ;

double lights[LIGHTS][PIGMSIZ] = {		/* standard light inten */
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  /* xenon */
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  /* spectrally uniform */
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0,

	/* These spectral light curves were calculated by "blackbody.n" */

0.8391, 0.8545,  0.869, 0.8827, 0.8956, 0.9076, 0.9189, 0.9293, /* sun */
0.9389, 0.9478, 0.9559, 0.9632, 0.9698, 0.9757, 0.9809, 0.9854, /* 5900 K */ 
0.9893, 0.9925, 0.9952, 0.9972, 0.9987, 0.9996,      1, 0.9999, 
0.9993, 0.9982, 0.9967, 0.9948, 0.9925, 0.9898, 0.9867, 0.9833, 
0.9795, 0.9754, 0.9711, 0.9664, 0.9615, 0.9563, 0.9509, 0.9453, 
0.9395, 0.9335, 0.9273,  0.921, 0.9145, 0.9078, 0.9011, 0.8942, 
0.8872, 0.8801, 0.8729, 0.8656, 0.8582, 0.8508, 0.8433, 0.8358, 
0.8282, 0.8206,  0.813, 0.8053, 0.7977,   0.79, 0.7823, 0.7746, 
0.7669, 0.7592, 0.7515, 0.7439, 0.7362, 0.7286,  0.721, 0.7134, 
0.7059, 0.6984, 0.6909, 0.6835, 0.6761, 0.6688, 0.6615, 0.6543, 
0.6471, 0.6399, 0.6328, 0.6258, 0.6188,

0.1186, 0.1284, 0.1386, 0.1492, 0.1602, 0.1715, 0.1833, 0.1953, /* tungsten */
0.2077, 0.2205, 0.2335, 0.2468, 0.2603, 0.2741, 0.2881, 0.3023, /* 3400 K */
0.3167, 0.3313,  0.346, 0.3608, 0.3757, 0.3907, 0.4058, 0.4208, 
 0.436, 0.4511, 0.4662, 0.4813, 0.4963, 0.5113, 0.5262,  0.541, 
0.5557, 0.5702, 0.5847, 0.5989, 0.6131,  0.627, 0.6408, 0.6543, 
0.6677, 0.6809, 0.6938, 0.7065,  0.719, 0.7312, 0.7432,  0.755, 
0.7664, 0.7777, 0.7886, 0.7993, 0.8097, 0.8198, 0.8297, 0.8393, 
0.8486, 0.8576, 0.8663, 0.8748,  0.883, 0.8909, 0.8985, 0.9058, 
0.9129, 0.9197, 0.9262, 0.9325, 0.9384, 0.9442, 0.9496, 0.9548, 
0.9597, 0.9644, 0.9688,  0.973, 0.9769, 0.9806,  0.984, 0.9873, 
0.9902,  0.993, 0.9955, 0.9978, 1, 

	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
	1.0, 1.0, 1.0, 1.0, 1.0,
	};
/*
  Black-body radiation:

  B = 2 * h * c^2 / ( l^5 * exp (h * c / ( k * l * T)) - 1)  W/(m^2*sr*um)

  c = 2.998e8
  l = wavelength (m)
  h = 6.626e-34
  k = 1.380e-23

Ref: http://www.star.ait.ac.th/star/rsnotes/CP1/T1-7-1.gif

*/




double luminos[LUMINS][PIGMSIZ] = {		/* standard luminosities */

	/* scotopic luminosity function (Wyszecki and Stiles, 1982) */

  5.89e-4, 1.108e-3, 2.209e-3, 4.53e-3, 9.29e-3, 1.852e-2, 3.484e-2, 6.04e-2,
  9.66e-2, 0.1436, 0.1998, 0.2625, 0.3281, 0.3931, 0.4550, 0.5130,
  0.5670, 0.6200, 0.6760, 0.7340, 0.7930, 0.8510, 0.9040, 0.9490,
  0.9820, 0.9980, 0.9970, 0.9750, 0.9350, 0.8800, 0.8110, 0.7330,
  0.6500, 0.5640, 0.4810, 0.4020, 0.3288, 0.2639, 0.2076, 0.1602,
  0.1212, 8.99e-2, 6.55e-2, 4.69e-2, 3.315e-2, 2.312e-2, 1.593e-2, 1.088e-2,
  7.37e-3, 4.97e-3, 3.335e-3, 2.235e-3, 1.497e-3, 1.005e-3, 6.77e-4, 4.59e-4,
  3.129e-4, 2.146e-4, 1.480e-4, 1.026e-4, 7.15e-5, 5.01e-5, 3.533e-5, 2.501e-5,
  1.780e-5, 1.273e-5, 9.14e-6, 6.60e-6, 4.78e-6, 3.482e-6, 2.546e-6, 1.870e-6,
  1.379e-6, 1.022e-6, 7.60e-7, 5.67e-7, 4.25e-7, 3.196e-7, 2.413e-7, 1.829e-7,
  1.390e-7, 1.01e-7, 7.7e-8, 5.7e-8, 4.3e-8,
	
	/* photopic luminosity function (Wyszecki and Stiles, 1982) */

	0.0000, 0.0000, 0.0001, 0.0003, 0.0004, 0.0009, 0.0012, 0.0028,
	0.0040, 0.0080, 0.0116, 0.0174, 0.0230, 0.0305, 0.0380, 0.0490,
	0.0600, 0.0755, 0.0910, 0.1150, 0.1390, 0.1735, 0.2080, 0.2655,
	0.3230, 0.4130, 0.5030, 0.6065, 0.7100, 0.7860, 0.8620, 0.9080,
	0.9540, 0.9745, 0.9950, 0.9950, 0.9950, 0.9735, 0.9520, 0.9110,
	0.8700, 0.8135, 0.7570, 0.6940, 0.6310, 0.5670, 0.5030, 0.4420,
	0.3810, 0.3230, 0.2650, 0.2200, 0.1750, 0.1410, 0.1070, 0.0840,
	0.0610, 0.0465, 0.0320, 0.0250, 0.0170, 0.0130, 0.0082, 0.0062,
	0.0041, 0.0031, 0.0021, 0.0016, 0.0010, 0.0008, 0.0005, 0.0004,
	.0003, 0.0002, 0.0001, 0.0001, 0.0001, 0.0000, 0.0000, 0.0000,
	0.0000, 0.0000, 0.0000, 0.0000, 0.0000,

	};

/*--------------------------------------------------------------*/

	/* Filter transmission tables: */

        /* transmission = exp   (-od * ln(10)) */
        /* transmission = exp10 (-od) */

double filts[FILTS][PIGMSIZ] = {

	/* 100% transmission filter */

	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0, 

	/* Macular pigment optical density */

	/* Macular densities from Wyszecki & Stiles (1982),Table II(2.4.6) */
	/* Check "odconv.n" to find calculations. */

	/* Extrapolated in the range 400 -> 380 nm (RGS)         */


	-0.025, -0.035, -0.045,  -0.06, -0.085,  -0.12,  -0.16, -0.225, 
	-0.3, -0.345, -0.365,  -0.38,   -0.4, -0.425,  -0.46,  -0.49, 
	-0.495,  -0.47, -0.445,  -0.41, -0.415,  -0.42,  -0.41,  -0.36, 
	-0.275, -0.195,  -0.13, -0.085,  -0.05, -0.025,  -0.01,     -0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0, 

	/* Ocular media optical density */
	/* Lens densities from Stockman, MacLeod and Johnson (1993),   */

	/* Proposed adjustment to the lens templates of Wyszecki & Stiles
	(1967) and van Norren & Vos (1974). The lens densities were
	optimized to transform cone fundamentals based on the CIE 1964
	10-deg color matching functions to cone photopigment curves.
	WIthout this small correction, the photopigment curves have a
	discontinuity in the violet. The discontinuity is dependent on
	the CMFs-- not on the choice of cone fundamentals (see Stockman,
	MacLeod & Johnson, 1993). 
	
	390-460 nm from Stockman, MacLeod & Johnson, 1993 (see Table 7) 
	465-730 nm from Wyszecki & Stiles, 1982 (see Table 1(2.4.6)). 
	
	These values are for a small pupil. For an open pupil, divide density
	by 1.16. 

	References:

	Stockman, A., MacLeod, D. I. A., & Johnson, N. E. (1993).
	Spectral sensitivities of human cones. Journal of the Optical
	Society of America A, 10, 2491-2521.

	van Norren, D. & Vos J.J. (1974) Spectral transmission of the
	human ocular media Vision Research, 14, 1237-1244.

	Wyszecki, G., & Stiles, W. S. (1967). Color Science New York: Wiley.

	Wyszecki, G., & Stiles, W. S. (1982). Color Science (2nd ed.).  
	New York: Wiley.

	---------------------------------------------------------------
	The above taken from Andrew Stockman's home page (Feb, 1998):
		http://www-cvrl.ucsd.edu
	*/

	/* Extrapolated in the range 390 -> 380 nm (RGS)         */

	-3.368, -2.862,   -2.4, -1.985,  -1.62, -1.308,  -1.05,  -0.84, 
	-0.675, -0.557, -0.468, -0.393, -0.335,  -0.29,  -0.26,  -0.24, 
	-0.225, -0.215, -0.203, -0.191,  -0.18, -0.168, -0.162, -0.151, 
	-0.145, -0.133, -0.128, -0.122, -0.116,  -0.11, -0.104, -0.099, 
	-0.093, -0.087, -0.081, -0.075,  -0.07, -0.064, -0.058, -0.052, 
	-0.046, -0.041, -0.036, -0.031, -0.028, -0.024, -0.021, -0.017, 
	-0.014, -0.012, -0.009, -0.007, -0.005, -0.003, -0.002, -0.001, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0, 

	/* Macular and lens filters in series: */
	/*  Check "odconv.n" for the calculations */

	-3.393, -2.897, -2.445, -2.045, -1.705, -1.428,  -1.21, -1.065, 
	-0.975, -0.902, -0.833, -0.773, -0.735, -0.715,  -0.72,  -0.73, 
	-0.72, -0.685, -0.648, -0.601, -0.595, -0.588, -0.572, -0.511, 
	-0.42, -0.328, -0.258, -0.207, -0.166, -0.135, -0.114, -0.099, 
	-0.093, -0.087, -0.081, -0.075,  -0.07, -0.064, -0.058, -0.052, 
	-0.046, -0.041, -0.036, -0.031, -0.028, -0.024, -0.021, -0.017, 
	-0.014, -0.012, -0.009, -0.007, -0.005, -0.003, -0.002, -0.001, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0,     0,     0,     0, 
	0,     0,     0,     0,     0

	};

/*------------------------------------*/

#define log10(arg) (log(arg)/LN10) 
#define exp10(arg) exp((arg)*LN10)

/*------------------------------------*/

int main(void)

/* This program prints 6 tables:

   1) the log of specific optical density for rods and the 3 cone types.
   2) the spectral amplitude of 4 light sources.
   3) the spectral transmittance of 2 filters (100% and macular pigment).
   4) the calibration sensitivity of 2 standard luminosity filters (scot, phot)
	with the 4 light sources.
   5) the average optical sensitivity of rods and cones with 4 light sources.
   6) the standard path lengths for rod and cones.

   All spectral sensitivities are calculated from 380 to 800 nm wavelength.

   The log10 of optical density is produced for output because
   it gives a better interpolation than the optical density.

   The original pigment sensitivity data is in the form of
   log10 sensitivity vs. wavelength.
   To get from sensitivity to optical density, use the following
   derivation:

   light absorbed = absolute sensitivity;
   light absorbed = incident light - transmitted light
   light absorbed, norm = (incident - transmitted) / incident 
   transmittance = transmitted light / incident light
   transmittance = 1 - light absorbed, norm
   optical density = - log10 ( transmittance )
   specific density * path length = optical density
   absolute sensitivity = rel sens * (1 - exp10 ( - recep axial density ))
   log10 (abs sens) = log10 (rel sens) + log10 (1-exp10(-recep axial dens))
 
   optical density = - log10 ( 1 - sensitivity )

   optical density = - log10 ( 1 - exp10 ( log10 ( sensitivity ) ) )

   specific density = -log10 ( 1 - exp10 ( log10(sensitivity) ) ) / path length

*/
   
{
   int i,j,w;
   double tot,val,specod;

				/* calculate pigment sens curves */
  calcwav(0);
  calcwav(1);
  calcwav(2);
  calcwav(3);

  calcwav(7);			/* photoreceptors without transduct cascade */
  calcwav(8);
  calcwav(9);

  calcwav(10);
  calcwav(11);
  calcwav(12);
  calcwav(13);
  calcwav(14);
  calcwav(15);

  calcwav(19);
  calcwav(20);
  calcwav(21);

  calcwav(22);

#ifdef TURTLECONES
  copywavt();			/* copy turtle cone data */
#endif

#ifdef GOLDFISHCONES
  copywavg();			/* copy goldfish cone data */
#endif

				/* calculate the calibrations */
				/*  for phot and scot luminosities */
				/*  with different light sources */
  for (i=0; i<LUMINS; i++) {	
    for (j=0; j<LIGHTS; j++) {	
      tot=0.0;
      for (w=0; w<PIGMSIZ; w++) {
        tot += luminos[i][w] * lights[j][w];	/* sum calib sensitivities */
      }
      lumlight[i][j] = tot/PIGMSIZ;		/* normalize sum */
    }
  }

				/* calculate the average sensitivity */
				/*  for standard light source */
  for (i=0; i<TOTREC; i++) {	
    for (j=0; j<LIGHTS; j++) {	
      tot = 0.0;
      for (w=0; w<PIGMSIZ; w++) {
        tot += exp10(pigmarr[i][w]) * lights[j][w];  /* sum sensitivities */
      }
      //pigmavg[i][j] = tot/PIGMSIZ * exp10(pigmconst[i]); /* normalize sum */
      pigmavg[i][j] = tot/PIGMSIZ * (1-exp10(-peakspecod[i]*pigmlen[i])); /* normalize sum */
    }
  }

				/* print the pigment density table */

  printf ("\n/* Log10 of pigment specific optical density spectrum: */\n\n");
  printf ("double specdens[TOTREC][PIGMSIZ] = {\n    ");
  p=0;
  for (i=0; i<TOTREC; i++) {	
    p=0;
    for (w=0; w<PIGMSIZ; w++) {
			/* abs sens of 1 um of pigm log10(1-exp10(-peakspecod)) */
			/* abs sens of 1 um of pigm log(1-exp(-peakspecod*LN10))/LN10 */
	val = pigmarr[i][w] + log10(1-exp10(-peakspecod[i]));	/* log10(abs sensitivity) */
        specod = -log10(1.0-exp10(val));
        // specod = -log10(1.0-exp10(val)) / pigmlen[i]; /* this for abs OD */
	wavout(log10(specod));
    }
    printf ("\n    ");
  }
  printf ("};\n");
						/* print the light table */
  printf ("\n/* Light source intensity spectrum: */\n\n");
  printf ("double lights[LIGHTS][PIGMSIZ] = {\n    ");
  p=0;
  for (i=0; i<LIGHTS; i++) {	
    p=0;
    for (w=0; w<PIGMSIZ; w++) {
	wavout(lights[i][w]);
    }
    printf ("\n    ");
  }
  printf ("};\n");

						/* print the filter table */
  printf ("\n/* Filter optical density spectrum: */\n\n");
  printf ("double filts[FILTS][PIGMSIZ] = {\n    ");
  p=0;
  for (i=0; i<FILTS; i++) {	
    p=0;
    for (w=0; w<PIGMSIZ; w++) {
	wavout(filts[i][w]);
    }
    printf ("\n    ");
  }
  printf ("};\n");

				/* print the luminosity calibration table */
  printf ("\n\n");

  printf ("double lumlight[LUMINS][LIGHTS] = {\n    ");
  for (i=0; i<LUMINS; i++) {	
    for (j=0; j<LIGHTS; j++) {
	val = lumlight[i][j];
        printf ("%8.5g,", val);
    }
    printf ("\n    ");
  }
  printf ("};\n");

				/* print the light source table */
  printf ("\n\n");

  printf ("double lightsens[TOTREC][LIGHTS] = {\n    ");
  for (i=0; i<TOTREC; i++) {	
    for (j=0; j<LIGHTS; j++) {
	val = pigmavg[i][j];
        printf ("%8.5g,", val);
    }
    printf ("\n    ");
  }
  printf ("};\n");


  printf ("double pigmlen[TOTREC] = {\n    ");
  for (i=0; i<TOTREC; i++) {	
	val = pigmlen[i];
        printf ("%8.5g,", val);
  }
  printf ("};\n");
  exit(0);
}

/*------------------------------------*/

double logwav(double w)
{
  return log(1/(w*.001 * LAMM/lamm)) / LN10;
}

/*------------------------------------*/

static double an[7] = {-5.2734, -87.403, 1228.4,
			  -3346.3, -5070.3, 30881.0, -31607.0};
void calcwav(int pigm)
           

/* Procedure to calculate spectral sensitivity
   curves for rod and three cone types.
   Spectral response curve taken from Baylor et al. 1987.
   Beyond 800 nm, the curve is linear with wavenumber,
   (1/w) (see Lewis, 1955 J.Physiol 130:45), and the curve is
   extended in this region using linear extrapolation
   as a funtion of 1/w, using the data for S= -3 to -5,
   using Baylor's slope measurements.  

   Wavelength increment is 5 nm, which gets quadratically
   interpolated when used in the simulator.

   The result of this routine is the log10 of sensitivity.
*/

{
    int i,n;
    double m, b, dy, dw, wold;
    double w, sum, x, xx, y, y0;

  switch (pigm) {

   case 0: 
          lamm = 500.0;			/* primate rod */
          b = 15000.0;
          m = 0;
          break;
	
   case 1: 
         lamm = 561.0;			/* primate red */
          b = 17200.0;
          m = 0;
          break;
		
   case 2: 
          lamm = 531.0;			/* primate green */
          b = 15900.0;
          m = 0;
          break;
		
   case 3: 
          lamm = 430.0;			/* primate blue */
          b = 12700.0;
          m = 0;
          break;
		
   case 7: 
          lamm = 561.0;			/* red without transduction */
          b = 17200.0;
          m = 0;
          break;
		
   case 8: 
          lamm = 531.0;			/* green without transduction */
          b = 15900.0;
          m = 0;
          break;

   case 9: 
          lamm = 430.0;			/* blue without transduction */
          b = 12700.0;
          m = 0;
          break;

   case 10: 
          lamm = 500.0;			/* rabbit rod */
          b = 15000.0;
          m = 0;
          break;

   case 11: 
          lamm = 527.0;			/* rabbit M cone */
          b = 15900.0;
          m = 0;
          break;

   case 12: 
          lamm = 430.0;			/* rabbit S cone */
          b = 12700.0;
          m = 0;
          break;

   case 13: 
          lamm = 500.0;			/* guinea pig rod */
          b = 15000.0;
          m = 0;
          break;

   case 14: 
          lamm = 527.0;			/* guinea pig M cone */
          b = 15900.0;
          m = 0;
          break;

   case 15: 
          lamm = 430.0;			/* guinea pig S cone */
          b = 12700.0;
          m = 0;
          break;

   case 19: 
         lamm = 561.0;			/* primate red, van hateren cone */
          b = 17200.0;
          m = 0;
          break;
		
   case 20: 
          lamm = 531.0;			/* primate green, van hateren cone */
          b = 15900.0;
          m = 0;
          break;
		
   case 21: 
          lamm = 430.0;			/* primate blue, van hateren cone */
          b = 12700.0;
          m = 0;
          break;

   case 22: 
          lamm = 500.0;			/* mouse rod, for invergo (2014) */
          b = 15000.0;
          m = 0;
          break;
  }	

  sum = 0;
  for (w=MINWAV; (sum > -5.5) && (w<=MAXWAV); w+= WAVINC) { 
    sum = 0;
    x = logwav(w); 
    for (n=0; n<=6; n++) {
       xx = 1.0;
       for (i=0; i<n; i++) {
	 xx *= x;	
       }
       sum += xx * an[n];
    } 
/*    printf ("%5g  %8.3g\n", w, sum);   	   /* log sensitivity */
/*    printf ("%5g  %8.3g\n", -1/w, sum);   	   /* log sensitivity */
/*    printf ("%5g  %8.3g\n", w, exp(LN10*sum)*.665);  /* sensitivity */
/*    printf ("%5g  %8.3g\n", w, -log(1.0-exp((sum-0.177)*LN10))/LN10 );  */
/*printf("%5g  %8.3g\n",-1/w,log(-log(1.0-exp((sum-0.177)*LN10))/LN10)/LN10); */
/*    wavout(sum); */
      arrfill(sum,pigm);
  }
  y0 = -27.75;
  y = sum;
  for (wold=w-5; w<=MAXWAV; wold=w, w+= WAVINC) {
    dw = 1/w - 1/wold;
    dy = m * 1/(1/lamm - 1/w) + b;
    y += dy * dw; 
/*    y = y0 + -(1/(w*.001 * LAMM/lamm)) * -17.82946;   */
/*    y = y0 + -(1/w) * -17829.46; */
/*    printf ("%5g  %8.3g\n", w, y);  /* */
/*    printf ("%5g  %8.3g\n", -1/w, y);  /* */
/*    printf ("%5g  %8.3g\n", w, exp(LN10*y)*.665);  /* */
/*    printf ("%5g  %8.3g\n", w, -log(1.0-exp((y-0.177)*LN10))/LN10 );  */
/*  printf("%5g  %8.3g\n",-1/w,log(-log(1.0-exp((y-0.177)*LN10))/LN10)/LN10);*/
/*    wavout(y); */
    arrfill(y,pigm);
  } 
}

/*------------------------------------*/


#ifdef TURTLECONES

/* Turtle sensitivity data taken from:
 
      Perlman, Itzhaki, Malik, and Alpern: 
      "The action spectra of cone photoreceptors in the turtle retina",
      Visual Neuroscience,(1994),Vol.11, page 247, Figs.4, 10.
 
    Data originally digitized by Amermuller et al, Oldenburg, DE.
    Peak shape modified for smoothness, RGS
    Extrapolated from 400 through 380, RGS
 
    Each curve has been extrapolated on the long wave end by
    the method of Lewis, 1955, (J.Physiol 130:45) in which the log(sens)
    of a visual pigment is supposed to be linear in the deep red when 
    plotted against wavenumbmer (1/wavelength).  The script "turtlesens.n"
    implements this method and was used to produce "turtle.h"
    which contains the log(sensitivity) functions extrapolated to
    800 nm.  In turn, "turtle.h" is included in "wave.cc" which
    generates log (specific OD) functions required in "wave.h"
    for "ncstim.cc".
*/

#include "turtle.h"

void copywavt(void) 

/* copy turtle log(sens) functions into pigmarr[] */

{
   int pigm,w;

  for (pigm=4; pigm<7; pigm++)
    for (w=0; w<PIGMSIZ; w++) {
      arrfill(turtlepigm[pigm-4][w],pigm);
    }
}
#endif

/*------------------------------------*/


#ifdef GOLDFISHCONES

/* Goldfish sensitivity data taken from:

    Carp fundamentals in van Dijk & Spekreijse (1984)
    fitted to equations in Govardovskii et al. (2000)
 
    Data fitted by Maarten Kamermans, in spreadsheet:
    goldfish_spectra.xlsx

    The script "goldfish.n" implements this method and was 
    used to produce "goldfish.h" which contains the log(sensitivity) 
    functions. In turn, "goldfish.h" is included in "wave.cc" which
    generates log (specific OD) functions required in "wave.h"
    for "ncstim.cc".
*/

#include "goldfish.h"

void copywavg(void) 

/* copy goldfish log(sens) functions into pigmarr[] */

{
   int pigm,w;

  for (pigm=16; pigm<19; pigm++)	// pigments for 16, 17, 18 (goldfish)
    for (w=0; w<PIGMSIZ; w++) {
      arrfill(gfpigm[pigm-16][w],pigm);
    }

  // for (pigm=19; pigm<22; pigm++)	// pigments for van hateren cones
  //  for (w=0; w<PIGMSIZ; w++) {
  //    arrfill(gfpigm[pigm-19][w],pigm);
  //  }
}
#endif

/*------------------------------------*/

void wavout (double val)
{
  printf ("%8.5g,", val);
  if (++p >= 8) {
     p = 0;
     printf ("\n    ");
  }    
}

/*------------------------------------*/

void arrfill (double val, int pigm)
             
           

/* fill pigment array */

{
  static int i,pg= -1;

  if (pg==pigm) {
    pigmarr[pigm][i++] = val;
  }
  else {
    pg = pigm;
    i = 0;
    pigmarr[pigm][i++] = val;
  } 
}
