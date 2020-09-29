/* Header file nconst.h for program "nc". */
/* Contains basic mathematical, physical and biological constants */

#define HMICSEC	0.0001			/* 100 microseconds in seconds */
#define MSEC	0.001			/* 1 msec in seconds */
#define HMICMS	0.1			/* 100 microseconds in msec */
#define MVOLT	1000.0			/* V to mV */
#define MSSEC	1000.0			/* msec to sec */
#define VTOMV   1000.0                  /* conversion from volts to mvolts */
#define DMUM    1e-5                    /* number of dm / um */
#define DMCM2   1e-2                    /* number of dm2 / cm2 */
#define UM2CM2  1e8                     /* number of um2 / cm2 */

#define BASETHH  6.3			/* temp for calc HH rates */
#define BASETC   22.0			/* temp for calc channel rates */
#define BASETSYN 22.0			/* temperature for calc synap rates */
#define BASETCA  22.0			/* temp for calc Ca pump/buf rates */

// #define RNDUP 1e-8
#define RNDUP 0

/* #define HHRATE (exp(log(3.0)*(BASETC-BASETHH)/10.0)) /* HH rate for Q10=3*/
#define    HHRATE 5.6115181
 

#define BASETDC 22.0			/* temperature for calc diffusion */
#define QDC   1.3			/* Q10 for for Ca diffusion constant */

#ifndef MPI
#define MPI   3.14159265358979323846264
#endif
#define LN10  2.30258509299404568402
#ifndef M_E
#define M_E   2.7182818284590452354	/* e */
#endif
// #define R     8.31441                   /* Gas constant J/mol/degK */
// #define Fd    96484.56                  /* Faraday's constant Coul/mol   */
#define R     8.3144621                 /* Gas constant J/mol/degK */
#define Fd    96485.3399		/* Faraday's constant Coul/mol    */
#define F2    5.1821762984668222083e-6  /* 1/2Fd converts amps Ca to moles/s */
#define RF    8.6173476875471067323e-5  /* R/F */
#define R2F   4.3086738437735533662e-5  /* R/2F */
#define FR    11604.498695638054414     /* F/R */
#define F2R   23208.997391276108829     /* 2F/R */
#define FFR   1.1196660e+09             /* (F*F/R) */
#define Bk    1.3823                    /* k Boltzman's const W sec/deg K */
#define KELVIN  273.16                  /* add to deg celsius for kelvin */
#define AF 2635.9917                    /* = .512 * (KELVIN+25)^1.5 */

#define EG    0.38                      /* Voltage barrier loc for caexch */
                                        /*  deSchutter & Smolen 1998 */
                                        /*  from Hille, 1992 */
                                        /*  see Gabbiani et al 1994) */ 
#define DCa   2e-6                      /* Diff const for ca (cm2/sec) */
                                        /*  slower in cytoplasm than in water */
                                        /*  see de Schutter&Smolen "Methods" */
                                        /*   (from Albritton et al, 1992) */

