/* Header file ncsub.h for ncsub.c */

#define NCSUBH 1

#ifndef MAXSHORT

#define _VALUES_H 1

#define BITSPERBYTE 8
#define BITS(type)  (BITSPERBYTE * (int)sizeof(type))
  /* Rename the above to __BITS so it doesn't intrude on users' namespace? */

#define CHARBITS    BITS(char)
#define SHORTBITS   BITS(short)
#define INTBITS     BITS(int)
#define LONGBITS    BITS(long)
#define PTRBITS     BITS(char*)
#define DOUBLEBITS  BITS(double)
#define FLOATBITS   BITS(float)

#define MINSHORT    ((short)(1 << (SHORTBITS - 1)))
#define MININT      (1 << (INTBITS - 1))
#define MINLONG     (1L << (LONGBITS - 1))

#define MAXSHORT    ((short)~MINSHORT)
#define MAXINT      (~MININT)
#define MAXLONG     (~MINLONG)

#define HIBITS	MINSHORT
#define HIBITL	MINLONG

#ifndef MAXDOUBLE
#define MAXDOUBLE   1.79769313486231570e+308
#define MINDOUBLE   4.94065645841246544e-324
#endif

#ifndef MAXFLOAT
#define MAXFLOAT    3.40282347e+38F
#define MINFLOAT    1.40129846e-45F
#endif

#elif defined(convex)
#define MAXDOUBLE   8.9884656743115785e+306
#define MAXFLOAT    ((float) 1.70141173e+38)
#define MINDOUBLE   5.5626846462680035e-308
#define MINFLOAT    ((float) 2.93873588e-39)
#endif

#ifndef MAXSHORT
#define MAXSHORT 32767			/* maximum value for short integer */
#endif

#ifndef MINSHORT
#define MINSHORT -32768			/* maximum value for short integer */
#endif

#define MAXNODE MAXSHORT		/* max value for node */
#define MINNODE MINSHORT		/* min value for node */

typedef int nodeint;			/* type of int used for node numbers */

#define MINVES 1e-6			/* minimum ves. release rate/step */
#define MAXVINT 10000			/* max interval for vesicle release */
#define MAXNODIM 4			/* max number of node dimensions */
#define RATESIZ 500			/* size of rate const lookup tables */
#define TRSCAL 1000.0                   /* conversion from 1000 to 1 (ligand) */
#define TRMAX  10000.0                  /* max value for trans in synapse */
#define TRMIN -10000.0                  /* min value for trans in synapse */
#define MAXPHOT 1e10                    /* maximum photons per .1 ms timestep */
#define BATTCOND 10.                    /* battery conductance (siemens) */
#define NUMREC   24                     /* number of pigments for receptors */
#define TOTREC   (NUMREC)               /* number of pigments for receptors */
#define LIGHTS   4                      /* number of light sources */
#define LUMINS   2                      /* number of standard luminosities */
#define FILTS    4                      /* number of filters */
#define MINWAV   380			/* min wavel for pigm table */
#define MAXWAV   800			/* max wavel for pigm table */
#define WAVINC   5			/* wavel inc for pigm act spec table */
#define PIGMSIZ ((MAXWAV-MINWAV)/WAVINC+1) /* size of pigment act spec. table */
#define SMALLCAP 1e-33			/* small capacitor for blank comp */
#define LARGERES 1e35			/* large capacitor for blank comp */
#define NAEXCHR 3.0			/* Na-Ca exchange rate */
#define CAZ   2.0			/* Ca++ valence */
#define CAZI  0.25			/* factor to equate Ca and K conc */

#define CONEVH1 19			/* new van hateren pigment number */
#define RODI 22				/* new invergo mouse rod pigment number */
#define ChR2 23				/* new channel rhodopsin pigment number */
#define NPHOS 7				/* number of phosphorylations in invergo mouse rod */
// #define NPHOS 1				/* number of phosphorylations in invergo mouse rod */

#define PORIG 0				/* original cone */
#define PVH  1				/* van hateren cone  */
#define PINVG 2				/* invergo rod */


#define MB 1048576.0			/* 1 MegaByte */
#define CASHD	0.1			/* calcium shell depth in um */
#define CACOMPRAD 5.0			/* calcium comp radius in um */

#define CNCOLOR   -1			/* display color from cell number */
#define VCOLOR    32768			/* display color from voltage */
#define LCOLOR    32767			/* display color from light inten */
#define RCOLOR    32766 		/* display color from cell region */
#define CACOLOR   32765			/* display color from light inten */
#define NAACOLOR  32764 		/* display color from Na activated state */
#define NAICOLOR  32763 		/* display color from Na inactivated state */
#define KACOLOR   32762 		/* display color from K state */
#define SGCOLOR   32761 		/* display color from syn cond */
#define SRCOLOR   32760 		/* display color from syn rate */

#define CASSTART 1.6e-6                 /* initialize IP3 [calcium] store to Cas2Start = 1.6e-6 M */
#define CAS2START 1.6e-6                /* initialize IP3 [calcium] store to Cas2Start = 1.6e-6 M */
#define IP3ISTART 0.4e-6                /* initialize [IP3] intracellular to Ip3iStart = 0.4e-6 M */
#define MIP3START 0                     /* IP3 store uptake m gate fraction initialize to 0 (default)   */
#define HIP3START 0.5                   /* IP3 store uptake h gate fraction initialize to 0.5 (default) */

#define DCAO 0.00115			/* default [cao] for ca voffset calc */
#define DMGO 0.0012			/* default [mgo] for ca voffset calc */

#define PNA   0				/* index for ions for channels */
#define PK    1
#define PCA   2
#define PCL   3
#define PCO3  4
#define PCS   5
#define PLI   6
#define PNH4  7
#define PTMA  8
#define NIONS 9				/* number of ions considered */

struct iontab {
	double ioni[NIONS];		/* internal ion concentration */
        double iono[NIONS];		/* external ion concentration */
        double ionp[NIONS];		/* relative ion permeability */
	};

#define GJOFF  0.015			/* gap junction offset voltage */

#define NRECA0 40                       /* filter before nt release */
#define NRECA1 41                       /* neurotransmitter in filter 1 */
#define NRECA2 42                       /* neurotransmitter in filter 2 */
#define NRECA3 43
#define NRECA4 44
#define NRECA8 45			/* transrate (actual vesicle rate */
#define NRECA9 46			/* transrate (actual vesicle rate */

#define NRECH0 47                       /* filter before nt release */
#define NRECH1 48                       /* neurotransmitter in filter 1 */
#define NRECH2 49                       /* neurotransmitter in filter 2 */
#define NRECH3 50                       /* neurotransmitter in filter 2 */
#define NRECH4 51                       /* neurotransmitter in filter 2 */

#define NRECB0 52                       /* filter after nt release */
#define NRECB1 53                       /* neurotransmitter in filter 1 */
#define NRECB2 54                       /* neurotransmitter in filter 2 */
#define NRECB3 55
#define NRECB4 56

#define NRECC0 57                       /* filter after saturation */
#define NRECC1 58                       /* neurotransmitter in filter 1 */
#define NRECC2 59                       /* neurotransmitter in filter 2 */
#define NRECC3 60
#define NRECC4 61
#define NRECC9 62

#define NRECG0  63
#define NRECG1  64
#define NRECG2  65
#define NRECG3  66

#define NRECG   67
#define NRECGI  68
#define NRECGV  69
#define NRECGC  70
#define NRECGM  71
#define NRECGH  72
#define NRECCA  73

#define CABLFRAC 10000			/* mult. for distance along cable */

					/* comp->miscfl values: */
#define VEXT 0x000f                     /* voltage clamp for compartment */
#define IEXT 0x00f0                     /* current clamp for compartment */
#define VBAT 0x0100                     /* voltage clamp is battery */
#define VBUF 0x0200                     /* voltage buffer (remote vclamp) */
#define CCA  0x0400			/* calcium integration in compartment */
