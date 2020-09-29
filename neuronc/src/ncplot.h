
#define PLOTNODSIZ 200			/* size of plotnod: max num of plots */
#define PNAMSIZ 30			/* length of name for frame, plname */

#define VREC 1                          /* record voltage for graphing node */
#define IREC 2                          /* record current for graphing node */
#define MREC 3                          /* record membrane current at node */
#define LREC 4                          /* record light for graphing node */
#define GREC 5                          /* number for graphing */
#define SREC 6                          /* plot symbol variable */
#define FREC 7                          /* plot function */
#define WREC 8                          /* record Vm (comp voltage - vext) for graphing */

#define CREC_GLU   9                     /* nt conc, cA/GMP for graphing */
#define CREC_AMPA 10			/* Keep these in order */
#define CREC_KAINATE 11			/* here and in nc.y */	
#define CREC_NMDA 12
#define CREC_CNQX 13
#define CREC_GABA 14
#define CREC_BIC  15
#define CREC_PTX  16
#define CREC_GLY  17
#define CREC_STRY 18
#define CREC_CAMP 19                     
#define CREC_CGMP 20	
#define CREC_CHRC 21	
#define CREC_PH   22	
#define CREC_ATP  23	

#define CREC_CACONC 24			/* Ca concentration outside */
#define CREC_CAV    25
#define CREC_CAIT   26
#define CREC_CAIP   27
#define CREC_CAIE   28
#define CREC_CAIPE  29
#define CREC_CABUF  30 	
#define CREC_CABUFB 31 	
#define CREC_CAS    32 	
#define CREC_CICR   33 	/* must change NREC in ncsub.h if this goes past 39 */

#define PLINES  'P'                     /* plots have only lines */
#define LINES   'L'                     /* plots have chars and lines */
#define NOLINES 'N'                     /* plots have chars but no lines */

struct lpfilt;

struct plotfr {                 	/* node numbers, vmin, vmax voltages */
        int cnod1;
        int cnod2;
        int cnod3;
        int cnod4;
        double pxmax;
        double pxmin;
        double pymax;
        double pymin;
        double xrange;
        double yrange;
        double oldx;
        double oldy;
        double olddx;
        double olddy;
        double yval;
        double plotsiz;		/* size of plot on screen */
        double plotval;		/* value to be passed into plot func */
        double plotval2;	/* second value to be passed into plot func */
        double plotval3;	/* third value to be passed into plot func */
        double plotval4;	/* fourth value to be passed into plot func */
        double plotval5;	/* fifth value to be passed into plot func */
        double plotval6;	/* sixth value to be passed into plot func */
        double plotval7;	/* seventh value to be passed into plot func */
        char *spnt;             /* pointer to symbol for plotting */
        double *var;		/* pointer to var for plot */
        double *arr;		/* array to plot into */
        int arrindex;		/* array index */
        int maxindex;		/* max array index */
        short int pmod;
        short int pmod2;
        short int pval;
        short int ppen;
        short int plotn;	/* plot on screen, can be assigned */
        int *func;		/* function for plot */
        int *funcc;		/* function for plot, C++ */
        int *vpen;		/* function for variable pen */
        int *vpenc;		/* function for variable pen, C++ */
	lpfilt *filt;		/* digital filter for plotting */
        short automax;		/* maxmin from automatic scaling */
        double csize;		/* size of char for plot */
        char charmode;		/* mode of char plot */
        char charfl;		/* char to be plotted */
        char plframe[PNAMSIZ];	/* frame to draw graph in */
        char plname[PNAMSIZ];	/* name of graph trace */
        short int newframe;	/* =1 -> first plot in new frame */
        };

// defined in makfilter.cc:

#define opt_be 0x00001  /* -Be          Bessel characteristic          */
#define opt_bu 0x00002  /* -Bu          Butterworth characteristic     */
#define opt_ch 0x00004  /* -Ch          Chebyshev characteristic       */
#define opt_re 0x00008  /* -Re          Resonator                      */
#define opt_pi 0x00010  /* -Pi          proportional-integral          */

#define opt_lp 0x00020  /* -Lp          lowpass                        */
#define opt_hp 0x00040  /* -Hp          highpass                       */
#define opt_bp 0x00080  /* -Bp          bandpass                       */
#define opt_bs 0x00100  /* -Bs          bandstop                       */
#define opt_ap 0x00200  /* -Ap          allpass                        */

#define opt_a  0x00400  /* -a           alpha value                    */
#define opt_l  0x00800  /* -l           just list filter parameters    */
#define opt_o  0x01000  /* -o           order of filter                */
#define opt_p  0x02000  /* -p           specified poles only           */
#define opt_w  0x04000  /* -w           don't pre-warp                 */
#define opt_z  0x08000  /* -z           use matched z-transform        */
#define opt_Z  0x10000  /* -Z           additional zero                */

