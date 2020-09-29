/* Header file ncomp.h for program "nc". */
/* Contains low-level data structures for runtime */

#include <stdint.h>

#define NCOMPH 1

#define NUMSTATE 16                     /* max number of states in channel */
#define NUMTRANS 6                      /* max transitions to other sts */
#define NCHANCONST 100                  /* max number of channel types */
//#define RNDSIZ  (3*sizeof(long))      /* size of random number state table for rnd_taus */
#define RNDSIZ  (4*sizeof(long))      /* size of random number state table for rnd_taus113 */
//#define RNDSIZ  (8*sizeof(uint32_t))    /* size of random number state table for rnd_tinymt */
//#define RNDSIZ  (625*sizeof(long))    /* size of random number state table for rnd_mt */
#define CHANRATE 1000.0			/* default chan rate = 1 msec tau */
#define NNTR 20				/* number of neurotransmitters */
#define NUMFILT 10			/* number of filters in lpfilt */
#define MAXQ 10				/* number of qrate values in chanparm */
#define MAXSTIMCHAN 200			/* max number of stimulus channels */

struct conn; 
struct dbuf;
struct chan;
struct hhchan; 
struct kcachan;
struct stconc;
struct sschan;
struct chanstate;
struct chantype;
struct lpfilt;
struct synap;
struct load;
struct photrec;
struct recpar;
struct recpari;
struct recparc;
struct recstim;
struct cacomp;
struct conlst;
struct dyadlst;
struct node;
struct ntcomp;

struct id {
	short int ctype;		/* type of object */
	short int stype;		/*  subtype */
	};

struct comp: id {			/* structure for compartment */
        short int miscfl;               /* misc flags: VEXT IEXT VBAT CA ... */
        int num;			/* comp number */
        int domain;			/* domain number (resistively coupled comps) */
        double v;                       /* voltage */
        double vext;                    /* voltage on outside of membrane */
        double oldvext;                 /* voltage on outside of membrane at prev time step */
        double vph;                     /* voltage shift for Ca gating from ext pH */
        double oldv;                    /* voltage at previous time step */
        double vest;                    /* volt estimate */
        double vesto;                   /* orig value of voltage est */
        double implf;                   /* implicit factor */
        double nodcur;                  /* node current varying w/vest */
        double nodc;                    /* precalc current not varying w/vest */
        double tcond;                   /* total static comp conductance */ 
        double tcondt;                  /* total cond not varying w/vest */ 
        double vrev;                    /* membrane reversal potential */
        double k;                       /* calculation constant */
        double tcondm;                   /* total cond through membrane */
        double verr;                     /* abs value of voltage error */
        double rm;                       /* membrane conductance in comp */
        double cap;                      /* capacitance in comp */
        double jnoise;                   /* membrane conductance in comp */
        double extv;                     /* voltage clamp voltage */
        double extv_old;                 /* previous voltage clamp voltage */
        double extvi;                    /* voltage clamp current */
        double exti;                     /* current clamp value */
        double exti_old;                 /* previous current clamp value */
        double cplam;                    /* local complam for condense */
        short int g;               	/* verr greater than zero */
        short int t;               	/* total verr calculated */
        short int virt;               	/* virt compartment */
	double relax;			/* relax for this comp only */
	char *jstate;			/* state array for Johnson noise */
	cacomp *capnt;			/* pointer to calcium compartment */
	ntcomp *ntpnt;			/* list of neurotr compartments */
        conlst *clst;           	/* list of conns to other comps */
        conlst *nodlst;         	/* (possible) nodes for this comp */
	comp *pvext;			/* pointer to outside comp (vext) */
        comp *hnext;             	/* pointer to next comp for hash */
        comp *next;             	/* pointer to next compartment */
        comp *last;             	/* pointer to last compartment */
        };

struct conn: id {                   	/* structure for comp connection */
	int num;			/* number of conn (only for synap, chan) */
        conn *next;             	/* pointer to next connection (gj) */
        conn *last;             	/* pointer to last connection (gj) */
        double conduct;                 /* conductance of connection */
        comp *comp1;                    /* first compartment connection */
        comp *comp2;                    /* second compartment connection */
        comp *compe;                    /* external compartment connection */
        comp *compp;                    /* external pH compartment connection */
	int num1;			/* comp # of first comp. */
	int num2;			/* comp # of second comp. */
        };

struct gj: conn {                   	/* structure for gap junction */
        double maxcond;                  /* maximum conductance */
        double vgain;               	/* mV/e-fold change (V sens) */
        double voff;               	/* volt. offset for expont. change */
        double taun;			/* activation tau */
        double rvgain;               	/* rev. mV/e-fold change (V sens) */
        double rvoff;               	/* reverse offset for expon. change */
        double rtaun;			/* activation tau */
        double n;			/* conductance fraction */
        double gnv;			/* ratio of min/max conductance */
	short int rev;			/* reverse parms set */
	short int rect;			/* gj is pure rectifier */
	short int modtyp;		/* CycA compartment */
	short int sign;			/* sign of modulation: open or close */
	node *nodec;			/* node to check for cycg modulation */
        };

struct pnx: gj {                   	/* structure for pannexin channel */
	double atp;			/* state variables for external ATP compartment */
	double atpx;			/* atp bound to atpase */
	double atpo;			/* atp minus h */
	double adp;
	double adpx;			/* adp bound to adpase */
	double adpo;			/* adp minus h */
	double amp;
	double ampx;
	double ad;
	double h;			/* protons */
	double pi;			/* phosphate */
	double pio;
	double atpase;
        double atp_decr;		/* decrement for atp per syn timestep (0-1) */
        double atp_gain;		/* atp gain from pnx current */
	double (*pnxfunc)(pnx*,double);	/* function to process pnx signal */
        };

struct ntcomp: id {			/* neurotrans compartment */
	short int n;			/* number of inputs */
	double val;			/* concentration */
	double tau;			/* time constant of decay */
	ntcomp *ntpnt;			/* pointer to next ntcomp */
	};

struct dbuf: conn {                   	/* structure for delayed buffer */
        double v;               	/* voltage output */
        double offset;               	/* voltage offset (volts) */
        double gain;               	/* voltage gain */
        double delay;               	/* size of delay buffer (in stiminc=1e-4sec)*/
	double *delbuf;			/* circular buffer pointer */
	double *delpnt;			/* input/output buffer pointer */
	lpfilt *filt;			/* low pass filter after delay */
        };

struct ndbuf: conn {                   	/* structure for ntrans buffer */
        double offset;               	/* voltage offset */
        double ntoffset;               	/* nt offset */
        double gain;               	/* voltage gain */
        int ntrans;               	/* ntrans output */
        };

struct chan: conn {            		/* basis stuct for membrane channel */
	double maxcond;			/* maximum total channel conductance */
        double rextern;			/* external resistance */
        double voffsm;			/* activation offset voltage */
        double voffsh;			/* inactivation offset voltage */
        double vrev;			/* reversal potential */
        double gvrev;			/* indiv slope vrev (if setvrev) */
        double gfrac;			/* indiv slope cond (if setvrev) */
	double *gvrevtab;		/* nonlinear vrev table */
	double *gfractab;		/* slope conductance table */
	chantype *chtyp;		/* channel constants */
        char *cstate;			/* channel random state */
        double cakd;			/* Ca dissociation constant (half max)*/
        double cahc;			/* Ca hill coefficient */
        double naratio;			/* nai/nao for channel */
        double kratio;			/* ki/ko for channel */
        double nchan;			/* total number of channels */
        int cno;			/* number of channels actually open */
        double arate;			/* forward rate (1/tau) */
        double brate;			/* reverse rate */
        double crate;			/* inactivation rate */
        double drate;			/* recovery from inactivation rate */
        double erate;			/* noise flicker rate */
        double frate;			/* noise flicker rate */
	short int setvrev;		/* = 1 -> vrev set by user */
	short int nocondens;		/* = 1 -> don't condense */
	short int numstate;		/* number of states used */
	iontab *ions;			/* ion concs, perms, defined in ncsub.h */ 
	stconc *conc;			/* state concentr (variable size arr) */
        };

struct hhchan: chan {                 	/* structure for membrane channel */
        double m;                       /* m,n activation, range 0 - 1 */
        double h;                  	/* h inactivation  range 0 - 1 */
        };

struct kcachan: hhchan {             	/* structure for KCa membrane channel */
        double d1;                  	/* voltage multiplier for alpha */
        double d2;                  	/* voltage multiplier for beta */
        double k1;                  	/* Ca multiplier for alpha */
        double k2;                  	/* Ca multiplier for beta */
        short initfl;                  	/* 1=> m initialized, 0=> m not init */
        };

struct ntchan: chan {                	/* struct for sequen-state synap chan */
	double mg;			/* magnesium conc. */
	double gly;			/* glycine conc. */
        };

struct chrchan: chan {                	/* struct for channel rhodopsin chan */
	double iflux;			/* intensity from photoreceptorc. */
	double p;			/* integrated intensity */
        };

struct dstconc {                 	/* structure for state concentration */
	double dcon;			/* change in state concentr estimate */
	double icon;			/* implicit change estimate */
	int   dchan;			/* change in channels in this state */
        };

struct stconc {                 	/* structure for state concentration */
	double cval;			/* state concentr */
	double cest;			/* state concentr */
	int   nchan;			/* number of channels in this state */
        };

struct chanstate {              	/* structure for chan state constants*/
        double cond;                     /* normalized conductance of state */
	short int numtrans;		/* number of transitions from state */
	char     trans[NUMTRANS];	/* four possible state transitions */
	double (*trate[NUMTRANS])
		  (chan *cpnt);		/* function returning rate for trans */
	double  ratemul[NUMTRANS];	/* rate multiplier for transition */
	char   rateo[NUMTRANS];		/* offset voltage for transition */
	};				/* rateo=0->voffsm; rateo=1->voffsh;*/

struct chanparm {                	/* structure for channel rate params */
	short int pn;			/* parm number m=0, h=1, etc. */
	short int nfval;		/* number of alpha, beta, etc. vals */
	short int nival;		/* number of implicit values */
	short int nq;			/* number of qrate values */
	double *fval;			/* alpha, beta, etc. vals for integ. */
	double *ival;			/* implicit values for integration */
        double *chanftab;		/* rate func tables */
        double *chanitab;		/* implicit tables if required */
        double (*chancalc)(double,int); /* default rate function */
        char *funcname;			/* user-defined rate function */
        double *qrate;			/* Q10 rate calc from dq, qt */
        double *dq;			/* Q10 factor for rate calc */
        double qt;			/* base temperature for rate funcs */
        double dqca;			/* Q10 for cao offset on chan voffset */
        double qcao;			/* base cao for chan offset: 1.15 mM */
        double qcavoff;			/* voltage offset calc from dqca,qcao */
        double voff;			/* voltage offset to modify rate funcs*/
        double oldvf;			/* Variables that need to be checked */
        double oldvi;
        double oldar;
        double oldbr;
        double oldtimincf;		/*  to see if table needs to be */
        double oldtiminci;
        double oldtempcel;		/*  regenerated. */
        double oldcao;
	};

struct chantype: id {                	/* master struct for chan constants */
        short int numstate;             /* number of states used */
	short int numparm;		/* number of params, set to 2 for HH */
	short int hh;			/* = 1 -> HH channel, otherwise SS */
	short int cabnd;		/* = 1 -> needs Ca binding */
	short int cgbnd;		/* = 1 -> needs cGMP binding & sat */
        double qc;			/* Q10 factor for conductance */
        double qt;			/* base temp for conductance */
        double qcond;			/* Q10 multiplier for conductance */
	double vrev;			/* default reversal potential */
	iontab *ions;			/* ion concs, perms, defined in ncsub.h */ 
	double unitary;			/* default unitary conductance */
	double trconc;			/* default tr conc multiplier */
	chanstate *state;		/* state definitions */
	chanparm *parm;			/* rate params */ 
	double *respamp;		/* response amplitudes */
	double *gvrevtab;		/* nonlinear vrev table */
	double *gfractab;		/* slope conductance table */
	double gamma;			/* weight for state O2 conductance */
	};


struct cacomp: id {			/* structure for calcium comp/chan */
        int num;                	/* compartment number */
        comp *comp1;                    /* voltage compartment     */
        int vcnum;                      /* voltage compartment num */
	double casf0;			/* shell factor for caflux */
	double casfn;			/* shell factor = 1/(dr*dr) */
	double casfc;			/* shell factor for core */
	double vrev;			/* calcium reversal potential */
	double *cais;			/* shells for inside ca diffusion */
	double casfno;			/* shell factor outside = 1/(dr*dr) */
	double casfco;			/* shell factor for outside core */
	double *caos;			/* shells for outside ca diffusion */
	double *cab;			/* buffer conc at each Ca shell */
	double cao;			/* calcium concentration outside */
	double cai;			/* starting [Ca]i */
	double tcai;			/* threshold [Ca]i for pump */
	double nao;			/* Na concentration outside */
	double nai;			/* Na conc inside */
	double ko;			/* K concentration outside */
	double ki;			/* K conc inside */
	double cavoff;			/* voltage offset due to cao */
	double cabti;			/* Ca buffer in first shell */
	double cabt;			/* Ca buffer in shell */
	double cabf;			/* forward rate to Ca.B /M/sec */
	double cabr;			/* reverse rate to Ca.B /sec */
	double cabnd;			/* frac Ca bound (not used with cabuf)*/
	double dr;			/* shell thickness (for prcomp() */
        double area;                    /* int area in dm */
	double ica;			/* current through calcium chan */
	double vmax;			/* calcium pump Vmax */
	double pkm;			/* calcium pump half-max */
	double kexi;			/* Na-calcium exchange rate const */
	double kexo;			/* Na-calcium exchange rate const */
	short int cashell;		/* number of ca shells */
	short int caoshell;		/* number of ca shells outside */
	double ego;			/* Voltage barrier location * F/R */
	double egi;			/* 1-Voltage barrier location * F/R */
	double ipump;			/* total current from pump */
	double iexch;			/* total current from exch */
	double slope;			/* "linearized" exchanger conductance */
	double eica;			/* exhanger y zero-crossing */
	double cas;			/* CICR [Ca] store */
	double vm2;			/* CICR vmax uptake */
	double vm3;			/* CICR vmax release */
	double ncicr;			/* CICR Hill coeff Ca uptake */
	double mcicr;			/* CICR Hill coeff Ca release */
	double pcicr;			/* CICR Ca release */
	double kacicr;			/* CICR thresh Ca release */
	double kfcicr;			/* CICR thresh Ca release */
	double krcicr;			/* CICR thresh Ca release */
	double k1cicr;                  /* CICR thresh Ca release */
	double k2cicr;			/* CICR thresh for CICR Ca release */
	double c1cicr;			/* ratio of ER to cytoplasmic vol for CICR */
	double cicrflux;		/* CICR Ca flux */
	double cas2;			/* IP3 [Ca] in IP3 store */
	double ip3i;			/* IP3 intracellular [IP3] */
	double bip3;			/* IP3 static frac rate Ca release */
	double vip3;			/* IP3 flux IP3 Ca release */
		                        /*     IP3 dynamics values */
	double b2ip3;                   /* IP3 backward constant for equilibrated Tau M value */
	double a2ip3;                   /* IP3 forward constant for equilibrated Tau M value */
	double v2ip3;                   /* IP3 flux constant for Ca release */
	double v4ip3;                   /* IP3 flux constant for Ca release */
	double v3ip3;                   /* IP3 flux constant for Ca uptake */
	double k3ip3;                   /* IP3 kd for Ca uptake */
	double d1ip3;                   /* IP3 constant for determining Q value for h gate uptake */
	double d2ip3;                   /* IP3 constant for determining Q value for h gate uptake */
	double d3ip3;                   /* IP3 constant for determining Q value for h gate uptake */
	double d4ip3;                   /* IP3 constant for determining m at inifinity */
	double hip3;                    /* IP3 h gate for uptake flux */
	double hip3temp;                /* IP3 h gate for uptake flux */
	double mip3;                    /* IP3 m gate for uptake flux */
	double a3ip3;                   /* IP3 contant for determining TauH */
	double oip3;                    /* IP3 hill coefficient for uptake */
	double mtypeip3;                /* dynamic or static m value */
	};

struct synap: conn {                  	/* structure for synaptic connection */

        double oldc1;                   /* old conduct,active */
        double oldc2;                   /* old conduct,active */
        double vrev;                    /* reversal potential */
        double transrate;               /* reversal potential */
        double thresh;                  /* offset (threshold) for lookup tab */
        double nkd;                      /* nt dissociation constant (half max)*/
	double npow;			/* Hill coeff for nt binding */
        double ngain;                    /* = 1/b (mvolts for efold increase) */
        double vgain;                    /* vesicle release gain (linear mult) */
        double caegain;                  /* calcium expon vesicle release gain */
        double ckd;                      /* cG dissociation constant (half max)*/
        double chc;                      /* cG Hill coeff for 2nd mesg chan */
        double cgain;                    /* gain (multiplier) after sat. */
        double coff;                     /* offset for 2nd mesg cascade subtr. */
        double maxcond;                  /* maximum conductance (num of chans) */
        double rextern;                  /* external resistance */
        double maxsrate;                 /* maximum sustained release rate */
        double rrpool;                   /* size of readily releasible pool */
        double rrpoolg;                  /* readily releasible pool gain mult */
        double maxrrpool;                /* max size of readily releasible pool*/
        double vsize;                    /* vesicle size */
        double vcov;                     /* vesicle size stdev/mean */
        double cov;             		/* Coeff of variation multiplier */
        char *vstate;           	/* vesicle random state */
        char *gstate;           	/* vesicle size random state */
        double vtime;              	/* ves. release time, comp to vrefr */
        double vtint;              	/* ves. threshold interval for reset */
	double trconc;			/* max neurotr conc at postsyn recept */
	double mesgconc;			/* max sec mesg conc */
        short int setvrev;             	/* =1 -> vrev set by user */
        short int vrefr;              	/* refractory period for ves. release */
        short int vflag;              	/* ves. release flag, after refr done */
        short int vsites;              	/* number of ves noise release sites */
        short int sens;               	/* presyn sensitivity, Ca, V (default)*/
        short int ntact;               	/* open or closed */
        short int secmsg;              	/* has secnd messenger pathway */
        short int mesg1;              	/* output to postsyn nt compartment */
        short int mesg2;              	/* output to postsyn nt compartment */
        short int curve;               	/* linear or exponential */
        ntchan *resp1;           	/* postsynaptic channel response */
        ntchan *resp2;           	/* postsyn 2nd messenger channel response */
        dyadlst *spost;            	/* pnt to dyad (correl. synapse) */
        dyadlst *spre;            	/* pnt to other presynaptic elems */
	lpfilt *filt1;			/* low pass filter before tr release */
	lpfilt *filt1h;			/* high pass filter before tr release */
	lpfilt *filt2;			/* low pass filter after tr release */
	lpfilt *filt3;			/* low pass filter after saturation */
        };

struct randstate {			/* random state indexed by number */
	int ngen;			/* number of generator */
	randstate *next;		/* pointer to next state */
	char rstate[RNDSIZ];		/* random state */
	};

struct lpfilt: id {			/* low pass filter for synapse */
        short int nfilt;                /* number of filters hi or lo-pass */
        double *lfilt;			/* array of filters to delay psp */
        double *ftau;			/* tau's for filters */
        double tfall;                   /* falling phase time constant */
        double offset;                  /* dc offset for high pass filter */
        double *xv;                     /* input values for bessel filter */
        double *yv;                     /* output values for bessel filter */
        double gain;                    /* gain for bessel filter */
	};

struct load: id {              		/* structure for synaptic connection */
        load *next;	            	/* pointer to next load */
        load *last;	            	/* pointer to last load */
        double conduct;                  /* conductance of channel */
        comp *comp1;                    /* first compartment connection */
        comp *comp2;                    /* second compartment connection */
        comp *compe;                    /* external compartment connection */
        comp *compp;                    /* external pH compartment connection, not used */
        int num1;			/* number of first compartment */
        double vrev;                    /* reversal potential */
	};

struct photrecv: load {             	/* structure for photoreceptor */
        double xloc;                    /* location of receptor */
        double yloc;
        nodeint recnm1;			/* receptor node num1 */ 
        nodeint recnm2;			/* receptor number 2 */ 
        nodeint recnm3;			/* receptor number 3 */ 
        nodeint recnm4;			/* receptor number 4 */ 
	double pathl;			/* path length through pigment */
	double attf;			/* attenuation factor */
	double tspeed;			/* time speed factor */
	double intenadj;		/* intensity adjust factor */
	double loopgain;		/* gain for ca feedback, sets stability */
	double dnois;			/* dark contin. noise (frac of phot)*/
	short int filt;			/* filter number (0=none,1=mac pigm)*/
	short int pnois;		/* poisson noise from photon flux */
	char *dstate;			/* state array for dark noise */
	char *pstate;			/* state array for photon noise */
        double area;                    /* photon collecting area of recept */
        double *aflux;                  /* average additive photon flux */
        double *mflux;                  /* average masking photon flux */
        double *mask;                   /* amount of masking 0=>none, 1=mask */
        double *chanw;                  /* stimulus channel weighting */
        int stimchan;                   /* private stimulus channel */
        double maxcond;                 /* multiplier to get true conductance */
        double iflux;                   /* instantaneous photon flux */
        chan *phchan;                   /* photoreceptor channel */
	};

   struct phvars {
        double cond;                    /* dark conductance normalized 0 to 1 */
        double rhod;                    /* rhodopsin */
        double rhodk;                   /* rhodopsin kinase */
        double rhod2;                   /* MR I */
        double rstar;                   /* MR II (R*)  */
        double gpr;                     /* G prot */
        double pde;                     /* phosphodiesterase */
        double gcyc;                    /* G cyclase enzyme activation level */
        double cycg;                    /* cyclic G level */
        double ca;                      /* calcium level */
        double cax;                     /* delayed calcium level */
        double cab;                     /* delayed calcium level */
	double ksens;			/* sens. for simple photoreceptor */

					// variables for van hateren cone type
	double resp_r;			// rhodopsin
	double resp_b;
	double resp_e;
	double resp_x;
	double resp_c;			// calcium
	double resp_os;			// response of outer segment
	double resp_is;			// response of inner segment
	double atten_i;
    };

struct photrec: photrecv {		/* orig type rod, cone photorec */
        recpar *chtyp;          	/* receptor constants */
	phvars vars;			/* state variables */
  };

	// Alternate rod phototransduction model from:
	//
	// Invergo BM, Dell'Orco D, Montanucci L, Koch KW, Bertranpetit J.
	// A comprehensive model of the phototransduction cascade in mouse rod cells.
	// Mol Biosyst. 2014 Jun;10(6):1481-9. doi: 10.1039/c3mb70584f
	// See supplementary electronic info file
	
    struct phvarsi{			//
        double cond;                    // dark conductance normalized 0 to 1
	double Rh;			// unphosphorylated rhodopsin
	double Gt;			//
	double R_Gt;			//
	double Rn[NPHOS];		// phosphorylated rhodopsin
	double RK;			// rhodopsin kinase
	double Rn_RKpre[NPHOS];		// rhodopsin kinase
	double Rn_RKpost[NPHOS];	// rhodopsin kinase
	double Arr;			// arrestin 
	double Rn_Arr[NPHOS];		// rhodopsin-arrestin complex
	double Ops;			// ligand-free receptor 
	double Arr_di;			// arrestin dimer
	double Arr_tetra;		// arrestin tetramer
	double Ops_Gt;			// receptor Gt complex
	double Ops_G;			// receptor G complex
	double Ops_Ggtp;		// receptor G complex
	double Ggtp;			// Ggtp 
	double Rn_Gt[NPHOS];		// receptor G complex
	double Rn_G[NPHOS];		// receptor G complex
	double Rn_Ggtp[NPHOS];		// receptor G complex
	double Gagtp;			// 
	double Gbg;			// 
	double PDE;			// 
	double PDE_Gagtp;		// PDE-Gagtp complex
	double PDEa_Gagtp;		// activated PDE-Gagtp complex
	double Gagtp_PDEa_Gagtp;	// activated PDE-Gagtp complex
	double Gagtp_aPDEa_Gagtp;	// activated PDE-Gagtp complex
	double RGS;			// 
	double RGS_PDEa_Gagtp;		// activated RGS-PDE-Gagtp complex
	double RGS_Gagtp_aPDEa_Gagtp;	// activated RGS-PDE-Gagtp complex
	double Gagdp;			// Gagdp 
	double Rect;			// 
	double Recr_Ca;			// 
	double Cafree;			// Ca free
	double Recr_Ca_RK;		// 
	double Cabuff;			// Ca buffered
	double cGMP;			// cyclic GMP
	double gain;			// gain for whole rod
	double gain_sp;			// gain for single photon response
     };

struct photreci: photrecv {		/* invergo mouse rod photorec */
        recpari *chtyp;          	/* receptor constants */
	phvarsi vars;			/* invergo state varuables */
    };

struct phvarsc{				/* for Williams et al (2013) channel rhodopsin */
	double cond;
	double SC1;
	double SO1;
	double SO2;
	double SC2;
	double p;
     };

struct photrecc: photrecv {		/* for Williams et al (2013) channel rhodopsin */
        recparc *chtyp;          	/* receptor constants */
	phvarsc vars;			/* invergo state varuables */
    };

struct recpar: id {                 	/* orig structure for receptor constants */
	short int pigm;			/* pigment type (0 - 23) */
	double pathl;			/* default pigment path length */
	double recslow;			/* slowness factor */
	double loopgain;			/* gain for ca feedback, sets stability */
	double qc;			/* Q10 value for transduction shape */
	double qcond;			/* Q10 value for transduction cond */
	double qt;			/* base temp for transduction Q10 */
	double oldtempcel;		/* old temp for transduction Q10 */
        double vrev;                    /* gated channel reversal potential */

        double cat;                      /* total calcium buffer */
        double lgain;                    /* rhodopsin gain */
        double rkgain;                   /* rhodopsin kinase gain */
        double rgain1;                   /* rhod gain */
        double rgain2;                   /* rhod2 gain */
        double rgain3;                   /* rstar gain */
        double ggain;                    /* G protein gain */
        double pdebase;                  /* dark level of PDE */
        double pdegain;
        double gcygain;
        double condgain;
        double capump;
        double gca;
        double gcab;			/* forward binding const for ca -> cab */
        double gcabr;			/* reverse binding const for cab -> ca */
        double cxgain;
        double kdgcyc;
        double kdrhodk;
        double decrstar;			/* decay constants */
        double pdedec;
        double deccond;
        double deccax;
	double powcag;			/* power cax -> cag */
	double powcak;			/* power cax -> cak */
	double powcycg;			/* power cycg -> cond */
        double totcond;                  /* absolute conductance (20 x dark) */

				     	/* alternate photoreceptor constants */
					// taken from:
					// van Hateren JH, Snippe HP (2006) 
					// Simulating human cones from mid-mesopic up to 
					// high-photopic luminances. Journal of Vision 7(4):1-11.

	double tau_r;			// R* lifetime (orig ms)	
	double tau_b;			// time constant of bleaching recovery
	double tau_c;			// time constant of Ca2+ extrusion
	double tau_e;			// E* lifetime
	double tau_vc;			// membrane time constant (for resp_is)
	double rk_b;			// parameter of bleaching recovery
	double cn;			// norm const for cone bleaching 
	double beta_e_max;		// parameter diffusion-limited cGMP hydrolysis
	double a_c;			// scaling const of GC activation
	double rnc;			// Hill coeff for GC activation
	double c_beta;			// dark PDE activity
	double rk_beta;			// E* dependence of PDE activity
	double rnx;			// apparent Hill coeff for CNG channels
	double a_is;
	double gamma_is;
	double tau_is;
	double dark_resp_os;		// dark response of outer segment
	double dark_resp_is;		// dark response of inner segment

	phvars vars;			// state variables from photrec above

	};

struct recpari: id {                 	/* orig structure for receptor constants */
	short int pigm;			/* pigment type (0 - 23) */
	double pathl;			/* default pigment path length */
	double recslow;			/* slowness factor */
	double loopgain;		/* gain for ca feedback, sets stability */
	double qc;			/* Q10 value for transduction shape */
	double qcond;			/* Q10 value for transduction cond */
	double qt;			/* base temp for transduction Q10 */
	double oldtempcel;		/* old temp for transduction Q10 */
        double vrev;                    /* gated channel reversal potential */

	/* constants */

	double Rtot;			/* total rhodopsin */
	double omGt;			/* expon rate of decay of Gt affinity for R* */
	double kG1_0;			/* binding rate of Gt to unphosph R* */
	double kG2;			/* dissoc rate of Rn_Gt complex */
	double kG3;			/* dissoc rate of Gdp from Rn_Gt complex */
	double kG4gdp;			/* assoc rate of Gdp to Rn_Gt complex */
	double kG5gtp;			/* assoc rate of Gtp to Rn_Gt complex */
	double kG6;			/* dissoc rate of Rn_Ggtp complex */
	double kG7;			/* dissoc rate of Ggtp into Gbg and Gagtp */
	double kOps;			/* assoc rate of Ops and Gt */
	double kP1;			/* binding rate of PDE to Gagtp */
	double kP1rev;			/* dissoc rate of PDE_Gagtp  */
	double kP2;			/* activ rate of first subunit of PDEg of PDE_Gagtp */
	double kP3;			/* binding rate of Gagtp to PDE_Gagtp */
	double kP4;			/* activ rate of 2nd PDEg to G_agtp_PDE_Gagtp */
	double kRK1_0;			/* binding rate of RK to unphosphor R* */
	double om;			/* expon rate of decay of RK affinity for R* with phos */
	double kRK2;			/* dissoc rate of R* from RK prior to phosphorylation */
	double kRK3atp;			/* binding rate of ATP to R*_RK */
	double kRK4;			/* dissoc rate of R* from R*_RK complex */
	double kArr;			/* binding rate of Arr to singly-phos R* */
	double kA2;			/* dissoc rate of R* from Arr_R* complex prior to R* inact */
	double mArr;			/* linear rate of increase Arr aff for R* with phosph */
	double kA3;			/* dissoc rate of R* from Arr_R* complex following R* inact */
	double kA4;			/* binding rate of Arr to form homo-oliomers */
	double kA5;			/* dissoc rate of Arr from homo-oliomers */
	double kRrecyc;			/* rate const for R regen from Ops */
	double ktherm;			/* thermal decay of R* */
	double kGrecyc;			/* binding rate for Gagdp to Gbg */
	double kGshutoff;		/* rate of Gagtp auto-catalytic GTPase activity */
	double kPDEshutoff;		/* rate of PDE-induced spon PDE_Gagtp shutoff */
	double kRGS1;			/* binding rate of RGS9_1 to PDE_Gagtp */
	double kRGS2;			/* rate of hydrolysis and dissoc of a PDE from PDE_Gagtp */
	double kRec1;			/* rate of Ca-triggered Rec conf change */
	double kRec2;			/* rate of Rec conf change */
	double kRec3;			/* binding rate of Ca_Rec to RK */
	double kRec4;			/* dissoc rate of RK from Rec_Ca */
	double Vcyto;			/* outer segment cytoplasmic volume */
	double Kc1;			/* EC50 for GCAP1-mediated Ca feedback on GC activity */
	double Kc2;			/* EC50 for GCAP2-mediated Ca feedback on GC activity */
	double m1;			/* Hill coeff for GCAP1 */
	double m2;			/* Hill coeff for GCAP2 */
	double Amax;			/* max rate of cGMP synthesis */
	double Bdark;			/* dark rate of cGMP hydrolysis */
	double Bsub;			/* rate const for one catalytic PDE subunit */
	double fCa;			/* frac of circulating current carried by Ca */
	double Jdark;			/* dark circulating current */
	double cGMPdark;		/* dark cGMP concentration */
	double ncg;			/* Hill coeff for opening cGMP-gated channels */
	double gCa;			/* rate of Ca extrusion Na/Ca K+ pump */
	double Ca_dark;			/* dark Ca concentration */
	double Ca_0;			/* minimum Ca concentration */
	double k1;			/* binding rate of Ca to cytoplasmic buffers */
	double k2;			/* dissoc rate of Ca from cytoplasmic buffers */
	double eT;			/* total Ca buffer molecules conc */
	double gain;			/* calibration for overall response */
	double gain_spr;		/* calibration for single photon resp */
	double gain_spr_mul;		/* mult factor for single photon resp, def 1.0 */
	double chan_eff;		/* loss factor for resp, leakage around pipette, def 0.75 */
	phvarsi vars;			// state variables from photreci above

	};

struct recparc: id {                 	/* orig structure for receptor constants */
	short int pigm;			/* pigment type (0 - 23) */
	double pathl;			/* default pigment path length */
	double recslow;			/* slowness factor */
	double loopgain;		/* gain for ca feedback, sets stability */
	double qc;			/* Q10 value for transduction shape */
	double qcond;			/* Q10 value for transduction cond */
	double qt;			/* base temp for transduction Q10 */
	double oldtempcel;		/* old temp for transduction Q10 */
        double vrev;                    /* gated channel reversal potential */

	double q10_gd1;			/* for Williams et al (2013) channel rhodopsin */
	double q10_gd2;
	double q10_e1;
	double q10_e2;
	double q10_gr;
	double q10_e12;
	double q10_e21;
	double tchr;			/* time constant of initial rise */
	double gamma;			/* conductance weight for state O2 */
	double gain_spr_mul;
	double chan_eff;

	phvarsc vars;			// state variables from photreci above
	};

struct recstim: id {                	/* structure for receptor stimulus */
        recstim *next;         		/* pointer to next rec stim */
        recstim *last;         		/* pointer to last rec stim */
        recstim *rnext;        		/* pointer to next stim in bin table */
        recstim *rlast;        		/* pointer to last stim in bin table */
        short int mode;                 /* absolute or delta */
        short int bgend;                /* beginning or end of delta */
        nodeint recnm1;                 /* receptor number num 1*/
        nodeint recnm2;                 /* receptor number 2 */
        nodeint recnm3;                 /* receptor number 3 */
        nodeint recnm4;                 /* receptor number 4 */
        double time;                    /* time of stimulus */
        double val;                     /* stimulus intensity */
        double mask;                    /* mask, 0 =>additive, 1=>mask */
        short int stimchan;             /* stimulus channel number, allows multiple masks */
        double wavel;                   /* stimulus wavelength */
        };

struct dyadlst {                        /* list of dyad-connected synapses */
        dyadlst *next;                  /* next in list */
        synap *sdyad;                   /* pointer to synapse */
        int sdyadnum;                   /* synapse element number */
        };

