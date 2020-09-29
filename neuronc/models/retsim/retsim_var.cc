
/* retsim_var.cc */
/* Variable defs for retsim.cc */

int xcone= XCONE;             /* Identity numbers for neurons */
int xrod = XROD;
int hbat = HBAT;
int ha   = HA;
int hb   = HB;
int rbp  = RBP;
int dbp1 = DBP1;              /* depolarizing cone bipolar */
int dbp2 = DBP2;              /* depolarizing cone bipolar */
int dbp3 = DBP3;              /* depolarizing cone bipolar */
int dbp4 = DBP4;              /* depolarizing cone bipolar */
int hbp1 = HBP1;              /* hyperpolarizing cone bipolar */
int hbp2 = HBP2;              /* hyperpolarizing cone bipolar */
int a17  = A17;
int aii  = AII;
int sbac = SBAC;
int am   = AM;
int am2  = AM2;
int am3  = AM3;
int am4  = AM4;
int amh  = AMH;
int amh2 = AMH2;
int ams  = AMS;
int amhs = AMHS;
int gca  = GCA;
int gcb  = GCB;
int dsgc = DSGC;
int gcaoff = GCAOFF;
int gcboff = GCBOFF;
int nceltypes = NCELTYPES;

int xglut    = XGLUT;		/* generic glutamate response */
int xampa    = XAMPA;		/* AMPA synaptic response */
int xampa1   = XAMPA;		/* AMPA type 2 synaptic response */
int xampa2   = XAMPA2;		/* AMPA type 2 synaptic response */
int xampa3   = XAMPA3;		/* AMPA type 3 synaptic response */
int xampa4   = XAMPA4;		/* AMPA type 4 synaptic response */
int xampa5   = XAMPA5;		/* AMPA type 5 synaptic response */
int xnmda    = XNMDA;		/* NMDA type 1 synaptic response */
int xnmda2   = XNMDA2;		/* NMDA type 2 synaptic response */
int xkainate = XKAINATE;	/* Kainate synaptic response */
int xmglur6  = XMGLUR6;		/* mGluR6 synaptic response */
int xgaba    = XGABA;		/* GABA (type 1) synaptic response */
int xgaba1   = XGABA1;		/* GABA type 1 synaptic response */
int xgaba2   = XGABA2;		/* GABA type 2 synaptic response */
int xgaba3   = XGABA3;		/* GABA type 3 synaptic response */
int xgaba4   = XGABA4;		/* GABA type 4 synaptic response */
int xgly     = XGLY;		/* Glycine synaptic response */
int xgapj    = XGAPJ;		/* gap junction synaptic response */
int xdyad    = XDYAD;		/* dyad synapse (uses other resp type) */
int nresptypes =NRESPTYPES;	/* number of synaptic types */

#include "nval_var.cc"		/* defs of nval parameter variables */

/* Regions: Column number 7 in morphology file, columns in density file */


int EGRAD     = -1;
int LGRAD     = -1;
int PSOMA     = -1;	/* Point soma for extended soma region, same as SOMA */
			/* PSOMA node has parent number set to one of its descendants */
			/* which is ignored */

int R1        = R_1;
int R2        = R_2;
int R3        = R_3;
int R4        = R_4;
int R5        = R_5;
int R6        = R_6;
int R7        = R_7;
int R8        = R_8;
int R9        = R_9;
int R10       = R_10;

int DEND      = R_11; 
int DENDD     = R_11;
int DEND_DIST = R_11;
int DENDM     = R_12;
int DEND_MED  = R_12;
int DENDP     = R_13;
int DEND_PROX = R_13;
int SOMA      = R_14;
int HCK       = R_15;
int HILLOCK   = R_15;
int AXONT     = R_16;
int AXON_THIN = R_16;
int AXON      = R_17;
int AXONP     = R_17;
int AXON_PROX = R_17;
int AXOND     = R_18;
int AXON_DIST = R_18;
int VARIC     = R_19;
int VARICOS   = R_19;

int NREGIONS  = R_NREGIONS; 
			 
/* Channel types: first column defines row numbers in "dens_xxx.n" */
/* These numbers are used as the chan index in celdens[ct][cn][chan][r] */

int _NA  = C_NA2;     /* NaV1.2, type 2, markov */
int _NA0 = C_NA0;     /* NaV1.2, no noise, type 0 */
int _NA1 = C_NA1;     /* NaV1.2, simple noise */
int _NA2 = C_NA2;     /* NaV1,2 same as type 1, but markov */
int _NA3 = C_NA3;     /* Na: type 3 */
int _NA4 = C_NA4;     /* Na: type 4 */
int _NA5 = C_NA5;     /* slowly inactivating Na: type 5 */
int _NA6 = C_NA6;     /* persistent Na: type 6 */
int _NA8 = C_NA8;     /* persistent Na: type 6 */

int _K0  = C_K0;      /* type 0 , no noise */
int _K1  = C_K1;      /* type 0 , noise: type 1 */
int _K2  = C_K2;      /* type 2 , noise: type 1 */
int _K3  = C_K3;      /* KA, type 3 */
int _K4  = C_K4;      /* Ih, type 4 */
int _K5  = C_K5;      /* Kir, Hz cell, type 5 */
int _K6  = C_K6;      /* Kv3, SBAC, type 6 */
int _K7  = C_K7;      /* Kv3, SBAC, type 7 */

// joesterle begin
int _HCN1 = C_HCN1;
int _HCN2 = C_HCN2;
int _HCN3 = C_HCN3;
int _HCN4 = C_HCN4;
// joesterle end

int _K   = C_K1;      /* type 1, noise */
int _KDR = C_K1;      /* type 1, noise */
int _KA  = C_K3;      /* KA, type 3 */
int _KH  = C_K4;      /* Ih, type 4 */
int _KIR = C_K5;      /* Kir, Hz cell, type 5 */
int _KHZ = C_K5;      /* Kir, Hz cell, type 5 */
int _KV3 = C_K6;      /* Kv3, SBAC, type 6 */
int _KV3B = C_K7;     /* Kv3b, SBAC, type 7 */

int _KCA0  = C_KCA0;    /* BK, HH-type */
int _KCA1  = C_KCA1;    /* SK, no voltage sensitivity, Markov */
int _KCA2  = C_KCA2;    /* BK, Markov */
int _KCA3  = C_KCA3;    /* BK, fast, Markov */
int _KCA4  = C_KCA4;    /* SK 1, moderately slow, Markov */
int _KCA5  = C_KCA5;    /* SK 2, slow, Markov */
int _KCA6  = C_KCA6;    /* SK 2, slow, Markov */

int _BKCA  = C_KCA3;    /* BK, fast, Markov */
int _SKCA1 = C_KCA4;    /* SK 1, no voltage sensitivity, Markov */
int _SKCA2 = C_KCA5;    /* SK 2, no voltage sensitivity, slow, Markov */

int _CLCA1 = C_CLCA1;   /* Calcium-actived chloride current */
int _CLCA2 = C_CLCA2;   /* Calcium-actived chloride current, core Ca */

int _CLCA  = C_CLCA1;   /* Calcium-actived chloride current */
int _CLCAC = C_CLCA2;   /* Calcium-actived chloride current, core Ca */

int _CA0  = C_CA0;     /* L-type Ca */
int _CA1  = C_CA1;     /* L-type Ca, Markov */
int _CA2  = C_CA2;     /* T-type Ca */
int _CA3  = C_CA3;     /* T-type Ca, Markov */
int _CA4  = C_CA4;     /* T-type Ca */
int _CA5  = C_CA5;     /* T-type Ca, Markov  */
int _CA6  = C_CA6;     /* T-type Ca, low thresh */
int _CA7  = C_CA7;     /* T-type Ca, low thresh */

int _CA   = C_CA0;     /* L-type Ca */
int _CA_L = C_CA0;     /* L-type Ca */
int _CA_T = C_CA6;     /* T-type Ca */

// joesterle begin
int _MGLUR = C_MGLUR;     /* MGLUR6*/
// joesterle end

int _NMDA1 = C_NMDA1;     /* NMDA simple */
int _NMDA2 = C_NMDA2;     /* NMDA complex */

int _AMPA1 = C_AMPA1;     /* AMPA type 1 */
int _AMPA2 = C_AMPA2;     /* AMPA type 2 */
int _AMPA3 = C_AMPA3;     /* AMPA type 3 */
int _AMPA4 = C_AMPA4;     /* AMPA type 4 */
int _AMPA5 = C_AMPA5;     /* AMPA type 5 */

int _GABA1 = C_GABA1;    /* GABA type 1 */
int _GABA2 = C_GABA2;    /* GABA type 2 */
int _GABA3 = C_GABA3;    /* GABA type 3 */
int _GABA4 = C_GABA4;    /* GABA type 4 */

int _GLY   = C_GLY;      /* GLY type 1 */

int _CGMP1  = C_CGMP1;      /* CGMP type 1 */
int _CGMP2  = C_CGMP2;      /* CGMP type 2 */
int _CGMP3  = C_CGMP3;      /* CGMP type 3 */
int _CGMP4  = C_CGMP4;      /* CGMP type 4 */
int _CGMP5  = C_CGMP5;      /* CGMP type 5 */
int _CGMP6  = C_CGMP6;      /* CGMP type 6 */
int _CGMP7  = C_CGMP7;      /* CGMP type 7 */
int _CGMP8  = C_CGMP8;      /* CGMP type 8 */
int _CGMP9  = C_CGMP9;      /* CGMP type 9 */
int _CGMP10 = C_CGMP10;     /* CGMP type 10 */
int _CGMP11 = C_CGMP11;     /* CGMP type 11 */

int _CAP    = C_CAPUMP;  /* capump vmax */
int _CAPUMP = C_CAPUMP;  /* capump factor */
int _CAPK   = C_CAPKM;   /* ca pump km */
int _CABV   = C_CABVMAX; /* cabuf vmax */
int _CABK   = C_CABKD;   /* cabuf kd */
int _CABT   = C_CABTOT;  /* cabuf btot */
int _CABI   = C_CABTOTI; /* cabuf btoti */
int _CASH   = C_CASHELL; /* ca nshells */
int _CAE    = C_CAEXCH;  /* ca exch factor */
int _CAEXCH = C_CAEXCH;  /* ca exch factor */

int _CAS    = C_CAS;		/* initial Ca in CICR store (M) */
int _VM2    = C_VM2;		/* max rate of Ca pump into ryanodine store */
int _VM3    = C_VM3;		/* max rate of Ca pump from ryanodine store */
int _KFC    = C_KFCICR;		/* passive leak rate from ryan store /s */
int _KFCICR = C_KFCICR;		/* passive leak rate from ryan store /s */
int _KAC    = C_KACICR;		/* thresh for ryan store release (M) */
int _KACICR = C_KACICR;		/* thresh for ryan store release (M) */
int _KRC    = C_KRCICR;		/* thresh for ryan store release (M) */
int _KRCICR = C_KRCICR;		/* thresh for ryan store release (M) */
int _K1C    = C_K1CICR;		/* thresh for ryan store uptake (M) */
int _K1CICR = C_K1CICR;		/* thresh for ryan store uptake (M) */
int _K2C    = C_K2CICR;		/* thresh for ryan store uptake (M) */
int _K2CICR = C_K2CICR;		/* thresh for ryan store uptake (M) */
int _NHC    = C_NCICR;		/* Hill Coeff, ryan store uptake */
int _NCICR  = C_NCICR;		/* Hill Coeff, ryan store uptake */
int _MHC    = C_MCICR;		/* Hill Coeff, ryan store release */
int _MCICR  = C_MCICR;		/* Hill Coeff, ryan store release */
int _PHC    = C_PCICR;		/* Hill Coeff, ryan store release */
int _PCICR  = C_PCICR;		/* Hill Coeff, ryan store release */
int _C1CICR = C_C1CICR;         /* volume fraction of calcium store */
int _CAS2   = C_CAS2;		/* initial Ca in IP3 store (M) */
int _IP3I   = C_IP3I;		/* initial [IP3] intrcellular */
int _BIP3   = C_BIP3;		/* init frac rate of IP3 release */
int _VIP3   = C_VIP3;		/* init IP3 store flux */
int _V2IP3  = C_V2IP3;		/* ip3 release rate */
int _V3IP3  = C_V3IP3;		/* ip3 uptake rate */
int _V4IP3  = C_V4IP3;		/* ip3 release rate */
int _D1IP3  = C_D1IP3;		/* IP3 store const for Q value for h gate, 0.13e-6 M */
int _D2IP3  = C_D2IP3;		/* IP3 store const for Q value for h gate, 1.049e-6 M  */
int _D3IP3  = C_D3IP3;		/* IP3 store const for Q value for h gate, 0.9434e-6 M */
int _D4IP3  = C_D4IP3;		/* IP3 store const for m gate at infinity, 0.823e-7 M */
int _A2IP3  = C_A2IP3;		/* IP3 store forward const for Tau M value, 4.2e8/ms */
int _A3IP3  = C_A3IP3;		/* IP3 store constant for IP3 Tau H, 2e5 1/MS */
int _B2IP3  = C_B2IP3;		/* IP3 store backward const for Tau M value, 4.1/ms */
int _K3IP3  = C_K3IP3;		/* IP3 store kd for Ca uptake, 0.1e-6 M) */

int _CVSTART= C_VSTART;  /* starting voltage */
int _VST    = C_VSTART;  /* starting voltage */
int _CVREV  = C_VREV;    /* membrane reversal potential */
int _VRV    = C_VREV;    /* membrane reversal potential */
int _RM     = C_RM;      /* membrane resistivity */
int _RI     = C_RI;      /* axial resistivity */
int _CM     = C_CM;      /* membrane capacitance */
int _DIA    = C_DDIA;    /* dia factor for region */
int _CPLAM  = C_CPLAM;   /* compartment lambda for region */
int _CPL    = C_CPLAM;   /* compartment lambda for region */
int _CMUL   = C_CMUL;    /* synaptic conductance mult for region */
int _COLOR  = C_COLOR;   /* color for region */
int _COL    = C_COLOR;   /* color for region */

/*--------------------------------------------------------------------*/

/* defs for chanparams */

#include "chanparams.cc"

/*--------------------------------------------------------------------*/

const char *cname[NCELTYPES];                   /* cell names */
const char *rname[NRESPTYPES];                  /* synaptic types */
double nval[NCELTYPES][NPARAMS];                /* cell build params */
extern const char *densfil[NCELTYPES];
int ncell_erased[NCELTYPES];			/* numbers of erased cells */

int cellconn [NCELTYPES][NCELTYPES][NSYNTYPES]; /* connections to other cell types */
int pickden [NCELTYPES] = {0};                  /* dendrites to allow connection to */
int celnode[NCELTYPES][MAXDEN] = {0};		/* numbers of dendrites */

// joesterle begin
const char* chname[] = { "Na0", "Na1", "Na2", "Na3", "Na4", "Na5", "Na6", "Na8", 
		"K0", "Kdr", "K2", "KA", "KH", "Kir", "Kv3", "Kv3b",
                "KCa0", "KCa1", "KCa2", "BKCa", "SKCa1", "SKCa2", "KCa6", "ClCa", "ClCaC",
		"Ca0", "Ca1", "Ca2", "Ca3", "Ca4", "Ca5", "Ca6", "Ca7", "MGLUR", "NMDA1", "NMDA2", 
		"AMPA", "AMPA2", "AMPA3", "AMPA4", "AMPA5", "GABA1", "GABA2", "GABA3", "GABAC", "GLY",
		"CGMP1", "CGMP2", "CGMP3", "CGMP4", "CGMP5", "CGMP6", "CGMP7", "CGMP8", "CGMP9", "CGMP10", "CGMP11",
    "HCN1", "HCN2", "HCN3", "HCN4",
		"Capump", "CapumpKm", "CabufVmax", "CabufKd", "Cabuftot", "Cabuftoti", "Cashell", "Caexch", 
		"Cas","Vmax2", "Vmax3", "KaCICR", "KfCICR", "KrCICR", "K1CICR", "K2CICR","NCICR","MCICR","PCICR",
		"C1CICR","Cas2", "IP3i", "BetaIP3","V1IP3","V3IP3","V2IP3","V4IP3","D1IP3","D2IP3","D3IP3",
		"D4IP3","A2IP3","A3IP3","B2IP3","K3IP3","Vstart", "Vrev", "Rm", "Ri", "Cm", "Dia", "cplam",
		"condmul", "color"};

const char* chnamea[] = { "NA0", "NA1", "NA2", "NA3", "NA4", "NA5", "NA6", "NA8", 
		"K0", "K1", "K2", "K3", "K4", "K5", "K6", "K7",
                "KCA0", "KCA1", "KCA2", "KCA3", "KCA4", "KCA5", "KCA6", "CLCA1", "CLCA2",
		"CA0", "CA1", "CA2", "CA3", "CA4", "CA5", "CA6", "CA7", "MGLUR", "NMDA1", "NMDA2",
		"AMPA1", "AMPA2", "AMPA3", "AMPA4", "AMPA5", "GABA1", "GABA2", "GABA3", "GABA4", "GLY",
		"CGMP1", "CGMP2", "CGMP3", "CGMP4", "CGMP5", "CGMP6", "CGMP7", "CGMP8", "CGMP9", "CGMP10", "CGMP11",
    "HCN1", "HCN2", "HCN3", "HCN4",
		"CAP", "CAPK", "CABV", "CABK", "CABT", "CABI", "CASH", "CAE", 
		"CAS", "VM2", "VM3", "KAC", "KFC", "KRC", "K1C", "K2C", "NHC", "MHC", "PHC","C1CICR",
		"CAS2", "IP3I", "BIP3", "VIP3", "V3IP3", "V2IP3", "V4IP3", "D1IP3", "D2IP3", "D3IP3",
		"D4IP3","A2IP3","A3IP3","B2IP3","K3IP3","VST", "VRV", "RM", "RI", "CM", "DIA", "CPL", "CMUL", "COL" };
// joesterle end
    
const char *regname[] =  { "R1","R2","R3","R4","R5","R6","R7","R8","R9", "R10", 
			"dend_dist","dend_med", "dend_prox", "soma","hillock",
			"axon_thin","axon_prox", "axon_dist", "varicos" };

const char *parmname[] =  { "offm","offh","taua","taub","tauc","taud"};

/* standard colors */
/*  see "nc/src/colors.h" */

int black=BLACK;
int blue=BLUE;
int green=GREEN;
int cyan=CYAN;
int red=RED;
int magenta=MAGENTA;
int brown=BROWN;
int white=WHITE;

int gray=GRAY;
int ltblue=LTBLUE;
int ltgreen=LTGREEN;
int ltcyan=LTCYAN;
int ltred=LTRED;
int ltmag=LTMAG;
int yellow=YELLOW;
int brtwht=BRTWHT;

const char *colornames[] = { "black","blue","green","cyan","red","magenta","brown","white",
			  "gray", "ltblue", "ltgreen", "ltcyan", "ltred", "ltmag", "yellow", "brtwht"};

const char *nvalfile;         	/* default neuron params */
const char *ha_file;	  	/* ha horiz cell anatomy */
const char *hb_file;	  	/* hb horiz cell anatomy */
const char *hbat_file;	  	/* hb horiz cell axon term anatomy */
const char *rbp_file;	  	/* rbp  bipolar cell anatomy */
const char *dbp1_file;	  	/* dbp1 bipolar cell anatomy */
const char *dbp2_file;	  	/* dbp2 bipolar cell anatomy */
const char *dbp3_file;	  	/* dbp3 bipolar cell anatomy */
const char *dbp4_file;	  	/* dbp4 bipolar cell anatomy */
const char *hbp1_file;	  	/* hbp1 bipolar cell anatomy */
const char *hbp2_file;	  	/* hbp2 bipolar cell anatomy */
const char *am_file;	  	/* amac cell anatomy */
const char *am2_file;	  	/* amac cell anatomy */
const char *am3_file;	  	/* amac cell anatomy */
const char *am4_file;	  	/* amac cell anatomy */
const char *amh_file;	  	/* amac cell anatomy */
const char *amh2_file;	  	/* amac cell anatomy */
const char *sbac_file; 		/* sbac cell anatomy */
const char *aii_file;  		/* aii  cell anatomy */
const char *a17_file;  		/* a17  cell anatomy */
const char *gca_file;     	/* gca  cell anatomy */
const char *gcb_file;     	/* gcb  cell anatomy */
const char *dsgc_file;     	/* dsgc cell anatomy */
const char *gcaoff_file;    	/* gcaoff cell anatomy */
const char *gcboff_file;    	/* gcboff cell anatomy */
const char *anatfiles[NCELTYPES]; /* default anatomy file names */

const char *expt;
double stfreq;
extern const char *plotlabel;
extern int nolabels;

/* params that control making cells */

int make_rods	= 0;
int make_cones	= 0;
int make_ha	= 0;
int make_hb	= 0; 
int make_hbat	= 0;
int make_dbp1	= 0;
int make_dbp2	= 0;
int make_dbp3	= 0;
int make_dbp4	= 0;
int make_hbp1	= 0;
int make_hbp2	= 0;
int make_rbp	= 0;
int make_am	= 0;
int make_am2	= 0;
int make_am3	= 0;
int make_am4	= 0;
int make_amh	= 0;
int make_amh2	= 0;
int make_ams	= 0;
int make_amhs	= 0;
int make_sbac	= 0;
int make_aii	= 0;
int make_a17	= 0;
int make_gca	= 0;
int make_gcb	= 0;
int make_dsgc	= 0;
int make_gcaoff	= 0;
int make_gcboff	= 0;

/* params that control making synapses */

int make_cone_cone	= 0;
int make_rod_rbp	= 1;
int make_cone_dbp1	= 1;
int make_cone_dbp2	= 1;
int make_cone_dbp3	= 1;
int make_cone_dbp4	= 1;
int make_cone_hbp1	= 1;
int make_cone_hbp2	= 1;
int make_rod_hbat	= 1;
int make_cone_ha	= 1;
int make_cone_hb	= 1;
int make_ha_ha		= 1;
int make_hb_hb		= 1;
int make_rbp_aii	= 1;
int make_aii_aii	= 0;
int make_dbp1_aii	= 0;
int make_dbp2_aii	= 0;
int make_dbp3_aii	= 0;
int make_dbp4_aii	= 0;
int make_aii_hbp1	= 0;
int make_aii_hbp2	= 0;
int make_rbp_a17	= 1;
int make_a17_rbp	= 1;
int make_dbp1_sbac	= 1;
int make_dbp2_sbac	= 0;
int make_dbp3_sbac	= 0;
int make_dbp4_sbac	= 0;
int make_dbp1_dbp1	= 0;
int make_dbp2_dbp2	= 0;
int make_dbp3_dbp3	= 0;
int make_dbp4_dbp4	= 0;
int make_hbp1_hbp1	= 0;
int make_hbp2_hbp2	= 0;
int make_dbp1_am	= 1;
int make_dbp1_am2	= 1;
int make_dbp1_am3	= 1;
int make_dbp1_am4	= 1;
int make_dbp1_amh	= 0;
int make_dbp1_amh2	= 0;
int make_dbp1_ams	= 1;
int make_dbp2_am	= 0;
int make_dbp2_am2	= 1;
int make_dbp2_am3	= 1;
int make_dbp2_am4	= 1;
int make_dbp2_amh	= 0;
int make_dbp2_amh2	= 0;
int make_dbp2_ams	= 1;
int make_dbp3_am	= 0;
int make_dbp3_am2	= 1;
int make_dbp3_am3	= 1;
int make_dbp3_am4	= 1;
int make_dbp3_amh	= 0;
int make_dbp3_amh2	= 0;
int make_dbp3_ams	= 1;
int make_dbp4_am	= 0;
int make_dbp4_am2	= 1;
int make_dbp4_am3	= 1;
int make_dbp4_am4	= 1;
int make_dbp4_amh	= 0;
int make_dbp4_amh2	= 0;
int make_dbp4_ams	= 1;
int make_hbp1_sbac	= 0;
int make_hbp1_amh	= 1;
int make_hbp1_amh2	= 1;
int make_hbp1_amhs	= 1;
int make_hbp2_sbac	= 0;
int make_hbp2_amh	= 1;
int make_hbp2_amh2	= 1;
int make_hbp2_amhs	= 1;
int make_dbp1_gca	= 1;
int make_dbp1_gcb	= 1;
int make_dbp2_gca	= 1;
int make_dbp3_gca	= 1;
int make_dbp4_gca	= 1;
int make_dbp2_gcb	= 1;
int make_dbp3_gcb	= 1;
int make_dbp4_gcb	= 1;
int make_sbac_dsgc	= 1;
int make_sbac_sbac	= 0;
int make_sbac_dbp1	= 0;
int make_sbac_dbp2	= 0;
int make_sbac_dbp3	= 0;
int make_sbac_dbp4	= 0;
int make_dbp1_dsgc	= 1;
int make_dbp2_dsgc	= 1;
int make_dbp3_dsgc	= 1;
int make_dbp4_dsgc	= 1;
int make_hbp1_dsgc	= 1;
int make_hbp2_dsgc	= 1;
int make_ams_dbp1	= 1;
int make_ams_dbp2	= 1;
int make_ams_dbp3	= 1;
int make_ams_dbp4	= 1;
int make_ams_dsgc	= 1;
int make_amhs_hbp1	= 1;
int make_amhs_dsgc	= 1;
int make_dsgc_dsgc	= 0;
int make_am_dbp1	= 1;
int make_am2_dbp1	= 1;
int make_am3_dbp1	= 0;
int make_am4_dbp1	= 0;
int make_am_dbp2	= 1;
int make_am2_dbp2	= 1;
int make_am3_dbp2	= 1;
int make_am4_dbp2	= 1;
int make_am_dbp3	= 1;
int make_am2_dbp3	= 1;
int make_am3_dbp3	= 1;
int make_am4_dbp3	= 1;
int make_am_dbp4	= 1;
int make_am2_dbp4	= 1;
int make_am3_dbp4	= 1;
int make_am4_dbp4	= 1;
int make_am_am		= 0;
int make_am2_am2	= 0;
int make_am3_am3	= 0;
int make_am4_am4	= 0;
int make_amh_amh	= 0;
int make_amh2_amh2	= 0;
int make_am_gca		= 1;
int make_am2_gca	= 1;
int make_am3_gca	= 1;
int make_am4_gca	= 1;
int make_am_gcaoff	= 1;
int make_am_dsgc	= 1;
int make_amh_hbp1	= 1;
int make_amh_hbp2	= 1;
int make_amh2_hbp1	= 0;
int make_amh2_hbp2	= 0;
int make_amh2_amh	= 0;
int make_amh_gcaoff	= 1;
int make_amh_gcboff	= 1;
int make_amh_dsgc	= 1;
int make_ams_gca	= 1;
int make_ams_gcb	= 1;
int make_ams_gcaoff	= 1;
int make_ams_gcboff	= 1;
int make_amhs_gcaoff	= 1;
int make_amhs_gcboff	= 1;
int make_hbp1_gcaoff	= 1;
int make_hbp2_gcaoff	= 1;
int make_hbp1_gcboff	= 1;
int make_hbp2_gcboff	= 1;
int make_ha_cone	= 1;
int make_ha_dbp1	= 0;
int make_ha_dbp2	= 1;
int make_ha_dbp3	= 1;
int make_ha_dbp4	= 1;
int make_ha_hbp1	= 1;
int make_ha_hbp2	= 1;
int make_hb_dbp1	= 1;
int make_hb_dbp2	= 1;
int make_hb_dbp3	= 1;
int make_hb_dbp4	= 1;
int make_hb_hbp1	= 1;
int make_hb_hbp2	= 1;

int made_gc_comps	= 0;

// int info = 1;    		/* level of information printout from binary code */
int ninfo = 1;    		/* level of information printout */
int info_chan = 0;		/* information about channels */
int info_cell = 1;		/* information about cells */
int info_syn  = 1;		/* information about synapses */
int info_disp = 0;		/* information about display */
int script = 0;			/* use when running higher level script */
extern int disp_ray;		/* =1 -> output povray format */

int gca_morph;
int gcb_morph;
int dsgc_morph;
int gcaoff_morph;
int gcboff_morph;
int sb_morph;
int am_morph;
int am2_morph;
int am3_morph;
int am4_morph;
int amh_morph;
int amh2_morph;
int rbp_morph;
int aii_morph;
int a17_morph;

int gca_biophys;
int gcb_biophys;
int dsgc_biophys;
int gcaoff_biophys;
int gcboff_biophys;
int am2_biophys;
int am3_biophys;
int am4_biophys;
int am_biophys;
int amh_biophys;
int amh2_biophys;
int rbp_biophys;
int aii_biophys;
int a17_biophys;
int sb_biophys;

double dbp1_thetax;
double dbp1_thetay;
double dbp1_thetaz;
int dbp1_flip;
double dbp2_thetax;
double dbp2_thetay;
double dbp2_thetaz;
int dbp2_flip;
double dbp3_thetax;
double dbp3_thetay;
double dbp3_thetaz;
int dbp3_flip;
double dbp4_thetax;
double dbp4_thetay;
double dbp4_thetaz;
int dbp4_flip;
double hbp1_thetax;
double hbp1_thetay;
double hbp1_thetaz;
int hbp1_flip;
double hbp2_thetax;
double hbp2_thetay;
double hbp2_thetaz;
int hbp2_flip;
double rbp_thetax;
double rbp_thetay;
double rbp_thetaz;
int rbp_flip;
double ha_thetax;
double ha_thetay;
double ha_thetaz;
int ha_flip;
double hb_thetax;
double hb_thetay;
double hb_thetaz;
int hb_flip;
double hbat_thetax;
double hbat_thetay;
double hbat_thetaz;
int hbat_flip;
double am_thetax;
double am_thetay;
double am_thetaz;
int am_flip;
double am2_thetax;
double am2_thetay;
double am2_thetaz;
int am2_flip;
double am3_thetax;
double am3_thetay;
double am3_thetaz;
int am3_flip;
double am4_thetax;
double am4_thetay;
double am4_thetaz;
int am4_flip;
double ams_thetax;
double ams_thetay;
double ams_thetaz;
int ams_flip;
double amh_thetax;
double amh_thetay;
double amh_thetaz;
int amh_flip;
double amh2_thetax;
double amh2_thetay;
double amh2_thetaz;
int amh2_flip;
double aii_thetax;
double aii_thetay;
double aii_thetaz;
int aii_flip;
double a17_thetax;
double a17_thetay;
double a17_thetaz;
int a17_flip;
double sbac_thetax;
double sbac_thetay;
double sbac_thetaz;
int sbac_flip;
double dsgc_thetax;
double dsgc_thetay;
double dsgc_thetaz;
int dsgc_flip;
double gc_thetax;
double gc_thetay;
double gc_thetaz;
int gc_flip;

int flip;			/* flip the cell in x */

extern double am_seglen;	/* definitions in sb_makfuncs.cc */
extern double am_den_seg;
extern double xradius;

double dvrev;			/* for dens_xxx.n files */
double dvst;

/* Use* Leave unset numbers = -1 for automatic algorithms below. */

int n_rods  = -1;
int n_rbp   = -1;
int n_cones = -1;
int n_ha    = -1;
int n_hb    = -1;
int n_hbat  = -1;
int n_dbp1  = -1;
int n_dbp2  = -1;
int n_dbp3  = -1;
int n_dbp4  = -1;
int n_hbp1  = -1;
int n_hbp2  = -1;
int n_aii   = -1;
int n_a17   = -1;
int n_am    = -1;
int n_am2   = -1;
int n_am3   = -1;
int n_am4   = -1;
int n_amh   = -1;
int n_amh2  = -1;
int n_ams   = -1;
int n_amhs  = -1;
int n_sbac  = -1;
int n_gca   = -1;
int n_gcb   = -1;
int n_dsgc  = -1;
int n_gcaoff = -1;
int n_gcboff = -1;

int nrods  = 0;
int nrbp   = 0;
int ncones = 0;
int nha    = 0;
int nhb    = 0;
int nhbat  = 0;
int ndbp1  = 0;
int ndbp2  = 0;
int ndbp3  = 0;
int ndbp4  = 0;
int nhbp1  = 0;
int nhbp2  = 0;
int naii   = 0;
int na17   = 0;
int nam    = 0;
int nam2   = 0;
int nam3   = 0;
int nam4   = 0;
int namh   = 0;
int namh2  = 0;
int nams   = 0;
int namhs  = 0;
int nsbac  = 0;
int ngca   = 0;
int ngcb   = 0;
int ndsgc  = 0;
int ngcaoff = 0;
int ngcboff = 0;

int na1nois;
int na2nois;
int na3nois;
int na4nois;
int na5nois;
int na6nois;
int na8nois;
int k1nois;
int k3nois;
int k4nois;
int k5nois;
int k6nois;
int k7nois;

// begin joesterle
int hcn1nois;
int hcn2nois;
int hcn3nois;
int hcn4nois;
// end joesterle


int kca1nois;
int kca3nois;
int kca4nois;
int kca5nois;
int kca6nois;
int ca1nois;
int ca3nois;
int ca5nois;
int ca6nois;
int ca7nois;
int clcanois;

int maxha;
int maxhb;
int maxhbat;
int maxrbp;
int maxdbp1;
int maxdbp2;
int maxdbp3;
int maxdbp4;
int maxhbp1;
int maxhbp2;
int maxa17;
int maxaii;
int maxsb;
int maxam;
int maxam2;
int maxam3;
int maxam4;
int maxamh;
int maxamh2;
int maxams;
int maxamhs;
int maxgca;
int maxgcb;
int maxdsgc;
int maxgcaoff;
int maxgcboff;

double node_scale;

double rod_nscale;
double cone_nscale;
double hbat_nscale;
double ha_nscale;
double hb_nscale;
double rbp_nscale;
double dbp1_nscale;
double dbp2_nscale;
double dbp3_nscale;
double dbp4_nscale;
double hbp1_nscale;
double hbp2_nscale;
double aii_nscale;
double sbac_nscale;
double am_nscale;
double am2_nscale;
double am3_nscale;
double am4_nscale;
double amh_nscale;
double amh2_nscale;
double ams_nscale;
double amhs_nscale;
double a17_nscale;
double gca_nscale;
double gcb_nscale;
double dsgc_nscale;
double gcaoff_nscale;
double gcboff_nscale;

double synap_scale = 1;
double gj_scale = 0.5;

double arrsiz;
double xarrsiz;
double yarrsiz;
double zarrsiz;
double garrsiz;
double gxarrsiz;
double gyarrsiz;
double amarrsiz;
double amxarrsiz;
double amyarrsiz;
double am2arrsiz;
double am2xarrsiz;
double am2yarrsiz;
double am3arrsiz;
double am3xarrsiz;
double am3yarrsiz;
double am4arrsiz;
double am4xarrsiz;
double am4yarrsiz;
double amharrsiz;
double amhxarrsiz;
double amhyarrsiz;
double amh2arrsiz;
double amh2xarrsiz;
double amh2yarrsiz;
double amsarrsiz;
double amsxarrsiz;
double amsyarrsiz;
double amhsarrsiz;
double amhsxarrsiz;
double amhsyarrsiz;
double aiiarrsiz;
double aiixarrsiz;
double aiiyarrsiz;
double dbp1arrsiz;
double dbp1xarrsiz;
double dbp1yarrsiz;
double dbp2arrsiz;
double dbp2xarrsiz;
double dbp2yarrsiz;
double dbp3arrsiz;
double dbp3xarrsiz;
double dbp3yarrsiz;
double dbp4arrsiz;
double dbp4xarrsiz;
double dbp4yarrsiz;
double hbp1arrsiz;
double hbp1xarrsiz;
double hbp1yarrsiz;
double hbp2arrsiz;
double hbp2xarrsiz;
double hbp2yarrsiz;
double rodarrsiz;
double rodxarrsiz;
double rodyarrsiz;
double conearrsiz;
double conexarrsiz;
double coneyarrsiz;
double rbparrsiz;
double rbpxarrsiz;
double rbpyarrsiz;
double sbarrsiz;
double sbxarrsiz;
double sbyarrsiz;
double arrcentx;
double arrcenty;
double arrcentz;
double rodarrcentx;
double rodarrcenty;
double conearrcentx;
double conearrcenty;
double dbp1arrcentx;
double dbp1arrcenty;
double dbp2arrcentx;
double dbp2arrcenty;
double dbp3arrcentx;
double dbp3arrcenty;
double dbp4arrcentx;
double dbp4arrcenty;
double hbp1arrcentx;
double hbp1arrcenty;
double hbp2arrcentx;
double hbp2arrcenty;

int dbp1_first_cent;
int dbp2_first_cent;
int dbp3_first_cent;
int dbp4_first_cent;
int hbp1_first_cent;
int hbp2_first_cent;

int setarrsiz;
int setxarrsiz;

double *amxarr=NULL,  *amyarr=NULL,  *amtharr=NULL;
double *am2xarr=NULL, *am2yarr=NULL, *am2tharr=NULL;
double *am3xarr=NULL, *am3yarr=NULL, *am3tharr=NULL;
double *am4xarr=NULL, *am4yarr=NULL, *am4tharr=NULL;
double *amhxarr=NULL,  *amhyarr=NULL , *amhtharr=NULL;
double *amh2xarr=NULL, *amh2yarr=NULL, *amh2tharr=NULL;
double *amsxarr=NULL, *amsyarr=NULL, *amstharr=NULL;
double *amhsxarr=NULL, *amhsyarr=NULL, *amhstharr=NULL;
double *aiixarr=NULL, *aiiyarr=NULL, *aiiytharr=NULL, *aiitharr=NULL;
double *sbxarr=NULL, *sbyarr=NULL, *sbytharr=NULL, *sbtharr=NULL;
double *a17xarr=NULL, *a17yarr=NULL, *a17tharr=NULL;
double *dbp1xarr=NULL, *dbp1yarr=NULL, *dbp1ytharr=NULL, *dbp1tharr=NULL;
double *dbp2xarr=NULL, *dbp2yarr=NULL, *dbp2ytharr=NULL, *dbp2tharr=NULL;
double *dbp3xarr=NULL, *dbp3yarr=NULL, *dbp3ytharr=NULL, *dbp3tharr=NULL;
double *dbp4xarr=NULL, *dbp4yarr=NULL, *dbp4ytharr=NULL, *dbp4tharr=NULL;
double *hbp1xarr=NULL, *hbp1yarr=NULL, *hbp1ytharr=NULL, *hbp1tharr=NULL;
double *hbp2xarr=NULL, *hbp2yarr=NULL, *hbp2ytharr=NULL, *hbp2tharr=NULL;
double *rodxarr=NULL, *rodyarr=NULL, *rodtharr=NULL;
double *conexarr=NULL, *coneyarr=NULL, *conetharr=NULL;
double *rbpxarr=NULL, *rbpyarr=NULL, *rbpytharr=NULL, *rbptharr=NULL;
double *gcxarr=NULL, *gcyarr=NULL, *gcytharr=NULL, *gctharr=NULL;
int *amhnarr=NULL, *amh2narr=NULL, *amsnarr=NULL, *amnarr=NULL, *am2narr=NULL, *am3narr=NULL, *am4narr=NULL, *aiinarr=NULL, *sbnarr=NULL, *a17narr=NULL, *dbp1narr=NULL, *dbp2narr=NULL, *dbp3narr=NULL, *dbp4narr=NULL, *hbp1narr=NULL, *hbp2narr=NULL, *rodnarr=NULL, *conenarr=NULL, *rbpnarr=NULL, *gcnarr=NULL;

double dispsize;
double maxsize;

double axarbdia;
double ath_dia;
double cbplam;
double cbplam2;

double dia_prox_rad;
double dia_prox_factor;
double gctheta;
double sbtheta;

int bp_dend1;

int hz_dend1;
int hz_dend2;
int hz_dend3;
int hz_dend4;
int hz_dend5;

int aii_dend1;
int aii_dend2;
int aii_dend3;
int aii_dend4;
int aii_dend5;

int sbac_dend1;
int sbac_dend2;
int sbac_dend3;
int sbac_dend4;
int sbac_dend5;

int am_dend1;
int am_dend2;
int am_dend3;
int am_dend4;
int am_dend5;

int gc_dend1;
int gc_dend2;
int gc_dend3;
int gc_dend4;
int gc_dend5;

double dnoise;
int pnoise, vnoise, cnoise;
int Chnoise;
double bg_inten;

int rod_rect_arr;
int cone_rect_arr;
int make_one_dbp1;
int make_one_dbp2;
int make_one_dbp3;
int make_one_dbp4;
int remove_nconns;
int print_conns;
int print_avg_conns;
int disp_ct;
int disp_cn;
int plot_freq;

double rodspac;
double conespac;

double cone_maxcond;
double rod_maxcond;
int cone_pigm;

double cone_timec;
double rod_timec;;

double cone_loopg;
double rod_loopg;;

int gcdistnod;

double xmax, xmin;
double ymax, ymin;
double zmax, zmin;
double set_zmax, set_zmin;

double disp_zmax, disp_zmin;
double disp_gc_zmax, disp_gc_zmin;
double disp_dsgc_zmax, disp_dsgc_zmin;

int xmaxnode, xminnode;
int ymaxnode, yminnode;
int zmaxnode, zminnode;

double disp_margin;
double disp_calib_len;
int limit_array;
int update_array;

double ttxbath;
double ttxsoma;
double ttxdend;
double tea;
double fourap; 
double zd7288; 
double ibtox; 

double dia_min_factor;
double dendd_dia_factor;
double dendi_dia_factor;
double dendp_dia_factor;
double dend_dia_factor;
double ax_dia_factor;
double cell_dia_factor;
double den_dia;
double taper;
double elec_resist1, elec_resist2;
double elec_capac1, elec_capac2;

double gc_dend_cplam;
double am_dend_cplam;
//double sbac_dend_cplam;  // defined in sb_makfuncs.cc

double gc_rm;
double gc_vr;
double gc_vs;
double dvs;
double mxrot;
double myrot;

int nshell;
int nshell_soma;

const char *confdir;
const char *chanparamsfile;
const char *def_densfile; 	/* def biophys */
const char *syn_savefile;
const char *syn_restorefile;

const char *cone_densfile;
const char *rod_densfile;
const char *hbat_densfile;
const char *ha_densfile;
const char *hb_densfile;
const char *rbp_densfile;
const char *dbp1_densfile;
const char *dbp2_densfile;
const char *dbp3_densfile;
const char *dbp4_densfile;
const char *hbp1_densfile;
const char *hbp2_densfile;
const char *a17_densfile;
const char *aii_densfile;
const char *sbac_densfile;
const char *am_densfile;
const char *am2_densfile;
const char *am3_densfile;
const char *am4_densfile;
const char *amh_densfile;
const char *amh2_densfile;
const char *ams_densfile;
const char *amhs_densfile;
const char *gca_densfile;
const char *gcb_densfile;
const char *dsgc_densfile;
const char *gcaoff_densfile;
const char *gcboff_densfile;

const char *cone_densfile2;
const char *rod_densfile2;
const char *hbat_densfile2;
const char *ha_densfile2;
const char *hb_densfile2;
const char *rbp_densfile2;
const char *dbp1_densfile2;
const char *dbp2_densfile2;
const char *dbp3_densfile2;
const char *dbp4_densfile2;
const char *hbp1_densfile2;
const char *hbp2_densfile2;
const char *a17_densfile2;
const char *aii_densfile2;
const char *sbac_densfile2;
const char *am_densfile2;
const char *am2_densfile2;
const char *am3_densfile2;
const char *am4_densfile2;
const char *amh_densfile2;
const char *amh2_densfile2;
const char *ams_densfile2;
const char *amhs_densfile2;
const char *gca_densfile2;
const char *gcb_densfile2;
const char *dsgc_densfile2;
const char *gcaoff_densfile2;
const char *gcboff_densfile2;


/*--------------------------------------------------------*/

#define SETPTRS 1

void setptrs(void)

{
   setptrn("rseed",      &rseed);
   setptrn("ninfo",      &ninfo);
   setptrn("info_chan",  &info_chan);
   setptrn("info_cell",  &info_cell);
   setptrn("info_syn",   &info_syn);
   setptrn("info_disp",  &info_disp);
   setptrn("script",     &script);

   setptrn("expt",	 &expt);
   setptrn("relax",	 &relax);
   setptrn("relincr",	 &relincr);

   setptrn("xcone",  &xcone);
   setptrn("xrod",   &xrod);
   setptrn("rbp",    &rbp);
   setptrn("ha",     &ha);
   setptrn("hb",     &hb);
   setptrn("hbat",   &hbat);
   setptrn("dbp1",   &dbp1);
   setptrn("dbp2",   &dbp2);
   setptrn("dbp3",   &dbp3);
   setptrn("dbp4",   &dbp4);
   setptrn("hbp1",   &hbp1);
   setptrn("hbp2",   &hbp2);
   setptrn("aii",    &aii);
   setptrn("a17",    &a17);
   setptrn("am",     &am);
   setptrn("am2",    &am2);
   setptrn("am3",    &am3);
   setptrn("am4",    &am4);
   setptrn("amh",    &amh);
   setptrn("amh2",   &amh2);
   setptrn("ams",    &ams);
   setptrn("amhs",   &amhs);
   setptrn("sbac",   &sbac);
   setptrn("gca",    &gca);
   setptrn("gcb",    &gcb);
   setptrn("dsgc",   &dsgc);
   setptrn("gcaoff", &gcaoff);
   setptrn("gcboff", &gcboff);

   setptrn("xglut",    &xglut);		/* generic glutamate response */
   setptrn("xampa",    &xampa);		/* AMPA (type 1) synaptic response */
   setptrn("xampa1",   &xampa1);	/* AMPA type 1 synaptic response */
   setptrn("xampa2",   &xampa2);	/* AMPA type 2 synaptic response */
   setptrn("xampa3",   &xampa3);	/* AMPA type 3 synaptic response */
   setptrn("xampa4",   &xampa4);	/* AMPA type 4 synaptic response */
   setptrn("xampa5",   &xampa5);	/* AMPA type 5 synaptic response */
   setptrn("xnmda",    &xnmda);		/* NMDA type 1 synaptic response */
   setptrn("xnmda2",   &xnmda2);	/* NMDA type 2 synaptic response */
   setptrn("xkainate", &xkainate);	/* Kainate synaptic response */
   setptrn("xmglur6",  &xmglur6);	/* mGluR6 synaptic response */
   setptrn("xgaba",    &xgaba);		/* GABA (type 1) synaptic response */
   setptrn("xgaba1",   &xgaba1);	/* GABA type 1 synaptic response */
   setptrn("xgaba2",   &xgaba2);	/* GABA type 2 synaptic response */
   setptrn("xgaba3",   &xgaba3);	/* GABA type 3 synaptic response */
   setptrn("xgaba4",   &xgaba4);	/* GABA type 4 synaptic response */
   setptrn("xgly",     &xgly);		/* Glycine synaptic response */
   setptrn("xgapj",    &xgapj);		/* gap junction synaptic response */
   setptrn("xdyad",    &xdyad);		/* dyad synapse (uses other resp type) */

   setptrn("black",	&black);	/* set colors for fread() */
   setptrn("blue",	&blue);
   setptrn("green",	&green);
   setptrn("cyan",	&cyan);
   setptrn("red",	&red);
   setptrn("magenta",	&magenta);
   setptrn("brown",	&brown);
   setptrn("white",	&white);

   setptrn("gray",	&gray);
   setptrn("ltblue",	&ltblue);
   setptrn("ltgreen",	&ltgreen);
   setptrn("ltcyan",	&ltcyan);
   setptrn("ltred",	&ltred);
   setptrn("ltmag",	&ltmag);
   setptrn("yellow",	&yellow);
   setptrn("brtwht",	&brtwht);


#include "nval_var_set.cc"

   setptrn("DEND",      &DEND); 
   setptrn("DENDP",     &DENDP); 
   setptrn("DEND_PROX", &DEND_PROX);
   setptrn("DENDD",     &DENDD);
   setptrn("DEND_DIST", &DEND_DIST);
   setptrn("PSOMA",     &PSOMA);
   setptrn("SOMA",      &SOMA);
   setptrn("HCK",    	&HCK);
   setptrn("HILLOCK",   &HILLOCK);
   setptrn("AXONT",     &AXONT);
   setptrn("AXON_THIN", &AXON_THIN);
   setptrn("AXON",      &AXON);
   setptrn("AXONP",     &AXONP);
   setptrn("AXON_PROX", &AXON_PROX);
   setptrn("AXOND",     &AXOND);
   setptrn("AXON_DIST", &AXON_DIST);
   setptrn("VARIC",     &VARIC);
   setptrn("VARICOS",   &VARICOS);

   setptrn("R1",      &R1);
   setptrn("R2",      &R2);
   setptrn("R3",      &R3);
   setptrn("R4",      &R4);
   setptrn("R5",      &R5);
   setptrn("R6",      &R6);
   setptrn("R7",      &R7);
   setptrn("R8",      &R8);
   setptrn("R9",      &R9);
   setptrn("R10",     &R10);
   setptrn("NREGIONS",  &NREGIONS);

   setptrn("EGRAD",   &EGRAD);
   setptrn("LGRAD",   &LGRAD);

   setptrn ("_NA0", &_NA0);
   setptrn ("_NA1", &_NA1);
   setptrn ("_NA2", &_NA2);
   setptrn ("_NA3", &_NA3);
   setptrn ("_NA4", &_NA4);
   setptrn ("_NA5", &_NA5);
   setptrn ("_NA6", &_NA6);
   setptrn ("_NA8", &_NA8);

   setptrn ("_NA",  &_NA);

   setptrn ("_K0",  &_K0);
   setptrn ("_K1",  &_K1);
   setptrn ("_K2",  &_K2);
   setptrn ("_K3",  &_K3);
   setptrn ("_K4",  &_K4);
   
   // joesterle begin
   setptrn ("_HCN1",  &_HCN1);
   setptrn ("_HCN2",  &_HCN2);
   setptrn ("_HCN3",  &_HCN3);
   setptrn ("_HCN4",  &_HCN4);
   // joesterle end
   
   setptrn ("_K5",  &_K5);
   setptrn ("_K6",  &_K6);
   setptrn ("_K7",  &_K7);

   setptrn ("_K",   &_K);
   setptrn ("_KA",  &_KA);
   setptrn ("_KH",  &_K4);

   setptrn ("_KDR", &_KDR);
   setptrn ("_KIR", &_KIR);
   setptrn ("_KHZ", &_KHZ);
   setptrn ("_KV3", &_KV3);

   setptrn ("_KCA0", &_KCA0);
   setptrn ("_KCA1", &_KCA1);
   setptrn ("_KCA2", &_KCA2);
   setptrn ("_KCA3", &_KCA3);
   setptrn ("_KCA4", &_KCA4);
   setptrn ("_KCA5", &_KCA5);
   setptrn ("_KCA6", &_KCA6);

   setptrn ("_BKCA",  &_BKCA);
   setptrn ("_SKCA1", &_SKCA1);
   setptrn ("_SKCA2", &_SKCA2);

   setptrn ("_CLCA1", &_CLCA1);
   setptrn ("_CLCA2", &_CLCA2);
   setptrn ("_CLCA",  &_CLCA);
   setptrn ("_CLCAC", &_CLCAC);

   setptrn ("_CA0", &_CA0);
   setptrn ("_CA1", &_CA1);
   setptrn ("_CA2", &_CA2);
   setptrn ("_CA3", &_CA3);
   setptrn ("_CA4", &_CA4);
   setptrn ("_CA5", &_CA5);
   setptrn ("_CA6", &_CA6);
   setptrn ("_CA7", &_CA7);

   // joesterle
   setptrn ("_MGLUR", &_MGLUR);
   // joesterle
   
   setptrn ("_NMDA1", &_NMDA1);
   setptrn ("_NMDA2", &_NMDA2);

   setptrn ("_CA", &_CA);
   setptrn ("_CA_L", &_CA_L);
   setptrn ("_CA_T", &_CA_T);

   setptrn ("_CAP",    &_CAP);
   setptrn ("_CAPK",   &_CAPK);
   setptrn ("_CABV",   &_CABV);
   setptrn ("_CABK",   &_CABK);
   setptrn ("_CABT",   &_CABT);
   setptrn ("_CABI",   &_CABI);
   setptrn ("_CASH",   &_CASH);
   setptrn ("_CAE",    &_CAE);
   setptrn ("_CAPUMP", &_CAPUMP);
   setptrn ("_CAEXCH", &_CAEXCH);

   setptrn ("_CAS",    &_CAS);
   setptrn ("_VM2",    &_VM2);
   setptrn ("_VM3",    &_VM3);
   setptrn ("_KFCICR", &_KFCICR);
   setptrn ("_KFC",    &_KFC);
   setptrn ("_KACICR", &_KACICR);
   setptrn ("_KAC",    &_KAC);
   setptrn ("_KRCICR", &_KRCICR);
   setptrn ("_K2C",    &_K2C);
   setptrn ("_K2CICR", &_K2CICR);
   setptrn ("_NHC",    &_NHC);
   setptrn ("_NCICR",  &_NCICR);
   setptrn ("_MHC",    &_MHC);
   setptrn ("_MCICR",  &_MCICR);
   setptrn ("_PHC",    &_PHC);
   setptrn ("_PCICR",  &_PCICR);

   setptrn ("_CAS2",   &_CAS2);
   setptrn ("_IP3I",   &_IP3I);
   setptrn ("_BIP3",   &_BIP3);
   setptrn ("_VIP3",   &_VIP3);
   setptrn ("_C1CICR", &_C1CICR);
   setptrn ("_V3IP3",  &_V3IP3);
   setptrn ("_V2IP3",  &_V2IP3);
   setptrn ("_V4IP3",  &_V4IP3);
   setptrn ("_D1IP3",  &_D1IP3);
   setptrn ("_D2IP3",  &_D2IP3);
   setptrn ("_D3IP3",  &_D3IP3);
   setptrn ("_D4IP3",  &_D4IP3);
   setptrn ("_K3IP3",  &_K3IP3);
   setptrn ("_A2IP3",  &_A2IP3);
   setptrn ("_A3IP3",  &_A3IP3);
   setptrn ("_B2IP3",  &_B2IP3);

   setptrn ("_VSTART", &_VSTART);
   setptrn ("_VST",    &_VST);
   setptrn ("_VRV",    &_VRV  );
   setptrn ("_VREV",   &_VREV );
   setptrn ("_RM",     &_RM   );
   setptrn ("_RI",     &_RI   );
   setptrn ("_CM",     &_CM   );
   setptrn ("_DIA",    &_DIA  );
   setptrn ("_CPLAM",  &_CPLAM);
   setptrn ("_CPL",    &_CPL);
   setptrn ("_CMUL",   &_CMUL);
   setptrn ("_COLOR",  &_COLOR);
   setptrn ("_COL",    &_COL);

   setchanptrs();

   setptr ("na1nois",  &na1nois);
   setptr ("na2nois",  &na2nois);
   setptr ("na3nois",  &na3nois);
   setptr ("na4nois",  &na4nois);
   setptr ("na5nois",  &na5nois);
   setptr ("na6nois",  &na6nois);
   setptr ("na8nois",  &na8nois);
   setptr ("k1nois",   &k1nois);
   setptr ("k3nois",   &k3nois);
   setptr ("k4nois",   &k4nois);
   setptr ("k5nois",   &k5nois);
   setptr ("k6nois",   &k6nois);
   setptr ("k7nois",   &k7nois);
   // begin joesterle
   setptr ("hcn1nois",   &hcn1nois);
   setptr ("hcn2nois",   &hcn2nois);
   setptr ("hcn3nois",   &hcn3nois);
   setptr ("hcn4nois",   &hcn4nois);
   // end joesterle
   setptr ("kca1nois", &kca1nois);
   setptr ("kca3nois", &kca3nois);
   setptr ("kca4nois", &kca4nois);
   setptr ("kca5nois", &kca5nois);
   setptr ("kca6nois", &kca6nois);
   setptr ("ca1nois",  &ca1nois);
   setptr ("ca3nois",  &ca3nois);
   setptr ("ca5nois",  &ca5nois);
   setptr ("ca6nois",  &ca6nois);
   setptr ("ca7nois",  &ca7nois);
   setptr ("clcanois", &clcanois);

   setptr("bp_dend1",   &bp_dend1);

   setptr("hz_dend1",   &hz_dend1);
   setptr("hz_dend2",   &hz_dend2);
   setptr("hz_dend3",   &hz_dend3);
   setptr("hz_dend4",   &hz_dend4);
   setptr("hz_dend5",   &hz_dend5);

   setptr("aii_dend1",   &aii_dend1);
   setptr("aii_dend2",   &aii_dend2);
   setptr("aii_dend3",   &aii_dend3);
   setptr("aii_dend4",   &aii_dend4);
   setptr("aii_dend5",   &aii_dend5);

   setptr("sbac_dend1",   &sbac_dend1);
   setptr("sbac_dend2",   &sbac_dend2);
   setptr("sbac_dend3",   &sbac_dend3);
   setptr("sbac_dend4",   &sbac_dend4);
   setptr("sbac_dend5",   &sbac_dend5);

   setptr("am_dend1",   &am_dend1);
   setptr("am_dend2",   &am_dend2);
   setptr("am_dend3",   &am_dend3);
   setptr("am_dend4",   &am_dend4);
   setptr("am_dend5",   &am_dend5);

   setptr("gc_dend1",   &gc_dend1);
   setptr("gc_dend2",   &gc_dend2);
   setptr("gc_dend3",   &gc_dend3);
   setptr("gc_dend4",   &gc_dend4);
   setptr("gc_dend5",   &gc_dend5);

   setptr("gctheta",   &gctheta);
   setptr("sbtheta",   &sbtheta);

   setptr("nvalfile",     &nvalfile);
   setptr("ha_file",      &ha_file);
   setptr("hb_file",      &hb_file);
   setptr("hbat_file",    &hbat_file);
   setptr("rbp_file",     &rbp_file);
   setptr("dbp1_file",    &dbp1_file);
   setptr("dbp2_file",    &dbp2_file);
   setptr("dbp3_file",    &dbp3_file);
   setptr("dbp4_file",    &dbp4_file);
   setptr("hbp1_file",    &hbp1_file);
   setptr("hbp2_file",    &hbp2_file);
   setptr("am_file",      &am_file);
   setptr("am2_file",     &am2_file);
   setptr("am3_file",     &am3_file);
   setptr("am4_file",     &am4_file);
   setptr("amh_file",     &amh_file);
   setptr("amh2_file",    &amh2_file);
   setptr("sbac_file",    &sbac_file);
   setptr("aii_file",     &aii_file);
   setptr("a17_file",     &a17_file);
   setptr("gca_file",     &gca_file);
   setptr("gcb_file",     &gcb_file);
   setptr("dsgc_file",    &dsgc_file);
   setptr("gcaoff_file",  &gcaoff_file);
   setptr("gcboff_file",  &gcboff_file);
   setptr("def_densfile", &def_densfile);

   setptr("confdir",     &confdir);
   setptr("gca_morph",   &gca_morph);
   setptr("gcb_morph",   &gcb_morph);
   setptr("dsgc_morph",  &dsgc_morph);
   setptr("gcaoff_morph",&gcaoff_morph);
   setptr("gcboff_morph",&gcboff_morph);
   setptr("sb_morph",    &sb_morph);
   setptr("rbp_morph",   &rbp_morph);
   setptr("aii_morph",   &aii_morph);
   setptr("a17_morph",   &a17_morph);
   setptr("am_morph",    &am_morph);
   setptr("am2_morph",   &am2_morph);
   setptr("am3_morph",   &am3_morph);
   setptr("am4_morph",   &am4_morph);
   setptr("amh_morph",   &amh_morph);
   setptr("amh2_morph",  &amh2_morph);

   setptr("gca_biophys",    &gca_biophys);
   setptr("gcb_biophys",    &gcb_biophys);
   setptr("dsgc_biophys",  &dsgc_biophys);
   setptr("gcaoff_biophys", &gcaoff_biophys);
   setptr("gcboff_biophys", &gcboff_biophys);
   setptr("am_biophys",    &am_biophys);
   setptr("am2_biophys",   &am2_biophys);
   setptr("am3_biophys",   &am3_biophys);
   setptr("am4_biophys",   &am4_biophys);
   setptr("amh_biophys",   &amh_biophys);
   setptr("amh2_biophys",  &amh2_biophys);
   setptr("rbp_biophys",   &rbp_biophys);
   setptr("aii_biophys",   &aii_biophys);
   setptr("a17_biophys",   &a17_biophys);
   setptr("sb_biophys",    &sb_biophys);
   setptr("Chnoise",       &Chnoise);

   setptr("dbp1_thetax",   &dbp1_thetax);
   setptr("dbp1_thetay",   &dbp1_thetay);
   setptr("dbp1_thetaz",   &dbp1_thetaz);
   setptr("dbp1_flip",     &dbp1_flip);
   setptr("dbp2_thetax",   &dbp2_thetax);
   setptr("dbp2_thetay",   &dbp2_thetay);
   setptr("dbp2_thetaz",   &dbp2_thetaz);
   setptr("dbp2_flip",     &dbp2_flip);
   setptr("dbp3_thetax",   &dbp3_thetax);
   setptr("dbp3_thetay",   &dbp3_thetay);
   setptr("dbp3_thetaz",   &dbp3_thetaz);
   setptr("dbp3_flip",     &dbp3_flip);
   setptr("dbp4_thetax",   &dbp4_thetax);
   setptr("dbp4_thetay",   &dbp4_thetay);
   setptr("dbp4_thetaz",   &dbp4_thetaz);
   setptr("dbp4_flip",     &dbp4_flip);
   setptr("hbp1_thetax",   &hbp1_thetax);
   setptr("hbp1_thetay",   &hbp1_thetay);
   setptr("hbp1_thetaz",   &hbp1_thetaz);
   setptr("hbp1_flip",     &hbp1_flip);
   setptr("hbp2_thetax",   &hbp2_thetax);
   setptr("hbp2_thetay",   &hbp2_thetay);
   setptr("hbp2_thetaz",   &hbp2_thetaz);
   setptr("hbp2_flip",     &hbp2_flip);
   setptr("rbp_thetax",    &rbp_thetax);
   setptr("rbp_thetay",    &rbp_thetay);
   setptr("rbp_thetaz",    &rbp_thetaz);
   setptr("rbp_flip",      &rbp_flip);
   setptr("ha_thetax",     &ha_thetax);
   setptr("ha_thetay",     &ha_thetay);
   setptr("ha_thetaz",     &ha_thetaz);
   setptr("ha_flip",       &ha_flip);
   setptr("hb_thetax",     &hb_thetax);
   setptr("hb_thetay",     &hb_thetay);
   setptr("hb_thetaz",     &hb_thetaz);
   setptr("hb_flip",       &hb_flip);
   setptr("hbat_thetax",   &hbat_thetax);
   setptr("hbat_thetay",   &hbat_thetay);
   setptr("hbat_thetaz",   &hbat_thetaz);
   setptr("hbat_flip",     &hbat_flip);
   setptr("am_thetax",     &am_thetax);
   setptr("am_thetay",     &am_thetay);
   setptr("am_thetaz",     &am_thetaz);
   setptr("am_flip",       &am_flip);
   setptr("am2_thetax",    &am2_thetax);
   setptr("am2_thetay",    &am2_thetay);
   setptr("am2_thetaz",    &am2_thetaz);
   setptr("am2_flip",      &am2_flip);
   setptr("am3_thetax",    &am3_thetax);
   setptr("am3_thetay",    &am3_thetay);
   setptr("am3_thetaz",    &am3_thetaz);
   setptr("am3_flip",      &am3_flip);
   setptr("am4_thetax",    &am4_thetax);
   setptr("am4_thetay",    &am4_thetay);
   setptr("am4_thetaz",    &am4_thetaz);
   setptr("am4_flip",      &am4_flip);
   setptr("ams_thetax",    &ams_thetax);
   setptr("ams_thetay",    &ams_thetay);
   setptr("ams_thetaz",    &ams_thetaz);
   setptr("ams_flip",      &ams_flip);
   setptr("amh_thetax",    &amh_thetax);
   setptr("amh_thetay",    &amh_thetay);
   setptr("amh_thetaz",    &amh_thetaz);
   setptr("amh_flip",      &amh_flip);
   setptr("amh2_thetax",   &amh2_thetax);
   setptr("amh2_thetay",   &amh2_thetay);
   setptr("amh2_thetaz",   &amh2_thetaz);
   setptr("amh2_flip",     &amh2_flip);
   setptr("aii_thetax",    &aii_thetax);
   setptr("aii_thetay",    &aii_thetay);
   setptr("aii_thetaz",    &aii_thetaz);
   setptr("aii_flip",      &aii_flip);
   setptr("a17_thetax",    &a17_thetax);
   setptr("a17_thetay",    &a17_thetay);
   setptr("a17_thetaz",    &a17_thetaz);
   setptr("a17_flip",      &a17_flip);
   setptr("sbac_thetax",   &sbac_thetax);
   setptr("sbac_thetay",   &sbac_thetay);
   setptr("sbac_thetaz",   &sbac_thetaz);
   setptr("sbac_flip",     &sbac_flip);
   setptr("dsgc_thetax",   &dsgc_thetax);
   setptr("dsgc_thetay",   &dsgc_thetay);
   setptr("dsgc_thetaz",   &dsgc_thetaz);
   setptr("dsgc_flip",     &dsgc_flip);
   setptr("gc_thetax",     &gc_thetax);
   setptr("gc_thetay",     &gc_thetay);
   setptr("gc_thetaz",     &gc_thetaz);
   setptr("gc_flip",       &gc_flip);

   setptr("flip",          &flip);

   setptr("dvrev", &dvrev);
   setptr("dvst", &dvst);

   setptrn("n_rods",   &n_rods);
   setptrn("n_cones",  &n_cones);
   setptrn("n_rbp",    &n_rbp);
   setptrn("n_ha",     &n_ha);
   setptrn("n_hb",     &n_hb);
   setptrn("n_hbat",   &n_hbat);
   setptrn("n_dbp1",   &n_dbp1);
   setptrn("n_dbp2",   &n_dbp2);
   setptrn("n_dbp3",   &n_dbp3);
   setptrn("n_dbp4",   &n_dbp4);
   setptrn("n_hbp1",   &n_hbp1);
   setptrn("n_hbp2",   &n_hbp2);
   setptrn("n_aii",    &n_aii);
   setptrn("n_a17",    &n_a17);
   setptrn("n_am",     &n_am);
   setptrn("n_am2",    &n_am2);
   setptrn("n_am3",    &n_am3);
   setptrn("n_am4",    &n_am4);
   setptrn("n_amh",    &n_amh);
   setptrn("n_amh2",   &n_amh2);
   setptrn("n_ams",    &n_ams);
   setptrn("n_amhs",   &n_amhs);
   setptrn("n_sbac",   &n_sbac);
   setptrn("n_gca",    &n_gca);
   setptrn("n_gcb",    &n_gcb);
   setptrn("n_dsgc",   &n_dsgc);
   setptrn("n_gcaoff",  &n_gcaoff);
   setptrn("n_gcboff",  &n_gcboff);

   setptr("maxha",   &maxha);
   setptr("maxhb",   &maxhb);
   setptr("maxhbat", &maxhbat);
   setptr("maxrbp",  &maxrbp);
   setptr("maxdbp1", &maxdbp1);
   setptr("maxdbp2", &maxdbp2);
   setptr("maxdbp3", &maxdbp3);
   setptr("maxdbp4", &maxdbp4);
   setptr("maxhbp1", &maxhbp1);
   setptr("maxhbp2", &maxhbp2);
   setptr("maxa17",  &maxa17);
   setptr("maxaii",  &maxaii);
   setptr("maxsb",   &maxsb);
   setptr("maxam",   &maxam);
   setptr("maxam2",  &maxam2);
   setptr("maxam3",  &maxam3);
   setptr("maxam4",  &maxam4);
   setptr("maxamh",  &maxamh);
   setptr("maxamh2", &maxamh2);
   setptr("maxams",  &maxams);
   setptr("maxamhs", &maxamhs);
   setptr("maxgca",  &maxgca);
   setptr("maxgcb",  &maxgcb);
   setptr("maxdsgc", &maxdsgc);
   setptr("maxgcaoff",&maxgcaoff);
   setptr("maxgcboff",&maxgcboff);

   setptr("node_scale",    &node_scale);

   setptr("rod_nscale",    &rod_nscale);
   setptr("cone_nscale",   &cone_nscale);
   setptr("hbat_nscale",   &hbat_nscale);
   setptr("ha_nscale",     &ha_nscale);
   setptr("hb_nscale",     &hb_nscale);
   setptr("rbp_nscale",    &rbp_nscale);
   setptr("dbp1_nscale",   &dbp1_nscale);
   setptr("dbp2_nscale",   &dbp2_nscale);
   setptr("dbp3_nscale",   &dbp3_nscale);
   setptr("dbp4_nscale",   &dbp4_nscale);
   setptr("hbp1_nscale",   &hbp1_nscale);
   setptr("hbp2_nscale",   &hbp2_nscale);
   setptr("aii_nscale",    &aii_nscale);
   setptr("sbac_nscale",   &sbac_nscale);
   setptr("am_nscale",     &am_nscale);
   setptr("am2_nscale",    &am2_nscale);
   setptr("am3_nscale",    &am3_nscale);
   setptr("am4_nscale",    &am4_nscale);
   setptr("amh_nscale",    &amh_nscale);
   setptr("amh2_nscale",   &amh2_nscale);
   setptr("ams_nscale",    &ams_nscale);
   setptr("amhs_nscale",   &amhs_nscale);
   setptr("a17_nscale",    &a17_nscale);
   setptr("gca_nscale",    &gca_nscale);
   setptr("gcb_nscale",    &gcb_nscale);
   setptr("dsgc_nscale",   &dsgc_nscale);
   setptr("gcaoff_nscale", &gcaoff_nscale);
   setptr("gcboff_nscale", &gcboff_nscale);

   setptrn("synap_scale",   &synap_scale);
   setptrn("gj_scale",      &gj_scale);

   setptrn("make_rods",     &make_rods);
   setptrn("make_cones",    &make_cones);
   setptrn("make_hbat",     &make_hbat);
   setptrn("make_ha",       &make_ha);
   setptrn("make_hb",       &make_hb);
   setptrn("make_rbp",      &make_rbp);
   setptrn("make_dbp1",     &make_dbp1);
   setptrn("make_dbp2",     &make_dbp2);
   setptrn("make_dbp3",     &make_dbp3);
   setptrn("make_dbp4",     &make_dbp4);
   setptrn("make_hbp1",     &make_hbp1);
   setptrn("make_hbp2",     &make_hbp2);
   setptrn("make_aii",      &make_aii);
   setptrn("make_sbac",     &make_sbac);
   setptrn("make_am",       &make_am);
   setptrn("make_am2",      &make_am2);
   setptrn("make_am3",      &make_am3);
   setptrn("make_am4",      &make_am4);
   setptrn("make_amh",      &make_amh);
   setptrn("make_amh2",     &make_amh2);
   setptrn("make_ams",      &make_ams);
   setptrn("make_amhs",     &make_amhs);
   setptrn("make_a17",      &make_a17);
   setptrn("make_gca",      &make_gca);
   setptrn("make_gcb",      &make_gcb);
   setptrn("make_dsgc",     &make_dsgc);
   setptrn("make_gcaoff",   &make_gcaoff);
   setptrn("make_gcboff",   &make_gcboff);

   setptrn("make_cone_cone", &make_cone_cone);
   setptrn("make_rod_rbp",   &make_rod_rbp);
   setptrn("make_cone_dbp1", &make_cone_dbp1);
   setptrn("make_cone_dbp2", &make_cone_dbp2);
   setptrn("make_cone_dbp3", &make_cone_dbp3);
   setptrn("make_cone_dbp4", &make_cone_dbp4);
   setptrn("make_cone_hbp1", &make_cone_hbp1);
   setptrn("make_cone_hbp2", &make_cone_hbp2);
   setptrn("make_rod_hbat",  &make_rod_hbat);
   setptrn("make_cone_ha",   &make_cone_ha);
   setptrn("make_cone_hb",   &make_cone_hb);
   setptrn("make_ha_ha",     &make_ha_ha);
   setptrn("make_hb_hb",     &make_hb_hb);
   setptrn("make_rbp_aii",   &make_rbp_aii);
   setptrn("make_aii_aii",   &make_aii_aii);
   setptrn("make_dbp1_aii",  &make_dbp1_aii);
   setptrn("make_dbp2_aii",  &make_dbp2_aii);
   setptrn("make_dbp3_aii",  &make_dbp3_aii);
   setptrn("make_dbp4_aii",  &make_dbp4_aii);
   setptrn("make_aii_hbp1",  &make_aii_hbp1);
   setptrn("make_aii_hbp2",  &make_aii_hbp2);
   setptrn("make_rbp_a17",   &make_rbp_a17);
   setptrn("make_a17_rbp",   &make_a17_rbp);
   setptrn("make_dbp1_sbac", &make_dbp1_sbac);
   setptrn("make_dbp2_sbac", &make_dbp2_sbac);
   setptrn("make_dbp3_sbac", &make_dbp3_sbac);
   setptrn("make_dbp4_sbac", &make_dbp4_sbac);
   setptrn("make_dbp1_dbp1", &make_dbp1_dbp1);
   setptrn("make_dbp2_dbp2", &make_dbp3_dbp3);
   setptrn("make_dbp4_dbp4", &make_dbp4_dbp4);
   setptrn("make_hbp1_hbp1", &make_hbp1_hbp1);
   setptrn("make_dbp1_am",   &make_dbp1_am);
   setptrn("make_dbp1_am2",  &make_dbp1_am2);
   setptrn("make_dbp1_am3",  &make_dbp1_am3);
   setptrn("make_dbp1_am4",  &make_dbp1_am4);
   setptrn("make_dbp1_amh",  &make_dbp1_amh);
   setptrn("make_dbp1_amh2", &make_dbp1_amh2);
   setptrn("make_dbp1_ams",  &make_dbp1_ams);
   setptrn("make_dbp2_am",   &make_dbp2_am);
   setptrn("make_dbp2_am2",  &make_dbp2_am2);
   setptrn("make_dbp2_am3",  &make_dbp2_am3);
   setptrn("make_dbp2_am4",  &make_dbp2_am4);
   setptrn("make_dbp2_amh",  &make_dbp2_amh);
   setptrn("make_dbp2_amh2", &make_dbp2_amh2);
   setptrn("make_dbp2_ams",  &make_dbp2_ams);
   setptrn("make_dbp3_am",   &make_dbp3_am);
   setptrn("make_dbp3_am2",  &make_dbp3_am2);
   setptrn("make_dbp3_am3",  &make_dbp3_am3);
   setptrn("make_dbp3_am4",  &make_dbp3_am4);
   setptrn("make_dbp3_amh",  &make_dbp3_amh);
   setptrn("make_dbp3_amh2", &make_dbp3_amh2);
   setptrn("make_dbp3_ams",  &make_dbp3_ams);
   setptrn("make_dbp4_am",   &make_dbp4_am);
   setptrn("make_dbp4_am2",  &make_dbp4_am2);
   setptrn("make_dbp4_am3",  &make_dbp4_am3);
   setptrn("make_dbp4_am4",  &make_dbp4_am4);
   setptrn("make_dbp4_amh",  &make_dbp4_amh);
   setptrn("make_dbp4_amh2", &make_dbp4_amh2);
   setptrn("make_dbp4_ams",  &make_dbp4_ams);
   setptrn("make_dbp1_gca",  &make_dbp1_gca);
   setptrn("make_dbp1_gcb",  &make_dbp1_gcb);
   setptrn("make_dbp2_gca",  &make_dbp2_gca);
   setptrn("make_dbp3_gca",  &make_dbp3_gca);
   setptrn("make_dbp4_gca",  &make_dbp4_gca);
   setptrn("make_dbp2_gcb",  &make_dbp2_gcb);
   setptrn("make_dbp3_gcb",  &make_dbp3_gcb);
   setptrn("make_dbp4_gcb",  &make_dbp4_gcb);
   setptrn("make_sbac_dsgc", &make_sbac_dsgc);
   setptrn("make_sbac_sbac", &make_sbac_sbac);
   setptrn("make_sbac_dbp1", &make_sbac_dbp1);
   setptrn("make_sbac_dbp2", &make_sbac_dbp2);
   setptrn("make_sbac_dbp3", &make_sbac_dbp3);
   setptrn("make_sbac_dbp4", &make_sbac_dbp4);
   setptrn("make_dbp1_dsgc", &make_dbp1_dsgc);
   setptrn("make_dbp2_dsgc", &make_dbp2_dsgc);
   setptrn("make_dbp3_dsgc", &make_dbp3_dsgc);
   setptrn("make_dbp4_dsgc", &make_dbp4_dsgc);
   setptrn("make_hbp1_sbac", &make_hbp1_sbac);
   setptrn("make_hbp1_dsgc", &make_hbp1_dsgc);
   setptrn("make_hbp2_dsgc", &make_hbp2_dsgc);
   setptrn("make_hbp1_amh",  &make_hbp1_amh);
   setptrn("make_hbp1_amh2", &make_hbp1_amh2);
   setptrn("make_hbp1_amhs", &make_hbp1_amhs);
   setptrn("make_hbp2_sbac", &make_hbp2_sbac);
   setptrn("make_hbp2_amh",  &make_hbp2_amh);
   setptrn("make_hbp2_amh2", &make_hbp2_amh2);
   setptrn("make_hbp2_amhs", &make_hbp2_amhs);
   setptrn("make_ams_dbp1",  &make_ams_dbp1);
   setptrn("make_ams_dbp2",  &make_ams_dbp2);
   setptrn("make_ams_dbp3",  &make_ams_dbp3);
   setptrn("make_ams_dbp4",  &make_ams_dbp4);
   setptrn("make_amhs_hbp1", &make_amhs_hbp1);
   setptrn("make_ams_dsgc",  &make_ams_dsgc);
   setptrn("make_amhs_dsgc", &make_amhs_dsgc);
   setptrn("make_dsgc_dsgc", &make_dsgc_dsgc);
   setptrn("make_am_dbp1",   &make_am_dbp1);
   setptrn("make_am2_dbp1",  &make_am2_dbp1);
   setptrn("make_am3_dbp1",  &make_am3_dbp1);
   setptrn("make_am4_dbp1",  &make_am4_dbp1);
   setptrn("make_am_dbp2",   &make_am_dbp2);
   setptrn("make_am2_dbp2",  &make_am2_dbp2);
   setptrn("make_am3_dbp2",  &make_am3_dbp2);
   setptrn("make_am4_dbp2",  &make_am4_dbp2);
   setptrn("make_am_dbp3",   &make_am_dbp3);
   setptrn("make_am2_dbp3",  &make_am2_dbp3);
   setptrn("make_am3_dbp3",  &make_am3_dbp3);
   setptrn("make_am4_dbp3",  &make_am4_dbp3);
   setptrn("make_am_dbp4",   &make_am_dbp4);
   setptrn("make_am2_dbp4",  &make_am2_dbp4);
   setptrn("make_am3_dbp4",  &make_am3_dbp4);
   setptrn("make_am4_dbp4",  &make_am4_dbp4);
   setptrn("make_am_gca",    &make_am_gca);
   setptrn("make_am2_gca",   &make_am2_gca);
   setptrn("make_am3_gca",   &make_am3_gca);
   setptrn("make_am4_gca",   &make_am4_gca);
   setptrn("make_am_am",     &make_am_am);
   setptrn("make_am2_am2",   &make_am2_am2);
   setptrn("make_am3_am3",   &make_am3_am3);
   setptrn("make_am4_am4",   &make_am4_am4);
   setptrn("make_am_gcaoff", &make_am_gcaoff);
   setptrn("make_amh_hbp1",  &make_amh_hbp1);
   setptrn("make_amh2_hbp1", &make_amh2_hbp1);
   setptrn("make_amh_hbp2",  &make_amh_hbp2);
   setptrn("make_amh2_hbp2", &make_amh2_hbp2);
   setptrn("make_amh2_amh",  &make_amh2_amh);
   setptrn("make_amh_dsgc",  &make_amh_dsgc);
   setptrn("make_amh_gcaoff",&make_amh_gcaoff);
   setptrn("make_amh_gcboff",&make_amh_gcboff);
   setptrn("make_amh_amh",   &make_amh_amh);
   setptrn("make_amh2_amh2", &make_amh2_amh2);
   setptrn("make_am_dsgc",   &make_am_dsgc);
   setptrn("make_ams_gca",   &make_ams_gca);
   setptrn("make_ams_gcb",   &make_ams_gcb);
   setptrn("make_hbp1_gcaoff",&make_hbp1_gcaoff);
   setptrn("make_hbp2_gcaoff",&make_hbp2_gcaoff);
   setptrn("make_amhs_gcaoff",&make_amhs_gcaoff);
   setptrn("make_hbp1_gcboff",&make_hbp1_gcboff);
   setptrn("make_hbp2_gcboff",&make_hbp2_gcboff);
   setptrn("make_amhs_gcboff",&make_amhs_gcboff);
   setptrn("make_ha_cone",   &make_ha_cone);
   setptrn("make_ha_dbp1",   &make_ha_dbp1);
   setptrn("make_ha_dbp2",   &make_ha_dbp2);
   setptrn("make_ha_dbp3",   &make_ha_dbp3);
   setptrn("make_ha_dbp4",   &make_ha_dbp4);
   setptrn("make_ha_hbp1",   &make_ha_hbp1);
   setptrn("make_ha_hbp2",   &make_ha_hbp2);
   setptrn("make_hb_dbp1",   &make_hb_dbp1);
   setptrn("make_hb_dbp2",   &make_hb_dbp2);
   setptrn("make_hb_dbp3",   &make_hb_dbp3);
   setptrn("make_hb_dbp4",   &make_hb_dbp4);
   setptrn("make_hb_hbp1",   &make_hb_hbp1);
   setptrn("make_hb_hbp2",   &make_hb_hbp2);

   setptr("arrsiz",&arrsiz);
   setptr("xarrsiz",&xarrsiz);
   setptr("yarrsiz",&yarrsiz);
   setptr("zarrsiz",&zarrsiz);
   setptr("garrsiz",&garrsiz);
   setptr("gxarrsiz",&gxarrsiz);
   setptr("gyarrsiz",&gyarrsiz);

   setptr("amharrsiz",&amharrsiz);
   setptr("amhxarrsiz",&amhxarrsiz);
   setptr("amhyarrsiz",&amhyarrsiz);
   setptr("amh2arrsiz",&amh2arrsiz);
   setptr("amh2xarrsiz",&amh2xarrsiz);
   setptr("amh2yarrsiz",&amh2yarrsiz);
   setptr("amsarrsiz",&amsarrsiz);
   setptr("amsxarrsiz",&amsxarrsiz);
   setptr("amsyarrsiz",&amsyarrsiz);
   setptr("amhsarrsiz",&amhsarrsiz);
   setptr("amhsxarrsiz",&amhsxarrsiz);
   setptr("amhsyarrsiz",&amhsyarrsiz);
   setptr("amarrsiz",&amarrsiz);
   setptr("amxarrsiz",&amxarrsiz);
   setptr("amyarrsiz",&amyarrsiz);
   setptr("am2arrsiz",&am2arrsiz);
   setptr("am2xarrsiz",&am2xarrsiz);
   setptr("am2yarrsiz",&am2yarrsiz);
   setptr("am3arrsiz",&am3arrsiz);
   setptr("am3xarrsiz",&am3xarrsiz);
   setptr("am3yarrsiz",&am3yarrsiz);
   setptr("am4arrsiz",&am4arrsiz);
   setptr("am4xarrsiz",&am4xarrsiz);
   setptr("am4yarrsiz",&am4yarrsiz);
   setptr("aiiarrsiz",&aiiarrsiz);
   setptr("aiixarrsiz",&aiixarrsiz);
   setptr("aiiyarrsiz",&aiiyarrsiz);
   setptr("sbarrsiz",&sbarrsiz);
   setptr("sbxarrsiz",&sbxarrsiz);
   setptr("sbyarrsiz",&sbyarrsiz);

   setptr("dbp1arrsiz", &dbp1arrsiz);
   setptr("dbp1xarrsiz",&dbp1xarrsiz);
   setptr("dbp1yarrsiz",&dbp1yarrsiz);
   setptr("dbp2arrsiz", &dbp2arrsiz);
   setptr("dbp2xarrsiz",&dbp2xarrsiz);
   setptr("dbp2yarrsiz",&dbp2yarrsiz);
   setptr("dbp3arrsiz", &dbp3arrsiz);
   setptr("dbp3xarrsiz",&dbp3xarrsiz);
   setptr("dbp3yarrsiz",&dbp3yarrsiz);
   setptr("dbp4arrsiz", &dbp4arrsiz);
   setptr("dbp4xarrsiz",&dbp4xarrsiz);
   setptr("dbp4yarrsiz",&dbp4yarrsiz);
   setptr("hbp1arrsiz", &hbp1arrsiz);
   setptr("hbp1xarrsiz",&hbp1xarrsiz);
   setptr("hbp1yarrsiz",&hbp1yarrsiz);
   setptr("hbp2arrsiz", &hbp2arrsiz);
   setptr("hbp2xarrsiz",&hbp2xarrsiz);
   setptr("hbp2yarrsiz",&hbp2yarrsiz);

   setptr("conearrsiz",&conearrsiz);
   setptr("conexarrsiz",&conexarrsiz);
   setptr("coneyarrsiz",&coneyarrsiz);
   setptr("rodarrsiz",&rodarrsiz);
   setptr("rodxarrsiz",&rodxarrsiz);
   setptr("rodyarrsiz",&rodyarrsiz);
   setptr("rbparrsiz",&rbparrsiz);
   setptr("rbpxarrsiz",&rbpxarrsiz);
   setptr("rbpyarrsiz",&rbpyarrsiz);
   setptr("arrcentx",&arrcentx);
   setptr("arrcenty",&arrcenty);
   setptr("arrcentz",&arrcentz);
   setptr("rodarrcentx",&rodarrcentx);
   setptr("rodarrcenty",&rodarrcenty);
   setptr("conearrcentx",&conearrcentx);
   setptr("conearrcenty",&conearrcenty);

   setptr("dbp1arrcentx",&dbp1arrcentx);
   setptr("dbp1arrcenty",&dbp1arrcenty);
   setptr("dbp2arrcentx",&dbp2arrcentx);
   setptr("dbp2arrcenty",&dbp2arrcenty);
   setptr("dbp3arrcentx",&dbp3arrcentx);
   setptr("dbp3arrcenty",&dbp3arrcenty);
   setptr("dbp4arrcentx",&dbp4arrcentx);
   setptr("dbp4arrcenty",&dbp4arrcenty);
   setptr("hbp1arrcentx",&hbp1arrcentx);
   setptr("hbp1arrcenty",&hbp1arrcenty);
   setptr("hbp2arrcentx",&hbp2arrcentx);
   setptr("hbp2arrcenty",&hbp2arrcenty);

   setptr("dbp1_first_cent",&dbp1_first_cent);
   setptr("dbp2_first_cent",&dbp2_first_cent);
   setptr("dbp3_first_cent",&dbp3_first_cent);
   setptr("dbp4_first_cent",&dbp4_first_cent);
   setptr("hbp1_first_cent",&hbp1_first_cent);
   setptr("hbp2_first_cent",&hbp2_first_cent);

   setptr("set_zmax",&set_zmax);			// for elimit 
   setptr("set_zmin",&set_zmin);

   setptr("disp_zmax",&disp_zmax);			// for dispcelltype() 
   setptr("disp_zmin",&disp_zmin);

   setptr("disp_gc_zmax",&disp_gc_zmax);		// for dispcelltype() 
   setptr("disp_gc_zmin",&disp_gc_zmin);

   setptr("disp_dsgc_zmax",&disp_dsgc_zmax);		// for dispcelltype() 
   setptr("disp_dsgc_zmin",&disp_dsgc_zmin);

   setptr("dispsize",&dispsize);
   setptr("disp_calib_len",&disp_calib_len);
   setptr("disp_margin",  &disp_margin);
   setptr("limit_array",  &limit_array);
   setptr("update_array", &update_array);

   setptr("conespac", &conespac);
   setptr("rodspac", &rodspac);

   setptr("cone_maxcond", &cone_maxcond);
   setptr("rod_maxcond", &rod_maxcond);
   setptr("cone_pigm",  &cone_pigm);

   setptr("cone_timec", &cone_timec);
   setptr("rod_timec",  &rod_timec);

   setptr("cone_loopg", &cone_loopg);
   setptr("rod_loopg",  &rod_loopg);

   setptr("plot_freq", &plot_freq);

   setptr("rod_rect_arr", &rod_rect_arr);
   setptr("make_one_dbp1", &make_one_dbp1);
   setptr("make_one_dbp2", &make_one_dbp2);
   setptr("make_one_dbp3", &make_one_dbp3);
   setptr("make_one_dbp4", &make_one_dbp4);
   setptr("gcdistnod",	  &gcdistnod);
   setptr("remove_nconns",&remove_nconns);
   setptr("print_conns",  &print_conns);
   setptr("print_avg_conns", &print_avg_conns);
   setptr("disp_cn", &disp_cn);
   setptr("disp_ct", &disp_ct);

   setptr("ttxbath", &ttxbath);
   setptr("ttxsoma", &ttxsoma);
   setptr("ttxdend", &ttxdend);
   setptr("tea",     &tea);
   setptr("fourap",  &fourap);
   setptr("zd7288",  &zd7288);
   setptr("ibtox",   &ibtox);

   setptr("dia_min_factor", &dia_min_factor);
   setptr("dendd_dia_factor", &dendd_dia_factor);
   setptr("dendi_dia_factor", &dendi_dia_factor);
   setptr("dendp_dia_factor", &dendp_dia_factor);
   setptr("dend_dia_factor", &dend_dia_factor);
   setptr("ax_dia_factor", &ax_dia_factor);
   setptr("cell_dia_factor", &cell_dia_factor);
   setptr("den_dia", &den_dia);
   setptr("taper", &taper);
   setptr("gc_dend_cplam", &gc_dend_cplam);
   setptr("am_dend_cplam", &am_dend_cplam);

   setptr("mxrot", &mxrot);
   setptr("myrot", &myrot);

   setptr("nshell", &nshell);
   setptr("nshell_soma", &nshell_soma);

   setptr("chanparamsfile",  &chanparamsfile);
   setptr("syn_savefile",    &syn_savefile);
   setptr("syn_restorefile", &syn_restorefile);

   setptr("cone_densfile",&cone_densfile);
   setptr("rod_densfile", &rod_densfile);
   setptr("hbat_densfile",&hbat_densfile);
   setptr("ha_densfile",  &ha_densfile);
   setptr("hb_densfile",  &hb_densfile);
   setptr("rbp_densfile", &rbp_densfile);
   setptr("dbp1_densfile",&dbp1_densfile);
   setptr("dbp2_densfile",&dbp2_densfile);
   setptr("dbp3_densfile",&dbp3_densfile);
   setptr("dbp4_densfile",&dbp4_densfile);
   setptr("hbp1_densfile",&hbp1_densfile);
   setptr("hbp2_densfile",&hbp2_densfile);
   setptr("a17_densfile", &a17_densfile);
   setptr("aii_densfile", &aii_densfile);
   setptr("sbac_densfile",&sbac_densfile);
   setptr("am_densfile",  &am_densfile);
   setptr("am2_densfile", &am2_densfile);
   setptr("am3_densfile", &am3_densfile);
   setptr("am4_densfile", &am4_densfile);
   setptr("amh_densfile", &amh_densfile);
   setptr("amh2_densfile",&amh2_densfile);
   setptr("ams_densfile", &ams_densfile);
   setptr("amhs_densfile",&amhs_densfile);
   setptr("gca_densfile", &gca_densfile);
   setptr("gcb_densfile", &gcb_densfile);
   setptr("dsgc_densfile",&dsgc_densfile);
   setptr("gcaoff_densfile",&gcaoff_densfile);
   setptr("gcboff_densfile",&gcboff_densfile);

   setptr("cone_densfile2",&cone_densfile2);
   setptr("rod_densfile2", &rod_densfile2);
   setptr("hbat_densfile2",&hbat_densfile2);
   setptr("ha_densfile2",  &ha_densfile2);
   setptr("hb_densfile2",  &hb_densfile2);
   setptr("rbp_densfile2", &rbp_densfile2);
   setptr("dbp1_densfile2",&dbp1_densfile2);
   setptr("dbp2_densfile2",&dbp2_densfile2);
   setptr("dbp3_densfile2",&dbp3_densfile2);
   setptr("dbp4_densfile2",&dbp4_densfile2);
   setptr("hbp1_densfile2",&hbp1_densfile2);
   setptr("hbp2_densfile2",&hbp2_densfile2);
   setptr("a17_densfile2", &a17_densfile2);
   setptr("aii_densfile2", &aii_densfile2);
   setptr("sbac_densfile2",&sbac_densfile2);
   setptr("am_densfile2",  &am_densfile2);
   setptr("am2_densfile2", &am2_densfile2);
   setptr("am3_densfile2", &am3_densfile2);
   setptr("am4_densfile2", &am4_densfile2);
   setptr("amh_densfile2", &amh_densfile2);
   setptr("amh2_densfile2",&amh2_densfile2);
   setptr("ams_densfile2", &ams_densfile2);
   setptr("amhs_densfile2",&amhs_densfile2);
   setptr("gca_densfile2", &gca_densfile2);
   setptr("gcb_densfile2", &gcb_densfile2);
   setptr("dsgc_densfile2",&dsgc_densfile2);
   setptr("gcaoff_densfile2",&gcaoff_densfile2);
   setptr("gcboff_densfile2",&gcboff_densfile2);

   setptr("bg_inten",&bg_inten);

   setptrn("pnoise",&pnoise);
   setptrn("dnoise",&dnoise);
   setptrn("vnoise",&vnoise);
   setptrn("cnoise",&cnoise);

   setptrn("drm",&drm);			// defs for density files, already init 
   setptrn("dcm",&dcm);
   setptrn("dri",&dri);
   setptrn("vk",&vk);
   setptrn("dvs",&dvs);
   setptrn("vcl",&vcl);
   setptrn("vna",&vna);
   setptrn("gc_rm", &gc_rm);		// must init 	
   setptrn("gc_vr", &gc_vr);		// must init 	
   setptrn("gc_vs", &gc_vs);		// must init 	

   setptr("cbplam", &cbplam);		// local complam for bipolar cell
   setptr("cbplam2", &cbplam2);		// local complam for bipolar cell, compact cell, densfile 2

   setptr("ath_dia",   &ath_dia);		// 
   setptrn("makestim", &makestim);		// =1 -> run stim 
   setptr("plotlabel", &plotlabel);		// 
}

