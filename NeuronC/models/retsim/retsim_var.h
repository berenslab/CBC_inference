/* retsim_var.h */
extern int xcone;		/* Identity numbers for neurons */
extern int xrod;
extern int hbat;
extern int ha;
extern int hb;
extern int rbp;
extern int dbp1;		/* depolarizing cone bipolar */
extern int dbp2;		/* depolarizing cone bipolar */
extern int dbp3;		/* depolarizing cone bipolar */
extern int dbp4;		/* depolarizing cone bipolar */
extern int hbp1;		/* hyperpolarizing cone bipolar */
extern int hbp2;		/* hyperpolarizing cone bipolar */
extern int a17;
extern int aii;
extern int sbac;
extern int am;
extern int am2;
extern int am3;
extern int am4;
extern int amh;
extern int amh2;
extern int ams;
extern int amhs;
extern int gca;
extern int gcb;
extern int dsgc;
extern int gcaoff;
extern int gcboff;
extern int nceltypes;


extern int xglut;
extern int xampa;
extern int xampa1;
extern int xampa2;
extern int xampa3;
extern int xampa4;
extern int xampa5;
extern int xnmda;
extern int xnmda2;
extern int xkainate;
extern int xmglur6;
extern int xgaba;
extern int xgaba1;
extern int xgaba2;
extern int xgaba3;
extern int xgaba4;
extern int xgly;
extern int xgapj;
extern int xdyad;

extern int R1;
extern int R2;
extern int R3;
extern int R4;
extern int R5;
extern int R6;
extern int R7;
extern int R8;
extern int R9;
extern int R10;
extern int EGRAD;
extern int LGRAD;

extern int DEND;
extern int DEND_DIST;
extern int DENDD;
extern int DEND_MED;
extern int DENDM;
extern int DEND_PROX;
extern int DENDP;
extern int SOMA;
extern int PSOMA;
extern int HILLOCK;
extern int HCK;
extern int AXON_THIN;
extern int AXONT;
extern int AXON;
extern int AXON_PROX;
extern int AXONP;
extern int AXON_DIST;
extern int AXOND;
extern int VARICOS;
extern int VARIC;
extern int NREGIONS;

#include "nval_var.h"

extern const char *cname[NCELTYPES];			/* cell names */
extern const char *rname[NRESPTYPES];			/* synaptic types */
extern const char *anatfiles[NCELTYPES];		/* cell anatomy file names */
extern double *cellmorphs[NCELTYPES];			/* morphology files read into arrays */
extern int cell_lines[NCELTYPES];			/* size of morphology arrays */

extern double nval[NCELTYPES][NPARAMS];		/* cell build params */
extern const char *densfil[NCELTYPES][NDENS];
extern int ncell_erased[NCELTYPES];		/* numbers of erased cells */

extern int cellconn [NCELTYPES][NCELTYPES][NSYNTYPES];	/* connections to other cell types */
extern int pickden [NCELTYPES];			/* dendrites to allow connection to */
extern int celnode[NCELTYPES][MAXDEN];		/* numbers of dendrites */

extern const char* chname[NCHANS];
extern const char* chnamea[NCHANS];

extern const char *regname [R_NREGIONS];
extern const char *regnamer[R_NREGIONS];

extern const char* parmname[NCHRATE];

/* standard colors */

extern int black;
extern int blue;
extern int green;
extern int cyan;
extern int red;
extern int magenta;
extern int brown;
extern int white;

extern int gray;
extern int ltblue;
extern int ltgreen;
extern int ltcyan;
extern int ltred;
extern int ltmag;
extern int yellow;
extern int brtwht;

extern char *colornames[];

extern const char *funcfile;
extern const char *segfile;
extern const char *makfile;
extern const char *colorfile;
extern const char *sbmakefile;
extern const char *nvalfile;
extern const char *ha_file;
extern const char *hb_file;
extern const char *hbat_file;
extern const char *rbp_file;
extern const char *dbp1_file;
extern const char *dbp2_file;
extern const char *dbp3_file;
extern const char *dbp4_file;
extern const char *hbp1_file;
extern const char *hbp2_file;
extern const char *am_file;
extern const char *am2_file;
extern const char *am3_file;
extern const char *am4_file;
extern const char *amh_file;
extern const char *amh2_file;
extern const char *sbac_file;
extern const char *aii_file;
extern const char *a17_file;
extern const char *gc_file;
extern const char *gca_file;
extern const char *gcb_file;
extern const char *dsgc_file;
extern const char *gcaoff_file;
extern const char *gcboff_file;
extern const char *anatfiles[];

extern const char *expt;
extern double stfreq;

/* params that control making cells */

extern int make_rods;
extern int make_cones;
extern int make_ha;
extern int make_hb;
extern int make_hbat;
extern int make_dbp1;
extern int make_dbp2;
extern int make_dbp3;
extern int make_dbp4;
extern int make_hbp1;
extern int make_hbp2;
extern int make_rbp;
extern int make_am;
extern int make_am2;
extern int make_am3;
extern int make_am4;
extern int make_amh;
extern int make_amh2;
extern int make_ams;
extern int make_amhs;
extern int make_sbac;
extern int make_aii;
extern int make_a17;
extern int make_gca;
extern int make_gcb;
extern int make_dsgc;
extern int make_gcaoff;
extern int make_gcboff;

/* params that control making synapses */

extern int make_cone_cone;
extern int make_rod_rbp;
extern int make_cone_dbp1;
extern int make_cone_dbp2;
extern int make_cone_dbp3;
extern int make_cone_dbp4;
extern int make_cone_hbp1;
extern int make_cone_hbp2;
extern int make_rod_hbat;
extern int make_cone_ha;
extern int make_cone_hb;
extern int make_ha_ha;
extern int make_hb_hb;
extern int make_rbp_aii;
extern int make_aii_aii;
extern int make_dbp1_aii;
extern int make_dbp2_aii;
extern int make_dbp3_aii;
extern int make_dbp4_aii;
extern int make_aii_hbp1;
extern int make_aii_hbp2;
extern int make_rbp_a17;
extern int make_a17_rbp;
extern int make_dbp1_sbac;
extern int make_dbp2_sbac;
extern int make_dbp3_sbac;
extern int make_dbp4_sbac;
extern int make_dbp1_dbp1;
extern int make_dbp2_dbp2;
extern int make_dbp3_dbp3;
extern int make_dbp4_dbp4;
extern int make_hbp1_hbp1;
extern int make_hbp2_hbp2;
extern int make_dbp1_am;
extern int make_dbp1_am2;
extern int make_dbp1_am3;
extern int make_dbp1_am4;
extern int make_dbp1_amh;
extern int make_dbp1_amh2;
extern int make_dbp1_ams;
extern int make_dbp2_am;
extern int make_dbp2_am2;
extern int make_dbp2_am3;
extern int make_dbp2_am4;
extern int make_dbp2_amh;
extern int make_dbp2_amh2;
extern int make_dbp2_ams;
extern int make_dbp3_am;
extern int make_dbp3_am2;
extern int make_dbp3_am3;
extern int make_dbp3_am4;
extern int make_dbp3_amh;
extern int make_dbp3_amh2;
extern int make_dbp3_ams;
extern int make_dbp4_am;
extern int make_dbp4_am2;
extern int make_dbp4_am3;
extern int make_dbp4_am4;
extern int make_dbp4_amh;
extern int make_dbp4_amh2;
extern int make_dbp4_ams;
extern int make_hbp1_sbac;
extern int make_hbp1_amh;
extern int make_hbp1_amh2;
extern int make_hbp1_amhs;
extern int make_hbp2_sbac;
extern int make_hbp2_amh;
extern int make_hbp2_amh2;
extern int make_hbp2_amhs;
extern int make_dbp1_gca;
extern int make_dbp1_gcb;
extern int make_dbp2_gca;
extern int make_dbp3_gca;
extern int make_dbp4_gcb;
extern int make_sbac_dsgc;
extern int make_sbac_sbac;
extern int make_sbac_dbp1;
extern int make_sbac_dbp2;
extern int make_sbac_dbp3;
extern int make_sbac_dbp4;
extern int make_dbp1_dsgc;
extern int make_dbp2_dsgc;
extern int make_dbp3_dsgc;
extern int make_dbp4_dsgc;
extern int make_hbp1_dsgc;
extern int make_hbp2_dsgc;
extern int make_ams_dbp1;
extern int make_ams_dbp2;
extern int make_ams_dbp3;
extern int make_ams_dbp4;
extern int make_ams_dsgc;
extern int make_amhs_hbp1;
extern int make_amhs_dsgc;
extern int make_dsgc_dsgc;
extern int make_am_dbp1;
extern int make_am2_dbp1;
extern int make_am3_dbp1;
extern int make_am4_dbp1;
extern int make_am_dbp2;
extern int make_am2_dbp2;
extern int make_am3_dbp2;
extern int make_am4_dbp2;
extern int make_am_dbp3;
extern int make_am2_dbp3;
extern int make_am3_dbp3;
extern int make_am4_dbp3;
extern int make_am_dbp4;
extern int make_am2_dbp4;
extern int make_am3_dbp4;
extern int make_am4_dbp4;
extern int make_am_am;
extern int make_am2_am2;
extern int make_am3_am3;
extern int make_am4_am4;
extern int make_amh_amh;
extern int make_amh2_amh2;
extern int make_am_gca;
extern int make_am2_gca;
extern int make_am3_gca;
extern int make_am4_gca;
extern int make_am_gcaoff;
extern int make_am_dsgc;
extern int make_amh_hbp1;
extern int make_amh_hbp2;
extern int make_amh2_hbp1;
extern int make_amh2_hbp2;
extern int make_amh2_amh ;
extern int make_amh_gcaoff;
extern int make_amh_gcboff;
extern int make_amh_dsgc;
extern int make_ams_gca;
extern int make_ams_gcb;
extern int make_ams_gcaoff;
extern int make_ams_gcboff;
extern int make_amhs_gcaoff;
extern int make_amhs_gcboff;
extern int make_hbp1_gcaoff;
extern int make_hbp2_gcaoff;
extern int make_hbp1_gcboff;
extern int make_hbp2_gcboff;
extern int make_ha_cone;
extern int make_ha_dbp1;
extern int make_ha_dbp2;
extern int make_ha_dbp3;
extern int make_ha_dbp4;
extern int make_ha_hbp1;
extern int make_ha_hbp2;
extern int make_hb_dbp1;
extern int make_hb_dbp2;
extern int make_hb_dbp3;
extern int make_hb_dbp4;
extern int make_hb_hbp1;
extern int make_hb_hbp2;

extern int made_gc_comps;

// extern int info;
extern int ninfo;
extern int info_chan;
extern int info_cell;
extern int info_syn;
extern int info_disp;
extern int script;
extern int disp_ray;            /* =1 -> output povray format */

extern int gca_morph;
extern int gcb_morph;
extern int dsgc_morph;
extern int gcoff_morph;
extern int sb_morph;
extern int am_morph;
extern int amh_morph;

extern int gca_biophys;
extern int gcb_biophys;
extern int dsgc_biophys;
extern int gcaoff_biophys;
extern int gcboff_biophys;
extern int am_biophys;
extern int amh_biophys;
extern int rbp_biophys;
extern int aii_biophys;
extern int a17_biophys;
extern int sb_biophys;

extern double dbp1_thetax;
extern double dbp1_thetay;
extern double dbp1_thetaz;
extern int dbp1_flip;
extern double dbp2_thetax;
extern double dbp2_thetay;
extern double dbp2_thetaz;
extern int dbp2_flip;
extern double dbp3_thetax;
extern double dbp3_thetay;
extern double dbp3_thetaz;
extern int dbp3_flip;
extern double dbp4_thetax;
extern double dbp4_thetay;
extern double dbp4_thetaz;
extern int dbp4_flip;
extern double hbp1_thetax;
extern double hbp1_thetay;
extern double hbp1_thetaz;
extern int hbp1_flip;
extern double hbp2_thetax;
extern double hbp2_thetay;
extern double hbp2_thetaz;
extern int hbp2_flip;
extern double rbp_thetax;
extern double rbp_thetay;
extern double rbp_thetaz;
extern int rbp_flip;
extern double ha_thetax;
extern double ha_thetay;
extern double ha_thetaz;
extern int ha_flip;
extern double hb_thetax;
extern double hb_thetay;
extern double hb_thetaz;
extern int hb_flip;
extern double hbat_thetax;
extern double hbat_thetay;
extern double hbat_thetaz;
extern int hbat_flip;
extern double aii_thetax;
extern double aii_thetay;
extern double aii_thetaz;
extern int aii_flip;
extern double am_thetax;
extern double am_thetay;
extern double am_thetaz;
extern int am_flip;
extern double am2_thetax;
extern double am2_thetay;
extern double am2_thetaz;
extern int am2_flip;
extern double am3_thetax;
extern double am3_thetay;
extern double am3_thetaz;
extern int am3_flip;
extern double am4_thetax;
extern double am4_thetay;
extern double am4_thetaz;
extern int am4_flip;
extern double ams_thetax;
extern double ams_thetay;
extern double ams_thetaz;
extern int ams_flip;
extern double amh_thetax;
extern double amh_thetay;
extern double amh_thetaz;
extern int amh_flip;
extern double amh2_thetax;
extern double amh2_thetay;
extern double amh2_thetaz;
extern int amh2_flip;
extern double a17_thetax;
extern double a17_thetay;
extern double a17_thetaz;
extern int a17_flip;
extern double sbac_thetax;
extern double sbac_thetay;
extern double sbac_thetaz;
extern int sbac_flip;
extern double dsgc_thetax;
extern double dsgc_thetay;
extern double dsgc_thetaz;
extern int dsgc_flip;
extern double gc_thetax;
extern double gc_thetay;
extern double gc_thetaz;
extern int gc_flip;
extern int flip;

extern double dvrev;		/* for dens_xxx.n files */
extern double dvst;

/* User may provide numbers from command line. */
/* Leave unset numbers = -1 for automatic algorithms below. */

extern int n_rods;
extern int n_rbp;
extern int n_cones;
extern int n_ha;
extern int n_hb;
extern int n_hbat;
extern int n_dbp1;
extern int n_dbp2;
extern int n_dbp3;
extern int n_dbp4;
extern int n_hbp1;
extern int n_hbp2;
extern int n_aii;
extern int n_a17;
extern int n_am;
extern int n_am2;
extern int n_am3;
extern int n_am4;
extern int n_amh;
extern int n_amh2;
extern int n_ams;
extern int n_amhs;
extern int n_sbac;
extern int n_gca;
extern int n_gcb;
extern int n_dsgc;
extern int n_gcaoff;
extern int n_gcboff;

extern int nrods;
extern int nrbp;
extern int ncones;
extern int nha;
extern int nhb;
extern int nhbat;
extern int ndbp1;
extern int ndbp2;
extern int ndbp3;
extern int ndbp4;
extern int nhbp1;
extern int nhbp2;
extern int naii;
extern int na17;
extern int nam;
extern int nam2;
extern int nam3;
extern int nam4;
extern int namh;
extern int namh2;
extern int nams;
extern int namhs;
extern int nsbac;
extern int ngca;
extern int ngcb;
extern int ndsgc;
extern int ngcaoff;
extern int ngcboff;


extern int _NA;
extern int _NA0;
extern int _NA1;
extern int _NA2;
extern int _NA3;
extern int _NA4;
extern int _NA5;
extern int _NA6;
extern int _NA8;

extern int _K0;
extern int _K1;
extern int _K2;
extern int _K3;
extern int _K4;
extern int _K5;
extern int _K6;
extern int _K7;

extern int _K;
extern int _KDR;
extern int _KA;
extern int _KH;
extern int _KIR;
extern int _KHZ;
extern int _KV3;

extern int _KCA0;
extern int _KCA1;
extern int _KCA2;
extern int _KCA3;
extern int _KCA4;
extern int _KCA5;

extern int _BKCA;
extern int _SKCA1;
extern int _SKCA2;

extern int _CLCA;
extern int _CLCAC;
extern int _CLCA1;
extern int _CLCA2;

extern int _CA0;
extern int _CA1;
extern int _CA2;
extern int _CA3;
extern int _CA4;
extern int _CA5;
extern int _CA6;
extern int _CA7;

extern int _NMDA1;
extern int _NMDA2;

extern int _AMPA1;
extern int _AMPA2;
extern int _AMPA3;
extern int _AMPA4;
extern int _AMPA5;

extern int _GABA1;
extern int _GABA2;
extern int _GABA3;
extern int _GABA4;

extern int _GLY;

extern int _CGMP1;
extern int _CGMP2;
extern int _CGMP3;
extern int _CGMP4;
extern int _CGMP5;
extern int _CGMP6;
extern int _CGMP7;
extern int _CGMP8;
extern int _CGMP9;
extern int _CGMP10;
extern int _CGMP11;

extern int _CA;
extern int _CA_L;
extern int _CA_T;

// joesterle begin
extern int _HCN1;
extern int _HCN2;
extern int _HCN3;
extern int _HCN4;
// joesterle end

extern int _CAP;
extern int _CAPUMP;
extern int _CAPK;
extern int _CABV;
extern int _CABK;
extern int _CABT;
extern int _CABI;
extern int _CASH;
extern int _CAE;
extern int _CAEXCH;

extern int _CAS;
extern int _VM2;
extern int _VM3;
extern int _KFC;
extern int _KFCICR;
extern int _KAC;
extern int _KACICR;
extern int _KRC;
extern int _KRCICR;
extern int _K1C;
extern int _K1CICR;
extern int _K2C;
extern int _K2CICR;
extern int _NHC;
extern int _NCICR;
extern int _MHC;
extern int _MCICR;
extern int _PHC;
extern int _PCICR;
extern int _CAS2;
extern int _IP3I;
extern int _BIP3;
extern int _VIP3;
extern int _V3IP3;
extern int _V2IP3;
extern int _V4IP3;
extern int _D1IP3;
extern int _D2IP3;
extern int _D3IP3;
extern int _D4IP3;
extern int _A2IP3;
extern int _A3IP3;
extern int _B2IP3;
extern int _K3IP3;
extern int _C1CICR;

extern int _VSTART;
extern int _VST;
extern int _VRV;
extern int _VREV;
extern int _RM;
extern int _RI;
extern int _CM;
extern int _DIA;
extern int _CPLAM;
extern int _CPL;
extern int _CMUL;
extern int _COLOR;
extern int _COL;

extern int na1nois;
extern int na2nois;
extern int na3nois;
extern int na4nois;
extern int na5nois;
extern int na6nois;
extern int na8nois;
extern int k1nois;
extern int k3nois;
extern int k4nois;
extern int k5nois;
extern int k6nois;
extern int k7nois;
// begin joesterle
extern int hcn1nois;
extern int hcn2nois;
extern int hcn3nois;
extern int hcn4nois;
// end joesterle
extern int kca1nois;
extern int kca3nois;
extern int kca4nois;
extern int kca5nois;
extern int kca6nois;
extern int ca1nois;
extern int ca3nois;
extern int ca5nois;
extern int ca6nois;
extern int ca7nois;
extern int clcanois;

extern int maxha;
extern int maxhb;
extern int maxhbat;
extern int maxrbp;
extern int maxdbp1;
extern int maxdbp2;
extern int maxdbp3;
extern int maxdbp4;
extern int maxhbp1;
extern int maxhbp2;
extern int maxa17;
extern int maxaii;
extern int maxsb;
extern int maxam;
extern int maxam2;
extern int maxam3;
extern int maxam4;
extern int maxamh;
extern int maxamh2;
extern int maxams;
extern int maxgc;
extern int maxdsgc;
extern int maxgcoff;

extern double node_scale;

extern double rod_nscale;
extern double cone_nscale;
extern double hbat_nscale;
extern double ha_nscale;
extern double hb_nscale;
extern double rbp_nscale;
extern double dbp1_nscale;
extern double dbp2_nscale;
extern double dbp3_nscale;
extern double dbp4_nscale;
extern double hbp1_nscale;
extern double hbp2_nscale;
extern double aii_nscale;
extern double sbac_nscale;
extern double am_nscale;
extern double am2_nscale;
extern double am3_nscale;
extern double am4_nscale;
extern double ams_nscale;
extern double a17_nscale;
extern double gca_nscale;
extern double gcb_nscale;
extern double dsgc_nscale;
extern double gcaoff_nscale;
extern double gcboff_nscale;

extern double synap_scale;
extern double gj_scale;

extern double arrsiz;
extern double xarrsiz;
extern double yarrsiz;
extern double zarrsiz;
extern double garrsiz;
extern double gxarrsiz;
extern double gyarrsiz;
extern double amarrsiz;
extern double amxarrsiz;
extern double amyarrsiz;
extern double am2arrsiz;
extern double am2xarrsiz;
extern double am2yarrsiz;
extern double am3arrsiz;
extern double am3xarrsiz;
extern double am3yarrsiz;
extern double am4arrsiz;
extern double am4xarrsiz;
extern double am4yarrsiz;
extern double amharrsiz;
extern double amhxarrsiz;
extern double amhyarrsiz;
extern double amh2arrsiz;
extern double amh2xarrsiz;
extern double amh2yarrsiz;
extern double amsarrsiz;
extern double amsxarrsiz;
extern double amsyarrsiz;
extern double aiiarrsiz;
extern double aiixarrsiz;
extern double aiiyarrsiz;
extern double sbarrsiz;
extern double sbxarrsiz;
extern double sbyarrsiz;
extern double dbp1arrsiz;
extern double dbp1xarrsiz;
extern double dbp1yarrsiz;
extern double dbp2arrsiz;
extern double dbp2xarrsiz;
extern double dbp2yarrsiz;
extern double dbp3arrsiz;
extern double dbp3xarrsiz;
extern double dbp3yarrsiz;
extern double dbp4arrsiz;
extern double dbp4xarrsiz;
extern double dbp4yarrsiz;
extern double hbp1arrsiz;
extern double hbp1xarrsiz;
extern double hbp1yarrsiz;
extern double hbp2arrsiz;
extern double hbp2xarrsiz;
extern double hbp2yarrsiz;
extern double rodarrsiz;
extern double rodxarrsiz;
extern double rodyarrsiz;
extern double conearrsiz;
extern double conexarrsiz;
extern double coneyarrsiz;
extern double arrcentx;
extern double arrcenty;
extern double arrcentz;
extern double rodarrcentx;
extern double rodarrcenty;
extern double conearrcentx;
extern double conearrcenty;
extern double dbp1arrcentx;
extern double dbp1arrcenty;
extern double dbp2arrcentx;
extern double dbp2arrcenty;
extern double dbp3arrcentx;
extern double dbp3arrcenty;
extern double dbp4arrcentx;
extern double dbp4arrcenty;
extern double hbp1arrcentx;
extern double hbp1arrcenty;
extern double hbp2arrcentx;
extern double hbp2arrcenty;

extern int dbp1_first_cent;
extern int dbp2_first_cent;
extern int dbp3_first_cent;
extern int dbp4_first_cent;
extern int hbp1_first_cent;
extern int hbp2_first_cent;

extern int setarrsiz;
extern int setxarrsiz;

extern double *amhxarr;
extern double *amhyarr;
extern double *amhtharr;
extern int *amhnarr;
extern double *amh2xarr;
extern double *amh2yarr;
extern double *amh2tharr;
extern int *amh2narr;
extern double *amsxarr;
extern double *amsyarr;
extern double *amstharr;
extern int *amsnarr;
extern double *amxarr;
extern double *amyarr;
extern double *amtharr;
extern int *amnarr;
extern double *am2xarr;
extern double *am2yarr;
extern double *am2tharr;
extern int *am2narr;
extern double *am3xarr;
extern double *am3yarr;
extern double *am3tharr;
extern int *am3narr;
extern double *am4xarr;
extern double *am4yarr;
extern double *am4tharr;
extern int *am4narr;
extern double *aiixarr;
extern double *aiiyarr;
extern double *aiitharr;
extern int *aiinarr;

extern double *rodxarr;
extern double *rodyarr;
extern double *rodtharr;
extern int *rodnarr;
extern double *conexarr;
extern double *coneyarr;
extern double *conetharr;
extern int *conenarr;
extern double *dbp1xarr;
extern double *dbp1yarr;
extern double *dbp1tharr;
extern double *dbp1ytharr;
extern int *dbp1narr;
extern double *dbp2xarr;
extern double *dbp2yarr;
extern double *dbp2tharr;
extern double *dbp2ytharr;
extern int *dbp2narr;
extern double *dbp3xarr;
extern double *dbp3yarr;
extern double *dbp3tharr;
extern double *dbp3ytharr;
extern int *dbp3narr;
extern double *dbp4xarr;
extern double *dbp4yarr;
extern double *dbp4tharr;
extern double *dbp4ytharr;
extern int *dbp4narr;
extern double *hbp1xarr;
extern double *hbp1yarr;
extern double *hbp1tharr;
extern double *hbp1ytharr;
extern int *hbp1narr;
extern double *hbp2xarr;
extern double *hbp2yarr;
extern double *hbp2tharr;
extern double *hbp2ytharr;
extern int *hbp2narr;
extern double *rbpxarr;
extern double *rbpyarr;
extern double *rbptharr;
extern double *rbpytharr;
extern int *rbpnarr;
extern double *sbxarr;
extern double *sbyarr;
extern double *sbytharr;
extern double *sbtharr;
extern int *sbnarr;
extern double *a17xarr;
extern double *a17yarr;
extern double *a17tharr;
extern int *a17narr;
extern double *gcxarr;
extern double *gcyarr;
extern double *gctharr;
extern double *gcytharr;
extern int *gcnarr;

extern double dispsize;
extern double maxsize;

extern double axarbdia;

extern double am_seglen;
extern double am_dia_prox_rad;
extern double am_dia_prox_factor;
extern int am_den_seg;
extern int n_amseg;
extern int makesbfile;
extern double xradius;

extern double axon_br_dia;
extern double cbplam;
extern double cbplam2;
extern double ath_dia;
extern double dia_prox_rad;
extern double dia_prox_factor;
extern double gctheta;

extern int bp_dend1;

extern int hz_dend1;
extern int hz_dend2;
extern int hz_dend3;
extern int hz_dend4;
extern int hz_dend5;

extern int aii_dend1;
extern int aii_dend2;
extern int aii_dend3;
extern int aii_dend4;
extern int aii_dend5;

extern int sbac_dend1;
extern int sbac_dend2;
extern int sbac_dend3;
extern int sbac_dend4;
extern int sbac_dend5;

extern int am_dend1;
extern int am_dend2;
extern int am_dend3;
extern int am_dend4;
extern int am_dend5;

extern int gc_dend1;
extern int gc_dend2;
extern int gc_dend3;
extern int gc_dend4;
extern int gc_dend5;

extern double sbtheta;

extern double dnoise;
extern int pnoise, vnoise, cnoise;
extern int Chnoise;
extern double bg_inten;

extern int rod_rect_arr;
extern int cone_rect_arr;
extern int make_one_dbp1;
extern int make_one_dbp2;
extern int make_one_dbp3;
extern int make_one_dbp4;
extern int remove_nconns;
extern int print_conns;
extern int print_avg_conns;
extern int disp_ct;
extern int disp_cn;
extern int plot_freq;

extern double rodspac;
extern double conespac;

extern double cone_maxcond;
extern double rod_maxcond;
extern int cone_pigm;

extern double cone_timec;
extern double rod_timec;

extern double cone_loopg;
extern double rod_loopg;

extern int gcdistnod;

extern double xmax, xmin;
extern double ymax, ymin;
extern double zmax, zmin;
extern double set_zmax, set_zmin;
extern double disp_zmax, disp_zmin;		// for dispcelltype()
extern double disp_gc_zmax, disp_gc_zmin;	// for dispcelltype()
extern double disp_dsgc_zmax, disp_dsgc_zmin;	// for dispcelltype()

extern int xmaxnode, xminnode;
extern int ymaxnode, yminnode;
extern int zmaxnode, zminnode;

extern double disp_margin;
extern double disp_calib_len;
extern int limit_array;
extern int update_array;

extern double ttxbath;
extern double ttxsoma;
extern double ttxdend;
extern double tea;
extern double fourap;
extern double zd7288;
extern double ibtox;

extern double dia_min_factor;
extern double dendd_dia_factor;
extern double dendi_dia_factor;
extern double dendp_dia_factor;
extern double dend_dia_factor;
extern double ax_dia_factor;
extern double cell_dia_factor;
extern double den_dia;
extern double taper;
extern double elec_resist1, elec_resist2;
extern double elec_capac1, elec_capac2;
extern double gc_dend_cplam;
extern double am_dend_cplam;

extern double gc_rm;
extern double gc_vr;
extern double gc_vs;
extern double dvs;
extern double mxrot;
extern double myrot;

extern int nshell;
extern int nshell_soma;

extern int nolabels;
extern const char *plotlabel;

extern const char *confdir;

extern const char *chanparamsfile;
extern const char *def_densfile;
extern const char *syn_savefile;
extern const char *syn_restorefile;

extern const char *cone_densfile;
extern const char *rod_densfile;
extern const char *hbat_densfile;
extern const char *ha_densfile;
extern const char *hb_densfile;
extern const char *rbp_densfile;
extern const char *dbp1_densfile;
extern const char *dbp2_densfile;
extern const char *dbp3_densfile;
extern const char *dbp4_densfile;
extern const char *hbp1_densfile;
extern const char *hbp2_densfile;
extern const char *a17_densfile;
extern const char *aii_densfile;
extern const char *sbac_densfile;
extern const char *am_densfile;
extern const char *am2_densfile;
extern const char *am3_densfile;
extern const char *am4_densfile;
extern const char *amh_densfile;
extern const char *amh2_densfile;
extern const char *ams_densfile;
extern const char *amhs_densfile;
extern const char *gca_densfile;
extern const char *gcb_densfile;
extern const char *dsgc_densfile;
extern const char *gcaoff_densfile;
extern const char *gcboff_densfile;

extern const char *cone_densfile2;
extern const char *rod_densfile2;
extern const char *hbat_densfile2;
extern const char *ha_densfile2;
extern const char *hb_densfile2;
extern const char *rbp_densfile2;
extern const char *dbp1_densfile2;
extern const char *dbp2_densfile2;
extern const char *dbp3_densfile2;
extern const char *dbp4_densfile2;
extern const char *hbp1_densfile2;
extern const char *hbp2_densfile2;
extern const char *a17_densfile2;
extern const char *aii_densfile2;
extern const char *sbac_densfile2;
extern const char *am_densfile2;
extern const char *am2_densfile2;
extern const char *am3_densfile2;
extern const char *am4_densfile2;
extern const char *amh_densfile2;
extern const char *amh2_densfile2;
extern const char *ams_densfile2;
extern const char *amhs_densfile2;
extern const char *gca_densfile2;
extern const char *gcb_densfile2;
extern const char *dsgc_densfile2;
extern const char *gcaoff_densfile2;
extern const char *gcboff_densfile2;

