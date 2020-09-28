/* retsim.h */
/* includes for retsim.cc */

#define abs(x)          ((x) < 0 ? -(x) : (x))
#define max(x, y)       (((x) < (y)) ? (y) : (x))
#define min(x, y)       (((x) < (y)) ? (x) : (y))

#define FILNAMSIZ 50	/* max size of file names for nval, density, chanparams */

#define PI M_PI
#define DEG 57.29577951308232087680 	/* deg/radian */

#include "nval.h"

#define soma  0		/* node number for soma */
#define axtrm 1		/* node number for axon terminals (bipolar cell, etc.) */
#define recpnt 10000	/* node number for sharp electrode recording point */
#define MAXDEN 2000	/* size of celnode[nceltypes][MAXDEN] */

#define NCELINFO 10	/* =0 -> cell count, = 1-9 -> cell numbers connecting */
#define NCELLS   0
#define CELN     1

#define NBRANCHED 0	/* not branched, separate dendrites for each input */
#define BRANCHED  1	/* branched dendritic tree */
#define HBRANCHED 2	/* highly branched */
#define SBRANCHED 3	/* starburst branching */

					/* defs for "disp" variable */
#define DISP  1                         /* allow "display" statement to run */
#define DCOMP 2                         /* display compartments */
#define DCONN 4                         /* display connections */
#define DNODE 8                         /* display nodes */
#define DSTIM 16                        /* display stimulus */
#define DMOVIE 32                       /* save neural elements for display */

#include "colors.h"			/* standard colors in nc/src/colors.h */

/*-------------------------------------------------*/

/*  Definition of columns for "anatfile" which defines cell morphology: */
/*

#    node  parent   dia     xbio     ybio     zbio     region   dendn
#
      0      0     15       -15.42   2.34     10       SOMA     0
      1      0     3.5      -12.12   5.4      15       AXON     0
      2      1     2.1      -13.61   6.71     15       DEND     1
      3      2     1.58     -12.25   4.76     15       DEND     2
      .
      .
      .
*/

#define C_NODE  0
#define C_PAR   1
#define C_DIA   2
#define C_AX    3
#define C_AY    4
#define C_AZ    5
#define C_REG   6
#define C_DENDN 7

/* Regions: Column number 1 */

// #define R_0	0  	// don't use region 0, because this is null value for synreg
#define R_1	1
#define R_2	2
#define R_3	3
#define R_4	4
#define R_5	5
#define R_6	6
#define R_7	7
#define R_8	8
#define R_9	9
#define R_10	10
#define R_11	11
#define R_12	12
#define R_13	13
#define R_14	14
#define R_15	15
#define R_16	16
#define R_17	17
#define R_18	18
#define R_19	19
// #define R_20	20

#define R_NREGIONS     20

/*-------------------------------------------------*/

#define LOCX 0
#define LOCY 1
#define DIR  2

#define MORPH_REAL  0
#define MORPH_A1    1
#define MORPH_A2    2
#define MORPH_A3    3
#define MORPH_A4    4
#define MORPH_SIMP  3
#define MORPH_SOMA  4

#define MORPH_T1    101
#define MORPH_T2    102
#define MORPH_T3    103
#define MORPH_T3A   104
#define MORPH_T3B   105
#define MORPH_T4    106
#define MORPH_T5    107
#define MORPH_T5A   108
#define MORPH_T5B   109
#define MORPH_T5R   110
#define MORPH_T5X   111
#define MORPH_TXBC  112
#define MORPH_T6    113
#define MORPH_T7    114
#define MORPH_T8    115
#define MORPH_T9    116


/*-------------------------------------------------*/

/* Channel types: Row numbers in "dens_xxx.n" */

#define C_NA0     0       /* NaV1.2, no noise, type 0 */
#define C_NA1     1       /* NaV1.2, simple noise */
#define C_NA2     2       /* NaV1,2 same as type 1, but markov */
#define C_NA3     3       /* Na: type 3 */
#define C_NA4     4       /* Na: type 4 */
#define C_NA5     5       /* slowly inactivating Na: type 5 */
#define C_NA6     6       /* persistent Na: type 6 */
#define C_NA8     7       /* persistent Na: type 6 */

#define C_K0      8       /* type 0, no noise */
#define C_K1      9       /* type 1, noise */
#define C_K2      10      /* type 2, noise */
#define C_K3      11      /* type 3, noise */
#define C_K4      12      /* type 4, Ih */
#define C_K5      13      /* type 5, Kir */
#define C_K6      14      /* type 6, Kv3 */
#define C_K7      15      /* type 6, Kv3 */

#define C_KCA0    16	  /* SK, like HH */
#define C_KCA1    17	  /* SK, Markov */
#define C_KCA2    18	  /* BK, Markov */
#define C_KCA3    19	  /* BK, Markov */
#define C_KCA4    20	  /* SK, Markov */
#define C_KCA5    21	  /* SK, slow, Markov */
#define C_KCA6    22	  /* BK, Markov */

#define C_CLCA1   23
#define C_CLCA2   24

#define C_CA0     25      /* L-type Ca */
#define C_CA1     26      /* L-type Ca, Markov */
#define C_CA2     27      /* T-type Ca */
#define C_CA3     28      /* T-type Ca, Markov */
#define C_CA4     29      /* T-type Ca */
#define C_CA5     30      /* T-type Ca, Markov  */
#define C_CA6     31      /* T-type Ca, low thresh */
#define C_CA7     32      /* T-type Ca, low thresh */

#define C_NMDA1   33      /* NMDA channel, simple */
#define C_NMDA2   34      /* NMDA channel, complex */

#define C_AMPA1   35      /* AMPA receptor, type 1 */
#define C_AMPA2   36      /* AMPA receptor, type 2 */
#define C_AMPA3   37      /* AMPA receptor, type 3 */
#define C_AMPA4   38      /* AMPA receptor, type 4 */
#define C_AMPA5   39      /* AMPA receptor, type 5 */

#define C_GABA1   40      /* GABA-A receptor, type 1 */
#define C_GABA2   41      /* GABA-A receptor, type 2 */
#define C_GABA3   42      /* GABA-A receptor, type 3 */
#define C_GABA4   43      /* GABA-A receptor, type 3 */

#define C_GLY     44      /* Glycine receptor, type 1 */

#define C_CGMP1   45      /* CGMP channel, type 1 */
#define C_CGMP2   46      /* CGMP channel, type 2 */
#define C_CGMP3   47      /* CGMP channel, type 3 */
#define C_CGMP4   48      /* CGMP channel, type 4 */
#define C_CGMP5   49      /* CGMP channel, type 5 */
#define C_CGMP6   50      /* CGMP channel, type 6 */
#define C_CGMP7   51      /* CGMP channel, type 7 */
#define C_CGMP8   52      /* CGMP channel, type 8 */
#define C_CGMP9   53      /* CGMP channel, type 9 */
#define C_CGMP10  54      /* CGMP channel, type 10 */
#define C_CGMP11  55      /* CGMP channel, type 11 */

// joesterle begin
#define C_HCN1    56
#define C_HCN2    57
#define C_HCN3    58
#define C_HCN4    59

#define C_MGLUR   60
// joesterle end

#define C_CAPUMP  61      /* capump vmax */
#define C_CAPKM   62      /* capump km */
#define C_CABVMAX 63   /* cabuf vmax */
#define C_CABKD   64   /* cabuf kd */
#define C_CABTOT  65   /* cabuf btot */
#define C_CABTOTI 66   /* cabuf btoti (first shell) */
#define C_CASHELL 67   /* ca nshells */
#define C_CAEXCH  68   /* ca exch rate A/mM4/cm2 */
#define C_CAS     69   /* initial Ca in CICR store (M) */
#define C_VM2     70   /* max rate of Ca pump into ryanodine store */
#define C_VM3     71   /* max rate of Ca pump from ryanodine store */
#define C_KFCICR  72   /* passive leak rate from ryan store /s */
#define C_KACICR  73   /* thresh for ryan store release (M) */
#define C_KRCICR  74   /* thresh for ryan store release (M) */
#define C_K1CICR  75   /* thresh for ryan store uptake (M) */
#define C_K2CICR  76   /* thresh for ryan store uptake (M) */
#define C_NCICR   77   /* Hill Coeff, ryan store uptake */
#define C_MCICR   78   /* Hill Coeff, ryan store release */
#define C_PCICR   79   /* Hill Coeff, ryan store release */
#define C_C1CICR  80   /* volume fraction of calcium store */
#define C_CAS2    81   /* initial Ca in IP3 store (M) */
#define C_IP3I    82   /* initial [IP3] intrcellular */
#define C_BIP3    83   /* frac rate of IP3 release */
#define C_VIP3    84   /* init IP3 store flux */
#define C_V2IP3   85   /* ip3 release rate */
#define C_V3IP3   86   /* ip3 uptake rate */
#define C_V4IP3   87   /* ip3 release rate */
#define C_D1IP3   88   /* IP3 store const for Q value for h gate */
#define C_D2IP3   89   /* IP3 store const for Q value for h gate */
#define C_D3IP3   90   /* IP3 store const for Q value for h gate */
#define C_D4IP3   91   /* IP3 store const for Q value for h gate */
#define C_A2IP3   92   /* IP3 store forward const for Tau M value */
#define C_A3IP3   93   /* IP3 store const for IP3 Tau H */
#define C_B2IP3   94   /* IP3 store backward const for Tau M value */
#define C_K3IP3   95   /* IP3 store kd for Ca uptake */
#define C_VSTART  96   /* starting voltage */
#define C_VREV    97   /* membrane reversal potential */
#define C_RM      98   /* membrane resistivity */
#define C_RI      99   /* axial resistivity */
#define C_CM      100  /* membrane capacitance */
#define C_DDIA    101  /* diameter factor for region */
#define C_CPLAM   102  /* compartment size for region */
#define C_CMUL    103  /* synaptic conductance mult for region */
#define C_COLOR   104  /* color for region */

#define NCHANS    105
#define ECHANS    60	  /* last channel type before Ca pumps, membr params */

/*-------------------------------------------------*/

/* defs for "chval[][]" */

#define CHOFFM  0             /* channel m offset */
#define CHOFFH  1             /* channel h offset */
#define CHTAUA  2             /* channel tau a */
#define CHTAUB  3             /* channel tau b */
#define CHTAUC  4             /* channel tau c */
#define CHTAUD  5             /* channel tau d */
#define NCHRATE 6
#define NCH NCHRATE

/*-------------------------------------------------*/

#define ALPHA 1
#define BETA  2

#define gs(sparam) getsv(prect,sparam,rcs)

#define NSYNTYPES 9

#define NDENS 2

/*-------------------------------------------------*/

extern int *ndens[NCELTYPES];	/* density array indexed by ct,cn */
extern double celdens[NCELTYPES][NDENS][NCHANS][R_NREGIONS];/* density data */

/*-------------------------------------------------*/

extern int setplsep;   /* in ncmain.cc */

/* in nc: ncm.cc */

void disperror(void);

/*-------------------------------------------------*/

void initneurvals(void );	/* initialize neuron params */

void initsynconn(void );	/* initialize synaptic connection table */
				/* after possible modifications in * setparams() */
/* functions in makcel.cc: */

void setmorphfiles(void);		/* set the filenames for realistic neuron morphologies */
int find_ct(const char *ct_name);
void make_ct(int ct);
void set_ncel(int ct, int n);
void   setn(int cel, int var, double val);
double getn(int cel, int var);
int anysetn (int var);
double getcv(int ctype, int var, int n);
void   setcv(int ctype, int var, double val);
double getsv(int ctype, int var, int n);
void   setsv(int ctype, int var, int n, double val);
void initsynconn(void);
int getconn (int a, int b);
int getconns(int a, int b, int i);
void enlarge_array(int ct, int cn);
int setupcells(int ct, int ncells);
int setupcells(int ct, int ncells, int first_cent);
int setupcells(int ct, int ncells, double lxarrsiz);
int setupcells(int ct, int ncells, double *sxarr, double *syarr);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *stharr);
int setupcells(int ct, int ncells, double *sxarr, double *syarr);
int setupcells(int ct, int ncells, double *sxarr, double *syarr, double *stharr);
int setupcells(int ct, int ncells, double *sxarr, double *syarr, double *stharr, int *snarr);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double arrcentx, double arrcenty);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double arrcentx, double arrcenty, int first_cent);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *stharr);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *stharr, int *snarr);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double arrcentx, double arrcenty, double *sxarr, double *syarr, double *stharr);
int setupcells(int ct, int ncells, double lxarrsiz, double lyarrsiz, double arrcentx, double arrcenty, double *sxarr, double *syarr, double *stharr, int *snarr);

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *sthyarr, double *stharr, int *snarr);
int setupcells(int ct,int ncells, double *sxarr, double *syarr, double *sthyarr, double *stharr, int *snarr);

void checkcellin (int ct);
void checkcellin (int ct, int cti);
void checkcellout (int ct);
void checkcellout (int ct, int ct2);
void checkcellout (int ct, int ct2, int ct3);
void checkcellout (int ct, int ct2, int ct3, int ct4);
void printfilenames(FILE *filout);
void makcell (int ctype, int n, double xpos, double ypos, double thetay, double thetaz, int flip);
void connect_types (int ct1, int ct2);
void connect_types (int ct1, int ct2, int make_rcr, int dyadc);
void dispcelltype(int ct, int start, int stop, double scal, double nscal);
void dispcelltype(int ct, int start, int stop, double scal, double nscal, double zrange1, double zrange2);
void draw_plotlabel(int rsize);

/* functions in celseg.cc: */

void findceldens(void);		/* read channel densities from density files */
void modchandens(void);		/* possibly modify channel densities */
void printchaninfo(void);	/* print channel information */
void make_celseg(int ct,int cn,int n1,int n2,double dia1, double dia2,int region,double cplam);
double getdens(int ctype, int dnum, int chantyp, int region); /* get channel density */
void set_chan_offsets(chattrib *a, int ct, int chvalIndex);
void set_chan_offset(int ct, int chvalIndex, int chrate, double val);

/* functions in celfuncs.cc: */

double rrange (double l, double h);
int ff(void);
int ffs(void);
int bptype(int ct);

double node_angle (node *npnt1, node *npnt2);
double node_angle (node *npnt);
double node_angle (int ct, int cn, int n);

double gauss(double r, double rad);
double atanx(double dx, double dy);
void makanatfile (int celltype,int cellnr);
int inrange(double a1, double a2, double t);
int findmid(int ct, double xoffset, double yoffset);
int findmida(int ct, double xoffset, double yoffset);
int findcell(int ct, double xoffset, double yoffset);
int findcella(int ct, double roffset, double theta);
int findnodloc(int ct, int cn, double xoffset, double yoffset);
int findnodloc(int ct, int cn, double xoffset, double yoffset, double maxdist);
int findnodlocr(int ct, int cn, double xoffset, double yoffset);
int findnodlocr(int ct, int cn, double xoffset, double yoffset, double maxdist);
int findnodlocra(int ct, int cn, double roffset, double theta);
int findnodlocraz(int ct, int cn, double roffset, double theta, double zmax, double zmin);
int findnodlocra(int ct, int cn, double roffset, double theta, double maxdist);
int findnodlocz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin);
int findnodlocrz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin);

int findsynlcpre (int ct, int cn, int ct2, int cn2);
int findsynlocp  (int ct, int cn, int ct2, int cn2);
int findsynlocpa (int ct, int cn, double xoffset, double yoffset, int ct2, int cn2);

synapse *findsyn(int ct, int cn, int nod);
synapse *findsyn(int ct, int cn, int nod, int ct2);
synapse *findsyn(int ct, int cn, int nod, int ct2, int cn2);
synapse *findsynloc(int ct, double xoffset, double yoffset);
synapse *findsynloc(int ct, int cn, double xoffset, double yoffset);
int      findsynlocn(int ct, int cn, double xoffset, double yoffset);
synapse *findsynloc(int ct, int cn, double xoffset, double yoffset, double vrev);
synapse *findsynloc(int ct, int cn, double xoffset, double yoffset, double vrev, double maxdist);
synapse *findsynloc(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev);
synapse *findsynloc(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev, double maxdist);
synapse *findsynloca(int ct, double roffset, double theta);
synapse *findsynloca(int ct, int cn, double roffset, double theta);
synapse *findsynloca(int ct, int cn, int ct2, int cn2, double roffset, double theta);
synapse *findsynloca(int ct, int ct2, int cn2, double roffset, double theta);
synapse *findsynlocr(int ct, int cn);
synapse *findsynlocr(int ct, int cn, double xoffset, double yoffset);
synapse *findsynlocr(int ct, int cn, double xoffset, double yoffset, double vrev);
synapse *findsynlocr(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset);
synapse *findsynlocr(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev);
synapse *findsynlocr(int ct, int cn, int ct2, int cn2, double xoffset, double yoffset, double vrev, double maxdist);
synapse *findsynlocra(int ct, int cn, double roffset, double theta);
synapse *findsynlocra(int ct,         double roffset, double theta);
synapse *findsynlocra(int ct, int cn, double roffset, double theta, double vrev);
synapse *findsynlocra(int ct, int cn, int ct2, int cn2, double roffset, double theta, double vrev);
synapse *findsynlocra(int ct, int cn, int ct2, int cn2, double roffset, double theta, double vrev, double maxdist);
synapse *findsynlocz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin);
synapse *findsynlocrz(int ct, int cn, double xoffset, double yoffset, double zmax, double zmin);
synapse *findsynlocaz(int ct, int cn, double roffset, double theta, double zmax, double zmin);
synapse *findsynlocarz(int ct, int cn, double roffset, double theta, double zmax, double zmin);
synapse *findsynlocaz(int ct, int cn, double roffset, double theta, double zmax, double zmin);
synapse *findsynlocaz(int ct,         double roffset, double theta, double zmax, double zmin);
synapse *findsynlocarz(int ct, 	      double roffset, double theta, double zmax, double zmin);


int findsynlocx(int ct, int cn, double xoffset, double yoffset);
int findmidc(int ct, int offset);
int findnext(int ct, int count);
int find_gtconn(int ct, int nconn);

double mod (double a,double b);
double round(double n,double p);
double modangl(double a);
double sindeg (double theta);
double cosdeg (double theta);

double rad_dist2 (int celtyp, int n, int tnod, int parent, int stflag);
double rad_dist (int celtyp, int n, int tnod, int parent, int stflag);
int rad_dir (int elnum, int c1, int c2, int c3);
double rad_diam (double dist,double diaspc,double diabs);
int taperdistden(int ct,int cellnr,int nod1,int nod2,double diabs,double diaspc,int region,int nden);
double taperdia (int n1a,int n1b,int n1c, int n2a,int n2b,int n2c,double diabs,double diaspc);
int taperden  (int ct,int cn,int n1,int n2, double cd1, double cd2,int nden);
double sigm(double xmin,double xmax,double ymin,double yrange,double x);
double comp_phase(double tfreq,double delaytime);
double sinewaves(double phase1,double phase2,double ampl1,double ampl2);

int mid(int siz);
int midrow(int siz);
void find_maxmin(int celltype,int cellnum);
double find_maxrad(int ct,int cn);
double qfact(double q10);
void makanatfile (int celltype,int cellnr) ;
int dendn_node (int ct, int dendn);
void make_gc_comps(int ct, int cn, int region, int ct2);

/* functions in synfuncs.cc: */

void print_connections (int ct);
void print_avg_connections (int ct);
int ncel_in (int to_celltype, int to_cellnum, int from_celltype);
int tot_ncel_in (int to_celltype, int to_cellnum);
int tot_ncel_ind (int to_celltype,   int to_cellnum);

int ncel_out(int from_celltype, int from_cellnum, int to_celltype);
int tot_ncel_out(int from_celltype, int from_cellnum);
int tot_ncel_outd(int from_celltype, int from_cellnum);
int connected (int from_celltype, int from_cellnum, int to_celltype, int to_cellnum);
int connected2 (int from_celltype, int from_cellnum, int to_ct2, int to_ct3);
int connected3 (int from_celltype, int from_cellnum, int to_ct2, int to_ct3, int to_ct4);
void save_synapses ();
void restore_synapses ();

/* Functions in sb_recfuncs.cc: */
void print_sb_out_syns(void);
void display_sb_out_syns(int i);
void mak_sbac(int ctype,int n, int aden, double xpos, double ypos, double thetax,
		double thetay, double thetaz, int flip, int morph);

/* Functions in plot_funcs.cc: */

void plot_l_nod(int ct,int cn, int n, double lmin, double lmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_v_lnod(int ct,int cn, int labl, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_vm_lnod(int ct,int cn, int labl, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_v_nod(int ct,int cn, int n, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_v_nod(node *npnt, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_vm_nod(int ct,int cn, int n, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_vm_nod(node *npnt, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_i_nod(int ct,int cn, int n, double vmin, double vmax,int pcolor,const char *label,
						int plotnum, double psize);
void plot_ca_nod(int ct, int cn, int n,double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_ca_nod(int ct, int cn, int n, int sh, double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_ca_nod(node *n, int sh, double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_ca_syn(synapse *s, int sh, double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_ca_syn(synapse *s, double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_cabufb_nod(int ct, int cn, int n,double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_cabufb_nod(int ct, int cn, int n, int sh, double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_cas_nod(int ct, int cn, int n,double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_cas_nod(int ct, int cn, int n, int sh, double maxca, int pcolor, const char *label,
						int plotnum, double psize);
void plot_ph_nod(node *npnt, double maxnt, double minnt, int pcolor, 
						const char *label, int plotnum, double psize);
void plot_ph_nod(int ct, int cn, int n, double maxnt, double minnt, int pcolor, 
						const char *label, int plotnum, double psize);

void plot_synrate(int ct, int cn, int nod, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_synrate(int ct, int cn, int nod, double rmin, double rmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_synrate(int ct, int cn, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_synrate(int ct, int cn, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_synrate(int ct, int cn, double rmin, double rmax, int pcolor, int plotnum, const char *plname, double plsize);

void plot_synves (int ct, int cn, int ves, double fmin, double mmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_synves (int ct, int cn, double fmin, double mmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_syncond(int ct, int cn, int nod, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_syncond(int ct, int cn, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);

void plot_synrate(synapse *s, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);

void plot_synrate(synapse *s, double rmin, double rmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_syncond(synapse *s, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_syncondp(synapse *s, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);
void plot_synves(synapse *s, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize);


void plot_synrate_out(int ct, int cn, double rmin, double rmax, int colr);
void plot_synrate_out(int ct, int cn, double rmin, double rmax, int colr, const char *plname);
void plot_synrate_out(int ct, int cn, double rmin, double rmax, int colr, double plsize);
void plot_synrate_out(int ct, int cn, double rmin, double rmax, int colr, const char *plname, double plsize);
void plot_synrate_out(int ct, int cn, double rmin, double rmax, double fmax, int colr, const char *plname, double plsize);

void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, int colr, int prate);
void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, int colr, int prate, const char *plname);
void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, int colr, int prate, double plsize);
void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, int colr, int prate, const char *plname, double plsize);

void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, double fmax, int colr, int prate);
void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, double fmax, int colr, int prate, const char *plname);
void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, double fmax, int colr, int prate, double plsize);
void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, double fmax, int colr, int prate, const char *plname, double plsize);

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, int colr,int prate);
void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, int colr,int prate, const char *plname);
void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, int colr,int prate, int pves, double plsize);
void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, int colr,int prate, int pves, const char *plname, double plsize);

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, double fmax, int colr,int prate);
void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, double fmax, int colr,int prate, const char *plname);
void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, double fmax, int colr,int prate, int pves, double plsize);
void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax, double fmax, int colr,int prate, int pves, const char *plname, double plsize);

void plot_spike_rate(int ct, int cn, int n, int pcolor, const char *label, int plotnum, double psize);

void plot_chan(int ct, int cn, int n, int ctype, int stype, int param, int pval, double pmax, double pmin);
void plot_chan(int ct, int cn, int n, int ctype, int stype, int param, double pmax, double pmin);
void plot_chan_current(int ct, int cn, int n, int ctype, int stype, double pmax, double pmin);
void plot_chan_current(int ct, int cn, int n, int ctype, int stype, double mult, double pmax, double pmin);
void plot_chan_cond(int ct, int cn, int n, int ctype, int stype, double pmax, double pmin);
void plot_chan_cond(int ct, int cn, int n, int ctype, int stype, double mult, double pmax, double pmin);
void plot_chan_states(int ct, int cn, int n, int ctype, int stype, int param, double pmax, double pmin,
		const char *plname, int plnum, double psize);

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2);
int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, double vrev);
int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, int connum);
int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, int connum, double disti, double disto);
int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, double disti, double disto);
int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, double disti, double disto, 
		double xloc, double yloc);
double rsyn_avg (double sl, double time);
double isyn_avg (double sl, double time);
double isyn_tot (double sl, double time);
double gsyn_avg (double sl, double time);
double gsyn_tot (double sl, double time);
double gnmda_syn_tot (double sl, double time);
double casyn_avg (double sl, double time);

chantype *getchantype (int ctype, int stype);

/* function in synfuncs.cc, allows use of nval for synaptic parameters: */

int connect_synapse (int prect, int precn, int prenode, int postct, int postcn, int postnode);
int connect_synapse (node *npre, node *npost);
void setcelconn (int from_celltype, int from_cellnum, int to_celltype, int to_cellnum);
void make_synapse_type (node *prenod, node *postnod, int syntype); 
void synfuncs_cleanup(void);
void printhash(void);
void rehashnodes(void);

/* Function in ncsetexpt.cc: */

void setexpt(void);
