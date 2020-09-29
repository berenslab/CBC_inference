/* ncelem.h */
/* neural element structures for "nc" */

struct elem;
struct attrib;
struct chattrib;
struct cattrib;
struct node;
struct conlst;
struct comp;
struct conn;

struct elem {

	short int ctype;		/* element type used to establish id */
	short int newelem;		/* =1-> new elem, must process at run time */
	short int nocondens;		/* =1-> don't condense because elem is saved by ename */
	int elnum;			/* the element number */
	const char *elabl;		/* element label */
	int region;			/* region of the cell defined in morph filel */
	elem *next;			/* pointer to next element */
	elem *last;			/* pointer to last element */
	elem *hnext;			/* pointer to next elnum hashed element */
	elem *hlast;			/* pointer to last elnum hashed element */
	elem *hnnext;			/* pointer to next node hashed element */
	elem *hnlast;			/* pointer to last node hashed element */
	nodeint node1a;			/* nodes at ends of branch */
	nodeint node1b;			/* second dim */
	nodeint node1c;			/* third dim */
	nodeint node1d;			/* fourth dim */
	nodeint node2a;
	nodeint node2b;
	nodeint node2c;
	nodeint node2d;
	node *nodp1;			/* pointer to node 1 */
	node *nodp2;			/* pointer to node 2 */
	attrib *attpnt;			/* pointer to attributes */
	double jnoise;			/* Johnson noise for Rm */
	int   rsd;			/* noise seed */
	int modif;			/* nonzero -> modification of element */
	elem *lptr;			/* pntr to low level struc for modif */
	};

struct sphere: elem {
	double dia;			/* diameter of sphere */
	double Rm;			/* Rm for sphere membrane  */
	double Cm;			/* Cm for sphere membrane  */
	double vrev; 			/* membrane reversal potential    */
	double vrest; 			/* membrane initial resting pot   */
	};

struct cable: elem {
	double dia;			/* diameter of cable first end */
	double Rm;			/* Rm for branch */
	double Cm;			/* Cm for branch */
	double vrev; 			/* membrane reversal potential    */
	double vrest; 			/* membrane initial resting pot   */
	double Ri;			/* Ri for branch */
	double dia2;			/* diameter of cable second end */
	double length;		
	double cplam;			/* space constant for compartments */
	};

struct synapse: elem {

	double thresh;			/* threshold for release  */
	double *timec1; 			/* time const for syn filter */
	double *timec1h; 		/* time const for syn filter */
	double *timec2;			/* time const for syn filter */
	double *timec3;			/* time const for syn filter */
	short int nfilt1; 		/* num of low pass filters */
	short int nfilt1h; 		/* num of high pass filters  */
	short int nfilt2; 		/* num of low pass filters */
	short int nfilt3; 		/* num of low pass filters */
	double filt1hg; 		/* gain for high pass filter */
	double filt1ho; 		/* dc offset for high pass filter */
	double tfall2; 			/* falling tau for synapse */
	double tfall3; 			/* falling tau for synapse */
	double maxcond; 		/* max conductance for synapse */
	double rextern; 		/* external resistance */
	double maxsrate; 		/* max sustained rate */
	double rrpool;			/* starting size of readily releasible pool*/
	double mrrpool;			/* max size of readily releasible pool*/
	double rrpoolg;			/* readily releasible pool gain mult */
	double trconc;			/* transmitter calib factor */
	double mesgconc;		/* transmitter calib factor for 2nd mesg chan */
	double vsize; 			/* vesicle size when no ves noise    */
	double vrev; 			/* membrane reversal potential    */
	double nkd; 			/* Kd for synapse */
	double npow; 			/* Binding Hill coeff for synapse */
	double ngain; 			/* value of xfer func */
	double vgain; 			/* vesicle release gain */
	double caegain; 		/* calcium expon vesicle release gain */
	double ckd; 			/* cG Kd for synapse */
	double chc; 			/* cG Hill coeff for 2nd mesg chan */
	double cgain; 			/* cycG gain for synapse */
	double coff; 			/* offset for cascade subraction */
	int spost;			/* postsyn element to connect to */
	short int connum;		/* postsynaptic connection number */
	short int sens;			/* presyn sensitivity: Ca, V (default)*/
	short int ntact;		/* neurotrans action: open, close */
	short int secmsg;		/* has second messenger pathway */
	short int curve;		/* transfer curve type: lin, expon */
	short int mesg1;		/* second messenger */
	short int mesg2;		/* second messenger */
	int dyadelem;			/* element for synaptic dyad */
	};

struct photorec: elem {
	double xpos;
	double ypos;
	double zpos;
	double dia;
	double maxcond;
	double pigm;
	double pathl;
	double attf;
	double timec1;
	double loopg;
	double darknoise;
	double qc;
	double qcond;
	double linit;
	double unitary;
	int chrseed;
	int dkrseed;
	int phrseed;
	short int channoise;
	short int photnoise;
	short int filt;
	short int stimchan;
	char  save;
	char  restore;
	};

struct gapjunc: elem {
	double area;
	double specres;
	double vgain			/* forward parameters */;
	double voff;
	double taun;
	double rvgain;			/* reverse parameters */
	double rvoff;
	double rtaun;
	double gmax;
	double gnv;			/* non-volt sens conductance */
	short int rect;			/* gj is pure rectifier */
	short int modtyp;		/* type of modulation cyca/g */
	short int sign;			/* modulated open or close by cyca */
	nodeint nodeca;			/* node that controls modulation */
	nodeint nodecb;
	nodeint nodecc;
	nodeint nodecd;
	};

struct pannex: gapjunc {
	double atp_decr;		/* decrement for atp per time step */
	double atp;			/* initial atp conc */
	double atp_gain;		/* atp gain from pnx curr */
	double adp;			/* initial adp conc */
	double amp;			/* initial amp conc */
	double h;			/* initial h conc */
	double (*pnxfunc)(pannex*,double); /* function to process atp signal */
	};

struct electrode: elem {
	double r;
	double vrev; 			/* membrane reversal potential    */
	double vrest; 			/* membrane starting potential    */
	double c;
	double ccomp;			/* cap added to comp at tip, can be neg */
	double dia;
	double length;
	};

struct resistor: elem {
	double r;			/* resistance */
	};

struct diode: elem {
	double r;			/* diode */
	};

struct loadelem: elem {
	double r;			/* resistance to gnd */
	double vrev; 			/* membrane reversal potential    */
	double vrest; 			/* membrane starting potential    */
	};

struct capac: elem {
	double c;			/* capacitance */
	double vrest;
	};

struct batt: elem {
	double v;			/* voltage */
	};

struct vbuf: elem {
	double delay;			/* delay in msec */
	double offset;			/* offset in volts */
	double gain;			/* voltage gain */
	double tau;			/* time const */
	short int lphp;			/* low pass or high pass */
	};

struct nbuf: elem {
	double gain;			/* voltage -> nt gain */
	double offset;			/* offset in volts */
	double ntoffset;		/* offset in ntrans */
	int ntrans;			/* neurotrans */
	};

struct attrib {				/* part of an element */
					/* base attribute */
	short int ctype;		/* type of attrib */
	short int stype;		/* channel sub-type */
	attrib *attpnt;			/* pointer to attributes */
	elem *lptr;			/* pntr to low level struc for modif */
	};

struct lattrib: attrib {		/* part of an element */
	double unitary;
	double tauf;
	double n;
	};

struct mattrib: lattrib {		/* part of an element */

	double vrev;
	double voffsm;
	double voffsh;
	double maxcond;
	double rextern;
	double density;
	double ndensity;
	double caperm;
	double cakd;
	double cahc;
	double taua;
	double taub;
	double tauc;
	double taud;
	double taue;
	};

struct chattrib: mattrib {		/* part of an element */
					/* channel attribute */
	double d1;			/* chtyp for K2 chan */
	double d2;
	double k1;
	double k2;
	double trconc;			/* transmitter calib factor */
	};

struct cattrib: chattrib {		/* part of an element */
					/* calcium attribute */
	short int exch;
	short int pump;
	short int cabuf;
	short int cicr;
	short int ip3;
	double cao;
	double cai;
	double tcai;
	double kex;
	double ekm;
	double vmax;
	double pkm;

	double cas;			/* CICR parameters */
	double cas2;
	double vm2;
	double vm3;
	double ncicr;
	double mcicr;
	double pcicr;
	double kacicr;
	double kfcicr;
	double krcicr;
	double k1cicr;
	double k2cicr;
        double c1cicr;

	double ip3i			/* IP3 parameters */;    
	double bip3;
	double b2ip3;
	double vip3;
	double v2ip3;
	double v3ip3; 
	double v4ip3;
	double mip3;
	double hip3;
	double oip3;
	double a2ip3;
	double a3ip3;
	double mtypeip3;
	double d1ip3;
	double d2ip3;
	double d3ip3;
	double d4ip3;
	double k3ip3;


	double bkd;
	double bmax;
	double btot;
	double btoti;
	double cabnd;
	double sarea;			/* area of shells. */
	double mg;			/* magnesium conc. */
	short int cashell;
	short int caoshell;
	short int caflg;		/* set if cacomp at synapse */
	};

struct nattrib: lattrib {		/* part of an element */
					/* noise attribute */
	double vsize;			/* vesicle size */
	double vcov;			/* vesicle size stdev / mean */
	double cov;			/* vesicle interval stdev / mean */
	double refr;			/* refractory period */
	unsigned int rseed;
	};

struct node {
	short int ctype;		/* type of node = NODE */
	nodeint nodenm1;		/* node dimension 1 */
	nodeint nodenm2;		/* node dimension 2 */
	nodeint nodenm3;		/* node dimension 3 */
	nodeint nodenm4;		/* node dimension 4 */
	double xloc;			/* location of node */
	double yloc;
	double zloc;
	conlst *elemlst;		/* element list for connections */
	comp *comptr;			/* compartment which represents node */
	int  compnum;			/* number of comp representing node */
	node *next;			/* pointer to next node for foreach */
	node *last;			/* pointer to last node for foreach */
	node *hnext;			/* pointer to next node in hash table */
	node *hlast;			/* pointer to last node in hash table */
	node *hcnext;			/* pointer to next node in hash table */
	node *hclast;			/* pointer to last node in hash table */
	short int region;		/* region of cell defined in morph file */
	short int dendn;		/* dendrite number defined in morph file */
	short int label;		/* display node as label when directed */
	const char *labeltext;		/* text for label */
	};

struct conlst {				/* list of connections */
	conlst *next;			/* next in list */
	conn *conpnt;			/* pointer to conn */
	};

#define NLABELS 20
#define NULLVAL (-32768)
#define NULND (-1)

#define RUN_SAVE	1		/* save neural elem for runtime */
#define RUN_SAVE_MASK	0xf		/* mask for neural elem runtime */
#define RUN_SET_MASK	0x10		/* mask for get channel runtime */

