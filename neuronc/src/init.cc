/* init.c */

#include "nc.h"
#include "y.tab.h"
#include "ncval.h"
#include "ncsub.h"
#include "ncelem.h"
#include "nconst.h"
#include "ncsetvar.h"
#include "controlx.h"
#include <string.h>
#include <time.h>

iontab *ions = {NULL};		/* ion concentrations,  permeabilities */
int    ionz[NIONS] = {0};	/* ion valences */

double Ii = 0;			/* total ionic strength inside */
double Io = 0;			/* total ionic strength outside */
double ktA = 0;			/* const for act() */

extern stype varval[VVSIZE];
extern int varset;
extern char *progname;
extern char *infile;
extern const char *plotlabel;
extern int cumcomp;
extern int setdisp;
extern int setprmap;
extern int setvid, unsetvid;
extern int pipefl;
extern int setplsep;
extern int nonodes;
extern int interp_only;
extern char *rootframe;
extern double setlamcrit;
extern double vidsiz;

extern datum	nclog(datum &x), nclog10(datum &x);
extern datum	ncsqrt(datum &x);
extern datum 	ncexp(datum &x), ncexp10(datum &x);
extern datum	ncsin(datum &x), ncasin(datum &x), nccos(datum &x); 
extern datum	ncacos(datum &x), nctan(datum &x), nctanh(datum&x), ncsinh(datum&x);
extern datum	ncatan(datum &x), ncatan2(datum &x, datum &y);
extern datum    ncfabs(datum &a), ncinteger(datum &x);
extern datum    amax(datum &a, datum &b), amin(datum&a, datum&b);
extern datum    ncmax(datum &a, datum &b), ncmin(datum&a, datum&b);
extern datum    transpose(datum &a), matinv(datum &a), matsolve(datum &a, datum &b);
extern datum	matmul(datum&a, datum&b);
extern datum	ncrand(void), ncrrand(datum &ngen),xgetpid(void);
extern datum	ncpow(datum &x, datum &y), ncabs(double &x);
extern datum    elap_timee(void), system_speedd(void), simwait(datum &d1);
extern datum	setvarval(void), xisnum(datum &d), xstrlen(datum &d), set_q10s(void);
extern datum	printsym(void), print_ncversion(datum &d);
extern datum	xstrtok(datum &d1,datum &d2), xstrcmp(datum &d1, datum &d2);
extern datum	xstrstr(datum &d1,datum &d2);
extern datum	xindex(datum &d1,datum &d2);
extern datum	xrindex(datum &d1,datum &d2);
extern datum	xsubstr(datum &d1,datum &d2,datum &d3);
extern datum	xccstr(datum &d1), xatof(datum &d1);
extern datum	xtransf(datum &x, datum &y, datum &z);
extern datum	ncgaus(datum &x, datum &r), ncpoisdev (datum &x);
extern datum	ncbinomdev(datum &x, datum &n), ncgamdev(datum &x);
extern datum	ncgasdev(void), ncrgasdev(datum &ngen);
extern datum	minf(datum &v, datum &elname);
extern datum	mtau(datum &v, datum &elname);
extern datum	hinf(datum &v, datum &elname);
extern datum	htau(datum &v, datum &elname);
extern datum	(*ninf) (datum &v, datum &elname);
extern datum	(*ntau) (datum &v, datum &elname);
extern datum	(*dinf) (datum &v, datum &elname);
extern datum	(*dtau) (datum &v, datum &elname);
extern datum	(*cinf) (datum &v, datum &elname);
extern datum	(*ctau) (datum &v, datum &elname);
extern datum	(*chinf) (datum &v,datum &elname);
extern datum	(*chtau) (datum &v,datum &elname);
extern datum	(ctaua) (datum &v,datum &elname);
extern datum	(ctaub) (datum &v,datum &elname);
extern datum	(ctauc) (datum &v,datum &elname);
extern datum	(ctaud) (datum &v,datum &elname);
extern datum	(ctaue) (datum &v,datum &elname);
extern datum	(ctauf) (datum &v,datum &elname);
extern datum	(copychan) (datum &sc, datum &ss, datum &dc, datum &ds, 
					datum &nstates);
extern datum	(setchan_ntrans)(datum &sc,datum &ss,datum &state, datum &val);
extern datum	(setchan_cond) (datum &sc,datum &ss,datum &state, datum &val);
extern datum	(setchan_trans)(datum &sc,datum &ss,datum &state, datum &trans, 
					datum &val);
extern datum	(setchan_trate)(datum &sc,datum &ss,datum &state, datum &trans, 
					datum &val);
extern datum	(setchan_mul) (datum &sc, datum &ss,datum &state, datum &trans, 
					datum &val);
extern datum	(setchan_rateo)(datum &sc,datum &ss,datum &state, datum &trans, 
					datum &val);

Symbol *timeptr;
Symbol *ncompptr;
Symbol *plotiptr;

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "stdplt.h"
#include "gr.h" 

int getpid(void);
// int time(int*);
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

#include "ncio.h"

double lookpv(const char *s);
Symbol *lookup(const char *s);
Symbol *lookups(const char *s);
void erasymtab();
void setval(Symbol *sym, double val, int vtyp);
void setval(Symbol *sym, int val, int vtyp);
void setval(Symbol *sym, char *cptr, int vtyp);
void setval(Symbol *sym, double *val, int vtyp);	// set pointers for vinit()
void setval(Symbol *sym, int *val, int vtyp);
void setval(Symbol *sym, char **cptr, int vtyp);
void dinit();
void vinit();
void varcopy();
void calcvar();
void initrec();
void setrand (int rseed);
void iinit();
double getval(const char *s);
void setvariable (stype *var);
char *emalloc(unsigned int n);

/*--------------------------------------------------------*/

static struct {		/* Keywords */
	const char *name;
	int kval;
} keywords[] = {
	"proc",		PROC,
	"func",		FUNC,
	"return",	RETURN,
	"exit",		EXIT,
	"debugf",	DEBUGF,
	"if",		IF,
	"else",		ELSE,
	"while",	WHILE,
	"print",	PRINT,
	"printf",	PRINTF,
	"sprintf",	SPRINTF,
	"fprintf",	FPRINTF,
	"fopen",	FOPEN,
	"fclose",	FCLOSE,
	"scanf",	SCANF,
	"fscanf",	FSCANF,
	"sscanf",	SSCANF,
	"fgets",	FGETS,
	"fgetc",	FGETC,
	"fputc",	FPUTC,
	"getflds",	GETFLDS,
	"read",		READ,
	"fread",	FREAD,
	"freads",	FREADS,
	"fwrite",	FWRITE,
	"unlink",	UNLINK,
	"notinit",	NOTINIT,
	"varnum",	VARNUM,
	"varstr",	VARSTR,
	"varchr",	VARCHR,
	"dim",		DIM,
	"sizeof",	SIZEOF,
	"dims",		DIMS,
	"fft",		FFT,
	"ifft",		IFFT,
	"pfft",		PFFT,
	"acov",		ACOV,
	"lmfit",	LMFIT,
	"lmfit2d",	LMFIT2D,
	"foreach",	FOREACH,
	"for",		FOR,
	"break",	BREAK,
	"continue",	CONTINUE,
	"edit",		EDIT,
	"system",	SYSTEM,
	"include",	INCLUDE,
	"local",	LOCAL,
	"graph",	GRAPH,
	"pen",		PEN,
	"restart",	RESTART,
	"init",		INIT,
	"X",		X,
	"Y",		Y,
	"Z",		Z,
	"max",		MAX,
	"min",		MIN,
	"char",		CHAR,
	"Char",		CCHAR,
	"vpen",		VPEN,
	"plname",	PLNAME,
	"plnum",	PLNUM,
	"plsize",	PLSIZE,
	"plval",	PLVAL,
	"plarr",	PLARR,
	"gmove",	GMOVE,
	"gdraw",	GDRAW,
	"grdraw",	GRDRAW,
	"grmove",	GRMOVE,
	"gpen",		GPEN,
	"gvpen",	GVPEN,
	"gcwidth",	GCWID,
	"gsize",	GSIZ,
	"gdash",	GDASH,
	"gwindow",	GWINDOW,
	"gframe",	GFRAME,
	"grot",		GROT,
	"gorigin",	GORIG,
	"gcrot",	GCROT,
	"gtext",	GTEXT,
	"glabel",	GLABEL,
	"gcirc",	GCIRC,
	"grect",	GRECT,
	"gpurge",	GPURGE,
	0,		0
};

static struct {		/* Keywords */
	const char *name;
	int kval;
} simwords[] = {
	"element",	ELEMENT,
	"matching",	MATCHING,
	"except",	EXCEPT,
	"comps",	COMPS,
	"range",	RANGE,
	"cable",	CABLE,
	"sphere",	SPHERE,
	"electrode",	ELECTRODE,
	"chan",		CHAN,
	"Na",		NA,
	"K",		K,
	"KCa",		KCa,
	"ClCa",		ClCa,
	"type",		TYPE,
	"ntype",	NTYPE,
	"stype",	STYPE,
	"nstate",	NSTATE,
	"density",	DENSITY,
	"ndensity",	NDENSITY,
	"tau",		TAU,
	"taum",		TAUM,
	"tauh",		TAUH,
	"taun",		TAUN,
	"taua",		TAUA,
	"taub",		TAUB,
	"tauc",		TAUC,
	"taud",		TAUD,
	"taue",		TAUE,
	"tauf",		TAUF,
	"Ca",		CA,
	"cacomp",	CACOMP,
	"cao",		CAO,
	"cai",		CAI,
	"tcai",		TCAI,
	"cbound",	CBOUND,
	"cabuf",	CABUF,
	"cabufb",	CABUFB,
	"btot",		BTOT,
	"btoti",	BTOTI,
	"cshell",	CASHELL,
	"caoshell",	CAOSHELL,
	"caperm",	CAPERM,
	"cicr",		CICR,
	"ip3",		IP3,
	"vmax",		VMAX,

	"cas",		CAS,
	"cas2",		CAS2,
	"ip3i",		IP3I,
	"vm2",		VM2,
	"vm3",		VM3,
	"ka",		KA,
	"kf",		KF,
	"kr",		KR,
	"bip3",		BIP3,
	"vip3",		VIP3,
	"c1",           C1,
	"oi",		OI,
	"a2",		A2,
	"a3",		A3,
	"b2",		B2,
	/* "d1",	D1, */
	/* "d2",	D2, */
	"d3",		D3,
	"d4",		D4,
	"k3",		K3,
	"v2",		V2,
	"v3",		V3,
	"v4",		V4,
	"mIP3",		MIP3,
	"hIP3",         HIP3,
	"mtype",	MTYPE,

	"km",		KM,
	"kex",		KEX,
	"capump",	CAPUMP,
	"caexch",	CAEXCH,
	"d1",		D1,
	"d2",		D2,
	"k1",		K1,
	"k2",		K2,
	"mg",		MG,
	"gj",		GJ,
	"pnx",		PNX,
	"atp",		ATP,
	"pH",		PH,
	"puff",		PUFF,
	"GLU",		GLU,
	"glu",		GLU,
	"AMPA",		AMPA,
	"ampa",		AMPA,
	"KAINATE",	KAINATE,
	"kainate",	KAINATE,
	"NMDA",		NMDA,
	"nmda",		NMDA,
	"CNQX",		CNQX,
	"cnqx",		CNQX,
	"GABA",		GABA,
	"gaba",		GABA,
	"BIC",		BIC,
	"bic",		BIC,
	"PTX",		PTX,
	"ptx",		PTX,
	"GLY",		GLY,
	"gly",		GLY,
	"STRY",		STRY,
	"stry",		STRY,
	"camp",		CAMP,
	"cAMP",		CAMP,
	"cgmp",		CGMP,
	"cGMP",		CGMP,
	"resp",		RESP,
	"syn2",		SYN2,
	"mesgin",	MESGIN,
	"mesgout",	MESGOUT,
	"synapse",	SYNAPSE,
	"load",		LOAD,
	"resistor",	RESISTOR,
	"diode",	DIODE,
	"transducer",	VTRANSDUCER,
	"itransducer",	ITRANSDUCER,
	"rod",		ROD,
	"cone",		CONE,
	"chr",		CHR,
	"linit",	LINIT,
	"spost",	SPOST,
	"dyad",		DYAD,
	"sdyad",	SDYAD,
	"scurve",	SCURVE,
	"numdyad",	NUMDYAD,
	"rsd",		RSD,
	"initrand",	INITRAND,
	"maxcond",	MAXCOND,
	"maxsrate",	MAXSRATE,
	"mrrpool",	MRRPOOL,
	"rrpool",	RRPOOL,
	"rrpoolg",	RRPOOLG,
	"length",	LENGTH,
	"area",		AREA,
	"dia",		DIA,
	"dia1",		DIA,
	"dia2",		DIA2,
	"radius",	RADIUS,
	"width",	WIDTH,
	"stim",		STIM,
	"bar",		BAR,
	"spot",		SPOT,
	"rect",		RECT,
	"sine",		SINE,
	"sineann",	SINEANN,
	"gabor",	GABOR,
	"Gabor",	GABOR,
	"windmill",	WINDMILL,
	"checkerboard",	CHECKERBOARD,
	"simage",	SIMAGE,
	"xenv",		XENV,
	"yenv",		YENV,
	"sq",		SQ,
	"xenon",	XENON,
	"sun",		SUN,
	"tungsten",	TUNGSTEN,
	"wavel",	WAVEL,
	"pigm",		PIGM,
	"pathl",	PATHL,
	"attf",		ATTF,
	"filt",		FILT,
	"save",		SAVE,
	"restore",	RESTORE,
	"gndbatt",	GNDBATT,
	"batt",		BATT,
	"cap",		CAP,
	"gndcap",	GNDCAP,
	"node",		NODE,
	"node1a",	NODE1A,
	"node1b",	NODE1B,
	"node1c",	NODE1C,
	"node1d",	NODE1D,
	"node2a",	NODE2A,
	"node2b",	NODE2B,
	"node2c",	NODE2C,
	"node2d",	NODE2D,
	"vclamp",	VCLAMP,
	"cclamp",	CCLAMP,
	"loc",		LOC,
	"xloc",		XLOC,
	"yloc",		YLOC,
	"zloc",		ZLOC,
	"center",	CENTER,
	"sscale",	SSCALE,
	"dur",		DUR,
	"tfreq",	TFREQ,
	"drift",	DRIFT,
	"orient",	ORIENT,
	"sphase",	SPHASE,
	"contrast",	CONTRAST,
	"start",	START,
	"inten",	INTEN,
	"mask",		MASK,
	"stimchan",	STIMCHAN,
	"backgr",	BACKGR,
	"blur",		BLUR,
	"scatter",	SCATTER,
	"file",		SFILE,
	"plot",		PLOT,
	"V",		V,
	"Vm",		VM,
	"I",		I,
	"Im",		IM,
	"L",		L,
	"P",		P,
	"S",		S,
	"FA0",		FA0,
	"FA1",		FA1,
	"FA2",		FA2,
	"FA3",		FA3,
	"FA4",		FA4,
	"FA9",		FA9,
	"FH0",		FH0,
	"FH1",		FH1,
	"FH2",		FH2,
	"FH3",		FH3,
	"FH4",		FH4,
	"FB0",		FB0,
	"FB1",		FB1,
	"FB2",		FB2,
	"FB3",		FB3,
	"FB4",		FB4,
	"FC0",		FC0,
	"FC1",		FC1,
	"FC2",		FC2,
	"FC3",		FC3,
	"FC4",		FC4,
	"FC9",		FC9,
	"G",		G,
	"M",		M,
	"H",		H,
	"G0",		G0,
	"G1",		G1,
	"G2",		G2,
	"IP",		IP,
	"IE",		IE,
	"IPE",		IPE,
	"display",	DISPLAY,
	"xrot",		XROT,
	"yrot",		YROT,
	"zrot",		ZROT,
	"size",		SIZE,
	"window",	WINDOW,
	"rmove",	RMOVE,
	"dscale",	DSCALE,
	"color",	COLOR,
	"cmap",		CMAP,
	"newpage",	NEWPAGE,
	"hide",		HIDE,
	"calibline",	CALIBLIN,
	"rm",		RM,
	"ri",		RI,
	"rs",		RS,
	"rg",		RG,
	"cm",		CM,
	"jnoise",	JNOISE,
	"vbuf",		BUF,
	"nbuf",		NBUF,
	"opamp",	OPAMP,
	"chnoise",	CCHNOISE,
	"photnoise",	PHOTNOISE,
	"channoise",	CHANNOISE,
	"darknoise",	DARKNOISE,
	"vesnoise",	VESNOISE,
	"loopg",	LOOPG,
	"N",		N,
	"unit",		UNIT,
	"vsize",	VSIZE,
	"vcov",		VCOV,
	"refr",		REFR,
	"CoV",		COV,
	"open",		OPEN,
	"close",	CLOSE,
	"sens",		SENS,
	"vrev",		VVREV,
	"vrest",	VREST,
	"cplam",	CPLAM,
	"thresh",	THRESH,
	"delay",	DELAY,
	"nfilt1",	NFILT1,
	"nfilt1h",	NFILT1H,
	"nfilt2",	NFILT2,
	"nfilt3",	NFILT3,
	"timec1",	TIMEC1,
	"timec1h",	TIMEC1H,
	"timec2",	TIMEC2,
	"tfall2",	TFALL2,
	"timec3",	TIMEC3,
	"tfall3",	TFALL3,
	"hgain",	HGAIN,
	"cgain",	CGAIN,
	"coff",		COFF,
	"ckd",		CKD,
	"chc",		CHC,
	"rev",		REV,
	"bgain",	BGAIN,
	"boffset",	BOFFSET,
	"vgain",	VGAIN,
	"gmax",		GMAX,
	"gnv",		GNV,
	"kd",		KD,
	"cakd",		CAKD,
	"cahc",		CAHC,
	"hcof",		HCOF,
	"linear",	LINEAR,
	"expon",	EXPON,
	"trconc",	TRCONC,
	"mesgconc",	MESGCONC,
	"exist",	EXIST,
	"numconn",	NUMCONN,
	"numsyn",	NUMSYN,
	"connect",	CONNECT,
	"conn",		CONNECT,
	"only",		ONLY,
	"at",		AT,
	"to",		TO,
	"edist",	EEDIST,
	"e3dist",	E3DIST,
	"e2dist",	E2DIST,
	"ezdist",	EZDIST,
	"efrac",	EFRAC,
	"ndist",	N3DIST,
	"n3dist",	N3DIST,
	"n2dist",	N2DIST,
	"nzdist",	NZDIST,
	"within2d",	WITHIN2D,
	"within3d",	WITHIN3D,
	"elabl",	ELABL,
	"ename",	ENAME,
	"modify",	MODIFY,
	"offset",	OFFSET,
	"offsetm",	OFFSETM,
	"offseth",	OFFSETH,
	"put",		PUT,
	"run",		RUN,	
	"step",		STEP,	
	"elimit",	ELIMIT,
	"erase",	ERASE,
	"model",	MODEL,
	"gausnn",	GAUSNN,
	"ginfo",	GINFO,
	"nnd",		NND,
	"nnstd",	NNSTD,
	"reg",		REG,
	0,		0
};
	
static struct {		/* control vars */
	const char *name;
} controls[] = {
	"vidmode",
	"debug",
	"debugz",
	"silent",		
	"info",		
	"rseed",
	"srseed",
	"chrseed",
	"dkrseed",
	"phrseed",
	"stderr",
	"stdin",
	"stdout",
	"nonodes",
	"version",
	"ncversion",
	"makestim",
	"infile",
	"progname",
	0,0
};

static struct {		/* control vars */
	const char *name;
} scontrols[] = {
	"vk",
	"vna",
	"vcl",
	"vglu",
	"dnai",
	"dnao",
	"dko",
	"dki",
	"dpnaca",
	"dpkca",
	"dpkna",
	"dpcana",
	"dpnak",
	"dpcak",
	"dpcaampa",
	"dpcacgmp",
	"dpcanmda",
	"dpcasyn2",
	"dcli",
	"dclo",
	"deleccap",
	"dcm",
	"dnadens",
	"dkdens",
	"dcadens",
	"dcavoff",
	"dcaspvrev",
	"ddca",
	"ddcao",
	"ddiacao",
	"dcoo",
	"dmgo",
	"dcao",
	"dcai",
	"dtcai",
	"dcabnd",
	"dcabt",
	"dcabti",
	"dcabf",
	"dcabr",
	"dcashell",
	"dcaoshell",
	"dcasarea",
	"dcapkm",
	"dcavmax",
	"dicafrac",
	"dcakex",
	"diexchfrac",
	"dmaxca", 
	"dcaoffs", 
	"dcataum",
	"dcatauh",
	"dkcataum",
	"dkcatauh",
	"dnatauf",
	"dktauf",
	"dcatauf",
	"dqm",
	"dqh",
	"dqn",
	"dqna",
	"dqnb",
	"dqca",
	"dqd",
	"dqkca",
	"dqsyn",
	"dqcavmax",
	"dqdc",
	"dqcab",
	"dqc",
	"dqri",
	"dqrm",
	"dqrec",
	"dqcrec",
	"drg",
	"dri",
	"drm",
	"drs",
	"dfta",
	"dftah",
	"dsfa",
	"dftb",
	"dsfb",
	"dscu",
	"dsvn",
	"dscavg",
	"dscaeg",
	"dvsz",
	"dvg",
	"dst",
	"dsc",
	"dmaxsyn", 
	"dskd",
	"dshc",
	"dstr",
	"dsmsgc",
	"dsmrrp",
	"dsrrg",
	"dsrrp",
	"dsms",
	"dcakd",
	"dcahc",
	"dckd",
	"dchc",
	"dsn", 
	"dsg", 
	"dmaxrod", 
	"dmaxcon", 
	"dmaxna", 
	"dmaxampa", 
	"dmaxnmda", 
	"dnmdamg", 
	"dmaxgaba", 
	"dnaoffsm", 
	"dnaoffsh", 
	"dnataum",
	"dnatauh",
	"dktaum",
	"dktauh",
	"dsyntau",
	"dsyntauf",
	"dgjoff", 
	"dgjtau", 
	"dgjatpdec", 
	"dgjnv", 
	"dbd1", 
	"dbd2", 
	"dbk1", 
	"dbk2", 
	"dsd1", 
	"dsd2", 
	"dsk1", 
	"dsk2", 
	"dmaxk", 
	"dmaxclca", 
	"dclcavs", 
	"dclcavsc", 
	"dsclcac", 
	"dnau",
	"dku", 
	"dkau", 
	"dkcabu", 
	"dkcasu", 
	"dclcasu", 
	"dkihu", 
	"dhcnu", 
	"dcatu", 
	"dcalu", 
	"dampau", 
	"dcgmpu", 
	"dgabau", 
	"dglyu", 
	"dpru", 
	"dkoffsn", 
	"dkoffsh", 
	"dbasetdc",
	"dbasetc",
	"dbasetca",
	"dbasetsyn",
	"dratehhm",
	"dratehhh",
	"dratehhna",
	"dratehhnb",
	"dqeff",
	"dsintinc",
	"dsintres",
	"tempcel",
	"silent",
	"label",
	"disp",	
	"prmap",	
	"stim_elem",	
	"complam",
	"lamcrit",
	"relincr",
	"timinc",	
	"stiminc",	
	"time",	
	"srtimestep",	
	"ncomps",	
	"endexp",	
	"ploti",
	"crit",		
	"plmax",
	"plmin",
	"relax",		
	"plsep",		
	"djnoise",		
	"rightyax",		
	"allowrev",		
	"dash",		
	"euler",		
	"implicit",		
	"implfk",
	"stimelem",		
	"vcolor",		
	"lcolor",		
	"rcolor",		
	"cacolor",		
	"nacolor",		
	"sgcolor",		
	"srcolor",		
	"nozerochan",		
	"oldsynapz",		
	"use_ghki",		
	"fread_expr",
	"synaptau",		
	"calcnernst",
	"setxmax",
	"setxmin",
	"setymax",
	"setymin",
	"numstimchan",
	"unmaskthr",
	"stimonh",
	"stimonl",
        "betaIP3",
        "v1IP3",
        "kfCICR",
        "vm2CICR",
        "vm3CICR",
	"k1CICR",
        "k2CICR",
        "krCICR",
        "kaCICR",
        "nCICR",
        "mCICR",
        "pCICR",
	"c1CICR",

	"oIP3",
	"mtypeIP3",
	"a2IP3",
	"a3IP3",
	"b2IP3",
	"d1IP3",
	"d2IP3",
	"d3IP3",
	"d4IP3",
	"k3IP3",
	"v4IP3",
	"v3IP3",
	"v2IP3",

	0,0
};


static struct {		/* Constants */
	const char *name;
	double cval;
} chtyp[] = {
	"PI", 		3.14159265358979323846,
	"E",		2.71828182845904523536,
	"GAMMA",	0.57721566490153286060,	/* Euler */
	"DEG",	       57.29577951308232087680,	/* deg/radian */
	"PHI",		1.61803398874989484820,	/* golden ratio */
	"UNINIT",	LARGENUM,
	(const char *)0,	0
};

struct builtin {		/* Built-ins */
	const char *name;
	Fpa  func;
};

static builtin builtins[] = {
	"log",	(Fpa )nclog,	/* checks range */
	"log10",(Fpa )nclog10,	/* checks range */
	"exp",	(Fpa )ncexp,	/* checks range */
	"exp10",(Fpa )ncexp10,	/* checks range */
	"pow",	(Fpa )ncpow,	/* checks range */
	"sin",	(Fpa )ncsin,	/* checks range */
	"asin",	(Fpa )ncasin,	/* checks range */
	"cos",	(Fpa )nccos,	/* checks range */
	"acos",	(Fpa )ncacos,	/* checks range */
	"tan",	(Fpa )nctan,	/* checks range */
	"tanh",	(Fpa )nctanh,	/* checks range */
	"sinh",	(Fpa )ncsinh,	/* checks range */
	"atan",	(Fpa )ncatan,	/* checks range */
	"atan2",(Fpa )ncatan2,	/* checks range */
	"sqrt",	(Fpa )ncsqrt,	/* checks range */
	"gauss",(Fpa )ncgaus,	/* checks range */
	"poisdev", (Fpa )ncpoisdev,
	"gamdev", (Fpa )ncgamdev, 
	"binomdev", (Fpa )ncbinomdev,/* checks range */
	"gasdev", (Fpa )ncgasdev,/* checks range */
	"rgasdev", (Fpa )ncrgasdev,/* checks range */
	"int",	(Fpa )ncinteger,
	"abs",	(Fpa )ncfabs,
	"amax",	(Fpa )amax,
	"amin",	(Fpa )amin,
	"fmax",(Fpa )ncmax,
	"fmin",(Fpa )ncmin,
	"transpose",(Fpa )transpose,
	"matmul",(Fpa )matmul,
	"matinv",(Fpa )matinv,
	"matsolve",(Fpa )matsolve,
	"rand",	(Fpa )ncrand,
	"rrand",(Fpa )ncrrand,
	"setvar", (Fpa )setvarval,
	"printsym", (Fpa )printsym,
	"getpid", (Fpa )xgetpid,
	"print_version", (Fpa )print_ncversion,
	"elap_time", (Fpa )elap_timee,
	"system_speed", (Fpa )system_speedd,
	"simwait", (Fpa )simwait,
	"isnum",  (Fpa )xisnum,
	"strlen", (Fpa )xstrlen,
	"strtok", (Fpa )xstrtok,
	"strcmp", (Fpa )xstrcmp,
	"strstr", (Fpa )xstrstr,
	"index", (Fpa )xindex,
	"rindex", (Fpa )xrindex,
	"substr", (Fpa )xsubstr,
	"ccstr", (Fpa )xccstr,
	"atof",  (Fpa )xatof,
	"set_q10s",(Fpa )set_q10s,
	(const char *)0, (Fpa )0
};

static builtin sbuiltins[] = {
	"minf", (Fpa )minf,
	"mtau", (Fpa )mtau,
	"hinf", (Fpa )hinf,
	"htau", (Fpa )htau,
	"ninf", (Fpa )ninf,
	"ntau", (Fpa )ntau,
	"dinf", (Fpa )dinf,
	"dtau", (Fpa )dtau,
	"cinf", (Fpa )cinf,
	"ctau", (Fpa )ctau,
	"chinf",(Fpa )chinf,
	"chtau",(Fpa )chtau,
	"ctaua",(Fpa )ctaua,
	"ctaub",(Fpa )ctaub,
	"ctauc",(Fpa )ctauc,
	"ctaud",(Fpa )ctaud,
	"ctaue",(Fpa )ctaue,
	"ctauf",(Fpa )ctauf,
	"copychan",(Fpa )copychan,
	"setchan_ntrans",(Fpa )setchan_ntrans,
	"setchan_cond",(Fpa )setchan_cond,
	"setchan_trans",(Fpa )setchan_trans,
	"setchan_trate",(Fpa )setchan_trate,
	"setchan_mul",(Fpa )setchan_mul,
	"setchan_rateo",(Fpa )setchan_rateo,
	"transf",(Fpa )xtransf,
	(char *)0, (Fpa )0
};

void init(void)	/* install constants and built-ins in table */
{
	int i;
	Symbol *s;
	int n, setvar(void);

        if (strcmp(progname,"nci")==0) { interp_only = 1; }
	erasymtab();			/* erase previous defs, in "symbol.c" */
	for (i=0; keywords[i].name; i++)
		s= install(keywords[i].name, keywords[i].kval, 0.0);
	for (i=0; controls[i].name; i++)
		s= install(controls[i].name, VAR, 0.0);
	for (i=0; chtyp[i].name; i++) {
		s = install(chtyp[i].name, CONST, chtyp[i].cval);
		setval (s,chtyp[i].cval,NUMBER);
	}
	for (i=0; builtins[i].name; i++) {
		s = install(builtins[i].name, BLTIN, 0.0);
		s->ptr = builtins[i].func;
	}
	if (!interp_only) {	/* simulator things */
	  for (i=0; simwords[i].name; i++)
		s= install(simwords[i].name, simwords[i].kval, 0.0);
	  for (i=0; scontrols[i].name; i++)
		s= install(scontrols[i].name, VAR, 0.0);
	  for (i=0; sbuiltins[i].name; i++) {
		s = install(sbuiltins[i].name, BLTIN, 0.0);
		s->ptr = sbuiltins[i].func;
	  }
	  vinit();		/* initialize pointers to variables in symbol table */
	  dinit();		/* initialize default values in symbol table and C variables */
	  varcopy();
	}
	iinit();		/* interpreter things */
	n = setvar();		/* set variables from command-line */
	if (!interp_only) {	/* simulator things */
	  initrec();		/* receptor constants */
	}
	frame (rootframe);
}

/****************************************/

int setint (double val) 

/* Set integer from double value, and maintain sign. */

{
    int retval;

     while (val > MAXINT) { val *= 0.5; val += 1.0; }
     while (val < MININT) { val *= 0.5; val -= 1.0; }
     retval   = int(val);
     return (retval);
}

/****************************************/

#define NR NUMBER

void setrseed(void)

/* run before any random sequence is generated */

{
    static int oldrseed = 0;

  rseed  = setint(getval ("rseed"));
  if (rseed==0) rseed = 1;
  if (rseed < 0) rseed = (getpid()*65536*+time((time_t*)NULL)) & 0x7fffffff;
   setval (lookup("rseed"),rseed,NR); /* random seed on comnd line*/
  if (rseed != oldrseed) {
	oldrseed = rseed;
	setrand (rseed);		/* set up random number seed */
  }
}

/****************************************/

int setvar(void)

/* set variables from command line.
    This overrides what values are set
    in a program.  You call this function
    "setvar()" in your C or nc script file.

*/

{
   int i;
   Symbol *s;

   /* Set values in interpreter symbol table, which are normally the default values */

   for (i=0; i<varset; i++) {
     if (varval[i].name) {
       if ((s=lookup(varval[i].name)) == 0)
           s=install(varval[i].name,VAR,varval[i].val);
       if (s->type==VAR) { 		/* set var from comnd line */
	   if (varval[i].type==VSTRING)
		 setval(lookup(varval[i].name),varval[i].cptr,STRING);
	   else  setval(lookup(varval[i].name),varval[i].val,NUMBER);
       }
     }
     else break;
   }

   /* Set variables in your C program (not for interpreter scripts) */
   /* see "ncsetvar.cc for examples */

   for (i=0; i<varset; i++) {
      setvariable(&varval[i]);	
   }

   setval (lookup("vidmode"), vidmode,NR);	/* video display flag */
   if (setdisp>=0) {
     disp = setdisp;	/* command line display flag */
     setval (lookup("disp"), disp,NR); /* display flag */
   }
   if (setprmap>=0) {
     prmap = setprmap;	/* command line print flag */
     setval (lookup("prmap"),prmap,NR); /* display flag */
   }

   if (setlamcrit>=0) {
     lamcrit = setlamcrit;   /* criterion for condensing comps */
     setval (lookup("lamcrit"),lamcrit,NR); 
   }

   if (setplsep) {
      plsep = setplsep;
      setval (lookup("plsep"),plsep,NR);  /* sep plots */
   }
   if (drseed) {				/* randomize if negative */
     setval (lookup("rseed"),drseed,NR); /* random seed on comnd line*/
     setrseed();
   }
   if (debug) {					/* -y debug  */
     setval (lookup("debug"),debug,NR);
   }
   if (debugz) {				/* -z debugz */
     setval (lookup("debugz"),debugz,NR); 
   }
   calcvar();

   return i;
}

/****************************************/

datum setvarval(void)

{
   static datum d={0};

  d.val = setvar();	/* call setvar() above to get command line vars */
  d.vtype = NUMBER;
  return d;
}

/****************************************/

double act (double Is) 
        
/* Calculate activity of monovalent ion, given total ionic strength (Is) */
/*  Use Davies Equation */

/* #define A 0.512		    /* in water at 25 deg C */
/* #define AF 2635.9917             /* = .512 * (KELVIN+25)^1.5 */
/* ktA   = AF * pow(ktemp,-1.5);    /* gives .512 @ 25 deg C, for act() */

{
     double g,sqrtI;

  sqrtI = sqrt(Is);
  g = exp (LN10 * -ktA * (sqrtI / (1 + sqrtI) - 0.2 * Is));
  return g;
}

/****************************************/

double act2 (double Is) 
        
/* Calculate activity of divalent ion, given total ionic strength (Is) */
/*  Use Davies Equation */

/* #define A 0.512		    /* in water at 25 deg C */
/* #define AF 2635.9917             /* = .512 * (KELVIN+25)^1.5 */
/* ktA   = AF * pow(ktemp,-1.5);    /* gives .512 @ 25 deg C, for act() */

{
     double g,sqrtI;

  sqrtI = sqrt(Is);
  g = exp (LN10 * -ktA * 4 * (sqrtI / (1 + sqrtI) - 0.2 * Is));
  return g;
}

/****************************************/

double nernstv (double cratio)

/* calculate Nernst potential for a monovalent cation */

{
   return (RF * ktemp * log (cratio * act(Io)/act(Ii)) ); 
}

/****************************************/

double nernstc (double vrev, double co, int z)

/* calculate internal concentration for a monovalent cation */

{
   return (co*act(Io) / exp(z*vrev * FR / ktemp)) / act(Ii); 
}

/****************************************/

double ghki (double v, int ion)

/* compute GHK current equation, z = valence of ion, for Cl, z = -1 */

{
    double df, i, io, p, z;

   io = ions->iono[ion];
   p = ions->ionp[ion];                      // matches conductance at large voltages
   z = ionz[ion];
   df = exp (-z * FFR * v / ktemp);
   i = p / io * z * z * (ions->ioni[ion] - io * df) / (1 - df);
   return i;
}

/****************************************/

iontab *init_ions(iontab *itab)
{
     int i;
     iontab *ions;

  if (itab!=NULL) return itab;
  if ((ions=(iontab *)emalloc(sizeof(iontab)))==NULL)
      ncfprintf (stderr,"makchan: can't alloc iontab\n");
  else { 
     for (i=0; i<NIONS; i++) {
         ions->ioni[i] = ions->iono[i] = ions->ionp[i] = 0.0;  
     }
  }
  return ions;
}

/****************************************/

void setion(void)

/* Set concentrations of internal ions. */
/* Compute total Ionic strength. */

{
    int i;
    double isumi, isumo, ionv2;
    double *ip;

 for (i=0; i<NIONS; i++) {
   ionz[i] = 1;
 }
 ions = init_ions(ions);

 ions->iono[PNA] = dnao;
 ions->ioni[PNA] = dnai;
 ions->ionp[PNA] = 1;
 ions->iono[PK]  = dko;
 ions->ioni[PK]  = dki;
 ions->ionp[PK]  = 1;  // for GHK curr eqn
 ions->iono[PCA] = dcao;
 ions->ioni[PCA] = dcai;
 ionz[PCA] = 2;
 ions->ionp[PCA] = 0;
 ions->iono[PCL] = dclo;
 ions->ioni[PCL] = dcli;
 ionz[PCL] = -1;
 ions->ionp[PCL] = 1;  // for GHK curr eqn

 isumi = isumo = 0;
 for (i=0; i<NIONS; i++) {
   ionv2 = ionz[i] * ionz[i];
   isumi += ions->ioni[i] * ionv2;
   isumo += ions->iono[i] * ionv2;
 }
 Io = 0.5 * isumo;		/* Ionic strength outside */
 /* Ii = 0.5 * isumi;		/* Ionic strength inside */
 Ii = .140;			/* assume Ionic strengths are equal */

  /* ncfprintf (stdout,"Io %g Ii %g\n",Io,Ii); /* */

}

/****************************************/

void setvrevs(void)

/* If Nernst potential for ions is not set, then create it from */
/*  their internal and external concentrations. */

/* Initially, if Nernst potential is set, then set the internal ion concs.  */
/*  from the external value and Nernst potential. */

/* If "calcnernst" = 1, when temperature or ion concentrations */
/* change, calculate new Nernst potential (vna, vk, and vcl). */

/* If "calcnernst" = 0, when temperature or Nernst potential values */
/* change, calculate new internal ion concentrations. */

/* If "calcnernst" is between 0 and 1, when temperature or ion
 * concentrations change, change Nernst potential
 * by the amount specified in calcnernst, and then make up the rest
 * of the change from ion concentration change */

{
 if ((getval("vna")==LARGENUM) || calcnernst) {
     double ovna=vna;
   vna = nernstv(dnao/dnai);
   if (ovna==LARGENUM) ovna = vna;
   if ((calcnernst > 0) && (calcnernst < 1)) {
      vna = ovna + (vna-ovna)*calcnernst;
      dnai = nernstc(vna,dnao,1);
      setval (lookup("dnai"),dnai,NR);  /* sodium internal concentration */
   }
 }
 else {
   dnai = nernstc(vna,dnao,1);
   setval (lookup("dnai"),dnai,NR);	/* sodium internal concentration */
 }
 setval (lookup("vna"),vna,NR);       /* sodium equilibrium potential    */


 if ((getval("vk")==LARGENUM) || calcnernst) {
     double ovk=vk;
   vk = nernstv(dko/dki);
   if (ovk==LARGENUM) ovk = vk;
   if ((calcnernst > 0) && (calcnernst < 1)) {
      vk = ovk + (vk-ovk)*calcnernst;
      dki = nernstc(vk,dko,1);
      setval (lookup("dki"),dki,NR);  /* sodium internal concentration */
   }
 }
 else {
   dki = nernstc(vk,dko,1);
   setval (lookup("dki"),dki,NR);	/* potassium external concentration */
 }
 setval (lookup("vk"),vk,NR);		/* potassium equilibrium potential */


 if ((getval("vcl")==LARGENUM) || calcnernst) {
     double ovcl=vcl;
   vcl = nernstv(dcli/dclo);
   if (ovcl==LARGENUM) ovcl = vcl;
   if ((calcnernst > 0) && (calcnernst < 1)) {
      vcl = ovcl + (vcl-ovcl)*calcnernst;
      dcli = nernstc(vcl,dclo,-1);	/* valence -1 */
      setval (lookup("dcli"),dcli,NR);  /* chloride internal concentration */
   }
 }
 else {
   dcli = nernstc(vcl,dclo,-1);		/* valence -1 */
   setval (lookup("dcli"),dcli,NR);	/* chloride internal concentration */
 }
 setval (lookup("vcl"),vcl,NR);	/* chloride equilibrium potential */

 if ((getval("vglu")==LARGENUM) || calcnernst) {
    vglu = 0;
 } 

}

/****************************************/

void calcvar(void)
{
     Symbol *plotsym;

  if (interp_only) return;

  ktemp = tempcel + KELVIN;
  ktA   = AF * pow(ktemp,-1.5);     /* gives .512 @ 25 deg C, for act() */
  frt   = FR / ktemp;
  f2rt  = F2R / ktemp;
  r2ft  = R2F * ktemp;

  setvrevs();
  setion();

  dpki  = dpkca * dki * CAZI;	     /* = perm K / perm Ca * [K]i / 4   */
				     /* see Hille 2nd ed, p 107 & p 344 */

  ktrb  = sqrt(4*Bk*ktemp/timinc); 		/* For Johnson noise */
		/* See Hille 1992 p 325 on Johnson noise */
  if (stiminc < timinc) {
        stiminc = timinc;
        setval (lookup("stiminc"),stiminc,NR);	
  }

	/* Set default mulitipliers for HH rate funcs. */

	/* The idea here is to retain the */
	/* original HH equations and consts */
	/* but to renormalize to 22 deg C */

	/* If user sets dqn, this overrides default */
	/* values for dqna and dqnb.  */

  if (dqn > 0) {
    dqna = dqn;
    dqnb = dqn;
  }
  if (implicit) { 
    setval (lookup("implfk"), 1.0,NR);	/* =1 -> implicit mode */
    implfk = getval ("implfk");
  }
  plotsym =  lookups ("plotlabel");
  if (plotsym!=NULL) plotlabel =  (const char*)plotsym->str;
}

/****************************************/

void varcopy(void)
{
     Symbol *plotsym;

  debug   = (int)getval ("debug");
  debugz  = (int)getval ("debugz");
  srseed   = (int)getval ("srseed");
  chrseed   = (int)getval ("chrseed");
  dkrseed   = (int)getval ("dkrseed");
  phrseed   = (int)getval ("phrseed");
  setrseed();
  nonodes  = (int)getval ("nonodes");
  info     = (int)getval ("info");
  if (silent) info = 0;
  if (interp_only) return;

  calcnernst = getval ("calcnernst");
  dnai    = getval ("dnai");
  dnao    = getval ("dnao");
  dko     = getval ("dko");
  dki     = getval ("dki");
  dpkna   = getval ("dpkna");
  dpcana  = getval ("dpcana");
  dpcaampa= getval ("dpcaampa");
  dpcacgmp= getval ("dpcacgmp");
  dpcanmda= getval ("dpcanmda");
  dpcasyn2= getval ("dpcasyn2");
  dpnak   = getval ("dpnak");
  dpcak   = getval ("dpcak");
  dpnaca  = getval ("dpnaca");
  dpkca   = getval ("dpkca");
  dcli    = getval ("dcli");
  dclo    = getval ("dclo");
  dcoo    = getval ("dcoo");
  dmgo    = getval ("dmgo");
  dcao    = getval ("dcao");
  dcai    = getval ("dcai");
  dtcai   = getval ("dtcai");
  dbasetdc= getval ("dbasetdc");
  dbasetc = getval ("dbasetc");
  dbasetca= getval ("dbasetca");
  dbasetsyn= getval ("dbasetsyn");
  tempcel = getval ("tempcel");

  betaIP3  = getval ("betaIP3");
  v1IP3    = getval ("v1IP3");
  kfCICR   = getval ("kfCICR");
  vm2CICR  = getval ("vm2CICR");
  vm3CICR  = getval ("vm3CICR");
  k1CICR   = getval ("k1CICR");
  k2CICR   = getval ("k2CICR");
  krCICR   = getval ("krCICR");
  kaCICR   = getval ("kaCICR");
  nCICR    = getval ("nCICR");
  mCICR    = getval ("mCICR");
  pCICR    = getval ("pCICR");
  c1CICR   = getval ("c1CICR");

  oIP3	    = getval ("oIP3");
  mtypeIP3  = getval ("mtypeIP3");
  a2IP3	    = getval ("a2IP3");
  a3IP3     = getval ("a3IP3");
  b2IP3     = getval ("b2IP3");
  d1IP3     = getval ("d1IP3");
  d2IP3     = getval ("d2IP3");
  d3IP3     = getval ("d3IP3");
  d4IP3     = getval ("d4IP3");
  k3IP3     = getval ("k3IP3");
  v2IP3     = getval ("v2IP3");
  v3IP3     = getval ("v3IP3");
  v4IP3     = getval ("v4IP3");

  deleccap= getval ("deleccap");
  dcm     = getval ("dcm");
  dnadens = getval ("dnadens");
  dkdens  = getval ("dkdens");
  dcadens = getval ("dcadens");
  dcavoff = getval ("dcavoff");
  dcaspvrev= getval ("dcaspvrev");
  ddca    = getval ("ddca");
  ddcao   = getval ("ddcao");
  ddiacao = getval ("ddiacao");
  dcabnd  = getval ("dcabnd");
  dcabt   = getval ("dcabt");
  dcabti  = getval ("dcabti");
  dcabf   = getval ("dcabf");
  dcabr   = getval ("dcabr");
  dcashell= (int)getval ("dcashell");
  dcaoshell= (int)getval ("dcaoshell");
  dcasarea= getval ("dcasarea");
  dcapkm  = getval ("dcapkm");
  dcavmax = getval ("dcavmax");
  dicafrac = getval ("dicafrac");
  dcakex  = getval ("dcakex");
  diexchfrac = getval ("diexchfrac");
  dmaxca  = getval ("dmaxca");
  dcaoffs  = getval ("dcaoffs");
  dcataum  = getval ("dcataum");
  dcatauh  = getval ("dcatauh");
  dkcataum = getval ("dkcataum");
  dkcatauh = getval ("dkcatauh");
  dnatauf  = getval ("dnatauf");
  dktauf   = getval ("dktauf");
  dcatauf  = getval ("dcatauf");
  dqm 	   = getval ("dqm");
  dqh 	   = getval ("dqh");
  dqn 	   = getval ("dqn");
  dqna 	   = getval ("dqna");
  dqnb 	   = getval ("dqnb");
  dqca 	   = getval ("dqca");
  dqd 	   = getval ("dqd");
  dqkca    = getval ("dqkca");
  dqsyn    = getval ("dqsyn");
  dqcavmax = getval ("dqcavmax");
  dqdc     = getval ("dqdc");
  dqcab    = getval ("dqcab");
  dqc 	   = getval ("dqc");
  dqri 	   = getval ("dqri");
  dqrm 	   = getval ("dqrm");
  dqrec    = getval ("dqrec");
  dqcrec   = getval ("dqcrec");
  dratehhm = getval ("dratehhm");
  dratehhh = getval ("dratehhh");
  dratehhna= getval ("dratehhna");
  dratehhnb= getval ("dratehhnb");
  dqeff    = getval ("dqeff");
  dsintinc = getval ("dsintinc");
  dsintres = getval ("dsintres");
  dri      = getval ("dri");
  drm      = getval ("drm");
  drs      = getval ("drs");
  drg      = getval ("drg");
  dsfa     = getval ("dsfa");
  dfta     = getval ("dfta");
  dftah    = getval ("dftah");
  dsfb     = getval ("dsfb");
  dftb     = getval ("dftb");
  dscu     = getval ("dscu");
  dsvn     = getval ("dsvn");
  dscavg   = getval ("dscavg");
  dscaeg   = getval ("dscaeg");
  dvg      = getval ("dvg");
  dvsz     = getval ("dvsz");
  dst      = getval ("dst");
  dsc      = getval ("dsc");
  dmaxsyn  = getval ("dmaxsyn");
  dskd     = getval ("dskd");
  dshc     = getval ("dshc");
  dstr     = getval ("dstr");
  dsmsgc   = getval ("dsmsgc");
  dsmrrp   = getval ("dsmrrp");
  dsrrg    = getval ("dsrrg");
  dsrrp    = getval ("dsrrp");
  dsms     = getval ("dsms");
  dcakd    = getval ("dcakd");
  dcahc    = getval ("dcahc");
  dckd     = getval ("dckd");
  dchc     = getval ("dchc");
  dsg      = getval ("dsg");
  dsn      = getval ("dsn");
  dmaxrod  = getval ("dmaxrod");
  dmaxcon  = getval ("dmaxcon");
  dmaxna   = getval ("dmaxna");
  dmaxampa = getval ("dmaxampa");
  dmaxnmda = getval ("dmaxnmda");
  dnmdamg  = getval ("dnmdamg");
  dnaoffsm = getval ("dnaoffsm");
  dnaoffsh = getval ("dnaoffsh");
  dnataum  = getval ("dnataum");
  dnatauh  = getval ("dnatauh");
  dktaum   = getval ("dktaum");
  dktauh   = getval ("dktauh");
  dsyntau  = getval ("dsyntau");
  dsyntauf = getval ("dsyntauf");
  dmaxk    = getval ("dmaxk");
  dmaxclca = getval ("dmaxclca");
  dclcavs  = getval ("dclcavs");
  dclcavsc = getval ("dclcavsc");
  dsclcac  = getval ("dsclcac");
  dnau     = getval ("dnau");
  dku      = getval ("dku");
  dkau     = getval ("dkau");
  dkcabu   = getval ("dkcabu");
  dkcasu   = getval ("dkcasu");
  dclcasu  = getval ("dclcasu");
  dkihu    = getval ("dkihu");
  dhcnu    = getval ("dhcnu");
  dcalu    = getval ("dcalu");
  dcatu    = getval ("dcatu");
  dcgmpu   = getval ("dcgmpu");
  dampau   = getval ("dampau");
  dgabau   = getval ("dgabau");
  dglyu    = getval ("dglyu");
  dpru     = getval ("dpru");
  dkoffsn  = getval ("dkoffsn");
  dkoffsh  = getval ("dkoffsh");
  dgjoff   = getval ("dgjoff");
  dgjtau   = getval ("dgjtau");
  dgjatpdec= getval ("dgjatpdec");
  dgjnv    = getval ("dgjnv");
  dbd1 	   = getval ("dbd1");
  dbd2 	   = getval ("dbd2");
  dbk1 	   = getval ("dbk1");
  dbk2 	   = getval ("dbk2");
  dsd1 	   = getval ("dsd1");
  dsd2 	   = getval ("dsd2");
  dsk1 	   = getval ("dsk1");
  dsk2 	   = getval ("dsk2");
  complam  = getval ("complam");
  lamcrit  = getval ("lamcrit");
  relincr  = getval ("relincr");
  timinc   = getval ("timinc");
  stiminc  = getval ("stiminc");
  simtime  = getval ("time");
  endexp   = getval ("endexp");
  ploti    = getval ("ploti");
  plmax    = getval ("plmax");
  plmin    = getval ("plmin");
  setxmax  = getval ("setxmax");
  setxmin  = getval ("setxmin");
  setymax  = getval ("setymax");
  setymin  = getval ("setymin");
  crit     = getval ("crit");
  relax    = getval ("relax");
  djnoise  = getval ("djnoise");
  stimelem = (int)getval ("stimelem");
  vcolor   = (int)getval ("vcolor");
  lcolor   = (int)getval ("lcolor");
  rcolor   = (int)getval ("rcolor");
  cacolor   = (int)getval ("cacolor");
  nacolor   = (int)getval ("nacolor");
  sgcolor   = (int)getval ("sgcolor");
  srcolor   = (int)getval ("srcolor");
  nozerochan= (int)getval ("nozerochan");
  oldsynapz = (int)getval ("oldsynapz");
  use_ghki = (int)getval ("use_ghki");
  fread_expr = (int)getval ("fread_expr");
  euler    = (int)getval ("euler");
  implicit = (int)getval ("implicit");
  implfk   = getval ("implfk");
  synaptau = getval ("synaptau");
  numstimchan  = getval ("numstimchan");
  unmaskthr  = getval ("unmaskthr");
  stimonh  = getval ("stimonh");
  stimonl  = getval ("stimonl");
  srtimestep  = getval ("srtimestep");
   
  vidmode  = (int)getval ("vidmode");
  plsep    = (int)getval ("plsep");
  allowrev = (int)getval ("allowrev");
  dashfl   = (int)getval ("dash");
  disp     = (int)getval ("disp");
  prmap    = (int)getval ("prmap");
  stim_elem = (int)getval ("stim_elem");
  calcvar();
}

/****************************************/

void varcopyu(void)

/* Set params from interpreter that might affect user-called functions. 
   Put into functions like "hinf", etc, which may be used without
   running a "step" or "run" which would initialize default
   variables. This is supposed to run faster than "varcopy" above. 
*/

{
  dcoo     = getval ("dcoo");
  dmgo     = getval ("dmgo");
  dcao     = getval ("dcao");
  dcavoff  = getval ("dcavoff");
  dnaoffsm = getval ("dnaoffsm");
  dnaoffsh = getval ("dnaoffsh");
  dkoffsn  = getval ("dkoffsn");
  dkoffsh  = getval ("dkoffsh");
  dcaoffs  = getval ("dcaoffs");
  timinc   = getval ("timinc");
  stiminc  = getval ("stiminc");
  dbasetc  = getval ("dbasetc");
  dbasetca = getval ("dbasetca");
  dbasetsyn= getval ("dbasetsyn");
  dratehhm = getval ("dratehhm");
  dratehhh = getval ("dratehhh");
  dratehhna= getval ("dratehhna");
  dratehhnb= getval ("dratehhnb");
  dqeff    = getval ("dqeff");
  ploti    = getval ("ploti");
  tempcel  = getval ("tempcel");
  dqm  = getval ("dqm");
  dqh  = getval ("dqh");
  dqn  = getval ("dqn");
  dqna  = getval ("dqna");
  dqnb  = getval ("dqnb");
  disp     = (int)getval ("disp");
  prmap    = (int)getval ("prmap");
  fread_expr = (int)getval ("fread_expr");
  if (stiminc < timinc) {
        stiminc = timinc;
        setval (lookup("stiminc"),stiminc,NR);	
  }
  implicit = (int) getval ("implicit");
  implfk = getval ("implfk");
  calcvar();
}

/****************************************/

void dinit()
{

#define NR NUMBER

  int len;

 Symbol *lookup(const char *str);

 setval (lookup("calcnernst"),0,NR);	/* =1 -> define vrevs from ion conc */
					/* =0 -> define ion conc from vrevs */
 setval (lookup("vna"),LARGENUM,NR);	/* sodium equilibrium potential    */
 setval (lookup("vk"),LARGENUM,NR);	/* potassium equilibrium potential */
 setval (lookup("vcl"),LARGENUM,NR);	/* chloride equilibrium potential  */
 setval (lookup("vglu"),LARGENUM,NR);	/* glutamate chan (nonselective) equil pot  */

 tempcel = 37.0;
 setval (lookup("tempcel"),tempcel,NR);	/* temperature for chan rate chtyp */
 ktemp = tempcel + KELVIN;
 ktA   = AF * pow(ktemp,-1.5);		/* gives .512 @ 25 deg C, for act() */

 dnao = 143e-3;				/* Ames medium, Ames & Nesbett 1981 */
 dko = 5.0e-3;				/* Ames medium = 3.6 mM */
					/* Set dko a little higher  */
 					/*  allows K vrev to be -.082 w/dpnak */
 dclo = 125.4e-3; 			/* Ames medium */

 dnai = nernstc(0.065,dnao,1);
 setval (lookup("dnao"),dnao,NR);	/* external [Na], used for Na vrev */
 setval (lookup("dnai"),dnai,NR);	/* internal [Na], used for Na vrev */
 setval (lookup("dpkna"),0.08,NR);	/* perm of K rel to Na, for Na vrev */
 setval (lookup("dpcana"),0.01,NR);	/* perm of Ca rel to Na, for Na vrev */
					/*  Baker et al, (1971) = 0.01       */
					/*  Meves and Vogel, (1973) = 0.1    */
 setval (lookup("dpcaampa"),0.1,NR);	/* rel perm of Ca for ampa synapse */
 setval (lookup("dpcacgmp"),0.15,NR);	/* rel perm of Ca for cgmp synapse */
 setval (lookup("dpcanmda"),0.2,NR);	/* rel perm of Ca for nmda synapse */
 setval (lookup("dpcasyn2"),0.1,NR);	/* rel perm of Ca for syn2 synapse */

 dki = nernstc(-0.0892,dko,1);		/* allows K vrev to be -.082 w/dpnak */
 setval (lookup("dko"),dko,NR);		/* external [K], used for K vrev */
 setval (lookup("dki"),dki,NR);		/* internal [K], used for K vrev */
 setval (lookup("dpnak"),0.01,NR);	/* perm of Na rel to K, for K vrev */
 setval (lookup("dpcak"),0.0001,NR);	/* perm of Ca rel to K, for K vrev */

					/* see Hille (1992) p 351-353 */
 setval (lookup("dpnaca"),1e-3,NR);	/* perm of Na rel to Ca, for Ca vrev */
 setval (lookup("dpkca"),3e-4,NR);	/* perm of K rel to Ca, for Ca vrev */
					/*  see Hess et al, (1986) */
					/*  and Tsien et al. (1987) */

 dcli = nernstc(-0.07,dclo,-1);
 setval (lookup("dclo"),dclo,NR);	/* external [Cl], used for Cl vrev */
 setval (lookup("dcli"),dcli,NR);	/* internal [Cl], used for Cl vrev */
 setion();			/* set outside ions and valences for a start */

 setval (lookup("deleccap"),1e-12,NR);	/* f electrode capacitance */
 setval (lookup("dcm"),1e-6,NR);	/* f/cm-cm     membrane capacitance */
 setval (lookup("dri"),200.,NR);	/* ohm-cm    axial resistance, 70-150 Stuart & Spruston, 1998 */
 					/* 140 ohm-cm @ 36 degC, Fohlmeister et al., 2010 */
 setval (lookup("drm"),40000.,NR);	/* ohm-cm-cm    membrane resistance */
 setval (lookup("drs"),10e6,NR);	/* ohm          electrode resistance */
 setval (lookup("drg"),5e6,NR);		/* ohm-um-um    gap junc resistance */
 setval (lookup("dfta"),0.2,NR);	/* synaptic time const (msec) */
 setval (lookup("dftah"),1,NR);		/* synaptic high-pass time c (msec) */
 setval (lookup("dsfa"),2.0,NR);	/* synaptic number of filters */
 setval (lookup("dftb"),.2,NR);		/* synaptic time const (msec) */
 setval (lookup("dsfb"),0.0,NR);	/* synaptic number of filters */
 setval (lookup("dscu"),25e-12,NR);	/* unitary conductance of chans */
 setval (lookup("dsvn"),5.0,NR);	/* number of syn vesicle rel sites */
 setval (lookup("dscavg"),1e6,NR);	/* synaptic ca vesicle release gain */
 setval (lookup("dscaeg"),1,NR);	/* synaptic ca power for release gain */
 setval (lookup("dvg"),1.0,NR);		/* vesicle release gain */
 setval (lookup("dvsz"),1.0,NR);	/* size of vesicles released */
 setval (lookup("dst"), -.04,NR);	/* synaptic threshold (volts) */
 setval (lookup("dsc"),(double)EXPON,NR);/* synaptic curve */
 setval (lookup("dmaxsyn"),200e-12,NR);	/* maximum cond for synapse (S) */
 setval (lookup("dskd"),1.0,NR);	/* half-max neurotr saturation (kd)*/
 setval (lookup("dshc"),1.0,NR);	/* neurotrans binding Hill coeff */
 setval (lookup("dstr"),100e-6,NR);	/* max neurotrans conc for synapse */
 setval (lookup("dsmsgc"),10e-6,NR);	/* max 2nd mesg conc for synapse */
 setval (lookup("dsmrrp"),500.0,NR);	/* max readily releasible pool */
 setval (lookup("dsrrg"),0.0,NR);	/* readily releasible pool gain mult */
 setval (lookup("dsrrp"),500.0,NR);	/* readily releasible pool */
 setval (lookup("dsms"),0.0,NR);	/* max sustained release rate (/sec) */
 setval (lookup("dcakd"),1e10,NR);	/* half-max Ca sat (kd) at chan */
 setval (lookup("dcahc"),1,NR);		/* hill coeff for Ca binding at chan */
 setval (lookup("dckd"),0.1,NR);	/* half-max cG saturation (kd) */
 setval (lookup("dchc"),1.0,NR);	/* hill coeff for cG binding (kd) */
 setval (lookup("dsn"),2.0,NR);		/* synaptic gain, lin or exp */
 setval (lookup("dsg"),1.0,NR);		/* synaptic cG gain */
 setval (lookup("dmaxrod"),44e-12,NR);/* dark cond of rod gated channel (S) */
 setval (lookup("dmaxcon"),1.4e-9,NR);	/* dark cond of cone gated chan (S)*/
 setval (lookup("dnadens"),0.25,NR);	/* density of Na V-gated chan(S/cm2)*/
 setval (lookup("dkdens"),0.0707,NR);	/* density of K V-gated chan (S/cm2)*/
 setval (lookup("dcadens"),0.005,NR);	/* density of Ca V-gated chan(S/cm2)*/
 setval (lookup("dcavoff"),.018,NR);	/* Q10 for voffset due to cao (V) (Hille) */
 setval (lookup("dcaspvrev"),.00324,NR);/* Ca surface potential (mult dcavoff) */
 setval (lookup("dcaphg"),0.45,NR);   /* Ca pH conductance (Barnes & Bui, 1991) */
 					/* Add to vrev.  (Hille)*/

					   /* Dan Emerson's CICR parameters */
					   /*   mostly taken from  Goldbeter, 1990 */

 setval (lookup("kfCICR"),1e-6,NR);	   /* passive leak rate const from ryanodine store, 1/s */
 setval (lookup("vm2CICR"),65e-6,NR);      /* max rate of Ca pump into ryanodine store, 65e-6/ms */
 setval (lookup("vm3CICR"),500e-6,NR);     /* max rate of Ca pump from ryanodine store,500e-6/ms */
 setval (lookup("k1CICR"),5.0e-6,NR);      /*assoc thresh pump const, ryan store uptake,5e-6 M */
 setval (lookup("k2CICR"),1e-6,NR);        /* assoc thresh pump const, ryan store uptake,1e-6 M */
 setval (lookup("krCICR"),2e-6,NR);        /* assoc thresh pump const, ryan store release,2e-6 M */
 setval (lookup("kaCICR"),0.9e-6,NR);      /* assoc thresh pump const, ryan store release,0.9e-6 M */
 setval (lookup("nCICR"),1,NR);   	   /* Hill coeff, coop pump binding ryan store uptake, 1 */ 
 setval (lookup("mCICR"),1,NR);		   /* Hill coeff, coop pump binding ryan store release, 1 */ 
 setval (lookup("pCICR"),4,NR);	           /* Hill coeff, coop pump binding ryan store release, 4 */ 

 setval (lookup ("c1CICR"),0.185,NR);      /* ratio of ER vol/Cytoplasmic vol, 0.185 (Rinzel and Li 1994) */ 
 
 setval (lookup("betaIP3"),0.31,NR);       /* init fractional rate of constant IP3 store release, default 0.31 */
 setval (lookup("v1IP3"),7.3e-6, NR);      /* init constant IP3 store flux to 7.3e-6/ms (Goldbeter 1990) */


					  /* constants for IP3 from  (Rinzel and * Li, 1994) */

 setval (lookup ("oIP3"),2,NR);          /* Hill coeff, coop pump binding ip3 store uptake, 2 */
 setval (lookup ("mtypeIP3"),0,NR);      /* 0 => static (def = use equilibrated inf m value) */
 setval (lookup ("a2IP3"),4.2e8,NR);     /* IP3 store forward const for Tau M value, 4.2e8/ms */
 setval (lookup ("a3IP3"),2e5,NR);       /* IP3 store constant for IP3 Tau H, 2e5 1/MS */
 setval (lookup ("b2IP3"),4.1,NR);       /* IP3 store backward const for Tau M value, 4.1/ms */
 setval (lookup ("d1IP3"),0.13e-6,NR);   /* IP3 store const for Q value for h gate, 0.13e-6 M */
 setval (lookup ("d2IP3"),1.049e-6,NR);  /* IP3 store const for Q value for h gate, 1.049e-6 M  */
 setval (lookup ("d3IP3"),0.9434e-6,NR); /* IP3 store const for Q value for h gate, 0.9434e-6 M */
 setval (lookup ("d4IP3"),0.823e-7,NR);  /* IP3 store const for m gate at infinity, 0.823e-7 M */
 setval (lookup ("k3IP3"),0.1e-6,NR);    /* IP3 store kd for Ca uptake, 0.1e-6 M) */
 setval (lookup ("v4IP3"),6,NR);         /* IP3 store rate for Ca release, 6 M S-1 */
 setval (lookup ("v2IP3"),0.11,NR);      /* IP3 store rate for Ca release, 0.11 S-1 */
 setval (lookup ("v3IP3"),0.9e-6,NR);    /* IP3 store constant for Ca uptake, 0.9e-6 M S-1 */


 setval (lookup("ddca"),DCa,NR);	/* Diffusion const for Ca = 2e-6 */
                                        /*  slower in cytoplasm than in H2O*/
					/* see de Schutter&Smolen "Methods" */
 setval (lookup("ddcao"),DCa,NR);	/* Diffusion const for ext Ca = 2e-6 */
                                        /*  slower in cytoplasm than in H2O*/
					/* see de Schutter&Smolen "Methods" */
					/*   (from Albritton et al, 1992) */
 setval (lookup("dcoo"),0.0,NR);	/* external Co++ (M) (for cavoff) */
 setval (lookup("dmgo"),DMGO,NR);	/* external Mg++ (M) (for cavoff) */
 setval (lookup("dcao"),DCAO,NR);	/* 1.15 mM = external Ca++ (M) */
 setval (lookup("dcai"),50e-9,NR);	/* conc. of internal Ca++ (M) */
 setval (lookup("dtcai"),10e-9,NR);	/* threshold internal Ca++ for pump (M)*/
 setval (lookup("dcabnd"),10.0,NR);	/* ratio of bound to free calcium */
 setval (lookup("dcabti"),30e-6,NR);	/* Ca buffer total in first shell */
 setval (lookup("dcabt"),3e-6,NR);	/* Ca buffer total in shells */
					/*  Yamada, Koch, Adams 1998) */
					/*  Simon, Llinas 1985 */
 setval (lookup("dcabf"),1e8,NR);	/* Ca buffer forw rate to Ca.B /M/sec */
 setval (lookup("dcabr"),100,NR);	/* Ca buffer reverse rate to Ca.B /sec */
 setval (lookup("dcashell"),10.0,NR);	/* Number of int Ca diffusion shells */
 setval (lookup("dcaoshell"),0.0,NR);	/* Number of ext Ca diffusion shells */
 setval (lookup("dcasarea"),0,NR);	/* Area of int shell for Ca flux, um2 */
 setval (lookup("ddiacao"),1e5,NR);	/* Diameter of ext Ca core */
 setval (lookup("dcapkm"),1e-6,NR);	/* 1/2 sat conc. of Ca++ pump (M) */
 setval (lookup("dcavmax"),2e-7,NR);	/* Vmax for Ca pump (ma/cm2)       */
 setval (lookup("dicafrac"),0,NR);	/* Fraction of Icapump added to comp Itot */
 setval (lookup("dcakex"),5e-11,NR);	/* Rate for Na/Ca exchanger(A/mM4/cm2)*/
 setval (lookup("diexcafrac"),0,NR);	/* Fraction of Icaexch added to comp Itot */
 setval (lookup("dmaxca"),1e-10,NR);	/* cond of Ca voltage gated channel */
 setval (lookup("dcaoffs"),0.0,NR);	/* activation offset for Ca chan  */
 setval (lookup("dcataum"),1.0,NR);	/* activation tau for Ca chan  */
 setval (lookup("dcatauh"),1.0,NR);	/* inactivation tau for Ca chan  */
 setval (lookup("dkcataum"),1.0,NR);	/* activation tau for Kca chan  */
 setval (lookup("dkcatauh"),1.0,NR);	/* inactivation tau for Kca chan  */
 setval (lookup("dnatauf"),1.0,NR);	/* tau for Na noise in 2-state modl */
 setval (lookup("dktauf"),1.0,NR);	/* tau for K noise in 2-state model */
 setval (lookup("dcatauf"),1.0,NR);	/* tau for Ca noise in 2-state model */

 setval (lookup("dqm"),2.3,NR);		/* Q10 for am, bm, (1.95@22-37deg, Fohlmeister et al. 2010) */
 setval (lookup("dqh"),2.3,NR);		/* Q10 for ah, bh, 1.95@22-37deg, 3.7@10-22deg "" */
 setval (lookup("dqn"),2.3,NR);		/* Q10 for an, bn (use dqna/b if 0) (1.9@22-37deg, 3.7@10-22deg) */
 setval (lookup("dqna"),2.3,NR);	/* Q10 for an, formerly 3.2 (1.9@22-37deg) */
 setval (lookup("dqnb"),2.3,NR);	/* Q10 for bn, formerly 2.8  (1.9@22-37deg)*/
 setval (lookup("dqca"),3.0,NR);	/* Q10 for ac, bc (Ca chan) */
 setval (lookup("dqd"),3.0,NR);		/* Q10 for ad, bd (KA h) */
 setval (lookup("dqkca"),3.0,NR);	/* Q10 for kca a, b */
 setval (lookup("dqsyn"),3.0,NR);	/* Q10 for synaptic chans (ampa,etc)*/
 setval (lookup("dqcavmax"),2.0,NR);	/* Q10 for Ca pump */
 setval (lookup("dqcab"),2.0,NR);	/* Q10 for Ca buffer */
 setval (lookup("dqcab"),2.0,NR);	/* Q10 for Ca buffer */
 setval (lookup("dqdc"),QDC,NR);	/* Q10 for diffusion constant, 1.3 */
 setval (lookup("dqc"),1.44,NR);	/* Q10 for channel conductance */
					/* 1.44, Rodriguez BM, et al. (1998) */
					/* J. Gen. Physiol.  112:223   */
 					/* 1.85@22-37 deg, Fohlmeister et * al., 2010) */
 					/* 3.0 @10-22 deg, Fohlmeister et * al., 2010) */
 setval (lookup("dqrm"),1.0,NR);	/* Q10 for Rm conductance, like dqc for channels */
 setval (lookup("dqri"),1.0,NR);	/* Q10 for Ri (axial conductance) (0.8, Fohlmeister et al, 2010) */
 setval (lookup("dqrec"),2.7,NR);	/* Q10 for photoreceptor conductance shape (not used) */
 setval (lookup("dqcrec"),1.8,NR);	/* Q10 for photoreceptor conductance (not used) */
					/* Baylor, Matthews & Yau, 1983) */

		/* the intention here is to preserve the */
		/* original HH rate functions, which were */
		/* specified for 6.3 deg C, but multiply them to */
 		/* normalize to the more common 22 deg C. */
		/* "dratehhm" multiplies alpham & betam, etc. */

		/* The variable "dbasetc" sets the temperature at */
		/* which voltage-gated channels are normalized. */

		/* Not used for other rate funcs which were originally */
		/* defined for 22 deg C */

 setval (lookup("dbasetc"),BASETC,NR);  /* base T for chan kinetics, 22 deg C */
 setval (lookup("dbasetca"),BASETCA,NR);/* base T for Ca pump kinet, 22 deg C */
 setval (lookup("dbasetsyn"),BASETSYN,NR);/* base T for synap kinet, 22 deg C */
 setval (lookup("dbasetdc"),BASETDC,NR); /* base T for Ca diffusion, 22 deg C */

 dbasetc = getval("dbasetc");

  /* Q10 values for 6.3 to 22 deg taken from Frankenhaeuser & Moore (1963) */

 dratehhm  = exp(log(3.0) * (dbasetc-BASETHH) / 10.0);
 dratehhh  = exp(log(3.0) * (dbasetc-BASETHH) / 10.0);

 dratehhna = exp(log(3.0) * (dbasetc-BASETHH) / 10.0); 	/* formerly 3.2 */
 dratehhnb = exp(log(3.0) * (dbasetc-BASETHH) / 10.0);  /* formerly 2.8 */

 setval (lookup("dratehhm"),dratehhm,NR); /* mult for HH m rate, 22 deg C */
 setval (lookup("dratehhh"),dratehhh,NR); /* mult for HH h rate, 22 deg C */
 setval (lookup("dratehhna"),dratehhna,NR); /* mult for HH alphan, 22 deg C */
 setval (lookup("dratehhnb"),dratehhnb,NR); /* mult for HH betan, 22 deg C */

 setval (lookup("dqeff"),0.67,NR);	/* quantum efficiency, R* per photon */
 setval (lookup("dsintinc"),0.001,NR);	/* min time incr for sine waves, sec */
 setval (lookup("dsintres"),0.002,NR);	/* time res for sine waves, per cyc */

 setval (lookup("dnmdamg"),2e-3,NR);	/* Mg++ extracellular conc */
 setval (lookup("dmaxna"),1e-9,NR);	/* cond of Na voltage gated channel(S)*/
 setval (lookup("dnaoffsm"),0.0,NR);	/* m offset voltage for Na chan */
 setval (lookup("dnaoffsh"),0.0,NR);	/* h offset voltage for Na chan */
 setval (lookup("dnataum"),1.0,NR);	/* activation tau for Na chan */
 setval (lookup("dnatauh"),1.0,NR);	/* inactivation tau for Na chan */
 setval (lookup("dgjoff"),GJOFF,NR);	/* gap junction voltage offset */
 setval (lookup("dgjtau"),1.0,NR);	/* activation tau for gap junction */
 setval (lookup("dgjatpdec"),0.999,NR);	/* decay per timestep for atp in synaptic cleft */
 setval (lookup("dgjnv"),0.2,NR);	/* non-volt-sens fraction of gj cond */
 setval (lookup("dbd1"),1.0,NR);	/* voltage multiplier for BK chan alph*/
 setval (lookup("dbd2"),1.0,NR);	/* voltage multiplier for BK chan bet */
 setval (lookup("dbk1"),1e-6,NR);	/* ca multiplier for BK chan alph */
 setval (lookup("dbk2"),1e-6,NR);	/* ca multiplier for BK chan bet  */
 setval (lookup("dsd1"),0.0,NR);	/* voltage multiplier for SK chan alph*/
 setval (lookup("dsd2"),0.0,NR);	/* voltage multiplier for SK chan bet */
 setval (lookup("dsk1"),1e-7,NR);	/* ca multiplier for SK chan alph */
 setval (lookup("dsk2"),1e-7,NR);	/* ca multiplier for SK chan bet  */
 setval (lookup("dmaxk"),2e-10,NR);	/* cond of K voltage gated channel(S)*/
 setval (lookup("dmaxclca"),2e-10,NR);	/* cond of Ca-gated Cl channel(S)*/
 setval (lookup("dclcavs"),2,NR);	/* ClCa chan voltage sens over 120 mv */
 setval (lookup("dclcavsc"),2,NR);	/* ClCa chan voltage sens over 120 mv, core */
 setval (lookup("dsclcac"),10,NR);	/* ClCa chan core Ca sens */
 setval (lookup("dnau"),22e-12,NR);	/* unitary cond of Na chan(S)=32pS@33 */
 setval (lookup("dku"), 11.5e-12,NR);	/* unitary cond of Kdr chan(S)=15pS@30 */
 setval (lookup("dkau"),22e-12,NR);	/* unitary cond of KA chan (S)=30ps@30 */
 setval (lookup("dkcasu"),14.2e-12,NR);	/* unitary cond of Kca SK ch=22 pS@35*/
 setval (lookup("dkcabu"),74.3e-12,NR);	/* unitary cond of Kca BK ch=115 pS@35*/
 setval (lookup("dclcasu"),14.2e-12,NR);/* unitary cond of Kca SK ch=22 pS@35*/
 setval (lookup("dkihu"),33e-12,NR);	/* unitary cond of Ih chan (S) */
 setval (lookup("dhcnu"),33e-12,NR);	/* unitary cond of Ih chan (S) */
 setval (lookup("dcalu"),8e-12,NR);	/* unit cond of Ca T chan=20 pS @35 */
 setval (lookup("dcatu"),3e-12,NR);	/* unit cond of Ca L chan=8 pS @35 */
 setval (lookup("dampau"),25e-12,NR);	/* unitary cond of AMPA chan (S) */
 setval (lookup("dgabau"),30e-12,NR);	/* unitary cond of GABA chan (S) */
 setval (lookup("dglyu"),30e-12,NR);	/* unitary cond of GLY chan (S) */
 setval (lookup("dpru"),2e-12,NR);	/* unitary cond of photoreceptor channel (S) */
 setval (lookup("dcgmpu"),25e-12,NR);	/* unitary cond of cGMP chan (S) */
 setval (lookup("dkoffsn"),0.0,NR);	/* offset voltage for K chan */
 setval (lookup("dkoffsh"),0.0,NR);	/* offset voltage for K chan */
 setval (lookup("dktaum"),1.0,NR);	/* activation tau for K chan */
 setval (lookup("dktauh"),1.0,NR);	/* inactivation tau for K chan */
 setval (lookup("dsyntau"),1.0,NR);	/* time const for synaptic chans */
 setval (lookup("dsyntauf"),1.0,NR);	/* time const for synaptic chan noise */
 setval (lookup("complam"), .1,NR);	/* fraction of lambda for compartment */
 setval (lookup("lamcrit"), .3,NR);	/* criterion for condensing comps */
 setval (lookup("timinc"), 1e-4,NR);	/* time incr for model (sec) */
 setval (lookup("stiminc"), 1e-4,NR);	/* synaptic/stim time incr (sec) */
 setval (lookup("ncomps"), 0.0,NR);	/* number of compartments */
 setval (lookup("endexp"), 50e-3,NR);	/* end of recording for model (sec) */
 setval (lookup("ploti"), 1e-3,NR);	/* plot increment (sec) */
 setval (lookup("plmax"), .04,NR);	/* max for voltage plots (volts) */
 setval (lookup("plmin"), -.08,NR);	/* min for voltage plots (volts) */
 setval (lookup("setxmax"),LARGENUM,NR);	/* override for x axis */
 setval (lookup("setxmin"),LARGENUM,NR);	/* override for x axis */
 setval (lookup("setymax"),LARGENUM,NR);	/* override for y axis */
 setval (lookup("setymin"),LARGENUM,NR);	/* override for y axis */
 setval (lookup("crit"), 1e-10,NR);	/* criterion for steps (volts) */
 setval (lookup("relax"), .1,NR);	/* over-relaxation constant */
 setval (lookup("relincr"), .0001,NR);	/* relax increment */
 setval (lookup("synaptau"), 1.0,NR);	/* synaptic time constant multiplier */
 setval (lookup("implfk"), 0.5,NR);	/* =0.5 -> C-N, =1 -> full implicit mode */
 setval (lookup("euler"), 0.0,NR);	/* =1 -> euler mode */
 setval (lookup("implicit"), 0.0,NR);	/* =1 -> implicit mode */
 setval (lookup("stimelem"), 0.0,NR);	/* "stim" understands neural elements */
 setval (lookup("vcolor"), VCOLOR,NR);	/* Use color from voltage */
 setval (lookup("lcolor"), LCOLOR,NR);	/* Use color from light flux */
 setval (lookup("rcolor"), RCOLOR,NR);	/* Use color from cell region */
 setval (lookup("cacolor"),CACOLOR,NR);	/* Use color from [Ca]i */
 setval (lookup("nacolor"),NAICOLOR,NR);/* Use color from Na chan */
 setval (lookup("sgcolor"),SGCOLOR,NR); /* Use color from syn cond */
 setval (lookup("srcolor"),SRCOLOR,NR); /* Use color from syn rate */
 setval (lookup("unmaskthr"),-1e20,NR); /* threshold for unmasking light */
 setval (lookup("numstimchan"),1,NR);   /* number of indep stimulus channels */
 setval (lookup("stimonh"),1.0,NR);	/* high thresh for transducer vclamp */
 setval (lookup("stimonl"),-1.0,NR);	/* low thresh for transducer vclamp */
 setval (lookup("nozerochan"),1.0,NR);	/* Don't make channels with maxcond=0*/
 setval (lookup("oldsynapz"),0.0,NR);	/* Revert to old synapse (getcurv after filt1) */
 setval (lookup("use_ghki"),0.0,NR);	/* Use GHK current eqn for channels */
 setval (lookup("fread_expr"),0.0,NR);	/* fread file expressions (no spaces) */
 setval (lookup("debug"), 0.0,NR);	/* turn on debugging features */
 setval (lookup("debugz"), 0.0,NR);	/* turn on debugging features */
 setval (lookup("dashfl"), 0.0,NR);	/* set dashes on plot */
 setval (lookup("allowrev"), 0.0,NR);	/* allow reverse point plotting */
 setval (lookup("djnoise"), 0.0,NR);	/* > 0 -> Johnson noise for all Rm's */
 setval (lookup("plsep"), 0.0,NR);	/* set multiple separate plotting */
 setval (lookup("disp"), 0.0,NR);	/* display control variable = "-d n" */
 setval (lookup("time"), 0.0,NR);	/* time at beginning */
 setval (lookup("srtimestep"), 1e-4,NR);/* time step for reading stimuli */
  len = strlen(progname);
  stdplt = stdout;
  if (len) {
#define VLINSIZ 40
       char runvidline[VLINSIZ];

    if (strcmp(progname,"nd")==0) disp |= DISP;
    if (strcmp(progname,"ndv")==0) { 
          FILE *ftemp;
    	disp |= DISP; 
	vidmode=1; 
	if (vidsiz <= 0) vidsiz = 1;
	sprintf (runvidline,"vid -w %g\n",vidsiz);
	if ((ftemp=popen(runvidline,"w"))==NULL) {
	    ncfprintf (stderr,"Error: nc can't open stdout\n");
	}
	else stdplt = ftemp;
	pipefl = 1;
    }
    if (strcmp(progname,"ncv")==0) {
          FILE *ftemp;
	vidmode=1;
	if (vidsiz <= 0) vidsiz = 1;
	sprintf (runvidline,"vid -w %g\n",vidsiz);
	if ((ftemp=popen(runvidline,"w"))==NULL) {
	    ncfprintf (stderr,"Error: nc can't open stdout\n");
	}
	else stdplt = ftemp;
	pipefl = 1;
    }
  }
 timeptr = lookup("time");		/* get pointer to time */
 plotiptr = lookup("ploti");		/* get pointer to ploti */
 ncompptr = lookup("ncomps");		/* get pointer to ncomps */

}

/*******************************/

void iinit()

/* initialize interpreter things */

{
    int len;

  setval (lookup("version"),   &ncversion,NR);/* nc version number */
  setval (lookup("ncversion"), &ncversion,NR);/* nc version number */
  setval (lookup("makestim"), &makestim,NR);  /* =1 -> running stim */
  setval (lookup("infile"),   &infile,STRING); /* input file */
  setval (lookup("progname"), &progname,STRING);/* nc commnd name*/
  len = strlen(progname);
  if (setvid) vidmode=1;
  else if (unsetvid) vidmode=0;

  setval (lookup("silent"),&silent,NR);		/* no printouts */
  info = 2;
  setval (lookup("info"),&info,NR);		/* level of info in printouts */
  if (silent) {
     info = 0;
  }

  rseed   = 123456.0; 	/* default random seed */
  srseed  = 0.0; 	/* stimulus random seed for physiology */
  chrseed = 5373142.0;	/* photoreceptor chan noise random seed */
  dkrseed = 2171828.0;	/* photoreceptor dark noise random seed */
  phrseed = 3141592.0;	/* photoreceptor photon noise random seed */
  
  setval (lookup("rseed"), &rseed,NR);		/* default random seed */
  setval (lookup("srseed"), &srseed,NR);   	/* stimulus random seed for physiology */
  setval (lookup("chrseed"), &chrseed,NR);	/* photoreceptor chan noise random seed */
  setval (lookup("dkrseed"), &dkrseed,NR);	/* photoreceptor dark noise random seed */
  setval (lookup("phrseed"), &phrseed,NR);	/* photoreceptor photon noise random seed */

  setval (lookup("stdin"), (double)(long int)stdin,NR);   /* standard input */
  setval (lookup("stdout"),(double)(long int)stdout,NR);  /* standard output */
  setval (lookup("stderr"),(double)(long int)stderr,NR);  /* standard err out */
}

#undef NR

/*******************************/

void vinit()
{

#define NR NUMBER

  int len;

 Symbol *lookup(const char *str);

 setval (lookup("calcnernst"),&calcnernst,NR);	/* =1 -> define vrevs from ion conc */
					/* =0 -> define ion conc from vrevs */
 setval (lookup("vna"),&vna,NR);	/* sodium equilibrium potential    */
 setval (lookup("vk"),&vk,NR);		/* potassium equilibrium potential */
 setval (lookup("vcl"),&vcl,NR);	/* chloride equilibrium potential  */
 setval (lookup("vglu"),&vglu,NR);	 /* glutamate (nonselective) equil pot  */

 setval (lookup("tempcel"),&tempcel,NR);	/* temperature for chan rate chtyp */
 setval (lookup("dnao"),&dnao,NR);	/* external [Na], used for Na vrev */
 setval (lookup("dnai"),&dnai,NR);	/* internal [Na], used for Na vrev */
 setval (lookup("dpkna"),&dpkna,NR);	/* perm of K rel to Na, for Na vrev */
 setval (lookup("dpcana"),&dpcana,NR);	/* perm of Ca rel to Na, for Na vrev */
					/*  Baker et al, (1971) = 0.01       */
					/*  Meves and Vogel, (1973) = 0.1    */
 setval (lookup("dpcaampa"),&dpcaampa,NR);/* rel perm of Ca for ampa synapse */
 setval (lookup("dpcacgmp"),&dpcacgmp,NR);/* rel perm of Ca for cgmp synapse */
 setval (lookup("dpcanmda"),&dpcanmda,NR);/* rel perm of Ca for nmda synapse */
 setval (lookup("dpcasyn2"),&dpcasyn2,NR);/* rel perm of Ca for syn2 synapse */

 setval (lookup("dko"),&dko,NR);	/* external [K], used for K vrev */
 setval (lookup("dki"),&dki,NR);	/* internal [K], used for K vrev */
 setval (lookup("dpnak"),&dpnak,NR);	/* perm of Na rel to K, for K vrev */
 setval (lookup("dpcak"),&dpcak,NR);	/* perm of Ca rel to K, for K vrev */

					/* see Hille (1992) p 351-353 */
 setval (lookup("dpnaca"),&dpnaca,NR);	/* perm of Na rel to Ca, for Ca vrev */
 setval (lookup("dpkca"),&dpkca,NR);	/* perm of K rel to Ca, for Ca vrev */
					/*  see Hess et al, (1986) */
					/*  and Tsien et al. (1987) */

 setval (lookup("dclo"),&dclo,NR);	/* external [Cl], used for Cl vrev */
 setval (lookup("dcli"),&dcli,NR);	/* internal [Cl], used for Cl vrev */

 setval (lookup("deleccap"),&deleccap,NR);/* f electrode capacitance */
 setval (lookup("dcm"),&dcm,NR);	/* f/cm-cm     membrane capacitance */
 setval (lookup("dri"),&dri,NR);	/* ohm-cm       axial resistance */
 setval (lookup("drm"),&drm,NR);	/* ohm-cm-cm    membrane resistance */
 setval (lookup("drs"),&drs,NR);	/* ohm          electrode resistance */
 setval (lookup("drg"),&drg,NR);	/* ohm-um-um    gap junc resistance */
 setval (lookup("dfta"),&dfta,NR);	/* synaptic time const (msec) */
 setval (lookup("dftah"),&dftah,NR);	/* synaptic high-pass time c (msec) */
 setval (lookup("dsfa"),&dsfa,NR);	/* synaptic number of filters */
 setval (lookup("dftb"),&dftb,NR);	/* synaptic time const (msec) */
 setval (lookup("dsfb"),&dsfb,NR);	/* synaptic number of filters */
 setval (lookup("dscu"),&dscu,NR);	/* unitary conductance of chans */
 setval (lookup("dsvn"),&dsvn,NR);	/* number of syn vesicle rel sites */
 setval (lookup("dscavg"),&dscavg,NR);	/* synaptic ca vesicle release gain */
 setval (lookup("dvg"),&dvg,NR);	/* vesicle release gain */
 setval (lookup("dvsz"),&dvsz,NR);	/* size of vesicles released */
 setval (lookup("dst"),&dst,NR);	/* synaptic threshold (volts) */
 setval (lookup("dsc"),&dsc,NR);	/* synaptic curve */
 setval (lookup("dmaxsyn"),&dmaxsyn,NR);/* maximum cond for synapse (S) */
 setval (lookup("dskd"),&dskd,NR);	/* half-max neurotr saturation (kd)*/
 setval (lookup("dshc"),&dshc,NR);	/* neurotrans binding Hill coeff */
 setval (lookup("dstr"),&dstr,NR);	/* max neurotrans conc for synapse */
 setval (lookup("dsmsgc"),&dsmsgc,NR);	/* max 2nd mesg conc for synapse */
 setval (lookup("dsmrrp"),&dsmrrp,NR);	/* max readily releasible pool */
 setval (lookup("dsrrg"),&dsrrg,NR);	/* readily releasible pool gain mult */
 setval (lookup("dsrrp"),&dsrrp,NR);	/* readily releasible pool */
 setval (lookup("dsms"),&dsms,NR);	/* max sustained release rate (/sec) */
 setval (lookup("dcakd"),&dcakd,NR);	/* half-max Ca sat (kd) at chan */
 setval (lookup("dcahc"),&dcahc,NR);	/* hill coeff for Ca binding at chan */
 setval (lookup("dckd"),&dckd,NR);	/* half-max cG saturation (kd) */
 setval (lookup("dchc"),&dchc,NR);	/* hill coeff for cG binding (kd) */
 setval (lookup("dsn"),&dsn,NR);	/* synaptic gain, lin or exp */
 setval (lookup("dsg"),&dsg,NR);	/* synaptic cG gain */
 setval (lookup("dmaxrod"),&dmaxrod,NR);/* dark cond of rod gated channel (S) */
 setval (lookup("dmaxcon"),&dmaxcon,NR);/* dark cond of cone gated chan (S)*/
 setval (lookup("dnadens"),&dnadens,NR);/* density of Na V-gated chan(S/cm2)*/
 setval (lookup("dkdens"),&dkdens,NR);	/* density of K V-gated chan (S/cm2)*/
 setval (lookup("dcadens"),&dcadens,NR);/* density of Ca V-gated chan(S/cm2)*/
 setval (lookup("dcavoff"),&dcavoff,NR);/* Q10 for voffset due to cao (V) */
 setval (lookup("dcaspvrev"),&dcaspvrev,NR); /* Ca surface potential (mult dcavoff) */
 setval (lookup("dcaphg"),&dcaphg,NR);   /* Ca pH conductance (Barnes & Bui, 1991) */

 						/* Add to vrev. */

					   /* Dan Emerson's CICR parameters */
 setval (lookup("kfCICR"),&kfCICR,NR);     /* passive leak rate const from ryanodine store */
 setval (lookup("vm2CICR"),&vm2CICR,NR);   /* max rate of Ca pump into ryanodine store */
 setval (lookup("vm3CICR"),&vm3CICR,NR);   /* max rate of Ca pump from ryanodine store */
 setval (lookup("k1CICR"),&k1CICR,NR);     /* assoc thresh pump const, ryan store release */
 setval (lookup("k2CICR"),&k2CICR,NR);     /* assoc thresh pump const, ryan store uptake */
 setval (lookup("krCICR"),&krCICR,NR);     /* assoc thresh pump const, ryan store release */
 setval (lookup("kaCICR"),&kaCICR,NR);     /* assoc thresh pump const, ryan store release */
 setval (lookup("nCICR"),&nCICR,NR);	   /* Hill coeff, coop pump binding ryan store uptake */
 setval (lookup("mCICR"),&mCICR,NR);	   /* Hill coeff, coop pump binding ryan store release */
 setval (lookup("pCICR"),&pCICR,NR);	   /* Hill coeff, coop pump binding ryan store release */
 setval (lookup("c1CICR"),&c1CICR,NR);     /* ratio ER vol/cytoplasmic vol */
 setval (lookup("betaIP3"),&betaIP3,NR);   /* init fractional rate of constant IP3 store release*/
 setval (lookup("v1IP3"),&v1IP3, NR);      /* init constant IP3 store flux */

					  /* constants for IP3 from  (Rinzel and * Li, 1994) */

 setval (lookup ("oIP3"),&oIP3,NR);           /* Hill coeff, coop pump binding ip3 store uptake, 2 */
 setval (lookup ("mtypeIP3"),&mtypeIP3,NR);   /* 0 => static (def = use equilibrated inf m value) */
 setval (lookup ("a2IP3"),&a2IP3,NR);         /* IP3 store forward const for Tau M value, */
 setval (lookup ("a3IP3"),&a3IP3,NR);         /* IP3 store constant for IP3 Tau H, 2e5 1/MS */
 setval (lookup ("b2IP3"),&b2IP3,NR);         /* IP3 store backward const for Tau M value */
 setval (lookup ("d1IP3"),&d3IP3,NR);         /* IP3 store const for Q value for h gate, 0.13e-6 M */
 setval (lookup ("d2IP3"),&d3IP3,NR);         /* IP3 store const for Q value for h gate, 1.049e-6 M  */
 setval (lookup ("d3IP3"),&d3IP3,NR);         /* IP3 store const for Q value for h gate, 0.9434e-6 M */
 setval (lookup ("d4IP3"),&d3IP3,NR);         /* IP3 store const for m gate at infinity, 0.82e-6 M */
 setval (lookup ("k3IP3"),&k3IP3,NR);         /* IP3 store kd for Ca uptake, 0.1 S-1) */
 setval (lookup ("v2IP3"),&v2IP3,NR);         /* IP3 store rate for Ca release, 0.11 S-1 */
 setval (lookup ("v3IP3"),&v3IP3,NR);         /* IP3 store constant for Ca uptake, 0.9 M S-1 */
 setval (lookup ("v4IP3"),&v4IP3,NR);         /* IP3 store constant for Ca uptake, 6 M S-1 */

 setval (lookup("ddca"),&ddca,NR);	/* Diffusion const for Ca = 2e-6 */
                                        /*  slower in cytoplasm than in H2O*/
					/* see de Schutter&Smolen "Methods" */
					/*   (from Albritton et al, 1992) */
 setval (lookup("ddcao"),&ddcao,NR);	/* Diffusion const for ext Ca = 2e-6 */
                                        /*  slower in cytoplasm than in H2O*/
					/* see de Schutter&Smolen "Methods" */
					/*   (from Albritton et al, 1992) */
 setval (lookup("ddiacao"),&ddiacao,NR);	/* dia of ext Ca core */
 setval (lookup("dcoo"),&dcoo,NR);	/* external Co++ (M) (for cavoff) */
 setval (lookup("dmgo"),&dmgo,NR);	/* external Mg++ (M) (for cavoff) */
 setval (lookup("dcao"),&dcao,NR);	/* 1.15 mM = external Ca++ (M) */
 setval (lookup("dcai"),&dcai,NR);	/* conc. of internal Ca++ (M) */
 setval (lookup("dtcai"),&dtcai,NR);	/* threshold internal Ca++ for pump (M)*/
 setval (lookup("dcabnd"),&dcabnd,NR);	/* ratio of bound to free calcium */
 setval (lookup("dcabti"),&dcabti,NR);	/* Ca buffer total in first shell */
 setval (lookup("dcabt"),&dcabt,NR);	/* Ca buffer total in shells */
					/*  Yamada, Koch, Adams 1998) */
					/*  Simon, Llinas 1985 */
 setval (lookup("dcabf"),&dcabf,NR);	/* Ca buffer forw rate to Ca.B /M/sec */
 setval (lookup("dcabr"),&dcabr,NR);	/* Ca buffer reverse rate to Ca.B /sec */
 setval (lookup("dcashell"),&dcashell,NR); /* Number of int Ca diffusion shells */
 setval (lookup("dcaoshell"),&dcaoshell,NR); /* Number of ext Ca diffusion shells */
 setval (lookup("dcasarea"),&dcasarea,NR); /* Area of shell for Ca flux */
 setval (lookup("dcapkm"),&dcapkm,NR);	/* 1/2 sat conc. of Ca++ pump (M) */
 setval (lookup("dcavmax"),&dcavmax,NR);/* Vmax for Ca pump (ma/cm2)       */
 setval (lookup("dicafrac"),&dicafrac,NR); /* Fraction of Icapump added to comp Itot   */
 setval (lookup("dcakex"),&dcakex,NR);	/* Rate for Na/Ca exchanger(A/mM4/cm2)*/
 setval (lookup("diexchfrac"),&diexchfrac,NR); /* Fraction of Icaexch added to comp Itot   */
 setval (lookup("dmaxca"),&dmaxca,NR);	/* cond of Ca voltage gated channel */
 setval (lookup("dcaoffs"),&dcaoffs,NR);/* activation offset for Ca chan  */
 setval (lookup("dcataum"),&dcataum,NR);/* activation tau for Ca chan  */
 setval (lookup("dcatauh"),&dcatauh,NR);/* inactivation tau for Ca chan  */
 setval (lookup("dkcataum"),&dcataum,NR);/* activation tau for Kca chan  */
 setval (lookup("dkcatauh"),&dcatauh,NR);/* inactivation tau for Kca chan  */
 setval (lookup("dnatauf"),&dnatauf,NR);/* tau for Na noise in 2-state modl */
 setval (lookup("dktauf"),&dktauf,NR);	/* tau for K noise in 2-state model */
 setval (lookup("dcatauf"),&dcatauf,NR);/* tau for Ca noise in 2-state model */

 setval (lookup("dqm"),&dqm,NR);	/* Q10 for am, bm */
 setval (lookup("dqh"),&dqh,NR);	/* Q10 for ah, bh */
 setval (lookup("dqn"),&dqn,NR);	/* Q10 for an, bn (use dqna/b if 0) */
 setval (lookup("dqna"),&dqna,NR);	/* Q10 for an, formerly 3.2 */
 setval (lookup("dqnb"),&dqnb,NR);	/* Q10 for bn, formerly 2.8 */
 setval (lookup("dqca"),&dqca,NR);	/* Q10 for ac, bc (Ca chan) */
 setval (lookup("dqd"),&dqd,NR);	/* Q10 for ad, bd (KA h) */
 setval (lookup("dqkca"),&dqkca,NR);	/* Q10 for kca a, b */
 setval (lookup("dqsyn"),&dqsyn,NR);	/* Q10 for synaptic chans (ampa,etc)*/
 setval (lookup("dqcavmax"),&dqcavmax,NR);/* Q10 for Ca pump */
 setval (lookup("dqcab"),&dqcab,NR);	/* Q10 for Ca buffer */
 setval (lookup("dqdc"),&dqdc,NR);	/* Q10 for diffusion constant, 1.3 */
 setval (lookup("dqc"),&dqc,NR);	/* Q10 for channel conductance */
					/* Rodriguez BM, et al. (1998) */
					/* J. Gen. Physiol.  112:223   */
 setval (lookup("dqrm"),&dqrm,NR);	/* Q10 for Rm conductance, like dqc for channels */
 setval (lookup("dqri"),&dqri,NR);	/* Q10 for Ri (axial conductance) (Moore et al,1978)*/
 setval (lookup("dqrec"),&dqrec,NR);	/* Q10 for photoreceptor transduction */
 setval (lookup("dqcrec"),&dqcrec,NR);	/* Q10 for photoreceptor transduction */
					/* Baylor, Matthews & Yau, 1983) */

		/* the intention here is to preserve the */
		/* original HH rate functions, which were */
		/* specified for 6.3 deg C, but multiply them to */
 		/* normalize to the more common 22 deg C. */
		/* "dratehhm" multiplies alpham & betam, etc. */

		/* The variable "dbasetc" sets the temperature at */
		/* which voltage-gated channels are normalized. */

		/* Not used for other rate funcs which were originally */
		/* defined for 22 deg C */

 setval (lookup("dbasetc"),&dbasetc,NR);  /* base T for chan kinetics, 22 deg C */
 setval (lookup("dbasetca"),&dbasetca,NR);/* base T for Ca pump kinet, 22 deg C */
 setval (lookup("dbasetsyn"),&dbasetsyn,NR);/* base T for synap kinet, 22 deg C */
 setval (lookup("dbasetdc"),&dbasetdc,NR); /* base T for Ca diffusion, 22 deg C */

  /* Q10 values for 6.3 to 22 deg taken from Frankenhaeuser & Moore (1963) */

 setval (lookup("dratehhm"),&dratehhm,NR); /* mult for HH m rate, 22 deg C */
 setval (lookup("dratehhh"),&dratehhh,NR); /* mult for HH h rate, 22 deg C */
 setval (lookup("dratehhna"),&dratehhna,NR); /* mult for HH alphan, 22 deg C */
 setval (lookup("dratehhnb"),&dratehhnb,NR); /* mult for HH betan, 22 deg C */

 setval (lookup("dqeff"),&dqeff,NR);	/* quantum efficiency, R* per photon */
 setval (lookup("dsintinc"),&dsintinc,NR);/* min time incr for sine waves, sec */
 setval (lookup("dsintres"),&dsintres,NR);/* time res for sine waves, per cyc */

 setval (lookup("dnmdamg"),&dnmdamg,NR);/* Mg++ extracellular conc */
 setval (lookup("dmaxna"),&dmaxna,NR);	/* cond of Na voltage gated channel(S)*/
 setval (lookup("dnaoffsm"),&dnaoffsm,NR);/* m offset voltage for Na chan */
 setval (lookup("dnaoffsh"),&dnaoffsh,NR);/* h offset voltage for Na chan */
 setval (lookup("dnataum"),&dnataum,NR);/* activation tau for Na chan */
 setval (lookup("dnatauh"),&dnatauh,NR);/* inactivation tau for Na chan */
 setval (lookup("dgjoff"),&dgjoff,NR);	/* gap junction voltage offset */
 setval (lookup("dgjtau"),&dgjtau,NR);	/* activation tau for gap junction */
 setval (lookup("dgjatpdec"),&dgjatpdec,NR);/* decrement for atp per time step */
 setval (lookup("dgjnv"),&dgjnv,NR);	/* non-volt-sens fraction of gj cond */
 setval (lookup("dbd1"),&dbd1,NR);	/* voltage multiplier for BK chan alph*/
 setval (lookup("dbd2"),&dbd2,NR);	/* voltage multiplier for BK chan bet */
 setval (lookup("dbk1"),&dbk1,NR);	/* ca multiplier for BK chan alph */
 setval (lookup("dbk2"),&dbk2,NR);	/* ca multiplier for BK chan bet  */
 setval (lookup("dsd1"),&dsd1,NR);	/* voltage multiplier for SK chan alph*/
 setval (lookup("dsd2"),&dsd2,NR);	/* voltage multiplier for SK chan bet */
 setval (lookup("dsk1"),&dsk1,NR);	/* ca multiplier for SK chan alph */
 setval (lookup("dsk2"),&dsk2,NR);	/* ca multiplier for SK chan bet  */
 setval (lookup("dmaxk"),&dmaxk,NR);	/* cond of K voltage gated channel(S)*/
 setval (lookup("dmaxclca"),&dmaxclca,NR);/* cond of K voltage gated channel(S)*/
 setval (lookup("dclcavs"),&dclcavs,NR); /* ClCa chan voltage sens over 120 mV */
 setval (lookup("dclcavsc"),&dclcavsc,NR);/* ClCa chan Ca core voltage sens over 120 mV */
 setval (lookup("dsclcac"),&dsclcac,NR); /* ClCa chan core Ca sens */
 setval (lookup("dnau"),&dnau,NR);	/* unitary cond of Na chan(S)=32pS@33 */
 setval (lookup("dku"),&dku,NR);		/* unitary cond of Kdr chan(S)=15pS@30 */
 setval (lookup("dkau"),&dkau,NR);	/* unitary cond of KA chan (S)=30ps@30 */
 setval (lookup("dkcasu"),&dkcasu,NR);	/* unitary cond of Kca SK ch=22 pS@35*/
 setval (lookup("dkcabu"),&dkcabu,NR);	/* unitary cond of Kca BK ch=115 pS@35*/
 setval (lookup("dclcasu"),&dclcasu,NR);	/* unitary cond of Kca SK ch=22 pS@35*/
 setval (lookup("dkihu"),&dkihu,NR);	/* unitary cond of Ih chan (S) */
 setval (lookup("dcalu"),&dcalu,NR);	/* unit cond of Ca T chan=20 pS @35 */
 setval (lookup("dcatu"),&dcatu,NR);	/* unit cond of Ca L chan=8 pS @35 */
 setval (lookup("dampau"),&dampau,NR);	/* unitary cond of AMPA chan (S) */
 setval (lookup("dgabau"),&dgabau,NR);	/* unitary cond of GABA chan (S) */
 setval (lookup("dglyu"),&dglyu,NR);	/* unitary cond of GLY chan (S) */
 setval (lookup("dpru"),&dpru,NR);	/* unitary cond of photoreceptor chan (S) */
 setval (lookup("dcgmpu"),&dcgmpu,NR);	/* unitary cond of cGMP chan (S) */
 setval (lookup("dkoffsn"),&dkoffsn,NR);/* offset voltage for K chan */
 setval (lookup("dkoffsh"),&dkoffsh,NR);/* offset voltage for K chan */
 setval (lookup("dktaum"),&dktaum,NR);	/* activation tau for K chan */
 setval (lookup("dktauh"),&dktauh,NR);	/* inactivation tau for K chan */
 setval (lookup("dsyntau"),&dsyntau,NR);/* time const for synaptic chans */
 setval (lookup("dsyntauf"),&dsyntauf,NR);/* time const for synaptic chan noise */
 setval (lookup("complam"),&complam,NR);/* fraction of lambda for compartment */
 setval (lookup("lamcrit"),&lamcrit,NR);/* criterion for condensing comps */
 setval (lookup("timinc"),&timinc,NR);	/* time incr for model (sec) */
 setval (lookup("stiminc"),&stiminc,NR);/* synaptic/stim time incr (sec) */
 setval (lookup("ncomps"),&cumcomp,NR);	/* number of compartments */
 setval (lookup("endexp"),&endexp,NR);	/* end of recording for model (sec) */
 setval (lookup("ploti"),&ploti,NR);	/* plot increment (sec) */
 setval (lookup("plmax"),&plmax,NR);	/* max for voltage plots (volts) */
 setval (lookup("plmin"),&plmin,NR);	/* min for voltage plots (volts) */
 setval (lookup("setxmax"),&setxmax,NR);/* override for x axis */
 setval (lookup("setxmin"),&setxmin,NR);/* override for x axis */
 setval (lookup("setymax"),&setymax,NR);/* override for y axis */
 setval (lookup("setymin"),&setymin,NR);/* override for y axis */
 setval (lookup("crit"),&crit,NR);	/* criterion for steps (volts) */
 setval (lookup("relax"),&relax,NR);	/* over-relaxation constant */
 setval (lookup("relincr"),&relincr,NR);/* relax increment */
 setval (lookup("synaptau"),&synaptau,NR);/* synaptic time constant multiplier */
 if (implicit) implfk = 1.0;
 setval (lookup("implfk"),&implfk,NR);	/* =0.5 -> C-N, =1 -> full implicit mode */
 setval (lookup("euler"), &euler,NR);	/* =1 -> euler mode */
 setval (lookup("implicit"),&implicit,NR);/* =1 -> implicit mode */
 setval (lookup("stimelem"),&stimelem,NR);/* "stim" understands neural elements */
 setval (lookup("vcolor"), &vcolor,NR);	/* Use color from voltage */
 setval (lookup("lcolor"), &lcolor,NR);	/* Use color from light flux */
 setval (lookup("lcolor"), &rcolor,NR);	/* Use color from cell region */
 setval (lookup("cacolor"),&cacolor,NR);/* Use color from [Ca]i */
 setval (lookup("nacolor"),&nacolor,NR); /* Use color from Na chan */
 setval (lookup("sgcolor"),&sgcolor,NR); /* Use color from syn cond */
 setval (lookup("srcolor"),&srcolor,NR); /* Use color from syn rate */
 setval (lookup("numstimchan"),&numstimchan,NR); /* threshold for unmasking light */
 setval (lookup("unmaskthr"),&unmaskthr,NR); /* threshold for unmasking light */
 setval (lookup("stimonh"),&stimonh,NR);/* high thresh for transducer vclamp */
 setval (lookup("stimonl"),&stimonl,NR);/* low thresh for transducer vclamp */
 setval (lookup("nozerochan"),&nozerochan,NR);	/* Don't make channels with maxcond=0*/
 setval (lookup("oldsynapz"),&oldsynapz,NR);	/* Revert to old synapse (getcurv after filt1*/
 setval (lookup("use_ghki"),&use_ghki,NR); /* Use GHK current equation for channels */
 setval (lookup("fread_expr"),&fread_expr,NR); /* fread file expressions (no spaces) */
 setval (lookup("debug"), &debug,NR);	/* turn on debugging features */
 setval (lookup("debugz"), &debugz,NR);	/* turn on debugging features */
 setval (lookup("nonodes"), &nonodes,NR); /* no node numbers in display */
 setval (lookup("dashfl"), &dashfl,NR);	/* set dashes on plot */
 setval (lookup("allowrev"), &allowrev,NR);/* allow reverse point plotting */
 setval (lookup("djnoise"), &djnoise,NR);/* > 0 -> Johnson noise for all Rm's */
 setval (lookup("plsep"), &plsep,NR);	/* set multiple separate plotting */
 setval (lookup("disp"), &disp,NR);	/* display control variable = "-d n" */
 setval (lookup("time"), &simtime,NR);	/* time at beginning */
 setval (lookup("srtimestep"), &srtimestep,NR);	/* time step for reading stimuli */
}

/*******************************/

double setq10s(void)

/* Set the Q10s for the simulation according to their variation
   with temperature. According to Fohlmeister Cohen & Newman (2010),
   the Q10s for channel gating and conductance change with temp.
   This overrides the original Q10s.
*/

{
  dqrm = 1.85;         /* Q10 for Rm, Fohlmeister et al. (2010) */
  dqri = 0.8;          /* Q10 for Ri, Trevelyan & Jack (2002) */
  dqrec = 1.0;         /* Q10 for photorec conductance shape */
  dqcrec = 1.8;        /* Q10 for photorec conductance cond */
  if (tempcel >= 22) {
     dqm = dqh = 1.95;
     dqca = 1.95;
     dqn = 1.9;
     dqd = 1.9;
     dqc = 1.85;
  }
  else { /* sim temp < 22 */
     dqm = dqh = 3.5;
     dqca = 3.5;
     dqn = 3.7;
     dqd = 3.5;
     dqc = 3.0;
  }
  vinit();
  return tempcel;
}

/*******************************/

datum set_q10s(void)
{
   datum a;

  varcopyu();
  a.val = setq10s();
  a.vtype = NUMBER;
  return a;
}

