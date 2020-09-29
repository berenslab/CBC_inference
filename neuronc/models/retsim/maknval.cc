/* maknval.n */

/* Script to create neuron parameter file.
   Run "maknval.n" > nval.n" to make.

   Note: 
      To add one column to an existing nval.n file, use nval_addcol
      To add one row to an existing nval.n file, use nval_addrow
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "colors.h"

/*-----------------------------------------*/

int xcone;            /* Identity numbers for neurons */
int xrod; 
int hbat;
int ha;
int hb;
int rbp;
int dbp1;             /* depolarizing cone bipolar */
int dbp2;             /* depolarizing cone bipolar */
int dbp3;             /* depolarizing cone bipolar */
int dbp4;             /* depolarizing cone bipolar */
int hbp1;             /* hyperpolarizing cone bipolar */
int hbp2;             /* hyperpolarizing cone bipolar */
int a17;
int aii;
int sbac;
int am;
int am2;
int am3;
int am4;
int amh;
int amh2;
int ams;
int amhs;
int gca;
int gcb;
int dsgc;
int gcaoff;
int gcboff;
int nceltypes;

/*-----------------------------------------*/

#define PNAMSIZ 550			/* neural parameters */
#define CNAMSIZ 30			/* presynaptic connections */
#define RNAMSIZ 30			/* response types */
#define SNAMSIZ 60			/* postsynaptic connections */
#define TNAMSIZ 50			/* cell types */

char *tname[TNAMSIZ][2];		/* cell type, name parameters */
char *pname[PNAMSIZ][2];		/* neural parameters */
char *rname[RNAMSIZ][2];		/* response types */
char *cname[CNAMSIZ][2];		/* presynaptic connections */
char *sname[SNAMSIZ][2];		/* postsynaptic connections */
int tindx = 0;
int pindx = 0;
int cindx = 0;
int rindx = 1;				/* responses start at 1  = xglut */
int sindx = 0;
int pind = 0;

void printparams(int stopval, char fch, char lch);

#define nval(ct,param) (*(nvalarr+ct*NPARAMS+param))

/*-----------------------------------------*/

char *strtolower (char *str)

{
  int i,len;
  static char lstr[20];

  len = strlen(str);
  for (i=0; i<len; i++) {
    lstr[i] = tolower(str[i]);
  } 
  lstr[i] = '\0';
  return lstr;
}

/*-----------------------------------------*/

/* print out param numbers */

void printparamnums(int stopval)
{
   int j;

  for (j=pind; j<stopval; j++) {
    printf ("#define %-10s %3d\t /* %s */\n",pname[j][0],j,pname[j][1]);
  };
  printf ("\n");
  pind = j;
}

/*-----------------------------------------*/

/* print out param variable defs for nval_var.h */

void printparamvardefs(int stopval)
{
   int j;

  for (j=pind; j<stopval; j++) {
    printf ("extern int _%s;\n",pname[j][0]);
  };
  printf ("\n");
  pind = j;
}

/*-----------------------------------------*/

/* print out param variables for nval_var.cc */

void printparamvars(int stopval)
{
   int j;

  for (j=pind; j<stopval; j++) {
    printf ("int _%-10s = %10s;\n",pname[j][0],pname[j][0]);
  };
  printf ("\n");
  pind = j;
}

/*-----------------------------------------*/

/* print out param variables set for nval_var_set.cc */

void printparamvarset(int stopval)
{
   int j;
   char var1[40];
   char var2[40];

  for (j=pind; j<stopval; j++) {
    sprintf(var1,"\"_%s\",",pname[j][0]);
    sprintf(var2,"_%s",pname[j][0]);
    printf ("setptrn(%-13s &%s);\n",var1,var2);
  };
  printf ("\n");
  pind = j;
}

/*-----------------------------------------*/

void printceltypes(char fch, char lch)

/* print a row of cell type names */

{
   int i;

  printf ("#\n");
  printf ("%c",fch);
  printf ("  ");
  for (i=0; i<nceltypes; i++ ) {
      printf ("%6s  ", strtolower(tname[i][0]));
  };
  if (lch!=' ') printf (" %c",lch);
  printf ("\n");
  printf ("#\n");
};


/*-----------------------------------------*/

int cksiz(int indx,int max,const char *labl)
{
   if (indx >= max) fprintf (stderr,"Overflow, please increase %s\n",labl);
   return indx;
}

/*-----------------------------------------*/

int makt (const char *var, const char *comment)

{
 tname[tindx][0] = (char *)var;
 tname[tindx][1] = (char *)comment;
 return cksiz(tindx++,TNAMSIZ,"TNAMSIZ");
}

int makn (const char *var, const char *comment)

{
 pname[pindx][0] = (char *)var;
 pname[pindx][1] = (char *)comment;
 return cksiz(pindx++,PNAMSIZ,"PNAMSIZ");
}

int makc (const char *var, const char *comment)

{
 cname[cindx][0] = (char *)var;
 cname[cindx][1] = (char *)comment;
 return cksiz(cindx++,CNAMSIZ,"CNAMSIZ");
}

int makr (const char *var, const char *comment)

{
 rname[rindx][0] = (char *)var;
 rname[rindx][1] = (char *)comment;
 return cksiz(rindx++,RNAMSIZ,"RNAMSIZ");
}

int maks (const char *var, const char *comment)

{
 sname[sindx][0] = (char *)var;
 sname[sindx][1] = (char *)comment;
 return cksiz(sindx++,SNAMSIZ,"SNAMSIZ");
}

/*-----------------------------------------*/

char *xsystem (const char *str)

{
    char c;
    int i;
    char *cbuf;
    FILE *cstdout;

#define CBUFSIZ 256

 if (!(cbuf= (char *)malloc(CBUFSIZ))) {
    fprintf (stderr,"maknval: can't allocate char buffer.\n");
    return (NULL);
 }
 if (cstdout=popen(str,"r")) {
   for (i=0; i<10; i++) cbuf[i] = 0;
   for (i=0;(c = getc(cstdout)) != EOF && i<CBUFSIZ; i++) {
     cbuf[i] = c;               /* copy output of pipe into buffer */
   }
   if ((i>0) && ((cbuf[i-1]=='\n') || (cbuf[i-1]=='\r'))) i--;
   cbuf[i++] = (char)NULL;
   pclose(cstdout);
  }
  return (cbuf);
}

/*-----------------------------------------*/

int MAKE;
int MAKE_DEND;
int MAKE_AXON;
int MAKE_DIST;
int NMADE;
int MAXNUM;
int NCOLOR;
int MAXCOV;
int MAXSYNI;
int MAXSYNO;
int DENS;
int REGU;
int MORPH;
int COMPLAM;
int BIOPHYS;
int CHNOISE;
int RATIOK;
int VSTART;
int VREV;
int NRM;
int NRI;

int SOMADIA;
int SOMAZ;
int DENDARB;
int DENDARBZ;
int DENZDIST;
int STRATDIA;
int DTIPDIA;
int DTREEDIA;
int ARBSCALE;
int DENDDIA;
int AXARBT;
int AXARBZ;
int AXTIPDIA;
int AXARBDIA;
int MAXSDIST;

int TAPERSPC;
int TAPERABS;
int NDENDR;
int GROWTHR;
int SEGLEN;

int CELPRE1;
int CONPRE1;
int CELCONV1;
int GROWPOST1;
int CELPRE2;
int CONPRE2;
int CELCONV2;
int GROWPOST2;
int CELPRE3;
int CONPRE3;
int CELCONV3;
int GROWPOST3;
int CELPRE4;
int CONPRE4;
int CELCONV4;
int GROWPOST4;
int CELPRE5;
int CONPRE5;
int CELCONV5;
int GROWPOST5;
int CELPRE6;
int CONPRE6;
int CELCONV6;
int GROWPOST6;
int CELPRE7;
int CONPRE7;
int CELCONV7;
int GROWPOST7;
int CELPRE8;
int CONPRE8;
int CELCONV8;
int GROWPOST8;
int CELPRE9;
int CONPRE9;
int CELCONV9;
int GROWPOST9;
int CELPRE10;
int CONPRE10;
int CELCONV10;
int GROWPOST10;

int CELPOST1;
int CONPOST1;
int CELDIV1;
int GROWPRE1;
int SYNREG1;
int SYNREGP1;
int SYNSPAC1;
int SYNANNI1;
int SYNANNO1;
int SYNANPI1;
int SYNANPO1;
int SYNANG1;
int SYNRNG1;
int USEDYAD1;
int DYADTYP1;
int AUTAPSE1;
int SYNNUM1;
int SENSCA1;
int SRRPOOL1;
int SRRPOOLG1;
int SMRRPOOL1;
int SMAXRATE1;
int SGAIN1;
int SVGAIN1;
int SDURH1;
int SNFILTH1;
int SHGAIN1;
int SHOFFS1;
int SVSIZ1;
int SCOND1;
int SCMUL1;
int SCGRAD1;
int SEGRAD1;
int STHRESH1;
int SVNOISE1;
int SCOV1;
int SDUR1;
int SFALL1;
int SNFILT1;
int STRCONC1;
int SRESP1;
int SPCA1;
int SCAVMAX1;
int SCAKM1;
int SCNFILT1;
int SCDUR1;
int SCGAIN1;
int SCOFF1;
int SCNOISE1;
int SNCHAN1;
int SUNIT1;
int SVREV1;

int CELPOST2;
int CONPOST2;
int CELDIV2;
int GROWPRE2;
int SYNREG2;
int SYNREGP2;
int SYNSPAC2;
int SYNANNI2;
int SYNANNO2;
int SYNANPI2;
int SYNANPO2;
int SYNANG2;
int SYNRNG2;
int USEDYAD2;
int DYADTYP2;
int AUTAPSE2;
int SYNNUM2;
int SENSCA2;
int SRRPOOL2;
int SRRPOOLG2;
int SMRRPOOL2;
int SMAXRATE2;
int SGAIN2;
int SVGAIN2;
int SDURH2;
int SNFILTH2;
int SHGAIN2;
int SHOFFS2;
int SVSIZ2;
int SCOND2;
int SCMUL2;
int SCGRAD2;
int SEGRAD2;
int STHRESH2;
int SVNOISE2;
int SCOV2;
int SDUR2;
int SFALL2;
int SNFILT2;
int STRCONC2;
int SRESP2;
int SPCA2;
int SCAVMAX2;
int SCAKM2;
int SCNFILT2;
int SCDUR2;
int SCGAIN2;
int SCOFF2;
int SCNOISE2;
int SNCHAN2;
int SUNIT2;
int SVREV2;

int CELPOST3;
int CONPOST3;
int CELDIV3;
int GROWPRE3;
int SYNREG3;
int SYNREGP3;
int SYNSPAC3;
int SYNANNI3;
int SYNANNO3;
int SYNANPI3;
int SYNANPO3;
int SYNANG3;
int SYNRNG3;
int USEDYAD3;
int DYADTYP3;
int AUTAPSE3;
int SYNNUM3;
int SENSCA3;
int SRRPOOL3;
int SRRPOOLG3;
int SMRRPOOL3;
int SMAXRATE3;
int SGAIN3;
int SVGAIN3;
int SDURH3;
int SNFILTH3;
int SHGAIN3;
int SHOFFS3;
int SVSIZ3;
int SCOND3;
int SCMUL3;
int SCGRAD3;
int SEGRAD3;
int STHRESH3;
int SVNOISE3;
int SCOV3;
int SDUR3;
int SFALL3;
int SNFILT3;
int STRCONC3;
int SRESP3;
int SPCA3;
int SCAVMAX3;
int SCAKM3;
int SCNFILT3;
int SCDUR3;
int SCGAIN3;
int SCOFF3;
int SCNOISE3;
int SNCHAN3;
int SUNIT3;
int SVREV3;
 
int CELPOST4;
int CONPOST4;
int CELDIV4;
int GROWPRE4;
int SYNREG4;
int SYNREGP4;
int SYNSPAC4;
int SYNANNI4;
int SYNANNO4;
int SYNANPI4;
int SYNANPO4;
int SYNANG4;
int SYNRNG4;
int USEDYAD4;
int DYADTYP4;
int AUTAPSE4;
int SYNNUM4;
int SENSCA4;
int SRRPOOL4;
int SRRPOOLG4;
int SMRRPOOL4;
int SMAXRATE4;
int SGAIN4;
int SVGAIN4;
int SDURH4;
int SNFILTH4;
int SHGAIN4;
int SHOFFS4;
int SVSIZ4;
int SCOND4;
int SCMUL4;
int SCGRAD4;
int SEGRAD4;
int STHRESH4;
int SVNOISE4;
int SCOV4;
int SDUR4;
int SFALL4;
int SNFILT4;
int STRCONC4;
int SRESP4;
int SPCA4;
int SCAVMAX4;
int SCAKM4;
int SCNFILT4;
int SCDUR4;
int SCGAIN4;
int SCOFF4;
int SCNOISE4;
int SNCHAN4;
int SUNIT4;
int SVREV4;
 
int CELPOST5;
int CONPOST5;
int CELDIV5;
int GROWPRE5;
int SYNREG5;
int SYNREGP5;
int SYNSPAC5;
int SYNANNI5;
int SYNANNO5;
int SYNANPI5;
int SYNANPO5;
int SYNANG5;
int SYNRNG5;
int USEDYAD5;
int DYADTYP5;
int AUTAPSE5;
int SYNNUM5;
int SENSCA5;
int SRRPOOL5;
int SRRPOOLG5;
int SMRRPOOL5;
int SMAXRATE5;
int SGAIN5;
int SVGAIN5;
int SDURH5;
int SNFILTH5;
int SHGAIN5;
int SHOFFS5;
int SVSIZ5;
int SCOND5;
int SCMUL5;
int SCGRAD5;
int SEGRAD5;
int STHRESH5;
int SVNOISE5;
int SCOV5;
int SDUR5;
int SFALL5;
int SNFILT5;
int STRCONC5;
int SRESP5;
int SPCA5;
int SCAVMAX5;
int SCAKM5;
int SCNFILT5;
int SCDUR5;
int SCGAIN5;
int SCOFF5;
int SCNOISE5;
int SNCHAN5;
int SUNIT5;
int SVREV5;

int CELPOST6;
int CONPOST6;
int CELDIV6;
int GROWPRE6;
int SYNREG6;
int SYNREGP6;
int SYNSPAC6;
int SYNANNI6;
int SYNANNO6;
int SYNANPI6;
int SYNANPO6;
int SYNANG6;
int SYNRNG6;
int USEDYAD6;
int DYADTYP6;
int AUTAPSE6;
int SYNNUM6;
int SENSCA6;
int SRRPOOL6;
int SRRPOOLG6;
int SMRRPOOL6;
int SMAXRATE6;
int SGAIN6;
int SVGAIN6;
int SDURH6;
int SNFILTH6;
int SHGAIN6;
int SHOFFS6;
int SVSIZ6;
int SCOND6;
int SCMUL6;
int SCGRAD6;
int SEGRAD6;
int STHRESH6;
int SVNOISE6;
int SCOV6;
int SDUR6;
int SFALL6;
int SNFILT6;
int STRCONC6;
int SRESP6;
int SPCA6;
int SCAVMAX6;
int SCAKM6;
int SCNFILT6;
int SCDUR6;
int SCGAIN6;
int SCOFF6;
int SCNOISE6;
int SNCHAN6;
int SUNIT6;
int SVREV6;

int CELPOST7;
int CONPOST7;
int CELDIV7;
int GROWPRE7;
int SYNREG7;
int SYNREGP7;
int SYNSPAC7;
int SYNANNI7;
int SYNANNO7;
int SYNANPI7;
int SYNANPO7;
int SYNANG7;
int SYNRNG7;
int USEDYAD7;
int DYADTYP7;
int AUTAPSE7;
int SYNNUM7;
int SENSCA7;
int SRRPOOL7;
int SRRPOOLG7;
int SMRRPOOL7;
int SMAXRATE7;
int SGAIN7;
int SVGAIN7;
int SDURH7;
int SNFILTH7;
int SHGAIN7;
int SHOFFS7;
int SVSIZ7;
int SCOND7;
int SCMUL7;
int SCGRAD7;
int SEGRAD7;
int STHRESH7;
int SVNOISE7;
int SCOV7;
int SDUR7;
int SFALL7;
int SNFILT7;
int STRCONC7;
int SRESP7;
int SPCA7;
int SCAVMAX7;
int SCAKM7;
int SCNFILT7;
int SCDUR7;
int SCGAIN7;
int SCOFF7;
int SCNOISE7;
int SNCHAN7;
int SUNIT7;
int SVREV7;

int CELPOST8;
int CONPOST8;
int CELDIV8;
int GROWPRE8;
int SYNREG8;
int SYNREGP8;
int SYNSPAC8;
int SYNANNI8;
int SYNANNO8;
int SYNANPI8;
int SYNANPO8;
int SYNANG8;
int SYNRNG8;
int USEDYAD8;
int DYADTYP8;
int AUTAPSE8;
int SYNNUM8;
int SENSCA8;
int SRRPOOL8;
int SRRPOOLG8;
int SMRRPOOL8;
int SMAXRATE8;
int SGAIN8;
int SVGAIN8;
int SDURH8;
int SNFILTH8;
int SHGAIN8;
int SHOFFS8;
int SVSIZ8;
int SCOND8;
int SCMUL8;
int SCGRAD8;
int SEGRAD8;
int STHRESH8;
int SVNOISE8;
int SCOV8;
int SDUR8;
int SFALL8;
int SNFILT8;
int STRCONC8;
int SRESP8;
int SPCA8;
int SCAVMAX8;
int SCAKM8;
int SCNFILT8;
int SCDUR8;
int SCGAIN8;
int SCOFF8;
int SCNOISE8;
int SNCHAN8;
int SUNIT8;
int SVREV8;

int CELPOST9;
int CONPOST9;
int CELDIV9;
int GROWPRE9;
int SYNREG9;
int SYNREGP9;
int SYNSPAC9;
int SYNANNI9;
int SYNANNO9;
int SYNANPI9;
int SYNANPO9;
int SYNANG9;
int SYNRNG9;
int USEDYAD9;
int DYADTYP9;
int AUTAPSE9;
int SYNNUM9;
int SENSCA9;
int SRRPOOL9;
int SRRPOOLG9;
int SMRRPOOL9;
int SMAXRATE9;
int SGAIN9;
int SVGAIN9;
int SDURH9;
int SNFILTH9;
int SHGAIN9;
int SHOFFS9;
int SVSIZ9;
int SCOND9;
int SCMUL9;
int SCGRAD9;
int SEGRAD9;
int STHRESH9;
int SVNOISE9;
int SCOV9;
int SDUR9;
int SFALL9;
int SNFILT9;
int STRCONC9;
int SRESP9;
int SPCA9;
int SCAVMAX9;
int SCAKM9;
int SCNFILT9;
int SCDUR9;
int SCGAIN9;
int SCOFF9;
int SCNOISE9;
int SNCHAN9;
int SUNIT9;
int SVREV9;

int NPARAMS; 

int xglut;
int xampa;
int xampa1;
int xampa2;
int xampa3;
int xampa4;
int xampa5;
int xnmda;
int xnmda2;
int xkainate;
int xmglur6;
int xgaba;
int xgaba1;
int xgaba2;
int xgaba3;
int xgaba4;
int xgly;
int xgapj;
int xdyad;
int nresptypes;

int CELPRE;
int CONPRE;
int CELCONV;
int GROWPOST;
int NCONNP;

int CELPOST;
int CONPOST;
int CELDIV;
int GROWPRE;
int SYNREG;
int SYNREGP;
int SYNSPAC;
int SYNANNI;
int SYNANNO;
int SYNANPI;
int SYNANPO;
int SYNANG;
int SYNRNG;
int USEDYAD;
int DYADTYP;
int AUTAPSE;
int SYNNUM;
int SENSCA;
int SRRPOOL;
int SRRPOOLG;
int SMRRPOOL;
int SMAXRATE;
int SGAIN;
int SVGAIN;
int SDURH;
int SNFILTH;
int SHGAIN;
int SHOFFS;
int SVSIZ;
int SCOND;
int SCMUL;
int SCGRAD;
int SEGRAD;
int STHRESH;
int SVNOISE;
int SCOV;
int SDUR;
int SFALL;
int SNFILT;
int STRCONC;
int SRESP;
int SPCA;
int SCAVMAX;
int SCAKM;
int SCNFILT;
int SCDUR;
int SCGAIN;
int SCOFF;
int SCNOISE;
int SNCHAN;
int SUNIT;
int SVREV;
int NSYNP;

double *nvalarr;

main (int argc, char **argv)

{
   int i,j,k,x,n,p;
 
// x = setvar();

for (j=0; j<TNAMSIZ; j++) {
    tname[j][0] = NULL;
    tname[j][1] = NULL;
};
for (j=0; j<PNAMSIZ; j++) {
    pname[j][0] = NULL;
    pname[j][1] = NULL;
};
for (j=0; j<CNAMSIZ; j++) {
    cname[j][0] = NULL;
    cname[j][1] = NULL; 
};
for (j=0; j<RNAMSIZ; j++) {
    rname[j][0] = NULL;
    rname[j][1] = NULL;
};
for (j=0; j<SNAMSIZ; j++) {
    sname[j][0] = NULL;
    sname[j][1] = NULL;
}


/*-----------------------------------------*/

/* cell types */

xcone     = makt("XCONE",	"/* cones */");
xrod      = makt("XROD",	"/* rods */");
hbat      = makt("HBAT",	"/* hbat */");
ha        = makt("HA",		"/* Type A horizontal cells */");
hb        = makt("HB",		"/* Type B horizontal cells */");
rbp       = makt("RBP",		"/* Rod bipolar cells */");
dbp1      = makt("DBP1",	"/* Depolarizing cone bipolar cell, type 1 */");
dbp2      = makt("DBP2",	"/* Depolarizing cone bipolar cell, type 2 */");
dbp3      = makt("DBP3",	"/* Depolarizing cone bipolar cell, type 3 */");
dbp4      = makt("DBP4",	"/* Depolarizing cone bipolar cell, type 4 */");
hbp1      = makt("HBP1",	"/* Hyperpolarizing bipolar cell, type 1 */");
hbp2      = makt("HBP2",	"/* Hyperpolarizing bipolar cell, type 2 */");
a17       = makt("A17",		"/* A17 amacrine cells, feedback to RBP */");
aii       = makt("AII",		"/* AII amacrine cells */");
sbac      = makt("SBAC",	"/* Starburst amacrine cells */");
am        = makt("AM",		"/* Amacrine cell type 1 */");
am2       = makt("AM2",		"/* Amacrine cell type 2 */");
am3       = makt("AM3",		"/* Amacrine cell type 3 */");
am4       = makt("AM4",		"/* Amacrine cell type 4 */");
amh       = makt("AMH",		"/* Amacrine cells, hyperpolarizing */");
amh2      = makt("AMH2",	"/* Amacrine cells, hyperpolarizing */");
ams       = makt("AMS",		"/* Amacrine cells, small-field */");
amhs      = makt("AMHS",	"/* Amacrine cells, small-field hyperpol */");
gca       = makt("GCA",		"/* Ganglion cells, On-type, alpha */");
gcb       = makt("GCB",		"/* Ganglion cells, On-type, beta */");
dsgc      = makt("DSGC",	"/* Direction-selective ganglion cells */");
gcaoff    = makt("GCAOFF",	"/* Ganglion cells, Off-type, alpha */");
gcboff    = makt("GCBOFF",	"/* Ganglion cells, Off-type, beta */");
nceltypes = makt("NCELTYPES",	"/* Number of cell types */");

 /* defs for "nval(,)" */

MAKE      = makn("MAKE",		"whether to make this cell type");
MAKE_DEND = makn("MAKE_DEND",		"whether to make dendrites");
MAKE_AXON = makn("MAKE_AXON",		"whether to make axon");
MAKE_DIST = makn("MAKE_DIST",		"whether to make axon distal");
NMADE     = makn("NMADE",		"number of cells made");
MAXNUM    = makn("MAXNUM",		"maximum number of cells of this type");
NCOLOR    = makn("NCOLOR",		"color of this cell type for display");
MAXCOV    = makn("MAXCOV",		"max coverage factor (for arrays)");
MAXSYNI   = makn("MAXSYNI",		"max number of syn input cells");
MAXSYNO   = makn("MAXSYNO",		"max number of syn output cells");
DENS      = makn("DENS",		"density of this type (per mm2)");
REGU      = makn("REGU",		"regularity (mean/stdev) of spacing");
MORPH     = makn("MORPH",		"morphology (=0 -> file, or artificial)");
COMPLAM   = makn("COMPLAM",		"compartment size (default=complam)");
BIOPHYS   = makn("BIOPHYS",		"add biophys properties (chan dens file)");
CHNOISE   = makn("CHNOISE",		"add membrane channel noise properties  ");
RATIOK    = makn("RATIOK",		"set K density values as ratio from Na");
VSTART    = makn("VSTART",		"initial resting potential");
VREV      = makn("VREV",		"membrane potential for Rm (VCl)");
NRM	  = makn("NRM",			"the cell's Rm, 0 => use default (drm)");
NRI	  = makn("NRI",			"the cell's Ri, 0 => use default (dri)");

SOMADIA   = makn("SOMADIA",		"Soma diameter");
SOMAZ     = makn("SOMAZ",		"Z location (x,y loc determ. by array)");
DENDARB   = makn("DENDARB",		"type of dendritic tree");
DENDARBZ  = makn("DENDARBZ",		"dendritic arborization level");
DENZDIST  = makn("DENZDIST",		"dendritic arborization z tolerance");
STRATDIA  = makn("STRATDIA",		"stratif. annulus dia (fract of treedia)");
DTIPDIA   = makn("DTIPDIA",		"diameter of dendritic tips");
DTREEDIA  = makn("DTREEDIA",		"diameter of dendritic tree");
ARBSCALE  = makn("ARBSCALE",		"scale for dia of real morph dend tree");
DENDDIA   = makn("DENDDIA",		"dend dia scale for real morph");
AXARBT    = makn("AXARBT",		"type of axonal tree");
AXARBZ    = makn("AXARBZ",		"axonal arborization level");
AXTIPDIA  = makn("AXTIPDIA",		"diameter of axonal tips");
AXARBDIA  = makn("AXARBDIA",		"diameter of axonal arbor");
MAXSDIST  = makn("MAXSDIST",		"maximum synaptic distance");

TAPERSPC  = makn("TAPERSPC",		"space constant of diameter taper");
TAPERABS  = makn("TAPERABS",		"abs diameter for taper");
NDENDR    = makn("NDENDR",		"number of first-order dendrites");
GROWTHR   = makn("GROWTHR",		"distance thresh for growth of dendrites");
SEGLEN    = makn("SEGLEN",		"length of dendrite segments");

CELPRE1   = makn("CELPRE1",		"cell type to connect to (neg, no conn)");
CONPRE1   = makn("CONPRE1",		"connection number of presyn cell");
CELCONV1  = makn("CELCONV1",		"number of presyn cells to connect to");
GROWPOST1 = makn("GROWPOST1",		"grow when making conn from presyn cell");
CELPRE2   = makn("CELPRE2",		"cell type to connect to (neg, no conn)");
CONPRE2   = makn("CONPRE2",		"connection number of presyn cell");
CELCONV2  = makn("CELCONV2",		"number of presyn cells to connect to");
GROWPOST2 = makn("GROWPOST2",		"grow when making conn from presyn cell");
CELPRE3   = makn("CELPRE3",		"cell type to connect to (neg, no conn)");
CONPRE3   = makn("CONPRE3",		"connection number of presyn cell");
CELCONV3  = makn("CELCONV3",		"number of presyn cells to connect to");
GROWPOST3 = makn("GROWPOST3",		"grow when making conn from presyn cell");
CELPRE4   = makn("CELPRE4",		"cell type to connect to (neg, no conn)");
CONPRE4   = makn("CONPRE4",		"connection number of presyn cell");
CELCONV4  = makn("CELCONV4",		"number of presyn cells to connect to");
GROWPOST4 = makn("GROWPOST4",		"grow when making conn from presyn cell");
CELPRE5   = makn("CELPRE5",		"cell type to connect to (neg, no conn)");
CONPRE5   = makn("CONPRE5",		"connection number of presyn cell");
CELCONV5  = makn("CELCONV5",		"number of presyn cells to connect to");
GROWPOST5 = makn("GROWPOST5",		"grow when making conn from presyn cell");
CELPRE6   = makn("CELPRE6",		"cell type to connect to (neg, no conn)");
CONPRE6   = makn("CONPRE6",		"connection number of presyn cell");
CELCONV6  = makn("CELCONV6",		"number of presyn cells to connect to");
GROWPOST6 = makn("GROWPOST6",		"grow when making conn from presyn cell");
CELPRE7   = makn("CELPRE7",		"cell type to connect to (neg, no conn)");
CONPRE7   = makn("CONPRE7",		"connection number of presyn cell");
CELCONV7  = makn("CELCONV7",		"number of presyn cells to connect to");
GROWPOST7 = makn("GROWPOST7",		"grow when making conn from presyn cell");
CELPRE8   = makn("CELPRE8",		"cell type to connect to (neg, no conn)");
CONPRE8   = makn("CONPRE8",		"connection number of presyn cell");
CELCONV8  = makn("CELCONV8",		"number of presyn cells to connect to");
GROWPOST8 = makn("GROWPOST8",		"grow when making conn from presyn cell");
CELPRE9   = makn("CELPRE9",		"cell type to connect to (neg, no conn)");
CONPRE9   = makn("CONPRE9",		"connection number of presyn cell");
CELCONV9  = makn("CELCONV9",		"number of presyn cells to connect to");
GROWPOST9 = makn("GROWPOST9",		"grow when making conn from presyn cell");
CELPRE10  = makn("CELPRE10",		"cell type to connect to (neg, no conn)");
CONPRE10  = makn("CONPRE10",		"connection number of presyn cell");
CELCONV10 = makn("CELCONV10",		"number of presyn cells to connect to");
GROWPOST10= makn("GROWPOST10",		"grow when making conn from presyn cell");

CELPOST1  = makn("CELPOST1",		"cell type to connect to (neg, no conn)");
CONPOST1  = makn("CONPOST1",		"connection number for postsyn cell");
CELDIV1   = makn("CELDIV1",		"number of postsyn cells to connect to");
GROWPRE1  = makn("GROWPRE1",		"grow when making conn to postsyn cell");
SYNREG1   = makn("SYNREG1",		"synaptic region in presyn dendritic tree");
SYNREGP1  = makn("SYNREGP1",		"synaptic region in postsyn dendritic tree");
SYNSPAC1  = makn("SYNSPAC1",		"synaptic spacing in presyn dendritic tree");
SYNANNI1  = makn("SYNANNI1",		"inner rad of annulus in presyn dendr tree");
SYNANNO1  = makn("SYNANNO1",		"outer rad of annulus in presyn dendr tree");
SYNANPI1  = makn("SYNANPI1",		"inner rad of annulus in postsyn dend tree");
SYNANPO1  = makn("SYNANPO1",		"outer rad of annulus in postsyn dend tree");
SYNANG1   = makn("SYNANG1",		"angle for postsynaptic cell");
SYNRNG1   = makn("SYNRNG1",		"range of angles for postsynaptic cell");
USEDYAD1  = makn("USEDYAD1",		"synapse is dyad using preexisting type");
DYADTYP1  = makn("DYADTYP1",		"type of dyad synapse to connect with");
AUTAPSE1  = makn("AUTAPSE1",		"synapse back to presynaptic node");
SYNNUM1   = makn("SYNNUM1",		"number of synapses per connection");
SENSCA1   = makn("SENSCA1",		"synaptic release sensitivity calcium");
SRRPOOL1  = makn("SRRPOOL1",		"synaptic readily releasable pool");
SRRPOOLG1 = makn("SRRPOOLG1",		"synaptic readily releasable pool gain");
SMRRPOOL1 = makn("SMRRPOOL1",		"synaptic max readily releasable pool");
SMAXRATE1 = makn("SMAXRATE1",		"maximum sustained synaptic release rate");
SGAIN1    = makn("SGAIN1",		"synaptic gain");
SVGAIN1   = makn("SVGAIN1",		"synaptic vgain");
SDURH1    = makn("SDURH1",		"synaptic high pass time const.");
SNFILTH1  = makn("SNFILTH1",		"synaptic high pass nfilt");
SHGAIN1   = makn("SHGAIN1",		"synaptic high pass gain");
SHOFFS1   = makn("SHOFFS1",		"synaptic high pass offset");
SVSIZ1    = makn("SVSIZ1",		"synaptic vesicle size");
SCOND1    = makn("SCOND1",		"synaptic conductance");
SCMUL1    = makn("SCMUL1",		"synaptic conductance mult for region");
SCGRAD1   = makn("SCGRAD1",		"synaptic conductance gradient from soma");
SEGRAD1   = makn("SEGRAD1",		"synaptic conductance expon grad fr soma");
STHRESH1  = makn("STHRESH1",		"synaptic threshold");
SVNOISE1  = makn("SVNOISE1",		"1 -> vesicle noise, override, vnoise=0");
SCOV1     = makn("SCOV1",		"1=Poisson, <1->more regular, gamma dist");
SDUR1     = makn("SDUR1",		"synaptic event time const.");
SFALL1    = makn("SFALL1",		"synaptic event fall time const.");
SNFILT1   = makn("SNFILT1",		"synaptic vesicle nfilt");
STRCONC1  = makn("STRCONC1",		"synaptic transmitter concentration.");
SRESP1    = makn("SRESP1",		"synaptic response (ampa,gaba,gj,etc.");
SPCA1     = makn("SPCA1",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX1  = makn("SCAVMAX1",		"synaptic postsyn Ca pump vmax.");
SCAKM1    = makn("SCAKM1",		"synaptic postsyn Ca pump Km.");
SCNFILT1  = makn("SCNFILT1",		"second mesng. nfilt");
SCDUR1    = makn("SCDUR1",		"second mesng. time const.");
SCGAIN1   = makn("SCGAIN1",		"synaptic second messenger gain");
SCOFF1    = makn("SCOFF1",		"synaptic second messenger offset");
SCNOISE1  = makn("SCNOISE1",		"1 -> channel noise, override cnoise=0");
SNCHAN1   = makn("SNCHAN1",		"number of channels");
SUNIT1    = makn("SUNIT1",		"synaptic channel unitary conductace");
SVREV1    = makn("SVREV1",		"synaptic reversal potential");

CELPOST2  = makn("CELPOST2",		"cell type to connect to (neg, no conn)");
CONPOST2  = makn("CONPOST2",		"connection number for postsyn cell");
CELDIV2   = makn("CELDIV2",		"number of postsyn cells to connect to");
GROWPRE2  = makn("GROWPRE2",		"grow when making conn to postsyn cell");
SYNREG2   = makn("SYNREG2",		"synaptic region in presyn dendritic tree");
SYNREGP2  = makn("SYNREGP2",		"synaptic region in postsyn dendritic tree");
SYNSPAC2  = makn("SYNSPAC2",		"synaptic spacing in presyn dendritic tree");
SYNANNI2  = makn("SYNANNI2",		"inner rad of annulus in presyn dendr tree");
SYNANNO2  = makn("SYNANNO2",		"outer rad of annulus in presyn dend tree");
SYNANPI2  = makn("SYNANPI2",		"inner rad of annulus in postsyn dend tree");
SYNANPO2  = makn("SYNANPO2",		"outer rad of annulus in postsyn dend tree");
SYNANG2   = makn("SYNANG2",		"angle for postsynaptic cell");
SYNRNG2   = makn("SYNRNG2",		"range of angles for postsynaptic cell");
USEDYAD2  = makn("USEDYAD2",		"synapse is dyad using preexisting type");
DYADTYP2  = makn("DYADTYP2",		"type of dyad synapse to connect with");
AUTAPSE2  = makn("AUTAPSE2",		"synapse back to presynaptic node");
SYNNUM2   = makn("SYNNUM2",		"number of synapses per connection");
SENSCA2   = makn("SENSCA2",		"synaptic release sensitivity calcium");
SRRPOOL2  = makn("SRRPOOL2",		"synaptic readily releasable pool");
SRRPOOLG2 = makn("SRRPOOLG2",		"synaptic readily releasable pool gain");
SMRRPOOL2 = makn("SMRRPOOL2",		"synaptic max readily releasable pool");
SMAXRATE2 = makn("SMAXRATE2",		"maximum sustained synaptic release rate");
SGAIN2    = makn("SGAIN2",		"synaptic gain");
SVGAIN2   = makn("SVGAIN2",		"synaptic vgain");
SDURH2    = makn("SDURH2",		"synaptic high pass time const.");
SNFILTH2  = makn("SNFILTH2",		"synaptic high pass nfilt");
SHGAIN2   = makn("SHGAIN2",		"synaptic high pass gain");
SHOFFS2   = makn("SHOFFS2",		"synaptic high pass offset");
SVSIZ2    = makn("SVSIZ2",		"synaptic vesicle size");
SCOND2    = makn("SCOND2",		"synaptic conductance");
SCMUL2    = makn("SCMUL2",		"synaptic conductance mult for region");
SCGRAD2   = makn("SCGRAD2",		"synaptic conductance grad from soma");
SEGRAD2   = makn("SEGRAD2",		"synaptic conductance expon grad fr soma");
STHRESH2  = makn("STHRESH2",		"synaptic threshold");
SVNOISE2  = makn("SVNOISE2",		"1 -> vesicle noise, override,vnoise=0");
SCOV2     = makn("SCOV2",		"1=Poisson, <1->more regular, gamma dist");
SDUR2     = makn("SDUR2",		"synaptic event time const.");
SFALL2    = makn("SFALL2",		"synaptic event fall time const.");
SNFILT2   = makn("SNFILT2",		"synaptic vesicle nfilt");
STRCONC2  = makn("STRCONC2",		"synaptic transmitter concentration.");
SRESP2    = makn("SRESP2",		"synaptic response (ampa,gaba,gj,etc.");
SPCA2     = makn("SPCA2",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX2  = makn("SCAVMAX2",		"synaptic postsyn Ca pump vmax.");
SCAKM2    = makn("SCAKM2",		"synaptic postsyn Ca pump Km.");
SCNFILT2  = makn("SCNFILT2",		"second mesng. nfilt");
SCDUR2    = makn("SCDUR2",		"second mesng. time const.");
SCGAIN2   = makn("SCGAIN2",		"synaptic second messenger gain");
SCOFF2    = makn("SCOFF2",		"synaptic second messenger offset");
SCNOISE2  = makn("SCNOISE2",		"1 -> channel noise, override,cnoise=0");
SNCHAN2   = makn("SNCHAN2",		"number of channels");
SUNIT2    = makn("SUNIT2",		"synaptic channel unitary conductace");
SVREV2    = makn("SVREV2",		"synaptic reversal potential");

CELPOST3  = makn("CELPOST3",		"cell type to connect to (neg, no conn)");
CONPOST3  = makn("CONPOST3",		"connection number for postsyn cell");
CELDIV3   = makn("CELDIV3",		"number of postsyn cells to connect to");
GROWPRE3  = makn("GROWPRE3",		"grow when making conn to postsyn cell");
SYNREG3   = makn("SYNREG3",		"synaptic region in presyn dendritic tree");
SYNREGP3  = makn("SYNREGP3",		"synaptic region in postsyn dendritic tree");
SYNSPAC3  = makn("SYNSPAC3",		"synaptic spacing in presyn dendritic tree");
SYNANNI3  = makn("SYNANNI3",		"inner rad of annulus in presyn dendr tree");
SYNANNO3  = makn("SYNANNO3",		"outer rad of annulus in presyn dendr tree");
SYNANPI3  = makn("SYNANPI3",		"inner rad of annulus in postsyn dend tree");
SYNANPO3  = makn("SYNANPO3",		"outer rad of annulus in postsyn dend tree");
SYNANG3   = makn("SYNANG3",		"angle for postsynaptic cell");
SYNRNG3   = makn("SYNRNG3",		"range of angles for postsynaptic cell");
USEDYAD3  = makn("USEDYAD3",		"synapse is dyad using preexisting type");
DYADTYP3  = makn("DYADTYP3",		"type of dyad synapse to connect with");
AUTAPSE3  = makn("AUTAPSE3",		"synapse back to presynaptic node");
SYNNUM3   = makn("SYNNUM3",		"number of synapses per connection");
SENSCA3   = makn("SENSCA3",		"synaptic release sensitivity calcium");
SRRPOOL3  = makn("SRRPOOL3",		"synaptic readily releasable pool");
SRRPOOLG3 = makn("SRRPOOLG3",		"synaptic readily releasable pool gain");
SMRRPOOL3 = makn("SMRRPOOL3",		"synaptic max readily releasable pool");
SMAXRATE3 = makn("SMAXRATE3",		"maximum sustained synaptic release rate");
SGAIN3    = makn("SGAIN3",		"synaptic gain");
SVGAIN3   = makn("SVGAIN3",		"synaptic vgain");
SDURH3    = makn("SDURH3",		"synaptic high pass time const.");
SNFILTH3  = makn("SNFILTH3",		"synaptic high pass nfilt");
SHGAIN3   = makn("SHGAIN3",		"synaptic high pass gain");
SHOFFS3   = makn("SHOFFS3",		"synaptic high pass offset");
SVSIZ3    = makn("SVSIZ3",		"synaptic vesicle size");
SCOND3    = makn("SCOND3",		"synaptic conductance");
SCMUL3    = makn("SCMUL3",		"synaptic conductance mult for region");
SCGRAD3   = makn("SCGRAD3",		"synaptic conductance gradient from soma");
SEGRAD3   = makn("SEGRAD3",		"synaptic conductance expon grad fr soma");
STHRESH3  = makn("STHRESH3",		"synaptic threshold");
SVNOISE3  = makn("SVNOISE3",		"1 -> vesicle noise, override,vnoise=0");
SCOV3     = makn("SCOV3",		"1=Poisson, <1->more regular, gamma dist");
SDUR3     = makn("SDUR3",		"synaptic event time const.");
SFALL3    = makn("SFALL3",		"synaptic event fall time const.");
SNFILT3   = makn("SNFILT3",		"synaptic vesicle nfilt");
STRCONC3  = makn("STRCONC3",		"synaptic transmitter concentration.");
SRESP3    = makn("SRESP3",		"synaptic response (ampa,gaba,gj,etc.");
SPCA3     = makn("SPCA3",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX3  = makn("SCAVMAX3",		"synaptic postsyn Ca pump vmax.");
SCAKM3    = makn("SCAKM3",		"synaptic postsyn Ca pump Km.");
SCNFILT3  = makn("SCNFILT3",		"second mesng. nfilt");
SCDUR3    = makn("SCDUR3",		"second mesng. time const.");
SCGAIN3   = makn("SCGAIN3",		"synaptic second messenger gain");
SCOFF3    = makn("SCOFF3",		"synaptic second messenger offset");
SCNOISE3  = makn("SCNOISE3",		"1 -> channel noise, override,cnoise=0");
SNCHAN3   = makn("SNCHAN3",		"number of channels");
SUNIT3    = makn("SUNIT3",		"synaptic channel unitary conductace");
SVREV3    = makn("SVREV3",		"synaptic reversal potential");
 
CELPOST4  = makn("CELPOST4",		"cell type to connect to (neg, no conn)");
CONPOST4  = makn("CONPOST4",		"connection number for postsyn cell");
CELDIV4   = makn("CELDIV4",		"number of postsyn cells to connect to");
GROWPRE4  = makn("GROWPRE4",		"grow when making conn to postsyn cell");
SYNREG4   = makn("SYNREG4",		"synaptic region in presyn dendritic tree");
SYNREGP4  = makn("SYNREGP4",		"synaptic region in postsyn dendritic tree");
SYNSPAC4  = makn("SYNSPAC4",		"synaptic spacing in presyn dendritic tree");
SYNANNI4  = makn("SYNANNI4",		"inner rad of annulus in presyn dendr tree");
SYNANNO4  = makn("SYNANNO4",		"outer rad of annulus in presyn dendr tree");
SYNANPI4  = makn("SYNANPI4",		"inner rad of annulus in postsyn dend tree");
SYNANPO4  = makn("SYNANPO4",		"outer rad of annulus in postsyn dend tree");
SYNANG4   = makn("SYNANG4",		"angle for postsynaptic cell");
SYNRNG4   = makn("SYNRNG4",		"range of angles for postsynaptic cell");
USEDYAD4  = makn("USEDYAD4",		"synapse is dyad using preexisting type");
DYADTYP4  = makn("DYADTYP4",		"type of dyad synapse to connect with");
AUTAPSE4  = makn("AUTAPSE4",		"synapse back to presynaptic node");
SYNNUM4   = makn("SYNNUM4",		"number of synapses per connection");
SENSCA4   = makn("SENSCA4",		"synaptic release sensitivity calcium");
SRRPOOL4  = makn("SRRPOOL4",		"synaptic readily releasable pool");
SRRPOOLG4 = makn("SRRPOOLG4",		"synaptic readily releasable pool gain");
SMRRPOOL4 = makn("SMRRPOOL4",		"synaptic max readily releasable pool");
SMAXRATE4 = makn("SMAXRATE4",		"maximum sustained synaptic release rate");
SGAIN4    = makn("SGAIN4",		"synaptic gain");
SVGAIN4   = makn("SVGAIN4",		"synaptic vgain");
SDURH4    = makn("SDURH4",		"synaptic high pass time const.");
SNFILTH4  = makn("SNFILTH4",		"synaptic high pass nfilt");
SHGAIN4   = makn("SHGAIN4",		"synaptic high pass gain");
SHOFFS4   = makn("SHOFFS4",		"synaptic high pass offset");
SVSIZ4    = makn("SVSIZ4",		"synaptic vesicle size");
SCOND4    = makn("SCOND4",		"synaptic conductance");
SCMUL4    = makn("SCMUL4",		"synaptic conductance mult for region");
SCGRAD4   = makn("SCGRAD4",		"synaptic conductance gradient from soma");
SEGRAD4   = makn("SEGRAD4",		"synaptic conductance expon grad fr soma");
STHRESH4  = makn("STHRESH4",		"synaptic threshold");
SVNOISE4  = makn("SVNOISE4",		"1 -> vesicle noise, override,vnoise=0");
SCOV4     = makn("SCOV4",		"1=Poisson, <1->more regular, gamma dist");
SDUR4     = makn("SDUR4",		"synaptic event time const.");
SFALL4    = makn("SFALL4",		"synaptic event fall time const.");
SNFILT4   = makn("SNFILT4",		"synaptic vesicle nfilt");
STRCONC4  = makn("STRCONC4",		"synaptic transmitter concentration.");
SRESP4    = makn("SRESP4",		"synaptic response (ampa,gaba,gj,etc.");
SPCA4     = makn("SPCA4",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX4  = makn("SCAVMAX4",		"synaptic postsyn Ca pump vmax.");
SCAKM4    = makn("SCAKM4",		"synaptic postsyn Ca pump Km.");
SCNFILT4  = makn("SCNFILT4",		"second mesng. nfilt");
SCDUR4    = makn("SCDUR4",		"second mesng. time const.");
SCGAIN4   = makn("SCGAIN4",		"synaptic second messenger gain");
SCOFF4    = makn("SCOFF4",		"synaptic second messenger offset");
SCNOISE4  = makn("SCNOISE4",		"1 -> channel noise, override,cnoise=0");
SNCHAN4   = makn("SNCHAN4",		"number of channels");
SUNIT4    = makn("SUNIT4",		"synaptic channel unitary conductace");
SVREV4    = makn("SVREV4",		"synaptic reversal potential");
 
CELPOST5  = makn("CELPOST5",		"cell type to connect to (neg, no conn)");
CONPOST5  = makn("CONPOST5",		"connection number for postsyn cell");
CELDIV5   = makn("CELDIV5",		"number of postsyn cells to connect to");
GROWPRE5  = makn("GROWPRE5",		"grow when making conn to postsyn cell");
SYNREG5   = makn("SYNREG5",		"synaptic region in presyn dendritic tree");
SYNREGP5  = makn("SYNREGP5",		"synaptic region in postsyn dendritic tree");
SYNSPAC5  = makn("SYNSPAC5",		"synaptic spacing in presyn dendritic tree");
SYNANNI5  = makn("SYNANNI5",		"inner rad of annulus in presyn dendr tree");
SYNANNO5  = makn("SYNANNO5",		"outer rad of annulus in presyn dendr tree");
SYNANPI5  = makn("SYNANPI5",		"inner rad of annulus in postsyn dend tree");
SYNANPO5  = makn("SYNANPO5",		"outer rad of annulus in postsyn dend tree");
SYNANG5   = makn("SYNANG5",		"angle for postsynaptic cell");
SYNRNG5   = makn("SYNRNG5",		"range of angles for postsynaptic cell");
USEDYAD5  = makn("USEDYAD5",		"synapse is dyad using preexisting type");
DYADTYP5  = makn("DYADTYP5",		"type of dyad synapse to connect with");
AUTAPSE5  = makn("AUTAPSE5",		"synapse back to presynaptic node");
SYNNUM5   = makn("SYNNUM5",		"number of synapses per connection");
SENSCA5   = makn("SENSCA5",		"synaptic release sensitivity calcium");
SRRPOOL5  = makn("SRRPOOL5",		"synaptic readily releasable pool");
SRRPOOLG5 = makn("SRRPOOLG5",		"synaptic readily releasable pool gain");
SMRRPOOL5 = makn("SMRRPOOL5",		"synaptic max readily releasable pool");
SMAXRATE5 = makn("SMAXRATE5",		"maximum sustained synaptic release rate");
SGAIN5    = makn("SGAIN5",		"synaptic gain");
SVGAIN5   = makn("SVGAIN5",		"synaptic vgain");
SDURH5    = makn("SDURH5",		"synaptic high pass time const.");
SNFILTH5  = makn("SNFILTH5",		"synaptic high pass nfilt");
SHGAIN5   = makn("SHGAIN5",		"synaptic high pass gain");
SHOFFS5   = makn("SHOFFS5",		"synaptic high pass offset");
SVSIZ5    = makn("SVSIZ5",		"synaptic vesicle size");
SCOND5    = makn("SCOND5",		"synaptic conductance");
SCMUL5    = makn("SCMUL5",		"synaptic conductance mult for region");
SCGRAD5   = makn("SCGRAD5",		"synaptic conductance gradient from soma");
SEGRAD5   = makn("SEGRAD5",		"synaptic conductance expon grad fr soma");
STHRESH5  = makn("STHRESH5",		"synaptic threshold");
SVNOISE5  = makn("SVNOISE5",		"1 -> vesicle noise, override,vnoise=0");
SCOV5     = makn("SCOV5",		"1=Poisson, <1->more regular, gamma dist");
SDUR5     = makn("SDUR5",		"synaptic event time const.");
SFALL5    = makn("SFALL5",		"synaptic event fall time const.");
SNFILT5   = makn("SNFILT5",		"synaptic vesicle nfilt");
STRCONC5  = makn("STRCONC5",		"synaptic transmitter concentration.");
SRESP5    = makn("SRESP5",		"synaptic response (ampa,gaba,gj,etc.");
SPCA5     = makn("SPCA5",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX5  = makn("SCAVMAX5",		"synaptic postsyn Ca pump vmax.");
SCAKM5    = makn("SCAKM5",		"synaptic postsyn Ca pump Km.");
SCNFILT5  = makn("SCNFILT5",		"second mesng. nfilt");
SCDUR5    = makn("SCDUR5",		"second mesng. time const.");
SCGAIN5   = makn("SCGAIN5",		"synaptic second messenger gain");
SCOFF5    = makn("SCOFF5",		"synaptic second messenger offset");
SCNOISE5  = makn("SCNOISE5",		"1 -> channel noise, override,cnoise=0");
SNCHAN5   = makn("SNCHAN5",		"number of channels");
SUNIT5    = makn("SUNIT5",		"synaptic channel unitary conductace");
SVREV5    = makn("SVREV5",		"synaptic reversal potential");
 
CELPOST6  = makn("CELPOST6",		"cell type to connect to (neg, no conn)");
CONPOST6  = makn("CONPOST6",		"connection number for postsyn cell");
CELDIV6   = makn("CELDIV6",		"number of postsyn cells to connect to");
GROWPRE6  = makn("GROWPRE6",		"grow when making conn to postsyn cell");
SYNREG6   = makn("SYNREG6",		"synaptic region in presyn dendritic tree");
SYNREGP6  = makn("SYNREGP6",		"synaptic region in postsyn dendritic tree");
SYNSPAC6  = makn("SYNSPAC6",		"synaptic spacing in presyn dendritic tree");
SYNANNI6  = makn("SYNANNI6",		"inner rad of annulus in presyn dendr tree");
SYNANNO6  = makn("SYNANNO6",		"outer rad of annulus in presyn dendr tree");
SYNANPI6  = makn("SYNANPI6",		"inner rad of annulus in postsyn dend tree");
SYNANPO6  = makn("SYNANPO6",		"outer rad of annulus in postsyn dend tree");
SYNANG6   = makn("SYNANG6",		"angle for postsynaptic cell");
SYNRNG6   = makn("SYNRNG6",		"range of angles for postsynaptic cell");
USEDYAD6  = makn("USEDYAD6",		"synapse is dyad using preexisting type");
DYADTYP6  = makn("DYADTYP6",		"type of dyad synapse to connect with");
AUTAPSE6  = makn("AUTAPSE6",		"synapse back to presynaptic node");
SYNNUM6   = makn("SYNNUM6",		"number of synapses per connection");
SENSCA6   = makn("SENSCA6",		"synaptic release sensitivity calcium");
SRRPOOL6  = makn("SRRPOOL6",		"synaptic readily releasable pool");
SRRPOOLG6 = makn("SRRPOOLG6",		"synaptic readily releasable pool gain");
SMRRPOOL6 = makn("SMRRPOOL6",		"synaptic max readily releasable pool");
SMAXRATE6 = makn("SMAXRATE6",		"maximum sustained synaptic release rate");
SGAIN6    = makn("SGAIN6",		"synaptic gain");
SVGAIN6   = makn("SVGAIN6",		"synaptic vgain");
SDURH6    = makn("SDURH6",		"synaptic high pass time const.");
SNFILTH6  = makn("SNFILTH6",		"synaptic high pass nfilt");
SHGAIN6   = makn("SHGAIN6",		"synaptic high pass gain");
SHOFFS6   = makn("SHOFFS6",		"synaptic high pass offset");
SVSIZ6    = makn("SVSIZ6",		"synaptic vesicle size");
SCOND6    = makn("SCOND6",		"synaptic conductance");
SCMUL6    = makn("SCMUL6",		"synaptic conductance mult for region");
SCGRAD6   = makn("SCGRAD6",		"synaptic conductance gradient from soma");
SEGRAD6   = makn("SEGRAD6",		"synaptic conductance expon grad fr soma");
STHRESH6  = makn("STHRESH6",		"synaptic threshold");
SVNOISE6  = makn("SVNOISE6",		"1 -> vesicle noise, override,vnoise=0");
SCOV6     = makn("SCOV6",		"1=Poisson, <1->more regular, gamma dist");
SDUR6     = makn("SDUR6",		"synaptic event time const.");
SFALL6    = makn("SFALL6",		"synaptic event fall time const.");
SNFILT6   = makn("SNFILT6",		"synaptic vesicle nfilt");
STRCONC6  = makn("STRCONC6",		"synaptic transmitter concentration.");
SRESP6    = makn("SRESP6",		"synaptic response (ampa,gaba,gj,etc.");
SPCA6     = makn("SPCA6",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX6  = makn("SCAVMAX6",		"synaptic postsyn Ca pump vmax.");
SCAKM6    = makn("SCAKM6",		"synaptic postsyn Ca pump Km.");
SCNFILT6  = makn("SCNFILT6",		"second mesng. nfilt");
SCDUR6    = makn("SCDUR6",		"second mesng. time const.");
SCGAIN6   = makn("SCGAIN6",		"synaptic second messenger gain");
SCOFF6    = makn("SCOFF6",		"synaptic second messenger offset");
SCNOISE6  = makn("SCNOISE6",		"1 -> channel noise, override,cnoise=0");
SNCHAN6   = makn("SNCHAN6",		"number of channels");
SUNIT6    = makn("SUNIT6",		"synaptic channel unitary conductace");
SVREV6    = makn("SVREV6",		"synaptic reversal potential");
 
CELPOST7  = makn("CELPOST7",		"cell type to connect to (neg, no conn)");
CONPOST7  = makn("CONPOST7",		"connection number for postsyn cell");
CELDIV7   = makn("CELDIV7",		"number of postsyn cells to connect to");
GROWPRE7  = makn("GROWPRE7",		"grow when making conn to postsyn cell");
SYNREG7   = makn("SYNREG7",		"synaptic region in presyn dendritic tree");
SYNREGP7  = makn("SYNREGP7",		"synaptic region in postsyn dendritic tree");
SYNSPAC7  = makn("SYNSPAC7",		"synaptic spacing in presyn dendritic tree");
SYNANNI7  = makn("SYNANNI7",		"inner rad of annulus in presyn dendr tree");
SYNANNO7  = makn("SYNANNO7",		"outer rad of annulus in presyn dendr tree");
SYNANPI7  = makn("SYNANPI7",		"inner rad of annulus in postsyn dend tree");
SYNANPO7  = makn("SYNANPO7",		"outer rad of annulus in postsyn dend tree");
SYNANG7   = makn("SYNANG7",		"angle for postsynaptic cell");
SYNRNG7   = makn("SYNRNG7",		"range of angles for postsynaptic cell");
USEDYAD7  = makn("USEDYAD7",		"synapse is dyad using preexisting type");
DYADTYP7  = makn("DYADTYP7",		"type of dyad synapse to connect with");
AUTAPSE7  = makn("AUTAPSE7",		"synapse back to presynaptic node");
SYNNUM7   = makn("SYNNUM7",		"number of synapses per connection");
SENSCA7   = makn("SENSCA7",		"synaptic release sensitivity calcium");
SRRPOOL7  = makn("SRRPOOL7",		"synaptic readily releasable pool");
SRRPOOLG7 = makn("SRRPOOLG7",		"synaptic readily releasable pool gain");
SMRRPOOL7 = makn("SMRRPOOL7",		"synaptic max readily releasable pool");
SMAXRATE7 = makn("SMAXRATE7",		"maximum sustained synaptic release rate");
SGAIN7    = makn("SGAIN7",		"synaptic gain");
SVGAIN7   = makn("SVGAIN7",		"synaptic vgain");
SDURH7    = makn("SDURH7",		"synaptic high pass time const.");
SNFILTH7  = makn("SNFILTH7",		"synaptic high pass nfilt");
SHGAIN7   = makn("SHGAIN7",		"synaptic high pass gain");
SHOFFS7   = makn("SHOFFS7",		"synaptic high pass offset");
SVSIZ7    = makn("SVSIZ7",		"synaptic vesicle size");
SCOND7    = makn("SCOND7",		"synaptic conductance");
SCMUL7    = makn("SCMUL7",		"synaptic conductance mult for region");
SCGRAD7   = makn("SCGRAD7",		"synaptic conductance gradient from soma");
SEGRAD7   = makn("SEGRAD7",		"synaptic conductance expon grad fr soma");
STHRESH7  = makn("STHRESH7",		"synaptic threshold");
SVNOISE7  = makn("SVNOISE7",		"1 -> vesicle noise, override,vnoise=0");
SCOV7     = makn("SCOV7",		"1=Poisson, <1->more regular, gamma dist");
SDUR7     = makn("SDUR7",		"synaptic event time const.");
SFALL7    = makn("SFALL7",		"synaptic event fall time const.");
SNFILT7   = makn("SNFILT7",		"synaptic vesicle nfilt");
STRCONC7  = makn("STRCONC7",		"synaptic transmitter concentration.");
SRESP7    = makn("SRESP7",		"synaptic response (ampa,gaba,gj,etc.");
SPCA7     = makn("SPCA7",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX7  = makn("SCAVMAX7",		"synaptic postsyn Ca pump vmax.");
SCAKM7    = makn("SCAKM7",		"synaptic postsyn Ca pump Km.");
SCNFILT7  = makn("SCNFILT7",		"second mesng. nfilt");
SCDUR7    = makn("SCDUR7",		"second mesng. time const.");
SCGAIN7   = makn("SCGAIN7",		"synaptic second messenger gain");
SCOFF7    = makn("SCOFF7",		"synaptic second messenger offset");
SCNOISE7  = makn("SCNOISE7",		"1 -> channel noise, override,cnoise=0");
SNCHAN7   = makn("SNCHAN7",		"number of channels");
SUNIT7    = makn("SUNIT7",		"synaptic channel unitary conductace");
SVREV7    = makn("SVREV7",		"synaptic reversal potential");
 
CELPOST8  = makn("CELPOST8",		"cell type to connect to (neg, no conn)");
CONPOST8  = makn("CONPOST8",		"connection number for postsyn cell");
CELDIV8   = makn("CELDIV8",		"number of postsyn cells to connect to");
GROWPRE8  = makn("GROWPRE8",		"grow when making conn to postsyn cell");
SYNREG8   = makn("SYNREG8",		"synaptic region in presyn dendritic tree");
SYNREGP8  = makn("SYNREGP8",		"synaptic region in postsyn dendritic tree");
SYNSPAC8  = makn("SYNSPAC8",		"synaptic spacing in presyn dendritic tree");
SYNANNI8  = makn("SYNANNI8",		"inner rad of annulus in presyn dendr tree");
SYNANNO8  = makn("SYNANNO8",		"outer rad of annulus in presyn dendr tree");
SYNANPI8  = makn("SYNANPI8",		"inner rad of annulus in postsyn dend tree");
SYNANPO8  = makn("SYNANPO8",		"outer rad of annulus in postsyn dend tree");
SYNANG8   = makn("SYNANG8",		"angle for postsynaptic cell");
SYNRNG8   = makn("SYNRNG8",		"range of angles for postsynaptic cell");
USEDYAD8  = makn("USEDYAD8",		"synapse is dyad using preexisting type");
DYADTYP8  = makn("DYADTYP8",		"type of dyad synapse to connect with");
AUTAPSE8  = makn("AUTAPSE8",		"synapse back to presynaptic node");
SYNNUM8   = makn("SYNNUM8",		"number of synapses per connection");
SENSCA8   = makn("SENSCA8",		"synaptic release sensitivity calcium");
SRRPOOL8  = makn("SRRPOOL8",		"synaptic readily releasable pool");
SRRPOOLG8 = makn("SRRPOOLG8",		"synaptic readily releasable pool gain");
SMRRPOOL8 = makn("SMRRPOOL8",		"synaptic max readily releasable pool");
SMAXRATE8 = makn("SMAXRATE8",		"maximum sustained synaptic release rate");
SGAIN8    = makn("SGAIN8",		"synaptic gain");
SVGAIN8   = makn("SVGAIN8",		"synaptic vgain");
SDURH8    = makn("SDURH8",		"synaptic high pass time const.");
SNFILTH8  = makn("SNFILTH8",		"synaptic high pass nfilt");
SHGAIN8   = makn("SHGAIN8",		"synaptic high pass gain");
SHOFFS8   = makn("SHOFFS8",		"synaptic high pass offset");
SVSIZ8    = makn("SVSIZ8",		"synaptic vesicle size");
SCOND8    = makn("SCOND8",		"synaptic conductance");
SCMUL8    = makn("SCMUL8",		"synaptic conductance mult for region");
SCGRAD8   = makn("SCGRAD8",		"synaptic conductance gradient from soma");
SEGRAD8   = makn("SEGRAD8",		"synaptic conductance expon grad fr soma");
STHRESH8  = makn("STHRESH8",		"synaptic threshold");
SVNOISE8  = makn("SVNOISE8",		"1 -> vesicle noise, override, vnoise=0");
SCOV8     = makn("SCOV8",		"1=Poisson, <1->more regular, gamma dist");
SDUR8     = makn("SDUR8",		"synaptic event time const.");
SFALL8    = makn("SFALL8",		"synaptic event fall time const.");
SNFILT8   = makn("SNFILT8",		"synaptic vesicle nfilt");
STRCONC8  = makn("STRCONC8",		"synaptic transmitter concentration.");
SRESP8    = makn("SRESP8",		"synaptic response (ampa,gaba,gj,etc.");
SPCA8     = makn("SPCA8",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX8  = makn("SCAVMAX8",		"synaptic postsyn Ca pump vmax.");
SCAKM8    = makn("SCAKM8",		"synaptic postsyn Ca pump Km.");
SCNFILT8  = makn("SCNFILT8",		"second mesng. nfilt");
SCDUR8    = makn("SCDUR8",		"second mesng. time const.");
SCGAIN8   = makn("SCGAIN8",		"synaptic second messenger gain");
SCOFF8    = makn("SCOFF8",		"synaptic second messenger offset");
SCNOISE8  = makn("SCNOISE8",		"1 -> channel noise, override,cnoise=0");
SNCHAN8   = makn("SNCHAN8",		"number of channels");
SUNIT8    = makn("SUNIT8",		"synaptic channel unitary conductace");
SVREV8    = makn("SVREV8",		"synaptic reversal potential");
 
CELPOST9  = makn("CELPOST9",		"cell type to connect to (neg, no conn)");
CONPOST9  = makn("CONPOST9",		"connection number for postsyn cell");
CELDIV9   = makn("CELDIV9",		"number of postsyn cells to connect to");
GROWPRE9  = makn("GROWPRE9",		"grow when making conn to postsyn cell");
SYNREG9   = makn("SYNREG9",		"synaptic region in presyn dendritic tree");
SYNREGP9  = makn("SYNREGP9",		"synaptic region in postsyn dendritic tree");
SYNSPAC9  = makn("SYNSPAC9",		"synaptic spacing in presyn dendritic tree");
SYNANNI9  = makn("SYNANNI9",		"inner rad of annulus in presyn dendr tree");
SYNANNO9  = makn("SYNANNO9",		"outer rad of annulus in presyn dendr tree");
SYNANPI9  = makn("SYNANPI9",		"inner rad of annulus in postsyn dend tree");
SYNANPO9  = makn("SYNANPO9",		"outer rad of annulus in postsyn dend tree");
SYNANG9   = makn("SYNANG9",		"angle for postsynaptic cell");
SYNRNG9   = makn("SYNRNG9",		"range of angles for postsynaptic cell");
USEDYAD9  = makn("USEDYAD9",		"synapse is dyad using preexisting type");
DYADTYP9  = makn("DYADTYP9",		"type of dyad synapse to connect with");
AUTAPSE9  = makn("AUTAPSE9",		"synapse back to presynaptic node");
SYNNUM9   = makn("SYNNUM9",		"number of synapses per connection");
SENSCA9   = makn("SENSCA9",		"synaptic release sensitivity calcium");
SRRPOOL9  = makn("SRRPOOL9",		"synaptic readily releasable pool");
SRRPOOLG9 = makn("SRRPOOLG9",		"synaptic readily releasable pool gain");
SMRRPOOL9 = makn("SMRRPOOL9",		"synaptic max readily releasable pool");
SMAXRATE9 = makn("SMAXRATE9",		"maximum sustained synaptic release rate");
SGAIN9    = makn("SGAIN9",		"synaptic gain");
SVGAIN9   = makn("SVGAIN9",		"synaptic vgain");
SDURH9    = makn("SDURH9",		"synaptic high pass time const.");
SNFILTH9  = makn("SNFILTH9",		"synaptic high pass nfilt");
SHGAIN9   = makn("SHGAIN9",		"synaptic high pass gain");
SHOFFS9   = makn("SHOFFS9",		"synaptic high pass offset");
SVSIZ9    = makn("SVSIZ9",		"synaptic vesicle size");
SCOND9    = makn("SCOND9",		"synaptic conductance");
SCMUL9    = makn("SCMUL9",		"synaptic conductance mult for region");
SCGRAD9   = makn("SCGRAD9",		"synaptic conductance gradient from soma");
SEGRAD9   = makn("SEGRAD9",		"synaptic conductance expon grad fr soma");
STHRESH9  = makn("STHRESH9",		"synaptic threshold");
SVNOISE9  = makn("SVNOISE9",		"1 -> vesicle noise, override, vnoise=0");
SCOV9     = makn("SCOV9",		"1=Poisson, <1->more regular, gamma dist");
SDUR9     = makn("SDUR9",		"synaptic event time const.");
SFALL9    = makn("SFALL9",		"synaptic event fall time const.");
SNFILT9   = makn("SNFILT9",		"synaptic vesicle nfilt");
STRCONC9  = makn("STRCONC9",		"synaptic transmitter concentration.");
SRESP9    = makn("SRESP9",		"synaptic response (ampa,gaba,gj,etc.");
SPCA9     = makn("SPCA9",		"synaptic postsyn Ca perm (ampa,nmda,etc.");
SCAVMAX9  = makn("SCAVMAX9",		"synaptic postsyn Ca pump vmax.");
SCAKM9    = makn("SCAKM9",		"synaptic postsyn Ca pump Km.");
SCNFILT9  = makn("SCNFILT9",		"second mesng. nfilt");
SCDUR9    = makn("SCDUR9",		"second mesng. time const.");
SCGAIN9   = makn("SCGAIN9",		"synaptic second messenger gain");
SCOFF9    = makn("SCOFF9",		"synaptic second messenger offset");
SCNOISE9  = makn("SCNOISE9",		"1 -> channel noise, override,cnoise=0");
SNCHAN9   = makn("SNCHAN9",		"number of channels");
SUNIT9    = makn("SUNIT9",		"synaptic channel unitary conductace");
SVREV9    = makn("SVREV9",		"synaptic reversal potential");
 
NPARAMS   = makn("NPARAMS",		"number of neuron parameters");


#define NCONNI 10			/* number of possible presynaptic cell types */
#define NCONNO 9			/* number of possible postsynaptic cell types */

// fprintf (stderr,"nparams %g\n", NPARAMS);

#define NBRANCHED 0			/* not branched, separate dendrites for each input */
#define BRANCHED  1			/* branched dendritic tree */
#define HBRANCHED 2			/* highly branched */
#define SBRANCHED 3			/* starburst branching */

/* defs for responses */

xglut     = makr("XGLUT",		"/* generic glutamate response */");
xampa     = makr("XAMPA",		"/* AMPA (type 1) synaptic response */");
xampa1    = makr("XAMPA1",		"/* AMPA type 1 synaptic response */");
xampa2    = makr("XAMPA2",		"/* AMPA type 2 synaptic response */");
xampa3    = makr("XAMPA3",		"/* AMPA type 3 synaptic response */");
xampa4    = makr("XAMPA4",		"/* AMPA type 4 synaptic response */");
xampa5    = makr("XAMPA5",		"/* AMPA type 5 synaptic response */");
xnmda     = makr("XNMDA",		"/* NMDA type 1 synaptic response */");
xnmda2    = makr("XNMDA2",		"/* NMDA type 2 synaptic response */");
xkainate  = makr("XKAINATE",		"/* Kainate synaptic response */");
xmglur6   = makr("XMGLUR6",		"/* mGluR6 synaptic response */");
xgaba     = makr("XGABA",		"/* GABA type 1 synaptic response */");
xgaba1    = makr("XGABA1",		"/* GABA type 1 synaptic response */");
xgaba2    = makr("XGABA2",		"/* GABA type 2 synaptic response */");
xgaba3    = makr("XGABA3",		"/* GABA type 3 synaptic response */");
xgaba4    = makr("XGABA4",		"/* GABA type 4 synaptic response */");
xgly      = makr("XGLY",		"/* Glycine synaptic response */");
xgapj     = makr("XGAPJ",		"/* gap junction synaptic response */");
xdyad     = makr("XDYAD",		"/* dyad synapse (uses other resp type) */");
nresptypes= makr("NRESPTYPES",		"/* number of synaptic types */");

/* defs for connections */
					/* change these here as well as above */
CELPRE   = makc("CELPRE",		"/* cell type to connect to (neg, no conn) */");
CONPRE   = makc("CONPRE",		"/* connection number of presyn cell */");
CELCONV  = makc("CELCONV",		"/* number of presyn cells to connect to */");
GROWPOST = makc("GROWPOST",		"/* grow when making conn from presyn cell */");
NCONNP   = makc("NCONNP",		"/* number of connection parameters */");

CELPOST  = maks("CELPOST",		"/* cell type to connect to (neg, no conn) */");
CONPOST  = maks("CONPOST",		"/* connection number for postsyn cell */");
CELDIV   = maks("CELDIV",		"/* number of postsyn cells to connect to */");
GROWPRE  = maks("GROWPRE",		"/* grow when making conn to postsyn cell */");
SYNREG   = maks("SYNREG",		"/* synaptic region in presyn dendritic tree */");
SYNREGP  = maks("SYNREGP",		"/* synaptic region in postsyn dendritic tree */");
SYNSPAC  = maks("SYNSPAC",		"/* synaptic spacing in presyn dendritic tree */");
SYNANNI  = maks("SYNANNI",		"/* inner dia of annulus in presyn dendr tree */");
SYNANNO  = maks("SYNANNO",		"/* outer dia of annulus in presyn dendr tree */");
SYNANPI  = maks("SYNANPI",		"/* inner dia of annulus in postsyn dend tree */");
SYNANPO  = maks("SYNANPO",		"/* outer dia of annulus in postsyn dend tree */");
SYNANG   = maks("SYNANG",		"/* angle for postsynaptic cell */");
SYNRNG   = maks("SYNRNG",		"/* range of angles for postsynaptic cell */");
USEDYAD  = maks("USEDYAD",		"/* synapse is dyad using preexisting type */");
DYADTYP  = maks("DYADTYP",		"/* type of dyad synapse to connect with */");
AUTAPSE  = maks("AUTAPSE",		"/* synapse back to presynaptic node */");
SYNNUM   = maks("SYNNUM",		"/* number of synapses per connection */");
SENSCA   = maks("SENSCA",		"/* synaptic release sensitivity calcium */");
SRRPOOL  = maks("SRRPOOL",		"/* synaptic readily releasable pool */");
SRRPOOLG = maks("SRRPOOLG",		"/* synaptic readily releasable pool gain */");
SMRRPOOL = maks("SMRRPOOL",		"/* synaptic max readily releasable pool */");
SMAXRATE = maks("SMAXRATE",		"/* maximum sustained synaptic release rate */");
SGAIN    = maks("SGAIN",		"/* synaptic gain */");
SVGAIN   = maks("SVGAIN",		"/* synaptic vgain */");
SDURH    = maks("SDURH",		"/* synaptic high pass time const. */");
SNFILTH  = maks("SNFILTH",		"/* synaptic high pass nfilt */");
SHGAIN   = maks("SHGAIN",		"/* synaptic high pass gain */");
SHOFFS   = maks("SHOFFS",		"/* synaptic high pass offset */");
SVSIZ    = maks("SVSIZ",		"/* synaptic vesicle size */");
SCOND    = maks("SCOND",		"/* synaptic conductance */");
SCMUL    = maks("SCMUL",		"/* synaptic conductance mult for region */");
SCGRAD   = maks("SCGRAD",		"/* synaptic conductance gradient from soma */");
SEGRAD   = maks("SEGRAD",		"/* synaptic conductance expon grad fr soma */");
STHRESH  = maks("STHRESH",		"/* synaptic threshold */");
SVNOISE  = maks("SVNOISE",		"/* 1 -> vesicle noise, override,vnoise=0 */");
SCOV     = maks("SCOV",			"/* 1=Poisson, <1->more regular, gamma dist */");
SDUR     = maks("SDUR",			"/* synaptic event time const. */");
SFALL    = maks("SFALL",		"/* synaptic event fall time const. */");
SNFILT   = maks("SNFILT",		"/* synaptic vesicle nfilt */");
STRCONC  = maks("STRCONC",		"/* synaptic transmitter concentration. */");
SRESP    = maks("SRESP",		"/* synaptic response (ampa,gaba,gj,etc. */");
SPCA     = maks("SPCA",	 		"/* synaptic postsyn Ca perm (ampa,nmda,etc. */");
SCAVMAX  = maks("SCAVMAX",		"/* synaptic postsyn Ca pump vmax. */");
SCAKM    = maks("SCAKM",		"/* synaptic postsyn Ca pump Km. */");
SCNFILT  = maks("SCNFILT",		"/* second mesng. nfilt */");
SCDUR    = maks("SCDUR",		"/* second mesng. time const. */");
SCGAIN   = maks("SCGAIN",		"/* synaptic second messenger gain */");
SCOFF    = maks("SCOFF",		"/* synaptic second messenger offset */");
SCNOISE  = maks("SCNOISE",		"/* 1 -> channel noise, override,cnoise=0 */");
SNCHAN   = maks("SNCHAN",		"/* number of channels */");
SUNIT    = maks("SUNIT",		"/* synaptic channel unitary conductace */");
SVREV    = maks("SVREV",		"/* synaptic reversal potential */");
NSYNP    = maks("NSYNP",		"/* number of synaptic parameters */");


/*-----------------------------------------*/


   nvalarr = (double *)malloc(NPARAMS*nceltypes*sizeof(double));  /* cell build params */

   for (n=0; n<nceltypes; n++) 	/* zero the array */
     for (p=0; p<NPARAMS; p++) 
        nval(n,p) = 0;

	/* Neurons */

  nval(xcone,MAKE_DEND) = 1;	/* Make dendrites */
  nval(xrod, MAKE_DEND) = 1;
  nval(ha,   MAKE_DEND) = 1;
  nval(hb,   MAKE_DEND) = 1;
  nval(hbat, MAKE_DEND) = 1;
  nval(rbp,  MAKE_DEND) = 1;
  nval(dbp1, MAKE_DEND) = 1;
  nval(dbp2, MAKE_DEND) = 1;
  nval(dbp3, MAKE_DEND) = 1;
  nval(dbp4, MAKE_DEND) = 1;
  nval(hbp1, MAKE_DEND) = 1;
  nval(hbp2, MAKE_DEND) = 1;
  nval(a17,  MAKE_DEND) = 1;
  nval(aii,  MAKE_DEND) = 1;
  nval(am,   MAKE_DEND) = 1;
  nval(am2,  MAKE_DEND) = 1;
  nval(am3,  MAKE_DEND) = 1;
  nval(am4,  MAKE_DEND) = 1;
  nval(amh,  MAKE_DEND) = 1;
  nval(amh2, MAKE_DEND) = 1;
  nval(ams,  MAKE_DEND) = 1;
  nval(amhs, MAKE_DEND) = 1;
  nval(sbac, MAKE_DEND) = 1;
  nval(gca,   MAKE_DEND) = 1;
  nval(gcb,   MAKE_DEND) = 1;
  nval(dsgc,  MAKE_DEND) = 1;
  nval(gcaoff,MAKE_DEND) = 1;
  nval(gcboff,MAKE_DEND) = 1;

  nval(xcone,MAKE_AXON) = 1;	/* Make axon */
  nval(xrod, MAKE_AXON) = 1;
  nval(dbp1, MAKE_AXON) = 1;
  nval(hbp1, MAKE_AXON) = 1;
  nval(rbp,  MAKE_AXON) = 1;
  nval(gca,   MAKE_AXON) = 1;
  nval(gcb,   MAKE_AXON) = 1;
  nval(dsgc,  MAKE_AXON) = 1;
  nval(gcaoff,MAKE_AXON) = 1;
  nval(gcboff,MAKE_AXON) = 1;

  nval(xcone,MAKE_DIST) = 0;	/* Axon distal */
  nval(xrod, MAKE_DIST) = 0;
  nval(dbp1, MAKE_DIST) = 0;
  nval(hbp1, MAKE_DIST) = 0;
  nval(rbp,  MAKE_DIST) = 0;
  nval(gca,   MAKE_DIST) = 0;
  nval(gcb,   MAKE_DIST) = 0;
  nval(dsgc,  MAKE_DIST) = 1;
  nval(gcaoff,MAKE_DIST) = 0;
  nval(gcboff,MAKE_DIST) = 0;

  nval(xcone,MAXNUM) = 20000;	/* maximum numbers of cells */
  nval(xrod, MAXNUM) = 50000;
  nval(ha,   MAXNUM) = 100;
  nval(hb,   MAXNUM) = 400;
  nval(hbat, MAXNUM) = 100;
  nval(rbp,  MAXNUM) = 2000;
  nval(dbp1, MAXNUM) = 10000;
  nval(dbp2, MAXNUM) = 10000;
  nval(dbp3, MAXNUM) = 10000;
  nval(dbp4, MAXNUM) = 10000;
  nval(hbp1, MAXNUM) = 10000;
  nval(hbp2, MAXNUM) = 10000;
  nval(a17,  MAXNUM) = 400;
  nval(aii,  MAXNUM) = 400;
  nval(sbac, MAXNUM) = 400;
  nval(am,   MAXNUM) = 400;
  nval(am2,  MAXNUM) = 400;
  nval(am3,  MAXNUM) = 400;
  nval(am4,  MAXNUM) = 400;
  nval(amh,  MAXNUM) = 400;
  nval(amh2, MAXNUM) = 400;
  nval(ams,  MAXNUM) = 2000;
  nval(amhs, MAXNUM) = 2000;
  nval(gca,   MAXNUM) = 100;
  nval(gcb,   MAXNUM) = 100;
  nval(dsgc,  MAXNUM) = 100;
  nval(gcaoff,MAXNUM) = 100;
  nval(gcboff,MAXNUM) = 100;

  nval(xcone,NCOLOR) = RED;	/* color */
  nval(xrod, NCOLOR) = MAGENTA;
  nval(ha,   NCOLOR) = YELLOW;
  nval(hb,   NCOLOR) = MAGENTA;
  nval(hbat, NCOLOR) = BLUE;
  nval(rbp,  NCOLOR) = GREEN;
  nval(dbp1, NCOLOR) = CYAN;
  nval(dbp2, NCOLOR) = CYAN;
  nval(dbp3, NCOLOR) = CYAN;
  nval(dbp4, NCOLOR) = CYAN;
  nval(hbp1, NCOLOR) = LTBLUE;
  nval(hbp2, NCOLOR) = LTBLUE;
  nval(a17,  NCOLOR) = LTCYAN;
  nval(aii,  NCOLOR) = YELLOW;
  nval(sbac, NCOLOR) = GREEN;
  nval(am,   NCOLOR) = GRAY;
  nval(am2,  NCOLOR) = GRAY;
  nval(am3,  NCOLOR) = GRAY;
  nval(am4,  NCOLOR) = GRAY;
  nval(amh,  NCOLOR) = BRTWHT;
  nval(amh2, NCOLOR) = BRTWHT;
  nval(ams,  NCOLOR) = LTCYAN;
  nval(amhs, NCOLOR) = LTRED;
  nval(gca,   NCOLOR) = -1;
  nval(gcb,   NCOLOR) = -1;
  nval(dsgc, NCOLOR)  = -1;
  nval(gcaoff,NCOLOR) = -1;
  nval(gcboff,NCOLOR) = -1;

  nval(xcone,MAXCOV) = 5;	
  nval(xrod, MAXCOV) = 5;
  nval(ha,   MAXCOV) = 10;	/* maximum coverage factors */
  nval(hb,   MAXCOV) = 10;
  nval(hbat, MAXCOV) = 10;
  nval(rbp,  MAXCOV) =  5;
  nval(dbp1, MAXCOV) =  3;
  nval(dbp2, MAXCOV) =  3;
  nval(dbp3, MAXCOV) =  3;
  nval(dbp4, MAXCOV) =  3;
  nval(hbp1, MAXCOV) = 10;
  nval(hbp2, MAXCOV) = 10;
  nval(a17,  MAXCOV) = 50;
  nval(aii,  MAXCOV) = 10;
  nval(sbac, MAXCOV) = 50;
  nval(am,   MAXCOV) = 30;
  nval(am2,  MAXCOV) = 30;
  nval(am3,  MAXCOV) = 30;
  nval(am4,  MAXCOV) = 30;
  nval(amh,  MAXCOV) = 30;
  nval(amh2, MAXCOV) = 30;
  nval(ams,  MAXCOV) = 10;
  nval(amhs, MAXCOV) = 10;
  nval(gca,   MAXCOV) = 10;
  nval(gcb,   MAXCOV) = 10;
  nval(dsgc,  MAXCOV) = 10;
  nval(gcaoff,MAXCOV) = 10;
  nval(gcboff,MAXCOV) = 10;
 
  nval(xcone,MAXSYNO) = 20;		/* maximum output synapses */
  nval(xrod, MAXSYNO) = 10;
  nval(ha,   MAXSYNO) = 800;
  nval(hb,   MAXSYNO) = 400;
  nval(hbat, MAXSYNO) = 1000;
  nval(dbp1, MAXSYNO) = 50;
  nval(dbp2, MAXSYNO) = 50;
  nval(dbp3, MAXSYNO) = 50;
  nval(dbp4, MAXSYNO) = 50;
  nval(hbp1, MAXSYNO) = 50;
  nval(hbp2, MAXSYNO) = 50;
  nval(rbp,  MAXSYNO) = 50;
  nval(a17,  MAXSYNO) = 500;
  nval(aii,  MAXSYNO) = 500;
  nval(sbac, MAXSYNO) = 500;
  nval(am,   MAXSYNO) = 500;
  nval(am2,  MAXSYNO) = 500;
  nval(am3,  MAXSYNO) = 500;
  nval(am4,  MAXSYNO) = 500;
  nval(amh,  MAXSYNO) = 500;
  nval(amh2, MAXSYNO) = 500;
  nval(ams,  MAXSYNO) = 50;
  nval(amhs, MAXSYNO) = 50;

  nval(xcone,MAXSYNI) = 20;		/* maximum input synapses */
  nval(xrod, MAXSYNI) = 10;
  nval(ha,   MAXSYNI) = 800;
  nval(hb,   MAXSYNI) = 400;
  nval(hbat, MAXSYNI) = 1000;
  nval(dbp1, MAXSYNI) = 50;
  nval(dbp2, MAXSYNI) = 50;
  nval(dbp3, MAXSYNI) = 50;
  nval(dbp4, MAXSYNI) = 50;
  nval(hbp1, MAXSYNI) = 50;
  nval(hbp2, MAXSYNI) = 50;
  nval(rbp,  MAXSYNI) = 50;
  nval(a17,  MAXSYNI) = 500;
  nval(aii,  MAXSYNI) = 500;
  nval(sbac, MAXSYNI) = 500;
  nval(am,   MAXSYNI) = 500;
  nval(am2,  MAXSYNI) = 500;
  nval(am3,  MAXSYNI) = 500;
  nval(am4,  MAXSYNI) = 500;
  nval(amh,  MAXSYNI) = 500;
  nval(amh2, MAXSYNI) = 500;
  nval(ams,  MAXSYNI) = 50;
  nval(amhs, MAXSYNI) = 50;
  nval(gca,   MAXSYNI) = 2000;
  nval(gcb,   MAXSYNI) = 1000;
  nval(dsgc,  MAXSYNI) = 1000;
  nval(gcaoff,MAXSYNI) = 2000;
  nval(gcboff,MAXSYNI) = 1000;

  nval(xcone,DENS) = 15000;		/* density, cells/mm2 */
  nval(xrod, DENS) = 250000;
  nval(ha,   DENS) = 220;
  nval(hb,   DENS) = 500;
  nval(hbat, DENS) = 1500;
  nval(rbp,  DENS) = 20000;
  nval(dbp1, DENS) = 3500;
  nval(dbp2, DENS) = 3500;
  nval(dbp3, DENS) = 3500;
  nval(dbp4, DENS) = 3500;
  nval(hbp1, DENS) = 2000;
  nval(hbp2, DENS) = 2000;
  nval(a17,  DENS) = 8000;
  nval(aii,  DENS) = 2500;
  nval(sbac, DENS) = 500;
  nval(am,   DENS) = 200;
  nval(am2,  DENS) = 200;
  nval(am3,  DENS) = 200;
  nval(am4,  DENS) = 200;
  nval(amh,  DENS) = 200;
  nval(amh2, DENS) = 200;
  nval(ams,  DENS) = 500;
  nval(amhs, DENS) = 500;
  nval(gca,   DENS) = 10;
  nval(gcb,   DENS) = 10;
  nval(dsgc,  DENS) = 50;
  nval(gcaoff,DENS) = 20;
  nval(gcboff,DENS) = 20;

  nval(xcone,REGU)= 10;		/* regularity of cell array */
  nval(xrod, REGU) = 20;
  nval(ha,   REGU) = 6;
  nval(hb,   REGU) = 6;
  nval(hbat, REGU) = 6;
  nval(rbp,  REGU) = 8;
  nval(dbp1, REGU) = 8;
  nval(dbp2, REGU) = 8;
  nval(dbp3, REGU) = 8;
  nval(dbp4, REGU) = 8;
  nval(hbp1, REGU) = 8;
  nval(hbp2, REGU) = 8;
  nval(a17,  REGU) = 6;
  nval(aii,  REGU) = 8;
  nval(sbac, REGU) = 8;
  nval(am,   REGU) = 6;
  nval(am2,  REGU) = 6;
  nval(am3,  REGU) = 6;
  nval(am4,  REGU) = 6;
  nval(amh,  REGU) = 6;
  nval(amh2, REGU) = 6;
  nval(ams,  REGU) = 6;
  nval(amhs, REGU) = 6;
  nval(gca,   REGU) = 6;
  nval(gcb,   REGU) = 6;
  nval(dsgc,  REGU) = 6;
  nval(gcaoff,REGU) = 8;
  nval(gcboff,REGU) = 8;

  nval(xcone,VSTART) = -0.0252;	/* initial resting potential */
  nval(xrod, VSTART) = -0.0351;
  nval(ha,   VSTART) = -0.04;
  nval(hb,   VSTART) = -0.04;
  nval(hbat, VSTART) = -0.04;
  nval(rbp,  VSTART) = -0.05;
  nval(dbp1, VSTART) = -0.043;
  nval(dbp2, VSTART) = -0.043;
  nval(dbp3, VSTART) = -0.043;
  nval(dbp4, VSTART) = -0.043;
  nval(hbp1, VSTART) = -0.044;
  nval(hbp2, VSTART) = -0.044;
  nval(a17,  VSTART) = -0.05;
  nval(aii,  VSTART) = -0.05;
  nval(sbac, VSTART) = -0.058;
  nval(am,   VSTART) = -0.05;
  nval(am2,  VSTART) = -0.05;
  nval(am3,  VSTART) = -0.05;
  nval(am4,  VSTART) = -0.05;
  nval(amh,  VSTART) = -0.05;
  nval(amh2, VSTART) = -0.05;
  nval(ams,  VSTART) = -0.05;
  nval(amhs, VSTART) = -0.05;
  nval(gca,   VSTART) = -0.06;
  nval(gcb,   VSTART) = -0.06;
  nval(dsgc,  VSTART) = -0.06;
  nval(gcaoff,VSTART) = -0.06;
  nval(gcboff,VSTART) = -0.06;

  nval(xcone,VREV) = -0.07;	/* membrane reversal potential */
  nval(xrod, VREV) = -0.0945;
  nval(ha,   VREV) = -0.07;
  nval(hb,   VREV) = -0.07;
  nval(hbat, VREV) = -0.07;
  nval(rbp,  VREV) = -0.048;
  nval(dbp1, VREV) = -0.050;
  nval(dbp2, VREV) = -0.050;
  nval(dbp3, VREV) = -0.050;
  nval(dbp4, VREV) = -0.050;
  nval(hbp1, VREV) = -0.052;
  nval(hbp2, VREV) = -0.052;
  nval(a17,  VREV) = -0.07;
  nval(aii,  VREV) = -0.05;
  nval(sbac, VREV) = -0.06;
  nval(am,   VREV) = -0.07;
  nval(am2,  VREV) = -0.07;
  nval(am3,  VREV) = -0.07;
  nval(am4,  VREV) = -0.07;
  nval(amh,  VREV) = -0.07;
  nval(amh2, VREV) = -0.07;
  nval(ams,  VREV) = -0.07;
  nval(amhs, VREV) = -0.07;
  nval(gca,   VREV) = -0.06;
  nval(gcb,   VREV) = -0.06;
  nval(dsgc,  VREV) = -0.06;
  nval(gcaoff,VREV) = -0.06;
  nval(gcboff,VREV) = -0.06;

  nval(xcone,NRM) = 1000;	/* Rm's */
  nval(xrod, NRM) = 1000;
  nval(ha,   NRM) = 10000;
  nval(hb,   NRM) = 10000;
  nval(hbat, NRM) = 20000;
  nval(dbp1, NRM) = 20000;
  nval(dbp2, NRM) = 20000;
  nval(dbp3, NRM) = 20000;
  nval(dbp4, NRM) = 20000;
  nval(hbp1, NRM) = 20000;
  nval(hbp2, NRM) = 20000;
  nval(rbp,  NRM) = 5000;
  nval(a17,  NRM) = 50000;
  nval(aii,  NRM) = 25000;
  nval(sbac, NRM) = 50000;
  nval(am,   NRM) = 50000;
  nval(am2,  NRM) = 50000;
  nval(am3,  NRM) = 50000;
  nval(am4,  NRM) = 50000;
  nval(amh,  NRM) = 50000;
  nval(amh2, NRM) = 50000;
  nval(ams,  NRM) = 50000;
  nval(amhs, NRM) = 50000;
  nval(gca,   NRM) = 50000;
  nval(gcb,   NRM) = 50000;
  nval(dsgc,  NRM) = 10000;
  nval(gcaoff,NRM) = 10000;
  nval(gcboff,NRM) = 10000;

  nval(xcone,NRI) = 0;		/* Ri's, defaults to dri */
  nval(xrod, NRI) = 0;
  nval(ha,   NRI) = 0;
  nval(hb,   NRI) = 0;
  nval(hbat, NRI) = 0;
  nval(dbp1, NRI) = 0;
  nval(dbp2, NRI) = 0;
  nval(dbp3, NRI) = 0;
  nval(dbp4, NRI) = 0;
  nval(hbp1, NRI) = 0;
  nval(hbp2, NRI) = 0;
  nval(rbp,  NRI) = 0;
  nval(a17,  NRI) = 0;
  nval(aii,  NRI) = 0;
  nval(sbac, NRI) = 0;
  nval(am,   NRI) = 0;
  nval(am2,  NRI) = 0;
  nval(am3,  NRI) = 0;
  nval(am4,  NRI) = 0;
  nval(amh,  NRI) = 0;
  nval(amh2, NRI) = 0;
  nval(ams,  NRI) = 0;
  nval(amhs, NRI) = 0;
  nval(gca,   NRI) = 0;
  nval(gcb,   NRI) = 0;
  nval(dsgc,  NRI) = 0;
  nval(gcaoff,NRI) = 0;
  nval(gcboff,NRI) = 0;

  nval(xcone,MORPH) = 1;	/* =1 -> artif morphology */
  nval(xrod, MORPH) = 1;
  nval(ha,   MORPH) = 1;
  nval(hb,   MORPH) = 1;
  nval(hbat, MORPH) = 1;
  nval(rbp,  MORPH) = 1;
  nval(dbp1, MORPH) = 1;
  nval(dbp2, MORPH) = 1;
  nval(dbp3, MORPH) = 1;
  nval(dbp4, MORPH) = 1;
  nval(hbp1, MORPH) = 1;
  nval(hbp2, MORPH) = 1;
  nval(a17,  MORPH) = 1;
  nval(aii,  MORPH) = 1;
  nval(sbac, MORPH) = 0;
  nval(am,   MORPH) = 1;
  nval(am2,  MORPH) = 1;
  nval(am3,  MORPH) = 1;
  nval(am4,  MORPH) = 1;
  nval(amh,  MORPH) = 1;
  nval(amh2, MORPH) = 1;
  nval(ams,  MORPH) = 1;
  nval(amhs, MORPH) = 1;
  nval(gca,   MORPH) = 0;
  nval(gcb,   MORPH) = 0;
  nval(dsgc,  MORPH) = 0;
  nval(gcaoff,MORPH) = 0;
  nval(gcboff,MORPH) = 0;

  nval(xcone,COMPLAM) = 0;	/* compartment size */
  nval(xrod, COMPLAM) = 0;
  nval(ha,   COMPLAM) = 0;
  nval(hb,   COMPLAM) = 0;
  nval(hbat, COMPLAM) = 0;
  nval(rbp,  COMPLAM) = 0;
  nval(dbp1, COMPLAM) = 0;
  nval(dbp2, COMPLAM) = 0;
  nval(dbp3, COMPLAM) = 0;
  nval(dbp4, COMPLAM) = 0;
  nval(hbp1, COMPLAM) = 0;
  nval(hbp2, COMPLAM) = 0;
  nval(a17,  COMPLAM) = 0;
  nval(aii,  COMPLAM) = 0;
  nval(sbac, COMPLAM) = 0;
  nval(am,   COMPLAM) = 0;
  nval(am2,  COMPLAM) = 0;
  nval(am3,  COMPLAM) = 0;
  nval(am4,  COMPLAM) = 0;
  nval(amh,  COMPLAM) = 0;
  nval(amh2, COMPLAM) = 0;
  nval(ams,  COMPLAM) = 0;
  nval(amhs, COMPLAM) = 0;
  nval(gca,  COMPLAM) = 0;
  nval(gcb,  COMPLAM) = 0;
  nval(dsgc, COMPLAM) = 0;
  nval(gcaoff,COMPLAM) = 0;
  nval(gcboff,COMPLAM) = 0;

  nval(xcone,BIOPHYS) = 0;	/* =1 -> add biophysical properties */
  nval(xrod, BIOPHYS) = 0;
  nval(ha,   BIOPHYS) = 0;
  nval(hb,   BIOPHYS) = 0;
  nval(hbat, BIOPHYS) = 0;
  nval(rbp,  BIOPHYS) = 0;
  nval(dbp1, BIOPHYS) = 0;
  nval(dbp2, BIOPHYS) = 0;
  nval(dbp3, BIOPHYS) = 0;
  nval(dbp4, BIOPHYS) = 0;
  nval(hbp1, BIOPHYS) = 0;
  nval(hbp2, BIOPHYS) = 0;
  nval(a17,  BIOPHYS) = 0;
  nval(aii,  BIOPHYS) = 1;
  nval(sbac, BIOPHYS) = 1;
  nval(am,   BIOPHYS) = 0;
  nval(am2,  BIOPHYS) = 0;
  nval(am3,  BIOPHYS) = 0;
  nval(am4,  BIOPHYS) = 0;
  nval(amh,  BIOPHYS) = 0;
  nval(amh2, BIOPHYS) = 0;
  nval(ams,  BIOPHYS) = 0;
  nval(amhs, BIOPHYS) = 0;
  nval(gca,   BIOPHYS) = 1;
  nval(gcb,   BIOPHYS) = 1;
  nval(dsgc,  BIOPHYS) = 1;
  nval(gcaoff,BIOPHYS) = 1;
  nval(gcboff,BIOPHYS) = 1;

  nval(xcone,CHNOISE) = 0;	/* =1 -> add membrane channel noise properties */
  nval(xrod, CHNOISE) = 0;
  nval(ha,   CHNOISE) = 0;
  nval(hb,   CHNOISE) = 0;
  nval(hbat, CHNOISE) = 0;
  nval(rbp,  CHNOISE) = 0;
  nval(dbp1, CHNOISE) = 0;
  nval(dbp2, CHNOISE) = 0;
  nval(dbp3, CHNOISE) = 0;
  nval(dbp4, CHNOISE) = 0;
  nval(hbp1, CHNOISE) = 0;
  nval(hbp2, CHNOISE) = 0;
  nval(a17,  CHNOISE) = 0;
  nval(aii,  CHNOISE) = 0;
  nval(am,   CHNOISE) = 0;
  nval(am2,  CHNOISE) = 0;
  nval(am3,  CHNOISE) = 0;
  nval(am4,  CHNOISE) = 0;
  nval(amh,  CHNOISE) = 0;
  nval(amh2, CHNOISE) = 0;
  nval(ams,  CHNOISE) = 0;
  nval(sbac, CHNOISE) = 0;
  nval(gca,   CHNOISE) = 0;
  nval(gcb,   CHNOISE) = 0;
  nval(dsgc,  CHNOISE) = 0;
  nval(gcaoff,CHNOISE) = 0;
  nval(gcboff,CHNOISE) = 0;

  nval(xcone,NDENDR) = 1;	/* Number of dendrites */
  nval(xrod, NDENDR) = 1;
  nval(ha,   NDENDR) = 4;
  nval(hb,   NDENDR) = 4;
  nval(hbat, NDENDR) = 8;
  nval(rbp,  NDENDR) = 25;
  nval(dbp1, NDENDR) = 2;
  nval(dbp2, NDENDR) = 2;
  nval(dbp3, NDENDR) = 2;
  nval(dbp4, NDENDR) = 2;
  nval(hbp1, NDENDR) = 5;
  nval(hbp2, NDENDR) = 5;
  nval(a17,  NDENDR) = 12;
  nval(aii,  NDENDR) = 5;
  nval(am,   NDENDR) = 5;
  nval(am2,  NDENDR) = 5;
  nval(am3,  NDENDR) = 5;
  nval(am4,  NDENDR) = 5;
  nval(amh,  NDENDR) = 5;
  nval(amh2, NDENDR) = 5;
  nval(ams,  NDENDR) = 5;
  nval(amhs, NDENDR) = 5;
  nval(sbac, NDENDR) = 5;
  nval(gca,   NDENDR) = 7;
  nval(gcb,   NDENDR) = 7;
  nval(dsgc,  NDENDR) = 8;
  nval(gcaoff,NDENDR) = 6;
  nval(gcboff,NDENDR) = 6;

  nval(xcone,SEGLEN) = 0;	/* Length of dendrite segments */	
  nval(xrod, SEGLEN) = 0;
  nval(ha,   SEGLEN) = 0;
  nval(hb,   SEGLEN) = 0;
  nval(hbat, SEGLEN) = 0;
  nval(rbp,  SEGLEN) = 0;
  nval(dbp1, SEGLEN) = 0;
  nval(dbp2, SEGLEN) = 0;
  nval(dbp3, SEGLEN) = 0;
  nval(dbp4, SEGLEN) = 0;
  nval(hbp1, SEGLEN) = 0;
  nval(hbp2, SEGLEN) = 0;
  nval(a17,  SEGLEN) = 0;
  nval(aii,  SEGLEN) = 0;
  nval(am,   SEGLEN) = 5;
  nval(am2,  SEGLEN) = 5;
  nval(am3,  SEGLEN) = 5;
  nval(am4,  SEGLEN) = 5;
  nval(amh,  SEGLEN) = 5;
  nval(amh2, SEGLEN) = 5;
  nval(ams,  SEGLEN) = 0;
  nval(amhs, SEGLEN) = 0;
  nval(sbac, SEGLEN) = 0;
  nval(gca,   SEGLEN) = 0;
  nval(gcb,   SEGLEN) = 0;
  nval(dsgc,  SEGLEN) = 0;
  nval(gcaoff,SEGLEN) = 0;
  nval(gcboff,SEGLEN) = 0;

  nval(gca,   GROWTHR) = 3;
  nval(dsgc, GROWTHR) = 3;
  nval(gcaoff,GROWTHR) = 3;

  nval(rbp,  AXTIPDIA) = 0.4;
  nval(dbp1, AXTIPDIA) = 0.4;
  nval(dbp2, AXTIPDIA) = 0.4;
  nval(dbp3, AXTIPDIA) = 0.4;
  nval(dbp4, AXTIPDIA) = 0.4;
  nval(hbp1, AXTIPDIA) = 0.4;
  nval(hbp2, AXTIPDIA) = 0.4;
  nval(ams,  AXTIPDIA) = 0.4;

  nval(xcone,AXARBDIA)  = 6.5;		/* axon arborization diameter */
  nval(xrod, AXARBDIA)  = 5;
  nval(ha,   AXARBDIA)  = 150;
  nval(hb,   AXARBDIA)  = 100;
  nval(rbp,  AXARBDIA)  = 15;
  nval(dbp1, AXARBDIA)  = 24;
  nval(dbp2, AXARBDIA)  = 24;
  nval(dbp3, AXARBDIA)  = 24;
  nval(dbp4, AXARBDIA)  = 24;
  nval(hbp1, AXARBDIA)  = 24;
  nval(hbp2, AXARBDIA)  = 24;
  nval(aii,  AXARBDIA)  = 0;
  nval(ams,  AXARBDIA)  = 24;

  nval(xcone,MAXSDIST) = 20;	/* maximum synaptic distance between cells */
  nval(xrod, MAXSDIST) = 10;
  nval(ha,   MAXSDIST) = 10;
  nval(hb,   MAXSDIST) = 10;
  nval(hbat, MAXSDIST) = 10;
  nval(rbp,  MAXSDIST) = 10;
  nval(dbp1, MAXSDIST) = 15;
  nval(dbp2, MAXSDIST) = 15;
  nval(dbp3, MAXSDIST) = 15;
  nval(dbp4, MAXSDIST) = 15;
  nval(hbp1, MAXSDIST) = 10;
  nval(hbp2, MAXSDIST) = 10;
  nval(a17,  MAXSDIST) = 10;
  nval(aii,  MAXSDIST) = 8;
  nval(sbac, MAXSDIST) = 10;
  nval(am,   MAXSDIST) = 10;
  nval(am2,  MAXSDIST) = 10;
  nval(am3,  MAXSDIST) = 10;
  nval(am4,  MAXSDIST) = 10;
  nval(amh,  MAXSDIST) = 10;
  nval(amh2, MAXSDIST) = 10;
  nval(ams,  MAXSDIST) = 10;
  nval(amhs, MAXSDIST) = 10;
  nval(gca,  MAXSDIST) = 10;
  nval(gcb,  MAXSDIST) = 10;
  nval(dsgc, MAXSDIST) = 10;
  nval(gcaoff,MAXSDIST) = 10;
  nval(gcboff,MAXSDIST) = 10;

  nval(xcone,TAPERSPC) = 1;	/* dendritic taper, space constant */
  nval(xrod, TAPERSPC) = 1;
  nval(ha,   TAPERSPC) = 70;
  nval(hb,   TAPERSPC) = 50;
  nval(hbat, TAPERSPC) = 80;
  nval(rbp,  TAPERSPC) = 50;
  nval(dbp1, TAPERSPC) = 10;
  nval(dbp2, TAPERSPC) = 10;
  nval(dbp3, TAPERSPC) = 10;
  nval(dbp4, TAPERSPC) = 10;
  nval(hbp1, TAPERSPC) = 10;
  nval(hbp2, TAPERSPC) = 10;
  nval(a17,  TAPERSPC) = 50;
  nval(aii,  TAPERSPC) = 50;
  nval(sbac, TAPERSPC) = 300;
  nval(am,   TAPERSPC) = 50;
  nval(am2,  TAPERSPC) = 50;
  nval(am3,  TAPERSPC) = 50;
  nval(am4,  TAPERSPC) = 50;
  nval(amh,  TAPERSPC) = 50;
  nval(amh2, TAPERSPC) = 50;
  nval(ams,  TAPERSPC) = 10;
  nval(amhs, TAPERSPC) = 10;
  nval(gca,   TAPERSPC) = 60;
  nval(gcb,   TAPERSPC) = 60;
  nval(dsgc, TAPERSPC) = 20;
  nval(gcaoff,TAPERSPC) = 60;
  nval(gcboff,TAPERSPC) = 60;

  nval(xcone,TAPERABS) = 1;	/* dendritic taper, absolute diameters */
  nval(xrod, TAPERABS) = 1;
  nval(ha,   TAPERABS) = 2;
  nval(hb,   TAPERABS) = 1;
  nval(hbat, TAPERABS) = 2;
  nval(rbp,  TAPERABS) = 0.8;
  nval(dbp1, TAPERABS) = 1;
  nval(dbp2, TAPERABS) = 1;
  nval(dbp3, TAPERABS) = 1;
  nval(dbp4, TAPERABS) = 1;
  nval(hbp1, TAPERABS) = 1;
  nval(hbp2, TAPERABS) = 1;
  nval(a17,  TAPERABS) = 0.2;
  nval(aii,  TAPERABS) = 0.5;
  nval(sbac, TAPERABS) = 0.2;
  nval(am,   TAPERABS) = 0.2;
  nval(am2,  TAPERABS) = 0.2;
  nval(am3,  TAPERABS) = 0.2;
  nval(am4,  TAPERABS) = 0.2;
  nval(amh,  TAPERABS) = 0.2;
  nval(amh2, TAPERABS) = 0.2;
  nval(ams,  TAPERABS) = 1;
  nval(amhs, TAPERABS) = 1;
  nval(gca,   TAPERABS) = 6;
  nval(gcb,   TAPERABS) = 6;
  nval(dsgc, TAPERABS) = 3;
  nval(gcaoff,TAPERABS) = 6;
  nval(gcboff,TAPERABS) = 6;

  nval(xcone,SOMADIA) = 3;	/* soma diameter */
  nval(xrod, SOMADIA) = 5;
  nval(ha,   SOMADIA) = 10;
  nval(hb,   SOMADIA) = 8;
  nval(hbat, SOMADIA) = 2;
  nval(rbp,  SOMADIA) = 7;
  nval(dbp1, SOMADIA) = 7;
  nval(dbp2, SOMADIA) = 7;
  nval(dbp3, SOMADIA) = 7;
  nval(dbp4, SOMADIA) = 7;
  nval(hbp1, SOMADIA) = 7;
  nval(hbp2, SOMADIA) = 7;
  nval(a17,  SOMADIA) = 9;
  nval(aii,  SOMADIA) = 8;
  nval(sbac, SOMADIA) = 8;
  nval(am,   SOMADIA) = 7;
  nval(am2,  SOMADIA) = 7;
  nval(am3,  SOMADIA) = 7;
  nval(am4,  SOMADIA) = 7;
  nval(amh,  SOMADIA) = 7;
  nval(amh2, SOMADIA) = 7;
  nval(ams,  SOMADIA) = 7;
  nval(amhs, SOMADIA) = 7;
  nval(gca,   SOMADIA) = 20;
  nval(gcb,   SOMADIA) = 15;
  nval(dsgc,  SOMADIA) = 15;
  nval(gcaoff,SOMADIA) = 20;
  nval(gcboff,SOMADIA) = 15;

  nval(xcone,SOMAZ) = 0;	/* soma Z position, layer within retina */
  nval(xrod, SOMAZ) = 0;
  nval(ha,   SOMAZ) = -8;
  nval(hb,   SOMAZ) = -6;
  nval(hbat, SOMAZ) = -4;
  nval(rbp,  SOMAZ) = -8;
  nval(dbp1, SOMAZ) = -8;
  nval(dbp2, SOMAZ) = -8;
  nval(dbp3, SOMAZ) = -8;
  nval(dbp4, SOMAZ) = -8;
  nval(hbp1, SOMAZ) = -8;
  nval(hbp2, SOMAZ) = -8;
  nval(a17,  SOMAZ) = -10;
  nval(aii,  SOMAZ) = -10;
  nval(sbac, SOMAZ) = -37;
  nval(am,   SOMAZ) = -10;
  nval(am2,  SOMAZ) = -10;
  nval(am3,  SOMAZ) = -10;
  nval(am4,  SOMAZ) = -10;
  nval(amh,  SOMAZ) = -10;
  nval(amh2, SOMAZ) = -10;
  nval(ams,  SOMAZ) = -10;
  nval(amhs, SOMAZ) = -10;
  nval(gca,   SOMAZ) = -45;
  nval(gcb,   SOMAZ) = -45;
  nval(dsgc,  SOMAZ) = -45;
  nval(gcaoff,SOMAZ) = -40;
  nval(gcboff,SOMAZ) = -40;

  nval(ha,   DENDARB) = HBRANCHED;	/* type of dendritic tree */
  nval(hb,   DENDARB) = HBRANCHED;
  nval(hbat, DENDARB) = HBRANCHED;
  nval(rbp,  DENDARB) = NBRANCHED;
  nval(dbp1, DENDARB) = NBRANCHED;
  nval(dbp2, DENDARB) = NBRANCHED;
  nval(dbp3, DENDARB) = NBRANCHED;
  nval(dbp4, DENDARB) = NBRANCHED;
  nval(hbp1, DENDARB) = NBRANCHED;
  nval(hbp2, DENDARB) = NBRANCHED;
  nval(a17,  DENDARB) = HBRANCHED;
  nval(aii,  DENDARB) = HBRANCHED;
  nval(sbac, DENDARB) = SBRANCHED;
  nval(am,   DENDARB) = BRANCHED;
  nval(am2,  DENDARB) = BRANCHED;
  nval(am3,  DENDARB) = BRANCHED;
  nval(am4,  DENDARB) = BRANCHED;
  nval(amh,  DENDARB) = BRANCHED;
  nval(amh2, DENDARB) = BRANCHED;
  nval(ams,  DENDARB) = NBRANCHED;
  nval(amhs, DENDARB) = NBRANCHED;
  nval(gca,   DENDARB) = BRANCHED;
  nval(gcb,   DENDARB) = BRANCHED;
  nval(dsgc,  DENDARB) = BRANCHED;
  nval(gcaoff,DENDARB) = BRANCHED;
  nval(gcboff,DENDARB) = BRANCHED;

  nval(xcone,DENDARBZ) = 0;	/* dendritic arborization layer within retina */
  nval(xrod, DENDARBZ) = 0;
  nval(ha,   DENDARBZ) = -6;
  nval(hb,   DENDARBZ) = -5;
  nval(hbat, DENDARBZ) = -4;
  nval(rbp,  DENDARBZ) = -1;
  nval(dbp1, DENDARBZ) = -1;
  nval(dbp2, DENDARBZ) = -1;
  nval(dbp3, DENDARBZ) = -1;
  nval(dbp4, DENDARBZ) = -1;
  nval(hbp1, DENDARBZ) = -1;
  nval(hbp2, DENDARBZ) = -1;
  nval(a17,  DENDARBZ) = -20;
  nval(aii,  DENDARBZ) = -25;
  nval(sbac, DENDARBZ) = -30;
  nval(am,   DENDARBZ) = -20;
  nval(am2,  DENDARBZ) = -20;
  nval(am3,  DENDARBZ) = -20;
  nval(am4,  DENDARBZ) = -20;
  nval(amh,  DENDARBZ) = -20;
  nval(amh2, DENDARBZ) = -20;
  nval(ams,  DENDARBZ) = -1;
  nval(amhs, DENDARBZ) = -1;
  nval(gca,   DENDARBZ) = -30;
  nval(gcb,   DENDARBZ) = -30;
  nval(dsgc,  DENDARBZ) = -30;
  nval(gcaoff,DENDARBZ) = -30;
  nval(gcboff,DENDARBZ) = -30;

  nval(xcone,DENZDIST) = 10;	/* dendritic z dist tolerance for synaptic input */
  nval(xrod, DENZDIST) = 10;
  nval(ha,   DENZDIST) = 10;
  nval(hb,   DENZDIST) = 10;
  nval(hbat, DENZDIST) = 10;
  nval(rbp,  DENZDIST) = 10;
  nval(dbp1, DENZDIST) = 10;
  nval(dbp2, DENZDIST) = 10;
  nval(dbp3, DENZDIST) = 10;
  nval(dbp4, DENZDIST) = 10;
  nval(hbp1, DENZDIST) = 10;
  nval(hbp2, DENZDIST) = 10;
  nval(a17,  DENZDIST) = 10;
  nval(aii,  DENZDIST) = 10;
  nval(sbac, DENZDIST) = 5;
  nval(am,   DENZDIST) = 10;
  nval(am2,  DENZDIST) = 10;
  nval(am3,  DENZDIST) = 10;
  nval(am4,  DENZDIST) = 10;
  nval(amh,  DENZDIST) = 10;
  nval(amh2, DENZDIST) = 10;
  nval(ams,  DENZDIST) = 10;
  nval(amhs, DENZDIST) = 10;
  nval(gca,  DENZDIST) = 15;
  nval(gcb,  DENZDIST) = 10;
  nval(dsgc, DENZDIST) = 10;
  nval(gcaoff,DENZDIST) = 15;
  nval(gcboff,DENZDIST) = 10;

  nval(rbp,   STRATDIA) = 0.3;
  nval(dbp1,  STRATDIA) = 0.3;
  nval(dbp2,  STRATDIA) = 0.3;
  nval(dbp3,  STRATDIA) = 0.3;
  nval(dbp4,  STRATDIA) = 0.3;
  nval(hbp1,  STRATDIA) = 0.3;
  nval(hbp2,  STRATDIA) = 0.3;
  nval(a17,   STRATDIA) = 0.2;
  nval(aii,   STRATDIA) = 0.2;
  nval(sbac,  STRATDIA) = 0.2;
  nval(am,    STRATDIA) = 0.2;
  nval(am2,   STRATDIA) = 0.2;
  nval(am3,   STRATDIA) = 0.2;
  nval(am4,   STRATDIA) = 0.2;
  nval(amh,   STRATDIA) = 0.2;
  nval(amh2,  STRATDIA) = 0.2;
  nval(ams,   STRATDIA) = 0.2;
  nval(amhs,  STRATDIA) = 0.2;
  nval(gca,   STRATDIA) = 0.25;
  nval(gcb,   STRATDIA) = 0.25;
  nval(dsgc,  STRATDIA) = 0.95;
  nval(gcaoff,STRATDIA) = 0.25;
  nval(gcboff,STRATDIA) = 0.25;

  nval(xrod,  DTIPDIA) = 0.1;
  nval(xcone, DTIPDIA) = 0.2;
  nval(ha,    DTIPDIA) = 0.2;	/* dendritic tip diameters */
  nval(hb,    DTIPDIA) = 0.2;
  nval(hbat,  DTIPDIA) = 0.2;
  nval(rbp,   DTIPDIA) = 0.2;
  nval(dbp1,  DTIPDIA) = 0.4;
  nval(dbp2,  DTIPDIA) = 0.4;
  nval(dbp3,  DTIPDIA) = 0.4;
  nval(dbp4,  DTIPDIA) = 0.4;
  nval(hbp1,  DTIPDIA) = 0.4;
  nval(hbp2,  DTIPDIA) = 0.4;
  nval(a17,   DTIPDIA) = 0.2;
  nval(aii,   DTIPDIA) = 0.2;
  nval(sbac,  DTIPDIA) = 0.5;
  nval(am,    DTIPDIA) = 0.2;
  nval(am2,   DTIPDIA) = 0.2;
  nval(am3,   DTIPDIA) = 0.2;
  nval(am4,   DTIPDIA) = 0.2;
  nval(amh,   DTIPDIA) = 0.2;
  nval(amh2,  DTIPDIA) = 0.2;
  nval(ams,   DTIPDIA) = 0.4;
  nval(amhs,  DTIPDIA) = 0.4;
  nval(gca,   DTIPDIA) = 0.2;
  nval(gcb,   DTIPDIA) = 0.2;
  nval(dsgc,  DTIPDIA) = 0.2;
  nval(gcaoff,DTIPDIA) = 0.2;
  nval(gcboff,DTIPDIA) = 0.2;

  nval(xcone, DTREEDIA) = 5;	/* dendritic tree diameters */
  nval(xrod,  DTREEDIA) = 10;
  nval(ha,    DTREEDIA) = 150;
  nval(hb,    DTREEDIA) = 100;
  nval(hbat,  DTREEDIA) = 50;
  nval(rbp,   DTREEDIA) = 15;
  nval(dbp1,  DTREEDIA) = 20;
  nval(dbp2,  DTREEDIA) = 20;
  nval(dbp3,  DTREEDIA) = 20;
  nval(dbp4,  DTREEDIA) = 20;
  nval(hbp1,  DTREEDIA) = 25;
  nval(hbp2,  DTREEDIA) = 25;
  nval(a17,   DTREEDIA) = 500;
  nval(aii,   DTREEDIA) = 40;
  nval(sbac,  DTREEDIA) = 500;
  nval(am,    DTREEDIA) = 400;
  nval(am2,   DTREEDIA) = 400;
  nval(am3,   DTREEDIA) = 400;
  nval(am4,   DTREEDIA) = 400;
  nval(amh,   DTREEDIA) = 500;
  nval(amh2,  DTREEDIA) = 500;
  nval(ams,   DTREEDIA) = 20;
  nval(amhs,  DTREEDIA) = 20;
  nval(gca,   DTREEDIA) = 600;
  nval(gcb,   DTREEDIA) = 200;
  nval(dsgc , DTREEDIA) = 600;
  nval(gcaoff,DTREEDIA) = 600;
  nval(gcboff,DTREEDIA) = 200;

  nval(rbp,   AXARBT) = NBRANCHED;
  nval(dbp1,  AXARBT) = NBRANCHED;
  nval(dbp2,  AXARBT) = NBRANCHED;
  nval(dbp3,  AXARBT) = NBRANCHED;
  nval(dbp4,  AXARBT) = NBRANCHED;
  nval(hbp1,  AXARBT) = NBRANCHED;
  nval(hbp2,  AXARBT) = NBRANCHED;
  nval(ams,   AXARBT) = NBRANCHED;

  nval(xcone, AXARBZ) =  0;	/* axonal arborization layer within the retina */
  nval(xrod,  AXARBZ) =  0;
  nval(rbp,   AXARBZ) = -32;
  nval(dbp1,  AXARBZ) = -27;
  nval(dbp2,  AXARBZ) = -27;
  nval(dbp3,  AXARBZ) = -27;
  nval(dbp4,  AXARBZ) = -27;
  nval(hbp1,  AXARBZ) = -22;
  nval(hbp2,  AXARBZ) = -22;
  nval(ams,   AXARBZ) = -27;
  nval(gca,   AXARBZ) = -60;
  nval(gcb,   AXARBZ) = -60;
  nval(dsgc,  AXARBZ) = -60;
  nval(gcaoff,AXARBZ) = -60;
  nval(gcboff,AXARBZ) = -60;


  /* Connections to presynaptic cell */

  for (n=0; n<nceltypes; n++) {     /* set synaptic params to defaults */
    for (k=j=0; j<NCONNI; j++,k+=NCONNP){ 
      nval(n,k+CELPRE1)  = -1;	    /* cell type to connect to (neg, no conn) */
      nval(n,k+CONPRE1)  = 1;	    /* connection number for presyn cell */
      nval(n,k+CELCONV1) = 1;	    /* number of postsyn cells to connect to */
      nval(n,k+GROWPOST1) = 0;	    /* =1->when postsyn, grow when making synapse */
    };
  };

  nval(ha,    GROWPOST1) = 1;		/* =1 -> when postsynaptic, grow to make synapse */
  nval(hb,    GROWPOST1) = 1;		/* =2 -> when postsynaptic, grow to make synapse */
  nval(hbat,  GROWPOST1) = 1;
  nval(rbp,   GROWPOST1) = 1;
  nval(dbp1,  GROWPOST1) = 1;
  nval(dbp2,  GROWPOST1) = 1;
  nval(dbp3,  GROWPOST1) = 1;
  nval(dbp4,  GROWPOST1) = 1;
  nval(hbp1,  GROWPOST1) = 1;
  nval(hbp2,  GROWPOST1) = 1;
  nval(a17,   GROWPOST1) = 0;
  nval(aii,   GROWPOST1) = 1;
  nval(sbac,  GROWPOST1) = 0;
  nval(am,    GROWPOST1) = 1;
  nval(am2,   GROWPOST1) = 1;
  nval(am3,   GROWPOST1) = 1;
  nval(am4,   GROWPOST1) = 1;
  nval(amh,   GROWPOST1) = 1;
  nval(amh2,  GROWPOST1) = 1;
  nval(ams,   GROWPOST1) = 1;
  nval(amhs,  GROWPOST1) = 1;
  nval(gca,   GROWPOST1) = 0;
  nval(gcb,   GROWPOST1) = 0;
  nval(dsgc,  GROWPOST1) = 0;
  nval(gcaoff,GROWPOST1) = 0;
  nval(gcboff,GROWPOST1) = 0;

  /* Connections to postsynaptic cell */

  for (n=0; n<nceltypes; n++) {      /* set synaptic params to defaults */
    for (k=j=0; j<NCONNO; j++,k+=NSYNP){ 
      nval(n,k+CELPOST1) = -1;	     /* cell type to connect to (neg, no conn) */
      nval(n,k+CONPOST1) = 0;	     /* connection number */
      nval(n,k+CELDIV1)  = 2;	     /* number of postsyn cells to connect to */
      nval(n,k+GROWPRE1) = 0;	     /* grow when making conn to postsyn cell */
      nval(n,k+SYNREG1)  = 0;	     /* synaptic region in presyn dendritic tree */
      nval(n,k+SYNREGP1) = 0;	     /* synaptic region in postsyn dendritic tree */
      nval(n,k+SYNSPAC1) = 0;	     /* synaptic spacing in presyn dendritic tree */
      nval(n,k+SYNANNI1) = 0;	     /* inner dia of annulus in presyn dendr tree */
      nval(n,k+SYNANNO1) = 0;	     /* outer dia of annulus in presyn dendr tree */
      nval(n,k+SYNANPI1) = 0;	     /* inner dia of annulus in postsyn dendr tree */
      nval(n,k+SYNANPO1) = 0;	     /* outer dia of annulus in postsyn dendr tree */
      nval(n,k+SYNANG1)  = 0;	     /* angle for postsynaptic cell */
      nval(n,k+SYNRNG1)  = 0;	     /* range of angles for postsynaptic cell */
      nval(n,k+USEDYAD1) = 0;	     /* synapse is dyad using preexisting type */
      nval(n,k+DYADTYP1) = -1;	     /* type of dyad synapse to connect with */
      nval(n,k+AUTAPSE1) = 0;	     /* synapse back to presynaptic node */
      nval(n,k+SYNNUM1)  = 1;	     /* number of synapses per connection */
      nval(n,k+SENSCA1)  = 0;	     /* synaptic sensitivity to calcium */
      nval(n,k+SRRPOOL1) = 0;	     /* synaptic readily releasable pool */
      nval(n,k+SRRPOOLG1)= 0;	     /* synaptic readily releasable pool gain */
      nval(n,k+SMRRPOOL1)= 0;	     /* synaptic max readily releasable pool */
      nval(n,k+SMAXRATE1)= 0;        /* maximum sustained synaptic release rate */
      nval(n,k+SGAIN1)   = 2;	     /* synaptic gain */
      nval(n,k+SVGAIN1)  = 1;	     /* synaptic vgain */
      nval(n,k+SDURH1)   = 0;	     /* synaptic high pass time const. */
      nval(n,k+SNFILTH1) = 0;	     /* synaptic high pass nfilt. */
      nval(n,k+SHGAIN1)  = 0;	     /* synaptic high pass gain. */
      nval(n,k+SHOFFS1)  = 0;	     /* synaptic high pass offset. */
      nval(n,k+SVSIZ1)   = 1.0;	     /* synaptic vesicle size */
      nval(n,k+SCOND1)   = 220e-12;  /* synaptic conductance */
      nval(n,k+SCMUL1)   = 0; 	     /* synaptic conductance mult for region */
      nval(n,k+SCGRAD1)  = 0; 	     /* synaptic conductance gradient from soma */
      nval(n,k+SEGRAD1)  = 0; 	     /* synaptic conductance expon gradient fr soma */
      nval(n,k+STHRESH1) = -0.05;    /* synaptic release threshold */
      nval(n,k+SVNOISE1) = 1;	     /* 1 -> synaptic vesicle noise */
      nval(n,k+SCOV1)    = 1; 	     /* 1=Poisson, <1->more regular, gamma dist */
      nval(n,k+SDUR1)    = 2;	     /* synaptic event time const. */
      nval(n,k+SFALL1)   = 0;	     /* synaptic event fall time const. */
      nval(n,k+SNFILT1)  = 1;	     /* number of synaptic vesicle filters. */
      nval(n,k+STRCONC1) = 100e-6;   /* synaptic transmitter concentration */
      nval(n,k+SRESP1)   = xampa;    /* synaptic response (ampa,gaba,gj,etc.)*/
      nval(n,k+SPCA1)    = 0; 	     /* synaptic postsyn Ca perm (ampa,nmda,etc.)*/
      nval(n,k+SCAVMAX1) = 0;	     /* synaptic postsyn Ca pump vmax. */
      nval(n,k+SCAKM1)   = 0;	     /* synaptic postsyn Ca pump Km. */
      nval(n,k+SCNFILT1) = 0;	     /* number of second messenger filters. */
      nval(n,k+SCDUR1)   = 1;	     /* second mesng. time const. */
      nval(n,k+SCGAIN1)  = 1;	     /* synaptic second messenger gain */
      nval(n,k+SCOFF1)   = 1;	     /* synaptic second messenger offset */
      nval(n,k+SCNOISE1) =  0;	     /* synaptic channel noise */
      nval(n,k+SNCHAN1)  = 20;	     /* number of channels */
      nval(n,k+SUNIT1)   = 22e-12;  /* unitary conductance */
      nval(n,k+SVREV1)   = 0;	     /* synaptic reversal potential */
    };
  };

  nval(xcone, GROWPRE1) = 0;		/* =1 -> when presynaptic, grow to make synapse */
  nval(xrod,  GROWPRE1) = 0;		/* =2 -> when presynaptic, grow when artif morph */
  nval(ha,    GROWPRE1) = 0;
  nval(hb,    GROWPRE1) = 0;
  nval(hbat,  GROWPRE1) = 0;
  nval(rbp,   GROWPRE1) = 0;
  nval(dbp1,  GROWPRE1) = 0;
  nval(dbp2,  GROWPRE1) = 0;
  nval(dbp3,  GROWPRE1) = 0;
  nval(dbp4,  GROWPRE1) = 0;
  nval(hbp1,  GROWPRE1) = 0;
  nval(hbp2,  GROWPRE1) = 0;
  nval(a17,   GROWPRE1) = 0;
  nval(aii,   GROWPRE1) = 0;
  nval(sbac,  GROWPRE1) = 0;
  nval(am,    GROWPRE1) = 0;
  nval(am2,   GROWPRE1) = 0;
  nval(am3,   GROWPRE1) = 0;
  nval(am4,   GROWPRE1) = 0;
  nval(amh,   GROWPRE1) = 0;
  nval(amh2,  GROWPRE1) = 0;
  nval(ams,   GROWPRE1) = 0;
  nval(amhs,  GROWPRE1) = 0;
  nval(gca,   GROWPRE1) = 0;
  nval(gcb,   GROWPRE1) = 0;
  nval(dsgc,  GROWPRE1) = 0;
  nval(gcaoff,GROWPRE1) = 0;
  nval(gcboff,GROWPRE1) = 0;

  nval(xcone, SYNSPAC1) = 0;		/* synaptic spacing in presyn dendritic tree */
  nval(xrod,  SYNSPAC1) = 0;
  nval(ha,    SYNSPAC1) = 0;
  nval(hb,    SYNSPAC1) = 0;
  nval(hbat,  SYNSPAC1) = 0;
  nval(rbp,   SYNSPAC1) = 0;
  nval(dbp1,  SYNSPAC1) = 0;
  nval(dbp2,  SYNSPAC1) = 0;
  nval(dbp3,  SYNSPAC1) = 0;
  nval(dbp4,  SYNSPAC1) = 0;
  nval(hbp1,  SYNSPAC1) = 0;
  nval(hbp2,  SYNSPAC1) = 0;
  nval(a17,   SYNSPAC1) = 0;
  nval(aii,   SYNSPAC1) = 0;
  nval(sbac,  SYNSPAC1) = 0;
  nval(am,    SYNSPAC1) = 0;
  nval(am2,   SYNSPAC1) = 0;
  nval(am3,   SYNSPAC1) = 0;
  nval(am4,   SYNSPAC1) = 0;
  nval(amh,   SYNSPAC1) = 0;
  nval(amh2,  SYNSPAC1) = 0;
  nval(ams,   SYNSPAC1) = 0;
  nval(amhs,  SYNSPAC1) = 0;
  nval(gca,   SYNSPAC1) = 0;
  nval(gcb,   SYNSPAC1) = 0;
  nval(dsgc,  SYNSPAC1) = 0;
  nval(gcaoff,SYNSPAC1) = 0;
  nval(gcboff,SYNSPAC1) = 0;

  nval(xcone, SYNANNI1) = 0;		/* inner dia of annulus in presyn dendr tree */
  nval(xrod,  SYNANNI1) = 0;
  nval(ha,    SYNANNI1) = 0;
  nval(hb,    SYNANNI1) = 0;
  nval(hbat,  SYNANNI1) = 0;
  nval(rbp,   SYNANNI1) = 0;
  nval(dbp1,  SYNANNI1) = 0;
  nval(dbp2,  SYNANNI1) = 0;
  nval(dbp3,  SYNANNI1) = 0;
  nval(dbp4,  SYNANNI1) = 0;
  nval(hbp1,  SYNANNI1) = 0;
  nval(hbp2,  SYNANNI1) = 0;
  nval(a17,   SYNANNI1) = 0;
  nval(aii,   SYNANNI1) = 0;
  nval(sbac,  SYNANNI1) = 0;
  nval(am,    SYNANNI1) = 0;
  nval(am2,   SYNANNI1) = 0;
  nval(am3,   SYNANNI1) = 0;
  nval(am4,   SYNANNI1) = 0;
  nval(amh,   SYNANNI1) = 0;
  nval(amh2,  SYNANNI1) = 0;
  nval(ams,   SYNANNI1) = 0;
  nval(amhs,  SYNANNI1) = 0;
  nval(gca,   SYNANNI1) = 0;
  nval(gcb,   SYNANNI1) = 0;
  nval(dsgc,  SYNANNI1) = 0;
  nval(gcaoff,SYNANNI1) = 0;
  nval(gcboff,SYNANNI1) = 0;

  nval(xcone, SYNANNO1) = 0;		/* outer dia of annulus in presyn dendr tree */
  nval(xrod,  SYNANNO1) = 0;
  nval(ha,    SYNANNO1) = 0;
  nval(hb,    SYNANNO1) = 0;
  nval(hbat,  SYNANNO1) = 0;
  nval(rbp,   SYNANNO1) = 0;
  nval(dbp1,  SYNANNO1) = 0;
  nval(dbp2,  SYNANNO1) = 0;
  nval(dbp3,  SYNANNO1) = 0;
  nval(dbp4,  SYNANNO1) = 0;
  nval(hbp1,  SYNANNO1) = 0;
  nval(hbp2,  SYNANNO1) = 0;
  nval(a17,   SYNANNO1) = 0;
  nval(aii,   SYNANNO1) = 0;
  nval(sbac,  SYNANNO1) = 0;
  nval(am,    SYNANNO1) = 0;
  nval(am2,   SYNANNO1) = 0;
  nval(am3,   SYNANNO1) = 0;
  nval(am4,   SYNANNO1) = 0;
  nval(amh,   SYNANNO1) = 0;
  nval(amh2,  SYNANNO1) = 0;
  nval(ams,   SYNANNO1) = 0;
  nval(amhs,  SYNANNO1) = 0;
  nval(gca,   SYNANNO1) = 0;
  nval(gcb,   SYNANNO1) = 0;
  nval(dsgc,  SYNANNO1) = 0;
  nval(gcaoff,SYNANNO1) = 0;
  nval(gcboff,SYNANNO1) = 0;

  nval(xcone, SYNANPI1) = 0;		/* inner dia of annulus in postsyn dendr tree */
  nval(xrod,  SYNANPI1) = 0;
  nval(ha,    SYNANPI1) = 0;
  nval(hb,    SYNANPI1) = 0;
  nval(hbat,  SYNANPI1) = 0;
  nval(rbp,   SYNANPI1) = 0;
  nval(dbp1,  SYNANPI1) = 0;
  nval(dbp2,  SYNANPI1) = 0;
  nval(dbp3,  SYNANPI1) = 0;
  nval(dbp4,  SYNANPI1) = 0;
  nval(hbp1,  SYNANPI1) = 0;
  nval(hbp2,  SYNANPI1) = 0;
  nval(a17,   SYNANPI1) = 0;
  nval(aii,   SYNANPI1) = 0;
  nval(sbac,  SYNANPI1) = 0;
  nval(am,    SYNANPI1) = 0;
  nval(am2,   SYNANPI1) = 0;
  nval(am3,   SYNANPI1) = 0;
  nval(am4,   SYNANPI1) = 0;
  nval(amh,   SYNANPI1) = 0;
  nval(amh2,  SYNANPI1) = 0;
  nval(ams,   SYNANPI1) = 0;
  nval(amhs,  SYNANPI1) = 0;
  nval(gca,   SYNANPI1) = 0;
  nval(gcb,   SYNANPI1) = 0;
  nval(dsgc,  SYNANPI1) = 0;
  nval(gcaoff,SYNANPI1) = 0;
  nval(gcboff,SYNANPI1) = 0;

  nval(xcone, SYNANPO1) = 0;		/* outer dia of annulus in postsyn dendr tree */
  nval(xrod,  SYNANPO1) = 0;
  nval(ha,    SYNANPO1) = 0;
  nval(hb,    SYNANPO1) = 0;
  nval(hbat,  SYNANPO1) = 0;
  nval(rbp,   SYNANPO1) = 0;
  nval(dbp1,  SYNANPO1) = 0;
  nval(dbp2,  SYNANPO1) = 0;
  nval(dbp3,  SYNANPO1) = 0;
  nval(dbp4,  SYNANPO1) = 0;
  nval(hbp1,  SYNANPO1) = 0;
  nval(hbp2,  SYNANPO1) = 0;
  nval(a17,   SYNANPO1) = 0;
  nval(aii,   SYNANPO1) = 0;
  nval(sbac,  SYNANPO1) = 0;
  nval(am,    SYNANPO1) = 0;
  nval(am2,   SYNANPO1) = 0;
  nval(am3,   SYNANPO1) = 0;
  nval(am4,   SYNANPO1) = 0;
  nval(amh,   SYNANPO1) = 0;
  nval(amh2,  SYNANPO1) = 0;
  nval(ams,   SYNANPO1) = 0;
  nval(amhs,  SYNANPO1) = 0;
  nval(gca,   SYNANPO1) = 0;
  nval(gcb,   SYNANPO1) = 0;
  nval(dsgc,  SYNANPO1) = 0;
  nval(gcaoff,SYNANPO1) = 0;
  nval(gcboff,SYNANPO1) = 0;

  nval(xcone, SYNANG1) = 0;		/* angle for postsynaptic cell */
  nval(xrod,  SYNANG1) = 0;
  nval(ha,    SYNANG1) = 0;
  nval(hb,    SYNANG1) = 0;
  nval(hbat,  SYNANG1) = 0;
  nval(rbp,   SYNANG1) = 0;
  nval(dbp1,  SYNANG1) = 0;
  nval(dbp2,  SYNANG1) = 0;
  nval(dbp3,  SYNANG1) = 0;
  nval(dbp4,  SYNANG1) = 0;
  nval(hbp1,  SYNANG1) = 0;
  nval(hbp2,  SYNANG1) = 0;
  nval(a17,   SYNANG1) = 0;
  nval(aii,   SYNANG1) = 0;
  nval(sbac,  SYNANG1) = 0;
  nval(am,    SYNANG1) = 0;
  nval(am2,   SYNANG1) = 0;
  nval(am3,   SYNANG1) = 0;
  nval(am4,   SYNANG1) = 0;
  nval(amh,   SYNANG1) = 0;
  nval(amh2,  SYNANG1) = 0;
  nval(ams,   SYNANG1) = 0;
  nval(amhs,  SYNANG1) = 0;
  nval(gca,   SYNANG1) = 0;
  nval(gcb,   SYNANG1) = 0;
  nval(dsgc,  SYNANG1) = 0;
  nval(gcaoff,SYNANG1) = 0;
  nval(gcboff,SYNANG1) = 0;

  nval(xcone, SYNRNG1) = 0;		/* range of angles for postsynaptic cell */
  nval(xrod,  SYNRNG1) = 0;
  nval(ha,    SYNRNG1) = 0;
  nval(hb,    SYNRNG1) = 0;
  nval(hbat,  SYNRNG1) = 0;
  nval(rbp,   SYNRNG1) = 0;
  nval(dbp1,  SYNRNG1) = 0;
  nval(dbp2,  SYNRNG1) = 0;
  nval(dbp3,  SYNRNG1) = 0;
  nval(dbp4,  SYNRNG1) = 0;
  nval(hbp1,  SYNRNG1) = 0;
  nval(hbp2,  SYNRNG1) = 0;
  nval(a17,   SYNRNG1) = 0;
  nval(aii,   SYNRNG1) = 0;
  nval(sbac,  SYNRNG1) = 0;
  nval(am,    SYNRNG1) = 0;
  nval(am2,   SYNRNG1) = 0;
  nval(am3,   SYNRNG1) = 0;
  nval(am4,   SYNRNG1) = 0;
  nval(amh,   SYNRNG1) = 0;
  nval(amh2,  SYNRNG1) = 0;
  nval(ams,   SYNRNG1) = 0;
  nval(amhs,  SYNRNG1) = 0;
  nval(gca,   SYNRNG1) = 0;
  nval(gcb,   SYNRNG1) = 0;
  nval(dsgc,  SYNRNG1) = 0;
  nval(gcaoff,SYNRNG1) = 0;
  nval(gcboff,SYNRNG1) = 0;

  nval(xcone,CELPOST1) = dbp1;	  	/* cell type to connect to */
  nval(xcone,CONPOST1) = 1;	  	/* connection number */
  nval(xcone,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xcone,SYNNUM1)  = 6;	  	/* number of synapses per conn */
  nval(xcone,SRESP1)   = xmglur6;	/* synaptic response */
  nval(xcone,SGAIN1)   = 3;	  	/* synaptic gain */
  nval(xcone,SVGAIN1)  = 1;	  	/* synaptic vgain */
  nval(xcone,SCGAIN1)  = 1.6;	  	/* synaptic second messenger gain */
  nval(xcone,SCOFF1)   = 1;	  	/* synaptic second messenger offset */
  nval(xcone,SCOND1)   = 2e-10;  	/* synaptic conductance */
  nval(xcone,SVREV1)   = 0;	  	/* synaptic reversal potential */
  nval(xcone,STHRESH1) = -0.0290;  	/* synaptic release threshold */
  nval(xcone,SNCHAN1)  = 20;	  	/* number of channels */
  nval(xcone,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SCDUR1)   = 2;	  	/* second mesng. time const. */
  nval(xcone,SCNFILT1) = 2;	  	/* second mesng. number of filters. */
  nval(xcone,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(xcone,STRCONC1) = 0.0009;	/* synaptic neurotransmitter conc */

  nval(dbp1,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE1)   = 1;	  	/* presynaptic connection number */
  nval(dbp1,CELCONV1)  = 8;		/* number of presyn cells to connect to */

  nval(xcone,CELPOST2) = dbp2;	  	/* cell type to connect to */
  nval(xcone,CONPOST2) = 1;	  	/* connection number */
  nval(xcone,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xcone,SYNNUM2)  = 6;	  	/* number of synapses per conn */
  nval(xcone,SRESP2)   = xmglur6;	/* synaptic response */
  nval(xcone,SGAIN2)   = 3;	  	/* synaptic gain */
  nval(xcone,SVGAIN2)  = 1;	  	/* synaptic vgain */
  nval(xcone,SCGAIN2)  = 1.6;	  	/* synaptic second messenger gain */
  nval(xcone,SCOFF2)   = 1;	  	/* synaptic second messenger offset */
  nval(xcone,SCOND2)   = 2e-10;  	/* synaptic conductance */
  nval(xcone,SVREV2)   = 0;	  	/* synaptic reversal potential */
  nval(xcone,STHRESH2) = -0.0290;  	/* synaptic release threshold */
  nval(xcone,SNCHAN2)  = 20;	  	/* number of channels */
  nval(xcone,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SCDUR2)   = 2;	  	/* second mesng. time const. */
  nval(xcone,SCNFILT2) = 2;	  	/* second mesng. number of filters. */
  nval(xcone,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(xcone,STRCONC2) = 0.0009;	/* synaptic neurotransmitter conc */

  nval(dbp2,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(dbp2,CONPRE1)   = 2;	  	/* presynaptic connection number */
  nval(dbp2,CELCONV1)  = 8;		/* number of presyn cells to connect to */

  nval(xcone,CELPOST3) = hbp1;	  	/* cell type to connect to */
  nval(xcone,CONPOST3) = 1;	  	/* connection number */
  nval(xcone,CELDIV3)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xcone,SYNNUM3)  = 6;	  	/* number of synapses per conn */
  nval(xcone,SRESP3)   = xampa;		/* synaptic response */
  nval(xcone,SGAIN3)   = 3;	  	/* synaptic gain */
  nval(xcone,SVGAIN3)  = 1;	  	/* synaptic vgain */
  nval(xcone,SCGAIN3)  = 1;	  	/* synaptic second messenger gain */
  nval(xcone,SCOND3)   = 1e-10;  	/* synaptic conductance */
  nval(xcone,STHRESH3) = -0.0290;  	/* synaptic release threshold */
  nval(xcone,SNCHAN3)  = 20;	  	/* number of channels */
  nval(xcone,SDUR3)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SVSIZ3)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT3)   = 22e-12;  	/* unitary conductance */
  nval(xcone,SVREV3)   = 0;	  	/* synaptic reversal potential */
  nval(xcone,STRCONC3) = 0.0001;	/* synaptic neurotransmitter conc */

  nval(hbp1,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(hbp1,CONPRE1)   = 3;	  	/* presynaptic connection number */
  nval(hbp1,CELCONV1)  = 8;		/* number of presyn cells to connect to */

  nval(xcone,CELPOST4) = hbp2;	  	/* cell type to connect to */
  nval(xcone,CONPOST4) = 1;	  	/* connection number */
  nval(xcone,CELDIV4)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xcone,SYNNUM4)  = 6;	  	/* number of synapses per conn */
  nval(xcone,SRESP4)   = xampa;		/* synaptic response */
  nval(xcone,SGAIN4)   = 3;	  	/* synaptic gain */
  nval(xcone,SVGAIN4)  = 1;	  	/* synaptic vgain */
  nval(xcone,SCGAIN4)  = 1;	  	/* synaptic second messenger gain */
  nval(xcone,SCOND4)   = 1e-10;  	/* synaptic conductance */
  nval(xcone,STHRESH4) = -0.0290;  	/* synaptic release threshold */
  nval(xcone,SNCHAN4)  = 20;	  	/* number of channels */
  nval(xcone,SDUR4)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SVSIZ4)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT4)   = 22e-12;  	/* unitary conductance */
  nval(xcone,SVREV4)   = 0;	  	/* synaptic reversal potential */
  nval(xcone,STRCONC4) = 0.0001;	/* synaptic neurotransmitter conc */

  nval(hbp2,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(hbp2,CONPRE1)   = 4;	  	/* presynaptic connection number */
  nval(hbp2,CELCONV1)  = 8;		/* number of presyn cells to connect to */

  nval(xcone,CELPOST5) = ha;	  	/* cell type to connect to */
  nval(xcone,CONPOST5) = 1;	  	/* connection number */
  nval(xcone,CELDIV5)  = 4;	  	/* number of postsyn cells to conn to */
  nval(xcone,USEDYAD5) = 1;		/* synapse is dyad using preexisting type */
  nval(xcone,DYADTYP5) = dbp1;		/* type of dyad synapse to connect with */
  nval(xcone,SYNNUM5)  = 1;	  	/* number of synapses per conn */
  nval(xcone,SGAIN5)   = 2;	  	/* synaptic gain */
  nval(xcone,SCGAIN5)  = 1;	  	/* synaptic second messenger gain */
  nval(xcone,SCOND5)   = 220e-12;  	/* synaptic conductance */
  nval(xcone,STHRESH5) = -0.050;  	/* synaptic release threshold */
  nval(xcone,SRESP5)   = xampa;		/* synaptic response */
  nval(xcone,SNCHAN5)  = 20;	  	/* number of channels */
  nval(xcone,SDUR5)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SCDUR5)   = 1;	  	/* second mesng. time const. */
  nval(xcone,SVSIZ5)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT5)   = 22e-12;  	/* unitary conductance */
  nval(xcone,SVREV5)   = 0;	  	/* synaptic reversal potential */

  nval(ha,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(ha,CONPRE1)   = 5;	  	/* presynaptic connection number */
  nval(ha,CELCONV1)  = 200;		/* number of presyn cells to connect to */

  nval(xcone,CELPOST6) = hb;	  	/* cell type to connect to */
  nval(xcone,CONPOST6) = 1;	  	/* connection number */
  nval(xcone,CELDIV6)  = 4;	  	/* number of postsyn cells to conn to */
  nval(xcone,USEDYAD6) = 1;		/* synapse is dyad using preexisting type */
  nval(xcone,DYADTYP6) = dbp1;		/* type of dyad synapse to connect with */
  nval(xcone,SYNNUM6)  = 1;	  	/* number of synapses per conn */
  nval(xcone,SRESP6)   = xampa;		/* synaptic response */
  nval(xcone,SGAIN6)   = 2;	  	/* synaptic gain */
  nval(xcone,SCGAIN6)  = 1;	  	/* synaptic second messenger gain */
  nval(xcone,SCOND6)   = 220e-12;  	/* synaptic conductance */
  nval(xcone,STHRESH6) = -0.050;  	/* synaptic release threshold */
  nval(xcone,SNCHAN6)  = 20;	  	/* number of channels */
  nval(xcone,SDUR6)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SCDUR6)   = 1;	  	/* second mesng. time const. */
  nval(xcone,SVSIZ6)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT6)   = 22e-12;  	/* unitary conductance */
  nval(xcone,SVREV6)   = 0;	  	/* synaptic reversal potential */

  nval(hb,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(hb,CONPRE1)   = 6;	  	/* presynaptic connection number */
  nval(hb,CELCONV1)  = 100;		/* number of presyn cells to connect to */

  nval(xcone,CELPOST7) = xcone;  	/* cell type to connect to */
  nval(xcone,CONPOST7) = 3;	  	/* connection number */
  nval(xcone,CELDIV7)  = 8;		/* number of postsyn cells to conn to */
  nval(xcone,GROWPRE7) = 0;		/* grow when making conn to postsyn cell */
  nval(xcone,SYNNUM7)  = 4;		/* number of synapses per conn */
  nval(xcone,SRESP7)   = xgapj;		/* synaptic response */
  nval(xcone,SGAIN7)   = 2;	  	/* synaptic gain */
  nval(xcone,SCGAIN7)  = 1;	  	/* synaptic second messenger gain */
  nval(xcone,SCOND7)   = 500e-12;  	/* synaptic conductance */
  nval(xcone,STHRESH7) = -0.050;  	/* synaptic release threshold */
  nval(xcone,SNCHAN7)  = 20;	  	/* number of channels */
  nval(xcone,SDUR7)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SCDUR7)   = 1;	  	/* second mesng. time const. */
  nval(xcone,SVSIZ7)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT7)   = 22e-12;  	/* unitary conductance */
  nval(xcone,SVREV7)   = 0;	  	/* synaptic reversal potential */

  nval(xcone,CELPRE3)   = xcone;	/* presynaptic cell type */
  nval(xcone,CONPRE3)   = 7;	  	/* presynaptic connection number */
  nval(xcone,GROWPOST3) = 0;	  	/* grow when making conn from presyn cell */
  nval(xcone,CELCONV3)  = 8;		/* number of postsyn cells to connect to */

  nval(xcone,CELPOST8) = dbp3;	  	/* cell type to connect to */
  nval(xcone,CONPOST8) = 1;	  	/* connection number */
  nval(xcone,CELDIV8)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xcone,SYNNUM8)  = 6;	  	/* number of synapses per conn */
  nval(xcone,SRESP8)   = xmglur6;	/* synaptic response */
  nval(xcone,SGAIN8)   = 3;	  	/* synaptic gain */
  nval(xcone,SVGAIN8)  = 1;	  	/* synaptic vgain */
  nval(xcone,SCGAIN8)  = 1.6;	  	/* synaptic second messenger gain */
  nval(xcone,SCOFF8)   = 1;	  	/* synaptic second messenger offset */
  nval(xcone,SCOND8)   = 2e-10;  	/* synaptic conductance */
  nval(xcone,SVREV8)   = 0;	  	/* synaptic reversal potential */
  nval(xcone,STHRESH8) = -0.0290;  	/* synaptic release threshold */
  nval(xcone,SNCHAN8)  = 20;	  	/* number of channels */
  nval(xcone,SDUR8)    = 2;	  	/* synaptic event time const. */
  nval(xcone,SCDUR8)   = 2;	  	/* second mesng. time const. */
  nval(xcone,SCNFILT8) = 2;	  	/* second mesng. number of filters. */
  nval(xcone,SVSIZ8)   = 1.0;	  	/* synaptic vesicle size */
  nval(xcone,SUNIT8)   = 22e-12;  	/* unitary conductance */
  nval(xcone,STRCONC8) = 0.0009;	/* synaptic neurotransmitter conc */

  nval(dbp3,CELPRE1)   = xcone;	  	/* presynaptic cell type */
  nval(dbp3,CONPRE1)   = 8;	  	/* presynaptic connection number */
  nval(dbp3,CELCONV1)  = 8;		/* number of presyn cells to connect to */

  nval(xrod,CELPOST1) = rbp;	  	/* cell type to connect to */
  nval(xrod,CONPOST1) = 1;	  	/* connection number */
  nval(xrod,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xrod,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(xrod,SRESP1)   = xmglur6;	/* synaptic response */
  nval(xrod,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(xrod,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(xrod,STHRESH1) = -0.050;  	/* synaptic release threshold */
  nval(xrod,SNCHAN1)  = 20;	  	/* number of channels */
  nval(xrod,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(xrod,SNFILT1)  = 1;	  	/* synaptic vesicle nfilt */
  nval(xrod,SCDUR1)   = 10;	  	/* second mesng. time const. */
  nval(xrod,SCGAIN1)  = 1.025;	  	/* synaptic second messenger gain */
  nval(xrod,SCNFILT1) = 2;	  	/* second mesng. nfilt  */
  nval(xrod,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(xrod,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(xrod,SVREV1)   = 0;	  	/* synaptic reversal potential */
  nval(xrod,STRCONC1) = 0.0001;		/* synaptic neurotransmitter conc */

  nval(rbp,CELPRE1)  = xrod;	  	/* presynaptic cell type */
  nval(rbp,CONPRE1)  = 1;	  	/* presynaptic connection number */
  nval(rbp,CELCONV1) = 25;		/* number of presyn cells to connect to */

  nval(xrod,CELPOST2) = hbat;	  	/* cell type to connect to */
  nval(xrod,CONPOST2) = 1;	  	/* connection number */
  nval(xrod,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(xrod,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(xrod,SRESP2)   = xampa;		/* synaptic response */
  nval(xrod,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(xrod,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(xrod,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(xrod,STHRESH2) = -0.050;  	/* synaptic release threshold */
  nval(xrod,SNCHAN2)  = 20;	  	/* number of channels */
  nval(xrod,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(xrod,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(xrod,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(xrod,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(xrod,SVREV2)   = 0;	  	/* synaptic reversal potential */

  nval(hbat,CELPRE1)  = xrod;	  	/* presynaptic cell type */
  nval(hbat,CONPRE1)  = 2;	  	/* presynaptic connection number */

  nval(hbat,CELPOST1) = xrod;	  	/* cell type to connect to */
  nval(hbat,CONPOST1) = 1;	  	/* connection number */
  nval(hbat,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(hbat,USEDYAD1) = 0;		/* synapse is dyad using preexisting type */
  nval(hbat,DYADTYP1) = rbp;		/* type of dyad synapse to connect with */
  nval(hbat,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(hbat,SRESP1)   = xgaba;		/* synaptic response */
  nval(hbat,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(hbat,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(hbat,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(hbat,STHRESH1) = -0.050;  	/* synaptic release threshold */
  nval(hbat,SNCHAN1)  = 20;	  	/* number of channels */
  nval(hbat,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(hbat,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(hbat,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(hbat,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(hbat,SVREV1)   = 0;	  	/* synaptic reversal potential */

  nval(xrod,CELPRE1)  = hbat;	  	/* presynaptic cell type */
  nval(xrod,CONPRE1)  = 1;	  	/* presynaptic connection number */

  nval(ha,CELPOST1)  = xcone;	  	/* cell type to connect to */
  nval(ha,CONPOST1)  = 1;	  	/* connection number */
  nval(ha,CELDIV1)   = 200;	  	/* number of postsyn cells to conn to */
  nval(ha,SYNNUM1)   = 1;	  	/* number of synapses per conn */
  nval(ha,SRESP1)    = xgaba;		/* synaptic response */
  nval(ha,SGAIN1)    = 2;	  	/* synaptic gain */
  nval(ha,SCGAIN1)   = 1;	  	/* synaptic second messenger gain */
  nval(ha,SCOND1)    = 220e-12;  	/* synaptic conductance */
  nval(ha,STHRESH1)  = -0.050;  		/* synaptic release threshold */
  nval(ha,SNCHAN1)   = 20;	  	/* number of channels */
  nval(ha,SDUR1)     = 2;	  	/* synaptic event time const. */
  nval(ha,SCDUR1)    = 1;	  	/* second mesng. time const. */
  nval(ha,SVSIZ1)    = 1.0;	  	/* synaptic vesicle size */
  nval(ha,SUNIT1)    = 22e-12;  	/* unitary conductance */
  nval(ha,SVREV1)    = -0.06;	  	/* synaptic reversal potential */

  nval(xcone,CELPRE1)  = ha;	  	/* presynaptic cell type */
  nval(xcone,CONPRE1)  = 1;	  	/* presynaptic connection number */

  nval(ha,CELPOST2)  = dbp1;	  	/* cell type to connect to */
  nval(ha,CONPOST2)  = 2;	  	/* connection number */
  nval(ha,CELDIV2)   = 70;	  	/* number of postsyn cells to conn to */
  nval(ha,SYNNUM2)   = 1;	  	/* number of synapses per conn */
  nval(ha,SRESP2)    = xgaba;		/* synaptic response */
  nval(ha,SGAIN2)    = 2;	  	/* synaptic gain */
  nval(ha,SCGAIN2)   = 1;	  	/* synaptic second messenger gain */
  nval(ha,SCOND2)    = 220e-12;  	/* synaptic conductance */
  nval(ha,STHRESH2)  = -0.050;  	/* synaptic release threshold */
  nval(ha,SNCHAN2)   = 20;	  	/* number of channels */
  nval(ha,SDUR2)     = 2;	  	/* synaptic event time const. */
  nval(ha,SCDUR2)    = 1;	  	/* second mesng. time const. */
  nval(ha,SVSIZ2)    = 1.0;	  	/* synaptic vesicle size */
  nval(ha,SUNIT2)    = 22e-12;  	/* unitary conductance */
  nval(ha,SVREV2)    = 0;	  	/* synaptic reversal potential */

  nval(dbp1,CELPRE2)  = ha;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE2)  = 2;	  	/* presynaptic connection number */
  nval(dbp1,CELCONV2) = 3;	  	/* number of presyn cells to connect to */

  nval(ha,CELPOST3)  = hbp1;	  	/* cell type to connect to */
  nval(ha,CONPOST3)  = 2;	  	/* connection number */
  nval(ha,CELDIV3)   = 2;	  	/* number of postsyn cells to conn to */
  nval(ha,SYNNUM3)   = 1;	  	/* number of synapses per conn */
  nval(ha,SRESP3)    = xgaba;		/* synaptic response */
  nval(ha,SGAIN3)    = 2;	  	/* synaptic gain */
  nval(ha,SCGAIN3)   = 1;	  	/* synaptic second messenger gain */
  nval(ha,SCOND3)    = 220e-12;  	/* synaptic conductance */
  nval(ha,STHRESH3)  = -0.050;  	/* synaptic release threshold */
  nval(ha,SNCHAN3)   = 20;	  	/* number of channels */
  nval(ha,SDUR3)     = 2;	  	/* synaptic event time const. */
  nval(ha,SCDUR3)    = 1;	  	/* second mesng. time const. */
  nval(ha,SVSIZ3)    = 1.0;	  	/* synaptic vesicle size */
  nval(ha,SUNIT3)    = 22e-12;  	/* unitary conductance */
  nval(ha,SVREV3)    = -0.06;	  	/* synaptic reversal potential */

  nval(hbp1,CELPRE2)  = ha;	  	/* presynaptic cell type */
  nval(hbp1,CONPRE2)  = 3;	  	/* presynaptic connection number */

  nval(ha,CELPOST4)  = ha;	  	/* cell type to connect to */
  nval(ha,CONPOST4)  = 2;	  	/* connection number */
  nval(ha,CELDIV4)   = 6;	  	/* number of postsyn cells to conn to */
  nval(ha,SYNNUM4)   = 1;	  	/* number of synapses per conn */
  nval(ha,SRESP4)    = xgapj;		/* synaptic response */
  nval(ha,SGAIN4)    = 2;	  	/* synaptic gain */
  nval(ha,SCGAIN4)   = 1;	  	/* synaptic second messenger gain */
  nval(ha,SCOND4)    = 220e-12;  	/* synaptic conductance */
  nval(ha,STHRESH4)  = -0.050;  	/* synaptic release threshold */
  nval(ha,SNCHAN4)   = 20;	  	/* number of channels */
  nval(ha,SDUR4)     = 2;	  	/* synaptic event time const. */
  nval(ha,SCDUR4)    = 1;	  	/* second mesng. time const. */
  nval(ha,SVSIZ4)    = 1.0;	  	/* synaptic vesicle size */
  nval(ha,SUNIT4)    = 22e-12;  	/* unitary conductance */
  nval(ha,SVREV4)    = 0;	  	/* synaptic reversal potential */

  nval(ha,CELPRE2)  = ha;	  	/* presynaptic cell type */
  nval(ha,CONPRE2)  = 4;	  	/* presynaptic connection number */
  nval(ha,CELCONV2) = 6;	  	/* number of presyn cells to connect to */

  nval(hb,CELPOST1) = xcone;	  	/* cell type to connect to */
  nval(hb,CONPOST1) = 2;	  	/* connection number */
  nval(hb,CELDIV1)  = 100;	  	/* number of postsyn cells to conn to */
  nval(hb,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(hb,SRESP1)   = xgaba;		/* synaptic response */
  nval(hb,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(hb,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(hb,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(hb,STHRESH1) = -0.050;  		/* synaptic release threshold */
  nval(hb,SNCHAN1)  = 20;	  	/* number of channels */
  nval(hb,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(hb,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(hb,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(hb,SUNIT1)   = 22e-12;  		/* unitary conductance */
  nval(hb,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  nval(xcone,CELPRE2)  = hb;	  	/* presynaptic cell type */
  nval(xcone,CONPRE2)  = 1;	  	/* presynaptic connection number */

  nval(hb,CELPOST2) = dbp1;	  	/* cell type to connect to */
  nval(hb,CONPOST2) = 3;	  	/* connection number */
  nval(hb,CELDIV2)  = 50;	  	/* number of postsyn cells to conn to */
  nval(hb,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(hb,SRESP2)   = xgaba;		/* synaptic response */
  nval(hb,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(hb,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(hb,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(hb,STHRESH2) = -0.050;  		/* synaptic release threshold */
  nval(hb,SNCHAN2)  = 20;	  	/* number of channels */
  nval(hb,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(hb,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(hb,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(hb,SUNIT2)   = 22e-12;  		/* unitary conductance */
  nval(hb,SVREV2)   = 0;	  	/* synaptic reversal potential */

  nval(dbp1,CELPRE3)  = hb;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE3)  = 2;	  	/* presynaptic connection number */
  nval(dbp1,CELCONV3) = 3;	  	/* number of presyn cells to connect to */

  nval(hb,CELPOST3)  = hbp1;	  	/* cell type to connect to */
  nval(hb,CONPOST3)  = 3;	  	/* connection number */
  nval(hb,CELDIV3)   = 2;	  	/* number of postsyn cells to conn to */
  nval(hb,SYNNUM3)   = 1;	  	/* number of synapses per conn */
  nval(hb,SRESP3)    = xgaba;		/* synaptic response */
  nval(hb,SGAIN3)    = 2;	  	/* synaptic gain */
  nval(hb,SCGAIN3)   = 1;	  	/* synaptic second messenger gain */
  nval(hb,SCOND3)    = 220e-12;  	/* synaptic conductance */
  nval(hb,STHRESH3)  = -0.050; 		/* synaptic release threshold */
  nval(hb,SNCHAN3)   = 20;	  	/* number of channels */
  nval(hb,SDUR3)     = 2;	  	/* synaptic event time const. */
  nval(hb,SCDUR3)    = 1;	  	/* second mesng. time const. */
  nval(hb,SVSIZ3)    = 1.0;	  	/* synaptic vesicle size */
  nval(hb,SUNIT3)    = 22e-12; 		/* unitary conductance */
  nval(hb,SVREV3)    = -0.06;	  	/* synaptic reversal potential */

  nval(hbp1,CELPRE3)  = hb;	  	/* presynaptic cell type */
  nval(hbp1,CONPRE3)  = 3;	  	/* presynaptic connection number */

  nval(hb,CELPOST4)  = hb;	  	/* cell type to connect to */
  nval(hb,CONPOST4)  = 2;	  	/* connection number */
  nval(hb,CELDIV4)   = 6;	  	/* number of postsyn cells to conn to */
  nval(hb,SYNNUM4)   = 1;	  	/* number of synapses per conn */
  nval(hb,SRESP4)    = xgapj;		/* synaptic response */
  nval(hb,SGAIN4)    = 2;	  	/* synaptic gain */
  nval(hb,SCGAIN4)   = 1;	  	/* synaptic second messenger gain */
  nval(hb,SCOND4)    = 22e-12;  	/* synaptic conductance */
  nval(hb,STHRESH4)  = -0.050; 		/* synaptic release threshold */
  nval(hb,SNCHAN4)   = 20;	  	/* number of channels */
  nval(hb,SDUR4)     = 2;	  	/* synaptic event time const. */
  nval(hb,SCDUR4)    = 1;	  	/* second mesng. time const. */
  nval(hb,SVSIZ4)    = 1.0;	  	/* synaptic vesicle size */
  nval(hb,SUNIT4)    = 22e-12; 		/* unitary conductance */
  nval(hb,SVREV4)    = 0;	  	/* synaptic reversal potential */

  nval(hb,CELPRE2)  = hb;	  	/* presynaptic cell type */
  nval(hb,CONPRE2)  = 4;	  	/* presynaptic connection number */
  nval(hb,CELCONV2) = 6;	  	/* number of presyn cells to connect to */

  nval(rbp,CELPOST1) = aii;	  	/* cell type to connect to */
  nval(rbp,CONPOST1) = 1;	  	/* connection number */
  nval(rbp,CELDIV1)  = 5;	  	/* number of postsyn cells to conn to */
  nval(rbp,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(rbp,SRRPOOL1) = 100;	  	/* synaptic readily releaseable pool */
  nval(rbp,SMRRPOOL1)= 100;	  	/* synaptic max readily releaseable pool */
  nval(rbp,SMAXRATE1)= 10;	  	/* synaptic maximum sustained release rate */
  nval(rbp,SRESP1)   = xampa;		/* synaptic response */
  nval(rbp,SGAIN1)   = 1.8;	  	/* synaptic gain */
  nval(rbp,SDURH1)   = 50;	  	/* synaptic high pass time const. */
  nval(rbp,SHGAIN1)  = 0.7;	  	/* synaptic high pass gain */ 
  nval(rbp,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(rbp,SCOND1)   = 620e-12;  	/* synaptic conductance */
  nval(rbp,STHRESH1) = -0.045;  	/* synaptic release threshold */
  nval(rbp,SNCHAN1)  = 20;	  	/* number of channels */
  nval(rbp,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(rbp,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(rbp,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(rbp,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(rbp,SVREV1)   = 0;	  	/* synaptic reversal potential */
 
  nval(aii,CELPRE1)  = rbp;	  	/* presynaptic cell type */
  nval(aii,CONPRE1)  = 1;	  	/* presynaptic connection number */
  nval(aii,CELCONV1) = 25;		/* number of presyn cells to connect to */

  nval(rbp,CELPOST2) = a17;	  	/* cell type to connect to */
  nval(rbp,CONPOST2) = 1;	  	/* connection number */
  nval(rbp,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(rbp,USEDYAD2) = 1;		/* synapse is dyad using preexisting type */
  nval(rbp,DYADTYP2) = rbp;		/* type of dyad synapse to connect with */
  nval(rbp,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(rbp,SRESP2)   = xampa;		/* synaptic response */
  nval(rbp,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(rbp,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(rbp,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(rbp,STHRESH2) = -0.050;  	/* synaptic release threshold */
  nval(rbp,SNCHAN2)  = 20;	  	/* number of channels */
  nval(rbp,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(rbp,SCDUR2)   = 2;	  	/* second mesng. time const. */
  nval(rbp,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(rbp,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(rbp,SVREV2)   = 0;	  	/* synaptic reversal potential */

  nval(a17,CELPRE1) = rbp;	  	/* presynaptic cell type */
  nval(a17,CONPRE1) = 2;	  	/* presynaptic connection number */

  nval(a17,CELPOST1) = rbp;	  	/* cell type to connect to */
  nval(a17,CONPOST1) = 2;	  	/* connection number */
  nval(a17,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(a17,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(a17,SRESP1)   = xgaba;		/* synaptic response */
  nval(a17,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(a17,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(a17,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(a17,STHRESH1) = -0.050;  	/* synaptic release threshold */
  nval(a17,SNCHAN1)  = 20;	  	/* number of channels */
  nval(a17,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(a17,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(a17,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(a17,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(a17,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  nval(rbp,CELPRE2)  = a17;	  	/* presynaptic cell type */
  nval(rbp,CONPRE2)  = 1;	  	/* presynaptic connection number */

  nval(dbp1,CELPOST1)  = gca;	  	/* cell type to connect to */
  nval(dbp1,CONPOST1) = 1;	  	/* connection number */
  nval(dbp1,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SRESP1)   = xampa5;		/* synaptic response */
  nval(dbp1,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(dbp1,SDURH1)   = 10;	  	/* synaptic high pass time const. */
  nval(dbp1,SNFILTH1) = 1;	  	/* high pass number of filters */
  nval(dbp1,SHGAIN1)  = 0.9;	  	/* high pass gain */
  nval(dbp1,SHOFFS1)  = -0.001;	        /* synaptic high pass offset. */
  nval(dbp1,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND1)   = 2e-10;  		/* synaptic conductance */
  nval(dbp1,STHRESH1) = -0.042;  	/* synaptic release threshold */
  nval(dbp1,SNCHAN1)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(dbp1,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(dbp1,SVREV1)   = 0;	  	/* synaptic reversal potential */
  nval(dbp1,STRCONC1) = 0.001;		/* synaptic neurotransmitter conc */

  nval(gca,CELPRE1)  = dbp1;	  	/* presynaptic cell type */
  nval(gca,CONPRE1)  = 1;	  	/* presynaptic connection number */
  nval(gca,CELCONV1) = 1e3;		/* number of presyn cells to connect to */

  nval(dbp1,CELPOST2) = dsgc;	  	/* cell type to connect to */
  nval(dbp1,CONPOST2) = 1;	  	/* connection number */
  nval(dbp1,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SGAIN2)   = 3;	  	/* synaptic gain */
  nval(dbp1,SDURH2)   = 40;	  	/* high pass time constant */
  nval(dbp1,SNFILTH2) = 1;	  	/* high pass number of filters */
  nval(dbp1,SHGAIN2)  = .7;	  	/* high pass gain */
  nval(dbp1,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(dbp1,STHRESH2) = -0.045;  	/* synaptic release threshold */
  nval(dbp1,SNCHAN2)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR2)    = 5;	  	/* synaptic event time const. */
  nval(dbp1,SRESP2)   = xampa5;		/* synaptic response */
  nval(dbp1,STRCONC2) = 200e-6;		/* synaptic neurotransmitter conc */
  nval(dbp1,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(dbp1,SVREV2)   = 0;	  	/* synaptic reversal potential */

  nval(dsgc,CELPRE1) = dbp1;	  	/* presynaptic cell type */
  nval(dsgc,CONPRE1) = 2;	  	/* presynaptic connection number */
  nval(dsgc,CELCONV1)= 1e3;		/* number of presyn cells to connect to */

  nval(dbp1,CELPOST3) = sbac;	  	/* cell type to connect to */
  nval(dbp1,CONPOST3) = 1;	  	/* connection number */
  nval(dbp1,CELDIV3)  = 6;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM3)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SRESP3)   = xampa5;		/* synaptic response */
  nval(dbp1,SGAIN3)   = 3;	  	/* synaptic gain */
  nval(dbp1,SCGAIN3)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND3)   = 220e-12;  	/* synaptic conductance */
  nval(dbp1,STHRESH3) = -0.045;  	/* synaptic release threshold */
  nval(dbp1,SNCHAN3)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR3)    = 5;	  	/* synaptic event time const. */
  nval(dbp1,SCDUR3)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ3)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT3)   = 22e-12;  	/* unitary conductance */
  nval(dbp1,SVREV3)   = 0;	  	/* synaptic reversal potential */

  nval(sbac,CELPRE1)  = dbp1;	  	/* presynaptic cell type */
  nval(sbac,CONPRE1)  = 3;	  	/* presynaptic connection number */
  nval(sbac,CELCONV1) = 400;		/* number of presyn cells to connect to */

  nval(dbp1,CELPOST4) = am;	  	/* cell type to connect to */
  nval(dbp1,CONPOST4) = 1;	  	/* connection number */
  nval(dbp1,CELDIV4)  = 10;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM4)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SRESP4)   = xampa;		/* synaptic response */
  nval(dbp1,SGAIN4)   = 2;	  	/* synaptic gain */
  nval(dbp1,SCGAIN4)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND4)   = 220e-12;  	/* synaptic conductance */
  nval(dbp1,STHRESH4) = -0.050;  	/* synaptic release threshold */
  nval(dbp1,SNCHAN4)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR4)    = 2;	  	/* synaptic event time const. */
  nval(dbp1,SCDUR4)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ4)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT4)   = 22e-12;  	/* unitary conductance */
  nval(dbp1,SVREV4)   = 0;	  	/* synaptic reversal potential */

  nval(am,CELPRE1)   = dbp1;	  	/* presynaptic cell type */
  nval(am,CONPRE1)   = 4;	  	/* presynaptic connection number */
  nval(am,CELCONV1) = 200;		/* number of presyn cells to connect to */

  nval(dbp1,CELPOST5) = dbp1;	  	/* cell type to connect to */
  nval(dbp1,CONPOST5) = 6;	  	/* connection number */
  nval(dbp1,CELDIV5)  = 8;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM5)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SRESP5)   = xgapj;		/* synaptic response */
  nval(dbp1,SGAIN5)   = 2;	  	/* synaptic gain */
  nval(dbp1,SCGAIN5)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND5)   = 8e-10;  		/* synaptic conductance */
  nval(dbp1,STHRESH5) = -0.050;  	/* synaptic release threshold */
  nval(dbp1,SNCHAN5)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR5)    = 2;	  	/* synaptic event time const. */
  nval(dbp1,SCDUR5)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ5)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT5)   = 22e-12;  	/* unitary conductance */
  nval(dbp1,SVREV5)   = 0;	  	/* synaptic reversal potential */

  nval(dbp1,CELPRE6)  = dbp1;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE6)  = 5;	  	/* presynaptic connection number */
  nval(dbp1,CELCONV6) = 8;		/* number of presyn cells to connect to */

  nval(dbp1,CELPOST6) = aii;	  	/* cell type to connect to */
  nval(dbp1,CONPOST6) = 3;	  	/* connection number */
  nval(dbp1,CELDIV6)  = 3;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM6)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SRESP6)   = xgapj;		/* synaptic response */
  nval(dbp1,SGAIN6)   = 2;	  	/* synaptic gain */
  nval(dbp1,SCGAIN6)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND6)   = 220e-12;  	/* synaptic conductance */
  nval(dbp1,STHRESH6) = -0.050;  		/* synaptic release threshold */
  nval(dbp1,SNCHAN6)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR6)    = 2;	  	/* synaptic event time const. */
  nval(dbp1,SNFILT6)  = 1;	  	/* synaptic event time const. */
  nval(dbp1,SCDUR6)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ6)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT6)   = 22e-12; 		/* unitary conductance */
  nval(dbp1,SVREV6)   = 0;	  	/* synaptic reversal potential */

  nval(aii,CELPRE3)  = dbp1;	  	/* presynaptic cell type */
  nval(aii,CONPRE3)  = 6;	  	/* presynaptic connection number */
  nval(aii,CELCONV3) = 20;	  	/* number of presyn cells to connect to */
  nval(aii,GROWPOST3)= 1;		/* grow when making conn from presyn cell */

  nval(dbp1,CELPOST7) = ams;	  	/* cell type to connect to */
  nval(dbp1,CONPOST7) = 1;	  	/* connection number */
  nval(dbp1,CELDIV7)  = 4;	  	/* number of postsyn cells to conn to */
  nval(dbp1,SYNNUM7)  = 1;	  	/* number of synapses per conn */
  nval(dbp1,SRESP7)   = xampa;		/* synaptic response */
  nval(dbp1,SGAIN7)   = 2;	  	/* synaptic gain */
  nval(dbp1,SCGAIN7)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp1,SCOND7)   = 220e-12;  	/* synaptic conductance */
  nval(dbp1,STHRESH7) = -0.050; 	/* synaptic release threshold */
  nval(dbp1,SNCHAN7)  = 20;	  	/* number of channels */
  nval(dbp1,SDUR7)    = 2;	  	/* synaptic event time const. */
  nval(dbp1,SNFILT7)  = 1;	  	/* synaptic event time const. */
  nval(dbp1,SCDUR7)   = 1;	  	/* second mesng. time const. */
  nval(dbp1,SVSIZ7)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp1,SUNIT7)   = 22e-12; 		/* unitary conductance */
  nval(dbp1,SVREV7)   = 0;	  	/* synaptic reversal potential */

  nval(ams,CELPRE1)  = dbp1;	  	/* presynaptic cell type */
  nval(ams,CONPRE1)  = 7;	  	/* presynaptic connection number */
  nval(ams,CELCONV1) = 50;	  	/* number of presyn cells to connect to */
  nval(ams,GROWPOST1)= 1;		/* grow when making conn from presyn cell */

  nval(dbp2,CELPOST1) = gcb;	  	/* cell type to connect to */
  nval(dbp2,CONPOST1) = 1;	  	/* connection number */
  nval(dbp2,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(dbp2,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(dbp2,SRESP1)   = xampa5;		/* synaptic response */
  nval(dbp2,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(dbp2,SDURH1)   = 10;	  	/* synaptic high pass time const. */
  nval(dbp2,SNFILTH1) = 1;	  	/* high pass number of filters */
  nval(dbp2,SHGAIN1)  = 0.9;	  	/* high pass gain */
  nval(dbp2,SHOFFS1)  = -0.001;	        /* synaptic high pass offset. */
  nval(dbp2,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(dbp2,SCOND1)   = 2e-10;  		/* synaptic conductance */
  nval(dbp2,STHRESH1) = -0.042;  	/* synaptic release threshold */
  nval(dbp2,SNCHAN1)  = 20;	  	/* number of channels */
  nval(dbp2,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(dbp2,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(dbp2,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(dbp2,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(dbp2,SVREV1)   = 0;	  	/* synaptic reversal potential */
  nval(dbp2,STRCONC1) = 0.001;		/* synaptic neurotransmitter conc */

  nval(gcb,CELPRE1)  = dbp2;	  	/* presynaptic cell type */
  nval(gcb,CONPRE1)  = 1;	  	/* presynaptic connection number */
  nval(gcb,CELCONV1) = 1e3;		/* number of presyn cells to connect to */

  nval(hbp1,CELPOST1) = gcaoff;	  	/* cell type to connect to */
  nval(hbp1,CONPOST1) = 1;	  	/* connection number */
  nval(hbp1,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(hbp1,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(hbp1,SRESP1)   = xampa5;		/* synaptic response */
  nval(hbp1,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(hbp1,SDURH1)   = 5;	  	/* synaptic high pass time const. */
  nval(hbp1,SNFILTH1) = 1;	  	/* high pass number of filters */
  nval(hbp1,SHGAIN1)  = 14;	  	/* high pass gain */
  nval(hbp1,SHOFFS1)  = -0.0045;	        /* synaptic high pass offset. */
  nval(hbp1,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(hbp1,SCOND1)   = 2e-10;		/* synaptic conductance */
  nval(hbp1,STHRESH1) = -0.045;  	/* synaptic release threshold */
  nval(hbp1,SNCHAN1)  = 20;	  	/* number of channels */
  nval(hbp1,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(hbp1,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(hbp1,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(hbp1,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(hbp1,SVREV1)   = 0;	  	/* synaptic reversal potential */
  nval(hbp1,STRCONC1) = 0.001;		/* synaptic neurotransmitter conc */
 
  nval(gcaoff,CELPRE1) = hbp1;	  	/* presynaptic cell type */
  nval(gcaoff,CONPRE1) = 1;	  	/* presynaptic connection number */
  nval(gcaoff,CELCONV1) = 1e3;		/* number of presyn cells to connect to */

  nval(hbp1,CELPOST2) = sbac;	  	/* cell type to connect to */
  nval(hbp1,CONPOST2) = 1;	  	/* connection number */
  nval(hbp1,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(hbp1,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(hbp1,SRESP2)   = xampa;		/* synaptic response */
  nval(hbp1,SGAIN2)   = 3;	  	/* synaptic gain */
  nval(hbp1,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(hbp1,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(hbp1,STHRESH2) = -0.030;  	/* synaptic release threshold */
  nval(hbp1,SNCHAN2)  = 5;	  	/* number of channels */
  nval(hbp1,SDUR2)    = .2;	  	/* synaptic event time const. */
  nval(hbp1,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(hbp1,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(hbp1,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(hbp1,SVREV2)   = 0;	  	/* synaptic reversal potential */

  nval(sbac,CELPRE2)  = hbp1;	  	/* presynaptic cell type */
  nval(sbac,CONPRE2)  = 2;	  	/* presynaptic connection number */

  nval(hbp1,CELPOST3) = amh;	  	/* cell type to connect to */
  nval(hbp1,CONPOST3) = 1;	  	/* connection number */
  nval(hbp1,CELDIV3)  = 2;	  	/* number of postsyn cells to conn to */
  nval(hbp1,SYNNUM3)  = 1;	  	/* number of synapses per conn */
  nval(hbp1,SRESP3)   = xampa;		/* synaptic response */
  nval(hbp1,SGAIN3)   = 2;	  	/* synaptic gain */
  nval(hbp1,SCGAIN3)  = 1;	  	/* synaptic second messenger gain */
  nval(hbp1,SCOND3)   = 220e-12;  	/* synaptic conductance */
  nval(hbp1,STHRESH3) = -0.050;  	/* synaptic release threshold */
  nval(hbp1,SNCHAN3)  = 20;	  	/* number of channels */
  nval(hbp1,SDUR3)    = 2;	  	/* synaptic event time const. */
  nval(hbp1,SCDUR3)   = 1;	  	/* second mesng. time const. */
  nval(hbp1,SVSIZ3)   = 1.0;	  	/* synaptic vesicle size */
  nval(hbp1,SUNIT3)   = 22e-12;  	/* unitary conductance */
  nval(hbp1,SVREV3)   = 0;	  	/* synaptic reversal potential */

  nval(amh,CELPRE1)  = hbp1;	  	/* presynaptic cell type */
  nval(amh,CONPRE1)  = 4;	  	/* presynaptic connection number */

  nval(hbp1,CELPOST4) = dsgc;	  	/* cell type to connect to */
  nval(hbp1,CONPOST4) = 3;	  	/* connection number */
  nval(hbp1,CELDIV4)  = 2;	  	/* number of postsyn cells to conn to */
  nval(hbp1,SYNNUM4)  = 1;	  	/* number of synapses per conn */
  nval(hbp1,SRESP4)   = xampa;		/* synaptic response */
  nval(hbp1,SGAIN4)   = 3;	  	/* synaptic gain */
  nval(hbp1,SCGAIN4)  = 1;	  	/* synaptic second messenger gain */
  nval(hbp1,SCOND4)   = 220e-12;  	/* synaptic conductance */
  nval(hbp1,STHRESH4) = -0.050;  	/* synaptic release threshold */
  nval(hbp1,SNCHAN4)  = 20;	  	/* number of channels */
  nval(hbp1,SDUR4)    = 5;	  	/* synaptic event time const. */
  nval(hbp1,SCDUR4)   = 1;	  	/* second mesng. time const. */
  nval(hbp1,SVSIZ4)   = 1.0;	  	/* synaptic vesicle size */
  nval(hbp1,SUNIT4)   = 22e-12;  	/* unitary conductance */
  nval(hbp1,SVREV4)   = 0;	  	/* synaptic reversal potential */

  nval(dsgc,CELPRE3) = hbp1;	  	/* presynaptic cell type */
  nval(dsgc,CONPRE3) = 4;	  	/* presynaptic connection number */
  nval(dsgc,CELCONV3) = 1e3;		/* number of presyn cells to connect to */

  nval(hbp2,CELPOST1) = gcboff;	  	/* cell type to connect to */
  nval(hbp2,CONPOST1) = 1;	  	/* connection number */
  nval(hbp2,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(hbp2,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(hbp2,SRESP1)   = xampa5;		/* synaptic response */
  nval(hbp2,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(hbp2,SDURH1)   = 5;	  	/* synaptic high pass time const. */
  nval(hbp2,SNFILTH1) = 1;	  	/* high pass number of filters */
  nval(hbp2,SHGAIN1)  = 14;	  	/* high pass gain */
  nval(hbp2,SHOFFS1)  = -0.0045;	        /* synaptic high pass offset. */
  nval(hbp2,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(hbp2,SCOND1)   = 2e-10;		/* synaptic conductance */
  nval(hbp2,STHRESH1) = -0.045;  	/* synaptic release threshold */
  nval(hbp2,SNCHAN1)  = 20;	  	/* number of channels */
  nval(hbp2,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(hbp2,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(hbp2,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(hbp2,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(hbp2,SVREV1)   = 0;	  	/* synaptic reversal potential */
  nval(hbp2,STRCONC1) = 0.001;		/* synaptic neurotransmitter conc */
 
  nval(gcboff,CELPRE1) = hbp2;	  	/* presynaptic cell type */
  nval(gcboff,CONPRE1) = 1;	  	/* presynaptic connection number */
  nval(gcboff,CELCONV1) = 1e3;		/* number of presyn cells to connect to */

  nval(aii,CELPOST1) = aii;	  	/* cell type to connect to */
  nval(aii,CONPOST1) = 2;	  	/* connection number */
  nval(aii,CELDIV1)  = 5;	  	/* number of postsyn cells to conn to */
  nval(aii,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(aii,SRESP1)   = xgapj;		/* synaptic response */
  nval(aii,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(aii,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(aii,SCOND1)   = 20e-12;  	/* synaptic conductance */
  nval(aii,STHRESH1) = -0.050; 		/* synaptic release threshold */
  nval(aii,SNCHAN1)  = 20;	  	/* number of channels */
  nval(aii,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(aii,SNFILT1)  = 1;	  	/* synaptic event time const. */
  nval(aii,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(aii,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(aii,SUNIT1)   = 22e-12; 		/* unitary conductance */
  nval(aii,SVREV1)   = 0;	  	/* synaptic reversal potential */

  nval(aii,CELPRE2)  = aii;	  	/* presynaptic cell type */
  nval(aii,CONPRE2)  = 1;	  	/* presynaptic connection number */
  nval(aii,CELCONV2) = 5;	  	/* number of presyn cells to connect to */

  nval(hbp1,CELPRE5)  = aii;	  	/* presynaptic cell type */
  nval(hbp1,CONPRE5)  = 3;	  	/* presynaptic connection number */
  nval(hbp1,CELCONV5) = 0;	  	/* number of presyn cells to connect to */

  nval(aii,CELPOST2) = dbp1;	  	/* cell type to connect to */
  nval(aii,CONPOST2) = 7;	  	/* connection number */
  nval(aii,CELDIV2)  = 0;	  	/* number of postsyn cells to conn to */
  nval(aii,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(aii,SRESP2)   = xgapj;		/* synaptic response */
  nval(aii,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(aii,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(aii,SCOND2)   = 0;  		/* synaptic conductance */
  nval(aii,STHRESH2) = 0; 		/* synaptic release threshold */
  nval(aii,SNCHAN2)  = 20;	  	/* number of channels */
  nval(aii,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(aii,SNFILT2)  = 1;	  	/* synaptic event time const. */
  nval(aii,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(aii,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(aii,SUNIT2)   = 22e-12; 		/* unitary conductance */
  nval(aii,SVREV2)   = 0;	  	/* synaptic reversal potential */
  nval(aii,STRCONC2) = 0;		/* synaptic neurotransmitter conc */

  nval(dbp1,CELPRE7)  = aii;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE7)  = 2;	  	/* presynaptic connection number */
  nval(dbp1,CELCONV7) = 0;	  	/* number of presyn cells to connect to */

  nval(aii,CELPOST3) = hbp1;	  	/* cell type to connect to */
  nval(aii,CONPOST3) = 5;	  	/* connection number */

  nval(sbac,CELPOST1) = dbp1;	  	/* cell type to connect to */
  nval(sbac,CONPOST1) = 4;	  	/* connection number */
  nval(sbac,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(sbac,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(sbac,SRESP1)   = xgaba;		/* synaptic response */
  nval(sbac,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(sbac,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(sbac,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(sbac,STHRESH1) = -0.050;  		/* synaptic release threshold */
  nval(sbac,SNCHAN1)  = 20;	  	/* number of channels */
  nval(sbac,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(sbac,SNFILT1)  = 1;	  	/* synaptic event time const. */
  nval(sbac,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(sbac,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(sbac,SUNIT1)   = 22e-12;  		/* unitary conductance */
  nval(sbac,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  nval(dbp1,CELPRE4) = sbac;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE4) = 1;	  	/* presynaptic connection number */

  nval(sbac,CELPOST2) = dsgc;	  	/* cell type to connect to */
  nval(sbac,CONPOST2) = 2;	  	/* connection number */
  nval(sbac,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(sbac,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(sbac,SRESP2)   = xgaba;		/* synaptic response */
  nval(sbac,SGAIN2)   = 3;	  	/* synaptic gain */
  nval(sbac,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(sbac,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(sbac,STHRESH2) = -0.050;  		/* synaptic release threshold */
  nval(sbac,SNCHAN2)  = 20;	  	/* number of channels */
  nval(sbac,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(sbac,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(sbac,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(sbac,SUNIT2)   = 22e-12;  		/* unitary conductance */
  nval(sbac,SVREV2)   = -0.06;	  	/* synaptic reversal potential */

  nval(dsgc,CELPRE2) = sbac;	  	/* presynaptic cell type */
  nval(dsgc,CONPRE2) = 2;	  	/* presynaptic connection number */

  nval(sbac,CELPOST3) = sbac;	  	/* cell type to connect to */
  nval(sbac,CONPOST3) = 3;	  	/* connection number */
  nval(sbac,CELDIV3)  = 2;	  	/* number of postsyn cells to conn to */
  nval(sbac,SYNNUM3)  = 1;	  	/* number of synapses per conn */
  nval(sbac,AUTAPSE3) = 0;	  	/* synapse back to presynaptic node */
  nval(sbac,SENSCA3)  = 1;	  	/* synaptic sensitivity to calcium */
  nval(sbac,SRESP3)   = xgaba;		/* synaptic response */
  nval(sbac,SGAIN3)   = 2;	  	/* synaptic gain */
  nval(sbac,SCGAIN3)  = 1;	  	/* synaptic second messenger gain */
  nval(sbac,SCOND3)   = 220e-12;  	/* synaptic conductance */
  nval(sbac,STHRESH3) = -0.050;  	/* synaptic release threshold */
  nval(sbac,SNCHAN3)  = 20;	  	/* number of channels */
  nval(sbac,SDUR3)    = 2;	  	/* synaptic event time const. */
  nval(sbac,SCDUR3)   = 1;	  	/* second mesng. time const. */
  nval(sbac,SVSIZ3)   = 1.0;	  	/* synaptic vesicle size */
  nval(sbac,SUNIT3)   = 22e-12;  	/* unitary conductance */
  nval(sbac,SVREV3)   = -0.06;	  	/* synaptic reversal potential */

  nval(sbac,CELPRE3) = sbac;	  	/* presynaptic cell type */
  nval(sbac,CONPRE3) = 3;	  	/* presynaptic connection number */

  nval(am,CELPOST1) = dbp1;	  	/* cell type to connect to */
  nval(am,CONPOST1) = 5;	  	/* connection number */
  nval(am,CELDIV1)  = 200;	  	/* number of postsyn cells to conn to */
  nval(am,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(am,SRESP1)   = xgaba;		/* synaptic response */
  nval(am,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(am,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(am,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(am,STHRESH1) = -0.050;  		/* synaptic release threshold */
  nval(am,SNCHAN1)  = 20;	  	/* number of channels */
  nval(am,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(am,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(am,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(am,SUNIT1)   = 22e-12;  		/* unitary conductance */
  nval(am,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  nval(dbp1,CELPRE5) = am;	  	/* presynaptic cell type */
  nval(dbp1,CONPRE5) = 1;	  	/* presynaptic connection number */
  nval(dbp1,CELCONV5) = 10;	  	/* number of presyn cells to connect to */

  nval(am,CELPOST2) = gca;	  	/* cell type to connect to */
  nval(am,CONPOST2) = 2;	  	/* connection number */
  nval(am,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(am,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(am,SRESP2)   = xgaba;		/* synaptic response */
  nval(am,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(am,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(am,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(am,STHRESH2) = -0.050;  		/* synaptic release threshold */
  nval(am,SNCHAN2)  = 20;	  	/* number of channels */
  nval(am,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(am,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(am,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(am,SUNIT2)   = 22e-12;  		/* unitary conductance */
  nval(am,SVREV2)   = -0.06;	  	/* synaptic reversal potential */

  nval(gca,CELPRE2)  = am;	  	/* presynaptic cell type */
  nval(gca,CONPRE2)  = 2;	  	/* presynaptic connection number */

  nval(am2,CELPOST1) = dbp2;	  	/* cell type to connect to */
  nval(am2,CONPOST1) = 5;	  	/* connection number */
  nval(am2,CELDIV1)  = 200;	  	/* number of postsyn cells to conn to */
  nval(am2,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(am2,SRESP1)   = xgaba;		/* synaptic response */
  nval(am2,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(am2,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(am2,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(am2,STHRESH1) = -0.050;  		/* synaptic release threshold */
  nval(am2,SNCHAN1)  = 20;	  	/* number of channels */
  nval(am2,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(am2,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(am2,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(am2,SUNIT1)   = 22e-12;  		/* unitary conductance */
  nval(am2,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  // nval(dbp2,CELPRE5) = am2;	  	/* presynaptic cell type */
  // nval(dbp2,CONPRE5) = 1;	  	/* presynaptic connection number */
  // nval(dbp2,CELCONV5) = 10;	  	/* number of presyn cells to connect to */

  nval(am2,CELPOST2) = gca;	  	/* cell type to connect to */
  nval(am2,CONPOST2) = 2;	  	/* connection number */
  nval(am2,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(am2,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(am2,SRESP2)   = xgaba;		/* synaptic response */
  nval(am2,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(am2,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(am2,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(am2,STHRESH2) = -0.050;  		/* synaptic release threshold */
  nval(am2,SNCHAN2)  = 20;	  	/* number of channels */
  nval(am2,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(am2,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(am2,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(am2,SUNIT2)   = 22e-12;  		/* unitary conductance */
  nval(am2,SVREV2)   = -0.06;	  	/* synaptic reversal potential */

//  nval(gca,CELPRE2)  = am;	  	/* presynaptic cell type */
//  nval(gca,CONPRE2)  = 2;	  	/* presynaptic connection number */

  nval(amh,CELPOST1) = hbp1;	  	/* cell type to connect to */
  nval(amh,CONPOST1) = 4;	  	/* connection number */
  nval(amh,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(amh,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(amh,SRESP1)   = xgaba;		/* synaptic response */
  nval(amh,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(amh,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(amh,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(amh,STHRESH1) = -0.050;  	/* synaptic release threshold */
  nval(amh,SNCHAN1)  = 20;	  	/* number of channels */
  nval(amh,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(amh,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(amh,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(amh,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(amh,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  nval(hbp1,CELPRE4)  = amh;	  	/* presynaptic cell type */
  nval(hbp1,CONPRE4)  = 1;	  	/* presynaptic connection number */

  nval(amh,CELPOST2) = gcaoff;	  	/* cell type to connect to */
  nval(amh,CONPOST2) = 2;	  	/* connection number */
  nval(amh,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(amh,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(amh,SRESP2)   = xgaba;		/* synaptic response */
  nval(amh,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(amh,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(amh,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(amh,STHRESH2) = -0.050;  		/* synaptic release threshold */
  nval(amh,SNCHAN2)  = 20;	  	/* number of channels */
  nval(amh,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(amh,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(amh,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(amh,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(amh,SVREV2)   = -0.06;	  	/* synaptic reversal potential */

  nval(gcaoff,CELPRE2) = amh;	  	/* presynaptic cell type */
  nval(gcaoff,CONPRE2) = 2;	  	/* presynaptic connection number */

  nval(amh2,CELPOST1) = hbp2;	  	/* cell type to connect to */
  nval(amh2,CONPOST1) = 4;	  	/* connection number */
  nval(amh2,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(amh2,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(amh2,SRESP1)   = xgaba;		/* synaptic response */
  nval(amh2,SGAIN1)   = 2;	  	/* synaptic gain */
  nval(amh2,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(amh2,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(amh2,STHRESH1) = -0.050;  	/* synaptic release threshold */
  nval(amh2,SNCHAN1)  = 20;	  	/* number of channels */
  nval(amh2,SDUR1)    = 2;	  	/* synaptic event time const. */
  nval(amh2,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(amh2,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(amh2,SUNIT1)   = 22e-12;  	/* unitary conductance */
  nval(amh2,SVREV1)   = -0.06;	  	/* synaptic reversal potential */

  // nval(hbp2,CELPRE4)  = amh2;	  	/* presynaptic cell type */
  // nval(hbp2,CONPRE4)  = 1;	  	/* presynaptic connection number */

  nval(amh2,CELPOST2) = gcaoff;	  	/* cell type to connect to */
  nval(amh2,CONPOST2) = 2;	  	/* connection number */
  nval(amh2,CELDIV2)  = 2;	  	/* number of postsyn cells to conn to */
  nval(amh2,SYNNUM2)  = 1;	  	/* number of synapses per conn */
  nval(amh2,SRESP2)   = xgaba;		/* synaptic response */
  nval(amh2,SGAIN2)   = 2;	  	/* synaptic gain */
  nval(amh2,SCGAIN2)  = 1;	  	/* synaptic second messenger gain */
  nval(amh2,SCOND2)   = 220e-12;  	/* synaptic conductance */
  nval(amh2,STHRESH2) = -0.050;  		/* synaptic release threshold */
  nval(amh2,SNCHAN2)  = 20;	  	/* number of channels */
  nval(amh2,SDUR2)    = 2;	  	/* synaptic event time const. */
  nval(amh2,SCDUR2)   = 1;	  	/* second mesng. time const. */
  nval(amh2,SVSIZ2)   = 1.0;	  	/* synaptic vesicle size */
  nval(amh2,SUNIT2)   = 22e-12;  	/* unitary conductance */
  nval(amh2,SVREV2)   = -0.06;	  	/* synaptic reversal potential */

  // nval(gcaoff,CELPRE2) = amh2;	/* presynaptic cell type */
  // nval(gcaoff,CONPRE2) = 2;	  	/* presynaptic connection number */

  nval(ams,CELPOST1) = dsgc;	  	/* cell type to connect to */
  nval(ams,CONPOST1) = 4;	  	/* connection number */
  nval(ams,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(ams,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(ams,SRESP1)   = xgaba;		/* synaptic response */
  nval(ams,SGAIN1)   = 3;	  	/* synaptic gain */
  nval(ams,SDURH1)   = 40;	  	/* high pass time constant */
  nval(ams,SNFILTH1) = 1;	  	/* high pass number of filters */
  nval(ams,SHGAIN1)  = .7;	  	/* high pass number gain */
  nval(ams,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(ams,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(ams,STHRESH1) = -0.045;  	/* synaptic release threshold */
  nval(ams,STRCONC1) = 200e-6;		/* synaptic neurotransmitter conc */
  nval(ams,SNCHAN1)  = 20;	  	/* number of channels */
  nval(ams,SDUR1)    = 15;	  	/* synaptic event time const. */
  nval(ams,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(ams,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(ams,SUNIT1)   = 22e-12; 		/* unitary conductance */
  nval(ams,SVREV1)   = -0.068;	  	/* synaptic reversal potential */

  nval(dsgc,CELPRE4) = ams;	  	/* presynaptic cell type */
  nval(dsgc,CONPRE4) = 1;	  	/* presynaptic connection number */
  nval(dsgc,CELCONV4) = 1e3;		/* number of presyn cells to connect to */

  nval(amhs,CELPOST1) = dsgc;	  	/* cell type to connect to */
  nval(amhs,CONPOST1) = 8;	  	/* connection number */
  nval(amhs,CELDIV1)  = 2;	  	/* number of postsyn cells to conn to */
  nval(amhs,SYNNUM1)  = 1;	  	/* number of synapses per conn */
  nval(amhs,SRESP1)   = xgaba;		/* synaptic response */
  nval(amhs,SGAIN1)   = 3;	  	/* synaptic gain */
  nval(amhs,SDURH1)   = 40;	  	/* high pass time constant */
  nval(amhs,SNFILTH1) = 1;	  	/* high pass number of filters */
  nval(amhs,SHGAIN1)  = .7;	  	/* high pass number gain */
  nval(amhs,SCGAIN1)  = 1;	  	/* synaptic second messenger gain */
  nval(amhs,SCOND1)   = 220e-12;  	/* synaptic conductance */
  nval(amhs,STHRESH1) = -0.045;  	/* synaptic release threshold */
  nval(amhs,STRCONC1) = 200e-6;		/* synaptic neurotransmitter conc */
  nval(amhs,SNCHAN1)  = 20;	  	/* number of channels */
  nval(amhs,SDUR1)    = 15;	  	/* synaptic event time const. */
  nval(amhs,SCDUR1)   = 1;	  	/* second mesng. time const. */
  nval(amhs,SVSIZ1)   = 1.0;	  	/* synaptic vesicle size */
  nval(amhs,SUNIT1)   = 22e-12; 		/* unitary conductance */
  nval(amhs,SVREV1)   = -0.068;	  	/* synaptic reversal potential */

  nval(dsgc,CELPRE8) = amhs;	  	/* presynaptic cell type */
  nval(dsgc,CONPRE8) = 1;	  	/* presynaptic connection number */
  nval(dsgc,CELCONV4) = 1e3;		/* number of presyn cells to connect to */

/*--------------------------------------------*/


  printf ("#\n");
  printf ("#  nval.n\n");
  printf ("#\n");
  printf ("#  Neuron parameters for retsim.cc simulation script\n");
  printf ("#  Generated  %s     \n", xsystem("date"));
  printf ("#\n");
  printf ("#    Modify this file to change parameter values.\n");
  printf ("#\n");
  printf ("#    Original created by \"maknval.n > nval.n\"\n");
  printf ("#    To add parameters, edit \"maknval.cc\", compile and run \"maknval > nval.n\",\n");
  printf ("#    then copy nval.n to \"nval.h\".  Remove the param defs from the end of nval.n.\n");
  printf ("#    Remove the nval.n table at the beginning of nval.h.\n");
  printf ("#    cp nval.h nval_var.h\n");
  printf ("#    cp nval.h nval_var.cc\n");
  printf ("#    cp nval.h nval_var_set.cc\n");
  printf ("#    and edit these files, then remove this content from nval.h.\n");
  printf ("#    Last, \"make clean\" and \"make retsim\".\n");

  pind = 0;
  printparams(CELPRE1,' ','0');
  printparams(CELPOST1,'#',' ');
  printparams(CELPOST2,'#',' ');
  printparams(CELPOST3,'#',' ');
  printparams(CELPOST4,'#',' ');
  printparams(CELPOST5,'#',' ');
  printparams(CELPOST6,'#',' ');
  printparams(CELPOST7,'#',' ');
  printparams(CELPOST8,'#',' ');
  printparams(CELPOST9,'#',' ');
  printparams(SVREV9+1,'#',' ');

  printf ("\n\n");

/*--------------------------------------------*/


  printf ("/*  Neuron parameters for makcel.cc simulation script */\n"); 
  printf ("\n");
  printf ("/* nval.h */\n");
  printf ("\n");
  printf ("/* To add parameters, edit \"maknval.cc\", compile and run \"maknval > nval.n\",\n");
  printf ("   then copy nval.n to \"nval.h\". Remove the param defs from the end of nval.n.\n");
  printf ("   Remove the nval.n table at the beginning of nval.h.\n");
  printf ("   Copy nval.h to \"nval_var.h\", \"nval_var.cc\", and \"nval_var_set.cc\"\n");
  printf ("   and edit these files, then remove this content from nval.h.\n");
  printf ("   Last, \"make clean\" and \"make retsim\".\n");
 
  printf ("*/\n");
  printf ("\n");

  for (j=0; j<=nceltypes; j++) {	/* print cell types */
    printf ("#define %-10s %3d\t%s\n",tname[j][0],j,tname[j][1]);
  };
  printf ("\n");

  pind = 0;
  printparamnums(CELPRE1);
  printparamnums(CELPOST1);
  printparamnums(CELPOST2);
  printparamnums(CELPOST3);
  printparamnums(CELPOST4);
  printparamnums(CELPOST5);
  printparamnums(CELPOST6);
  printparamnums(CELPOST7);
  printparamnums(CELPOST8);
  printparamnums(CELPOST9);
  printparamnums(NPARAMS);
  printparamnums(NPARAMS+1);


  for (j=0; j<=NCONNP; j++) {	/* print presynaptic connection params */
    printf ("#define %-10s %3d\t%s\n",cname[j][0],j,cname[j][1]);
  };
  printf ("\n");

  for (j=0; j<=NSYNP; j++) {	/* print postsynaptic conn/syn params */
    printf ("#define %-10s %3d\t%s\n",sname[j][0],j,sname[j][1]);
  }

  printf ("\n");
  printf ("#define %-10s %3d\t%s\n","NCONNI",NCONNI,"/* number of input connection cell types */");
  printf ("#define %-10s %3d\t%s\n","NCONNO",NCONNO,"/* number of output connection cell types  */");

  printf ("\n");
  for (j=1; j<=nresptypes; j++) {	/* print response types */
    printf ("#define %-10s %3d\t%s\n",rname[j][0],j,rname[j][1]);
  };
  printf ("\n");

  printf ("// nval_var.h\n\n");
  pind = 0;
  printparamvardefs(CELPRE1);
  printparamvardefs(CELPOST1);
  printparamvardefs(CELPOST2);
  printparamvardefs(CELPOST3);
  printparamvardefs(CELPOST4);
  printparamvardefs(CELPOST5);
  printparamvardefs(CELPOST6);
  printparamvardefs(CELPOST7);
  printparamvardefs(CELPOST8);
  printparamvardefs(CELPOST9);
  printparamvardefs(NPARAMS);

  printf ("// nval_var.cc\n\n");
  pind = 0;
  printparamvars(CELPRE1);
  printparamvars(CELPOST1);
  printparamvars(CELPOST2);
  printparamvars(CELPOST3);
  printparamvars(CELPOST4);
  printparamvars(CELPOST5);
  printparamvars(CELPOST6);
  printparamvars(CELPOST7);
  printparamvars(CELPOST8);
  printparamvars(CELPOST9);
  printparamvars(NPARAMS);

  printf ("// nval_var_set.cc\n\n");
  pind = 0;
  printparamvarset(CELPRE1);
  printparamvarset(CELPOST1);
  printparamvarset(CELPOST2);
  printparamvarset(CELPOST3);
  printparamvarset(CELPOST4);
  printparamvarset(CELPOST5);
  printparamvarset(CELPOST6);
  printparamvarset(CELPOST7);
  printparamvarset(CELPOST8);
  printparamvarset(CELPOST9);
  printparamvarset(NPARAMS);
}

/*-----------------------------------------*/

void printparams(int stopval, char fch, char lch)
{
   int i,j;

 if (pind<=stopval) printceltypes(fch,lch);

  /* Print out numeric parameter values, */ 
  /*  but string names for four rows: CELPOST, SRESP, DYADTYP, and CELPRE */

  for (j=pind; j<stopval; j++) {
    printf ("  ");
    for (i=0; i<nceltypes; i++ ) {
      if (j >= CELPOST1 && (j <=(CELPOST1+NCONNO*NSYNP)) && ((j-CELPOST1)%NSYNP)==0) {
        if (nval(i,j) >= 0)
	  printf ("%7s ", strtolower(tname[(int)nval(i,j)][0]));
        else
	  printf ("%7d ", -1);
      } else
      if (j >= SRESP1 && (j <= (SRESP1+NCONNO*NSYNP)) && ((j-SRESP1)%NSYNP)==0) {
        if (nval(i,j) >= 0)
	  printf ("%7s ", strtolower(rname[(int)nval(i,j)][0]));
        else
	  printf ("%7d ", -1);
      } else
      if (j >= DYADTYP1 && (j <= (DYADTYP1+NCONNO*NSYNP)) && ((j-DYADTYP1)%NSYNP)==0) {
        if (nval(i,j) >= 0)
	  printf ("%7s ", strtolower(tname[(int)nval(i,j)][0]));
        else
	  printf ("%7d ", -1);
      } else
      if (j >= CELPRE1 && (j <= (CELPRE1+NCONNI*NCONNP)) && ((j-CELPRE1)%NCONNP)==0) {
        if (nval(i,j) >= 0)
	  printf ("%7s ", strtolower(tname[(int)nval(i,j)][0]));
        else
	  printf ("%7d ", -1);
      }
      else 
	printf ("%7.4g ", nval(i,j));
    };
    printf ("  _%-9s # %s\n",pname[j][0],pname[j][1]);
    //printf (" # %3d %s\n",j,pname[j][1]);
    //printf (" # %3d %-9s %s\n",j,pname[j][0],pname[j][1]); /* the old way */
    //printf (" #  %3d %s\n",j,pname[j][0]);
  }
  pind = j;
}

