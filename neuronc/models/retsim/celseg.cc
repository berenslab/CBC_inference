/* module celseg.cc */

/* Biophysical properties for makcel.n.  */
/* Creates cable segments and membrane channels. */

extern "C" {
#include <math.h>
#include <stdio.h>
#include <string.h>
}
#include "ncio.h"
#include "namedparams.h"

#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "celsegca.h"

extern int malloc(int);

int *ndens[NCELTYPES] = {NULL};		/* density file indexed by [ct][cn] */
int *cellnums[NCELTYPES] = {NULL};	/* cell numbers indexed by [ct][n], [0] is num  */
const char *densfil[NCELTYPES][NDENS] = {NULL}; /* density file filenames */
double celdens[NCELTYPES][NDENS][NCHANS][R_NREGIONS];/* density data */
double zerodens[NCHANS][R_NREGIONS]; 	    /* zero density file, first col = label (int) */
double chval[NCELTYPES][NCHANS][NCHRATE];   /* channel offsets, rates */
double chanunit[NCHANS];   		    /* channel unitary conductances */
int chanelem[NCELTYPES][NCHANS] = {0};      /* channel element enames (soma only) */
double *cdens = NULL;			    /* temp array for density data */

/* suffixes to append when reading the channel properties file */
//const char *chanfile_props[]  = {string(".offm"), string(".offh"), string(".taua"),
//		   	    string(".taub"), string(".tauc"), string(".taud")};
const char *chanfile_props[]  = {".offm", ".offh", ".taua",".taub", ".tauc", ".taud"};

double ratio_kdr = 0;
double ratio_ka  = 0;
double ratio_sk1 = 0;
double ratio_sk2 = 0;
double ratio_bk  = 0;
double ratio_k   = 0;
extern int R5,_KDR;

/*---------------------------------------------------------*/


void findchanparms(char *chparmfil)

/* read channel parameters in chparm array */

{
   int ch, ct, r;
   int nrows, ncols;
   int chantyp, chindex, paramtyp, param;
   double *chparms = NULL;                       /* temp array for channel params data */
   static const char delim[] = {":"};

  chparms = fread(chparmfil, &nrows, &ncols);  /* read density file into temp arr */

  if (chparms == NULL) {
	fprintf(stderr, "Cannot open channel parameter file '%s'!\n", chparmfil);

  } else {	

     // Set up chanparams array
     // Start by copying default values

   for (ct=0; ct<NCELTYPES; ct++) {             
      for (ch=0; ch<NCHANS; ch++) {            
         chval[ct][ch][CHOFFM] = 0;
         chval[ct][ch][CHOFFH] = 0;
         chval[ct][ch][CHTAUA] = 1;
         chval[ct][ch][CHTAUB] = 1;
         chval[ct][ch][CHTAUC] = 1;
         chval[ct][ch][CHTAUD] = 1;
       }
    }

        // Chanparams file has variable format (like density file), 
	//   nrows, ncols are set by index labels at top or left column.
        // Skip values with indexes out of range

    if (ninfo > 3) {
       fprintf(stderr, "#\n");
       fprintf(stderr, "# Channel params read from chanparams file '%s'\n",chparmfil);
       fprintf(stderr, "#\n");
       fprintf(stderr, "       0%s     ",delim);
       for (r=1; r<ncols; r++) {                /* skip first col = label */
           ct = *(chparms+r);                       /* cell index is first row */
           fprintf(stderr, "%-6s  ",cname[ct]);
       }
       fprintf(stderr, "\n");
       fprintf(stderr, "#\n");
    }
    for (ch=1; ch<nrows; ch++) {              /* copy chanparam data by row  */
         chindex = ch*ncols;                  /* channel index is first col (0) */
         paramtyp = *(chparms + chindex);     /* get channel parameter (col 0) */
         chantyp = paramtyp / NCHRATE;	      /* channel type */ 
	 param = paramtyp - (chantyp*NCHRATE);		  
         if (paramtyp<0 || paramtyp>NCHANS*NCHRATE) continue;  /* ignore if channel parameter isn't set */
         for (r=1; r<ncols; r++) {                /* skip first col = label */
            ct = *(chparms+r);                       /* cell index is first row */
            if (ct<0 || ct>NCELTYPES) continue;  /* ignore if cell isn't set */
            chval[ct][chantyp][param] = *(chparms+chindex+r);
         }
	 if (ninfo > 3) {
		 char chnambuf[30];
	     sprintf(chnambuf, "%s.%s%s  ", chnamea[chantyp],parmname[param],delim);
	     fprintf(stderr, "%-10.10s    ", chnambuf);
             for (r=1; r<ncols; r++) {                /* skip first col = label */
                 ct = *(chparms+r);                       /* cell index is first row */
		 fprintf(stderr, "%-6g  ", chval[ct][chantyp][param]);
             }
	     fprintf(stderr,"\n");
             if (param==NCHRATE-1) fprintf(stderr,"#\n");
         }
    }
    efree (chparms);
  }
}

/*----------------------------------------------------------*/

void celseg_init(void)
{
	int dn, n, ch, r, rt, ncells;
        char *chparmfil;

	if (notinit(ttxbath))		ttxbath = 0;
	if (notinit(ttxsoma))		ttxsoma = 0;
	if (notinit(ttxdend))		ttxdend = 0;
	if (notinit(tea))		tea     = 0;
	if (notinit(fourap))		fourap  = 0;
	if (notinit(zd7288))		zd7288  = 0;
	if (notinit(ibtox))		ibtox   = 0;

	if (notinit(Chnoise))           Chnoise = 0;
	
	for (n=0; n<nceltypes; n++) { 		       /* set default offsets, rates */
	  for (ch=0; ch<NCHANS; ch++) {
	    chanelem[n][ch] = MININT;
	    for (rt=0; rt<NCHRATE; rt++) {
	      chval[n][ch][CHOFFM] = 0;
	      chval[n][ch][CHOFFH] = 0;
	      chval[n][ch][CHTAUA] = 1;
	      chval[n][ch][CHTAUB] = 1;
	      chval[n][ch][CHTAUC] = 1;
	      chval[n][ch][CHTAUD] = 1;
	    }
	  }
	}
	
	for (ch=0; ch<NCHANS; ch++) { 
	  for (r=0; r<R_NREGIONS; r++) { 
	    zerodens[ch][r] = 0;  		/* zero array */
	  }
	}
	for (r=0; r<R_NREGIONS; r++) {     /* set reasonable defaults */
	    zerodens[C_VSTART][r] = vcl;    /* for membrane properties */
	    zerodens[C_VREV][r]   = vcl;
	    zerodens[C_RM][r]     = drm;
	    zerodens[C_CM][r]     = dcm;
	    zerodens[C_RI][r]     = dri;
	    zerodens[C_COLOR][r]  = 1;	    /* default color blue */
	}
	for (n=0; n<nceltypes; n++) {	     /* zero density array that will be filled */
	  for (dn=0; dn<NDENS; dn++) { 
	    for (ch=0; ch<NCHANS; ch++) {
	      for (r=0; r<R_NREGIONS; r++) { 
              celdens[n][dn][ch][r] = 0;
	      }
	    }
	  }
	}
	for (n=0; n<nceltypes; n++) {	     /* set up ndens and cellnums to be indexed by [ct][cn] */
	  if (getn(n,MAKE)>0) {
		ncells = getn(n,MAXNUM)+1;			  // get max number of cells 
		ndens[n] = (int *)emalloc(ncells*sizeof(int));    // start filled with zeros 
		cellnums[n] = (int *)emalloc(ncells*sizeof(int)); 
	  } else {
		ndens[n] = NULL;
		cellnums[n] = NULL;
          }
	}	

	/* Format of param files for channels in a cell's different regions (S/cm2): */
	
	/* (Note that predefined variables can be set in file) */
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	# Example membrane properties (density) file
	#
	# The number of rows and columns is variable.
	#  The values for a missing column or row are set 
	#  to zero in "celdens[][][][]".
	#
	# Note that the first column is the channel label. 
	#   which is a variable set in retsim_var.cc
	#
	# ch dend dend_prox dend_dist  soma  hillock  thin   	axon   axon_dist ax_node varicos
	#
	#     R1	R2	R3	R4	R5	R6	R7	R8     R9     R10
        #
	0     DENDD     DENDP   SOMA    HILLCK  AXONT   AXON    AXOND   VARIC  R9
	#
	_NA   40e-3	0e-3	93e-3	85e-3	80e-3	55e-3	0e-3    0e-3   0e-3   # Na
	_NA5  0e-3      0e-3    0e-3 	0e-3   	0e-3    0e-3  	0e-3    0e-3   0e-3   # Na5
	_NA6  0e-3	0e-3  	0e-3    0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # Na6
	_NA8  0e-3	0e-3  	0e-3    0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # Na8
	_KDR  10e-3	0e-3	10e-3   10e-3	10e-3	10e-3   0e-3    0e-3   0e-3   # Kdr
	_KA   5e-3 	0e-3	10e-3	10e-3	10e-3	0e-3 	0e-3    0e-3   0e-3   # KA
	_K4   0e-3	0e-3	0e-3	0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # KIh
	_K5  0e-3	0e-3	0e-3	0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # Khz
	_K6  0e-3	0e-3	0e-3	0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # Kv3
	_K7  0e-3	0e-3	0e-3	0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # Kv3b
	_sKCA1 0.03e-3 0e-3  	0.03e-3 0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # sKCa1
	_sKCA2 0.05e-3 0e-3  	0.05e-3 0e-3	0e-3	0e-3	0e-3    0e-3   0e-3   # sKCa2
	_BKCA  0.7e-3	0e-6	0.7e-3	0e-6	0e-6	0e-6	0e-6    0e-3   0e-3   # bKCa
	_CLCA1 0e-3	0e-3	0e-3    0e-3	0e-3    0e-3	0e-3    0e-3   0e-3   # ClCa
	_CLCA2 0e-3	0e-3	0e-3    0e-3	0e-3    0e-3	0e-3    0e-3   0e-3   # ClCaC
	_CA0   0.2e-3	0e-3	0.2e-3  0.2e-3	0e-3    0e-3	0e-3    0e-3   0e-3   # Ca
	_CA5   0e-3	0e-3	0e-3    0e-3	0e-3    0e-3	0e-3    0e-3   0e-3   # Ca5
	_CAP   5e-6      0e-3    5e-6    5e-6	1e-12  	0e-3    0e-3    0e-3   0e-3   # CAPUMP
	_CAE   0e-6      0e-3    0e-6    0e-6	0e-12  	0e-3    0e-3    0e-3   0e-3   # CAEXCH
	_VST  dvs       dvs     dvs     dvs     dvs     dvs     dvs     dvs    dvs    # vstart
	_VRV  vcl       vcl     vcl     vcl     vcl     vcl     vcl     vcl    vcl    # vrev
	_RM   drm       drm     drm     drm     drm     drm     drm     drm    drm    # Rm
	_CM   dcm       dcm     dcm     dcm     dcm     dcm     dcm     dcm    dcm    # Cm
	_DIA  ddia      ddia    ddia    ddia    ddia    ddia    ddia    ddia   ddia   # Dia
	_CPL  cplam     cplam   cplam   cplam   cplam   cplam   cplam   cplam  cplam  # cplam
	_CMUL  1     	1   	1   	1   	1   	1   	1   	1  	1     # cmul
	_COL  green     blue    red     blue    blue    magenta yellow  gray   ltred  # color

	#
	# */
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	/* set channel unitary currents */
	
	chanunit[C_NA0]  = dnau;
	chanunit[C_NA1]  = dnau;
	chanunit[C_NA2]  = dnau;
	chanunit[C_NA3]  = dnau;
	chanunit[C_NA4]  = dnau;
	chanunit[C_NA5]  = dnau;
	chanunit[C_NA6]  = dnau;
	chanunit[C_NA8]  = dnau;
	chanunit[C_CA0]  = dcalu;
	chanunit[C_CA1]  = dcalu;
	chanunit[C_CA2]  = dcalu;
	chanunit[C_CA3]  = dcalu;
	chanunit[C_CA4]  = dcalu;
	chanunit[C_CA5]  = dcalu;
	chanunit[C_CA6]  = dcalu;
	chanunit[C_CA7]  = dcalu;
	chanunit[C_K0]   = dku;
	chanunit[C_K1]   = dku;
	chanunit[C_K2]   = dku;
	chanunit[C_K3]   = dkau;
	chanunit[C_K4]   = dkau;
  
  // joesterle begin
  chanunit[C_HCN1]   = dkau;
  chanunit[C_HCN2]   = dkau;
  chanunit[C_HCN3]   = dkau;
  chanunit[C_HCN4]   = dkau;
  // joesterle end
  
	chanunit[C_K5]   = dkau;
	chanunit[C_K6]   = dku;
	chanunit[C_K7]   = dku;
  
	chanunit[C_KCA0] = dkcasu;
	chanunit[C_KCA1] = dkcasu;
	chanunit[C_KCA2] = dkcabu;
	chanunit[C_KCA3] = dkcabu;	// BK
	chanunit[C_KCA4] = dkcasu;
	chanunit[C_KCA5] = dkcasu;
	chanunit[C_KCA6] = dkcabu;	// BK
	chanunit[C_CLCA1] = dku;
	chanunit[C_CLCA2]= dku;
	chanunit[C_AMPA1]= dampau;
	chanunit[C_AMPA2]= dampau;
	chanunit[C_AMPA3]= dampau;
	chanunit[C_AMPA4]= dampau;
	chanunit[C_AMPA5]= dampau;
	chanunit[C_GABA1]= dgabau;
	chanunit[C_GABA2]= dgabau;
	chanunit[C_GABA3]= dgabau;
	chanunit[C_GABA4]= dgabau;
	chanunit[C_GLY]  = dglyu;
	
	/* Note that if left zero, Kdr and KA densities and
	    dend_prox and axon_dist densities are automatically filled in.  */
	
	ratio_k=1;	/* if K curr-dens = 0, use ratio to Na to compute curr dens */
	
	dpnak = 0;	/* set Na-permeability of K channel to zero */
	dpkna = 0;	/* set K-permeability of Na channel to zero */
	
	
	if (notinit(ratio_kdr)) ratio_kdr = 0.25;
	if (notinit(ratio_ka))  ratio_ka  = 0.4;
	if (notinit(ratio_sk1)) ratio_sk1 = 0.001;
	if (notinit(ratio_sk2)) ratio_sk2 = 0.002;
	if (notinit(ratio_bk))  ratio_bk  = 0.005;
	
	if (notinit(nshell)) nshell           = 1;
	if (notinit(nshell_soma)) nshell_soma = 10;
	
	//ddca = 1e-10;
	//ddca = 2e-6;
	
	//try to open the chan params file

        chparmfil = emalloc(min(FILNAMSIZ,strlen(confdir)) + min(FILNAMSIZ,strlen(chanparamsfile)) + 2);

        if (streq(confdir,"")) sprintf (chparmfil,"%.50s",chanparamsfile);
	else 		       sprintf (chparmfil,"%.50s/%.50s",confdir,chanparamsfile);

         findchanparms(chparmfil);

/* original way to read chanparams, using "namedparams.cc" */


/*	FILE *fd = fopen(chparmfil, "r");
	if (fd == NULL) {
		fprintf(stderr, "Cannot open channel parameter file '%s'!\n", chparmfil);
	} else {		
	   //read chan params file and set channel parameters in chval array accordingly
	   fclose(fd);
	   namedparams *np = new namedparams(chparmfil);				
		
	   if (ninfo > 3) fprintf(stderr, "Channel params read from %s:\n",chparmfil);
	   for (n = 0; n < NCELTYPES; n++) {			
	       for (ch = 0; ch < NCHANS; ch++) {
		   string cellname = string(cname[n]);
		   string channame = string(chnamea[ch]);
		   string soffm = (channame + chanfile_props[0]);
		   string soffh = (channame + chanfile_props[1]);
		   string staua = (channame + chanfile_props[2]);
		   string staub = (channame + chanfile_props[3]);
		   string stauc = (channame + chanfile_props[4]);
		   string staud = (channame + chanfile_props[5]);
				
		   chval[n][ch][CHOFFM] = np->get(cellname, soffm, 0.000);
		   chval[n][ch][CHOFFH] = np->get(cellname, soffh, 0.000);
		   chval[n][ch][CHTAUA] = np->get(cellname, staua, 1);
		   chval[n][ch][CHTAUB] = np->get(cellname, staub, 1);
		   chval[n][ch][CHTAUC] = np->get(cellname, stauc, 1);
		   chval[n][ch][CHTAUD] = np->get(cellname, staud, 1);
		   if (ninfo > 3 && getn(n,MAKE) && getn(n,BIOPHYS)) {
		fprintf(stderr, "cellname=%s, chname=%s, offm=%g, offh=%g, taua=%g, taub=%g, tauc=%g, taud=%g\n",
				cellname.c_str(), channame.c_str(), chval[n][ch][CHOFFM], chval[n][ch][CHOFFH],
				chval[n][ch][CHTAUA], chval[n][ch][CHTAUB], chval[n][ch][CHTAUC], chval[n][ch][CHTAUD]);
		   }
	       }
	   }
	   delete np;
	} /* orig chanparams method */


/*	if (ninfo > 3) {
		int i,j,k;

	   for (j=0; j<NCHANS; j++) {
	     for (k=0; k<NCHRATE; k++) {
	       for (i=0; i<NCELTYPES; i++) {
	 	 fprintf (stderr,"%g ", chval[i][j][k]);
               }
	       fprintf (stderr,"\n");
             }
	     fprintf (stderr,"\n");
	   }
	 } /* */

/*	
	if (notinit(na1nois))  na1nois = 0;
	if (notinit(na2nois))  na2nois = 0;
	if (notinit(na3nois))  na3nois = 0;
	if (notinit(na4nois))  na4nois = 0;
	if (notinit(na5nois))  na5nois = 0;
	if (notinit(na6nois))  na6nois = 0;
	if (notinit(na8nois))  na8nois = 0;
	if (notinit(k1nois))   k1nois = 0;
	if (notinit(k3nois))   k3nois  = 0;
	if (notinit(k4nois))   k4nois  = 0;
	if (notinit(k5nois))   k5nois  = 0;
	if (notinit(k6nois))   k6nois  = 0;
	if (notinit(k7nois))   k7nois  = 0;
	if (notinit(hcn1nois))   hcn1nois  = 0;
	if (notinit(hcn2nois))   hcn2nois  = 0;
	if (notinit(hcn3nois))   hcn3nois  = 0;
	if (notinit(hcn4nois))   hcn4nois  = 0;
	if (notinit(kca1nois)) kca1nois = 0;
	if (notinit(kca3nois)) kca3nois = 0;
	if (notinit(kca4nois)) kca4nois = 0;
	if (notinit(kca5nois)) kca5nois = 0;
	if (notinit(kca6nois)) kca6nois = 0;
	if (notinit(ca1nois))  ca1nois  = 0;
	if (notinit(ca3nois))  ca3nois  = 0;
	if (notinit(ca5nois))  ca5nois  = 0;
	if (notinit(ca6nois))  ca6nois  = 0;
	if (notinit(ca7nois))  ca7nois  = 0;
	if (notinit(clcanois)) clcanois = 0;
*/

	/*----------------------------------------------------------------*/
	
	if (notinit(dia_min_factor))
	  dia_min_factor = 0; 		/* factor for minimum diam */
	
	if (notinit(dend_dia_factor))
	  dend_dia_factor = 1; 		/* dia factor for dends */
	
	if (notinit(dendd_dia_factor))
	  dendd_dia_factor = 1; 	/* dia factor for distal dends */
	
	if (notinit(dendi_dia_factor))
	  dendi_dia_factor = 1; 	/* dia factor for intermediate dends */
	
	if (notinit(dendp_dia_factor))
	  dendp_dia_factor = 1; 	/* dia factor for proximal dends */
	
	if (notinit(ax_dia_factor))
	  ax_dia_factor = 1; 		/* factor for diam of axons */
	
	if (notinit(cell_dia_factor))
	  cell_dia_factor = 1; 		/* factor for diam of whole cell */
	

}  /* celseg_init() */

/*--------------------------------------------------------*/

double qcond (double cond)
/* convert from "density" S/cm2 to "ndensity" N/um2 at 22 deg C. */
/*  Unitary conductances are defined at 22 deg C. */
{
   return (1e8 * cond);
}

/*--------------------------------------------------------*/

const char *setdensfile(int ctype, int cd, const char *filenam)

/* Set a channel density file name */

{
  return densfil[ctype][cd] = filenam;
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double getdens(int ctype, int dn, int chantyp, int region)

/* Get a channel density value */

{
  return celdens[ctype][dn][chantyp][region];
};

/*-------------------------------------------------*/

void print_chan_offset(int ct)

{
     int ch,chrate;

   if (getn(ct,MAKE) && getn(ct,BIOPHYS)) {
      printf ("#   Channel rates for '%s' cell type, from chanparams array %s after setdens()\n",
		      		cname[ct], chanparamsfile);
      printf ("#\n");
      printf ("#      ");
      for (chrate=0; chrate<NCHRATE; chrate++) {
         printf ("%-8s ",chanfile_props[chrate]);
      }
      printf ("\n");
      for (ch = 0; ch <= ECHANS; ch++) {
        printf ("# %-7s ",chnamea[ch]);
        for (chrate=0; chrate<NCHRATE; chrate++) {
            printf ("%-8.4g ",chval[ct][ch][chrate]);

        }
        printf ("\n");
      }
   }
}

/*-------------------------------------------------*/

void printdens (int ct, int dn, const char *label, double uden)

{
  int ch;
  int region;
  double mult;
  const char *densfile;

  printf ("# \n");
  if ((densfile=densfil[ct][dn]) == NULL) densfile = "zerodens";
  printf ("#   Channel densities for '%s' cell type, density array %d, (%s) from '%s'\n",cname[ct],dn+1,label, densfile);
  printf ("# \n");

  printf ("#          ");
  for (ch=0; ch<(R_NREGIONS-1); ch++)
    printf ("%-8s ",regname[ch]);
  printf ("\n");
  
 // printf ("#          ");
 // for (ch=0; ch<R_NREGIONS; ch++)
 //   printf ("%-7s  ",regnamer[ch]);
 // printf ("\n");

  printf ("# \n");

  for (ch=0; ch<NCHANS; ch++) {
    printf ("#   ");
    printf ("%-6s ",chnamea[ch]);
    for (region=R_1; region<R_NREGIONS; region++) {
      if (ch<=ECHANS) {   /* ECHANS is last channel in table before pumps, membr params */
	if (uden) mult = 1;
	else      mult = 1000 * qcond(chanunit[ch]);
	printf ("%-8.4g ",celdens[ct][dn][ch][region] * mult);
      }
      else
	printf ("%-8.4g ",celdens[ct][dn][ch][region]);
    };
    // printf (" %s\n",chnamea[ch]);
    printf ("\n");
  };
  printf ("# \n");
}

/*---------------------------------------------------------*/

/* Find all the membrane property files. */
/* If any have not been defined, try to use the default file. */
/* If that is not found, use "zerodens[][]" which has zero   */
/* for all the channel densities. */

void findceldens() 

{
#define FNSIZ 120
   int ct, dn, ch, r;
   char filnam[FNSIZ];
   double *cdenssav = NULL;

 celseg_init();

 for (ct=0; ct<nceltypes; ct++) {
   for (dn=0; dn<NDENS; dn++) {
    const char *fname = NULL;
    char *cfname = NULL;
    FILE *x = NULL;

  /* density files from command line */

    fname = NULL;
    if (dn==0) {
    if      ((ct==xcone)&& (!notinit(cone_densfile))) fname = cone_densfile;
    else if ((ct==xrod) && (!notinit(rod_densfile)))  fname = rod_densfile;
    else if ((ct==hbat) && (!notinit(hbat_densfile))) fname = hbat_densfile;
    else if ((ct==ha)   && (!notinit(ha_densfile)))   fname = ha_densfile;
    else if ((ct==hb)   && (!notinit(hb_densfile)))   fname = hb_densfile;
    else if ((ct==rbp)  && (!notinit(rbp_densfile)))  fname = rbp_densfile;
    else if ((ct==dbp1) && (!notinit(dbp1_densfile))) fname = dbp1_densfile;
    else if ((ct==dbp2) && (!notinit(dbp2_densfile))) fname = dbp2_densfile;
    else if ((ct==dbp3) && (!notinit(dbp3_densfile))) fname = dbp3_densfile;
    else if ((ct==dbp4) && (!notinit(dbp4_densfile))) fname = dbp4_densfile;
    else if ((ct==hbp1) && (!notinit(hbp1_densfile))) fname = hbp1_densfile;
    else if ((ct==hbp2) && (!notinit(hbp2_densfile))) fname = hbp2_densfile;
    else if ((ct==a17)  && (!notinit(a17_densfile)))  fname = a17_densfile;
    else if ((ct==aii)  && (!notinit(aii_densfile)))  fname = aii_densfile;
    else if ((ct==sbac) && (!notinit(sbac_densfile))) fname = sbac_densfile;
    else if ((ct==am)   && (!notinit(am_densfile)))   fname = am_densfile;
    else if ((ct==am2)  && (!notinit(am2_densfile)))  fname = am2_densfile;
    else if ((ct==am3)  && (!notinit(am3_densfile)))  fname = am3_densfile;
    else if ((ct==am4)  && (!notinit(am4_densfile)))  fname = am4_densfile;
    else if ((ct==amh)  && (!notinit(amh_densfile)))  fname = amh_densfile;
    else if ((ct==amh2) && (!notinit(amh2_densfile))) fname = amh2_densfile;
    else if ((ct==ams)  && (!notinit(ams_densfile)))  fname = ams_densfile;
    else if ((ct==amhs) && (!notinit(amhs_densfile))) fname = amhs_densfile;
    else if ((ct==gca)  && (!notinit(gca_densfile)))  fname = gca_densfile;
    else if ((ct==gcb)  && (!notinit(gcb_densfile)))  fname = gcb_densfile;
    else if ((ct==dsgc) && (!notinit(dsgc_densfile))) fname = dsgc_densfile;
    else if ((ct==gcaoff)&& (!notinit(gcaoff_densfile)))fname = gcaoff_densfile;
    else if ((ct==gcboff)&& (!notinit(gcboff_densfile)))fname = gcboff_densfile;
    }
    else if (dn==1) {
    if      ((ct==xcone)&& (!notinit(cone_densfile2))) fname = cone_densfile2;
    else if ((ct==xrod) && (!notinit(rod_densfile2)))  fname = rod_densfile2;
    else if ((ct==hbat) && (!notinit(hbat_densfile2))) fname = hbat_densfile2;
    else if ((ct==ha)   && (!notinit(ha_densfile2)))   fname = ha_densfile2;
    else if ((ct==hb)   && (!notinit(hb_densfile2)))   fname = hb_densfile2;
    else if ((ct==rbp)  && (!notinit(rbp_densfile2)))  fname = rbp_densfile2;
    else if ((ct==dbp1) && (!notinit(dbp1_densfile2))) fname = dbp1_densfile2;
    else if ((ct==dbp2) && (!notinit(dbp2_densfile2))) fname = dbp2_densfile2;
    else if ((ct==dbp3) && (!notinit(dbp3_densfile2))) fname = dbp3_densfile2;
    else if ((ct==dbp4) && (!notinit(dbp4_densfile2))) fname = dbp4_densfile2;
    else if ((ct==hbp1) && (!notinit(hbp1_densfile2))) fname = hbp1_densfile2;
    else if ((ct==hbp2) && (!notinit(hbp2_densfile2))) fname = hbp2_densfile2;
    else if ((ct==a17)  && (!notinit(a17_densfile2)))  fname = a17_densfile2;
    else if ((ct==aii)  && (!notinit(aii_densfile2)))  fname = aii_densfile2;
    else if ((ct==sbac) && (!notinit(sbac_densfile2))) fname = sbac_densfile2;
    else if ((ct==am)   && (!notinit(am_densfile2)))   fname = am_densfile2;
    else if ((ct==am2)  && (!notinit(am2_densfile2)))  fname = am2_densfile2;
    else if ((ct==am3)  && (!notinit(am3_densfile2)))  fname = am3_densfile2;
    else if ((ct==am4)  && (!notinit(am4_densfile2)))  fname = am4_densfile2;
    else if ((ct==amh)  && (!notinit(amh_densfile2)))  fname = amh_densfile2;
    else if ((ct==amh2) && (!notinit(amh2_densfile2))) fname = amh2_densfile2;
    else if ((ct==ams)  && (!notinit(ams_densfile2)))  fname = ams_densfile2;
    else if ((ct==amhs) && (!notinit(amhs_densfile2))) fname = amhs_densfile2;
    else if ((ct==gca)  && (!notinit(gca_densfile2)))  fname = gca_densfile2;
    else if ((ct==gcb)  && (!notinit(gcb_densfile2)))  fname = gcb_densfile2;
    else if ((ct==dsgc) && (!notinit(dsgc_densfile2))) fname = dsgc_densfile2;
    else if ((ct==gcaoff)&& (!notinit(gcaoff_densfile2)))fname = gcaoff_densfile2;
    else if ((ct==gcboff)&& (!notinit(gcboff_densfile2)))fname = gcboff_densfile2;
    }
   if (getn(ct,MAKE) && getn(ct,BIOPHYS)) {
       char *tfilnam = NULL;

    sprintf(filnam, "dens_%s.n", cname[ct]);
    cfname = emalloc(min(FNSIZ,strlen(filnam)) + 1);
    if (fname == NULL) {			// if user fname isn't given, use cell type
      if (dn==1) fname = densfil[ct][0];	// use densfile if densfile2 not defined
      else fname = filnam;
    }
    strcpy(cfname,filnam);
    tfilnam = emalloc(min(FNSIZ,strlen(fname)) + 1);
    sprintf(tfilnam, "%.50s", fname);
    sprintf(filnam, "%.50s/%.50s", confdir, tfilnam);
    if ((x=fopen(filnam,"r"))!=NULL) {             /* check for file with density data */
      fclose(x);
      densfil[ct][dn] = tfilnam;
    } else {				   /* check for default density file */
      if (ninfo>=2) {
	   ncfprintf(stderr,"# findceldens(): can't find density file %s specified for %s, ", filnam,cname[ct]);
           ncfprintf(stderr,"using standard density file %s\n",cfname);
      }
      if (tfilnam!=NULL) efree(tfilnam);
      tfilnam = emalloc(min(FILNAMSIZ,strlen(cfname)) + 1);
      sprintf(tfilnam, "%.50s", cfname);
      sprintf(filnam, "%.50s/%.50s", confdir, tfilnam);
      if ((x=fopen(filnam,"r"))!=NULL) {             /* check for file with density data */
        fclose(x);
        densfil[ct][dn] = tfilnam;
      } else {				   /* check for default density file */
        if (ninfo>=2) {
          ncfprintf(stderr,"# findceldens(): can't find standard density file %s, ", filnam);
          ncfprintf(stderr,"Using default density file %s\n",def_densfile);
	}
        if (tfilnam!=NULL) efree(tfilnam);
        tfilnam = emalloc(min(FILNAMSIZ,strlen(def_densfile)) + 1);
        sprintf(tfilnam, "%.50s", def_densfile);
        sprintf(filnam, "%.50s/%.50s", confdir, tfilnam);
        if ((x = fopen(filnam, "r"))!=NULL) { 	  /* check for default file with density data */
	  fclose(x);
	  densfil[ct][dn] = tfilnam;
        } 
        else {
          if (ninfo>=2) {
            ncfprintf(stderr,"# findceldens(): can't find standard density file %s, ", filnam);
            ncfprintf(stderr,"Using zero density file %s\n",zerodens);
	  }
	  densfil[ct][dn] = NULL;			 /* otherwise use zerodens */
          if (tfilnam!=NULL) efree(tfilnam);
        }
      }
    }
    fname = NULL;
    if (ninfo>=2) ncfprintf(stderr,"#\n");
   }
  }  /* for (dn;;) */
 }  /* for (ct;;) */


 for (ct=0; ct<nceltypes; ct++) { 	/* read density files */
   for (dn=0; dn<NDENS; dn++) {
    if (getn(ct,MAKE) && getn(ct,BIOPHYS)) {
        int nrows;
        int ncols;
        int chindex, chantyp, regindex;

      if (ninfo>=3) {
	if (densfil[ct][dn])
	  ncfprintf (stderr,"#   celseg.cc: %5s dens arr %d channel densities from '%s'.\n",
		  cname[ct], dn, densfil[ct][dn]);
	else ncfprintf (stderr,
	"#   celseg.cc: %5s biophysics specified, but channel density file not present.\n",cname[ct]);
      }
      if (!densfil[ct][dn]) {
        cdens = NULL;		// no cell density file
        if (ninfo>=3) {
          ncfprintf (stderr,"#   celseg.cc, densities for '%s' dens arr %d not specified, using '%s'.\n",
	       cname[ct],dn,"zerodens");
        }
        nrows = NCHANS;
        ncols = R_NREGIONS;
      }
      else {
        sprintf(filnam, "%.50s/%.50s", confdir, densfil[ct][dn]);
 	cdens = fread(filnam, &nrows, &ncols);	/* read density file into temp arr */
      } 
        // Set up density array
        // Start by copying zerodens as default values

      for (ch=0; ch<NCHANS; ch++) {		/* copy density data from zerodens as default  */
        for (r=0; r<NREGIONS; r++) {		/*   in case row or column is skipped below */
          celdens[ct][dn][ch][r] = zerodens[ch][r];
        }
      }

	// Density file has variable format, nrows, ncols are set
	//    by index labels at top or left column
        // Skip values with indexes out of range

       if (cdens!=NULL) {
        cdenssav = cdens;
        for (ch=1; ch<nrows; ch++) {		/* copy density data  */
         chindex = ch*ncols;			/* channel index is first col (0) */
         chantyp = *(cdens + chindex);		/* get channel type (col 0) */
         if (chantyp<0 || chantyp>NCHANS) continue;  /* ignore if channel type isn't set */
         for (r=1; r<ncols; r++) {		/* skip first col = label */
           regindex = *(cdens+r);			/* region index is first row */
           if (regindex<0 || regindex>NREGIONS) continue;  /* ignore if region isn't set */
           celdens[ct][dn][chantyp][regindex] = *(cdens+chindex+r);
         }
        }
        if (cdenssav!=NULL && cdenssav != &zerodens[0][0]) {
	     efree (cdenssav);
        }
       }
     } /* if (MAKE && BIOPHYS) */
   } /* for (dn;;) */
 } /* for (ct;;) */
 
 if (ninfo>=3) ncfprintf (stderr,"\n");
}

/*---------------------------------------------------------*/

void modchandens()

/* modify channel densities */

{
   int ct, ch, dn;
   int region;

 for (ct=0; ct<nceltypes; ct++) {
   for (dn=0; dn<NDENS; dn++) {

   /* possibly modify channel densities */

  /* if RATIOK, set undetermined values */

   if (getn(ct,RATIOK)==1) {

    /* use ratios for Ka and Kdr where not already specified */
     for (region=0; region<R_NREGIONS; region ++) {
       
       if (celdens[ct][dn][_KDR][region]==0)
	   celdens[ct][dn][_KDR][region] = ratio_kdr * celdens[ct][dn][_NA] [region];
       if (celdens[ct][dn][_KA][region]==0)
	   celdens[ct][dn][_KA][region]  = ratio_ka  * celdens[ct][dn][_NA][region];
       if ((region != AXON) && (region != AXON_DIST)) {
	 if (celdens[ct][dn][_SKCA1][region]==0)
	     celdens[ct][dn][_SKCA1][region]  = ratio_sk1  * celdens[ct][dn][_NA][region];
	 if (celdens[ct][dn][_SKCA2][region]==0)
	     celdens[ct][dn][_SKCA2][region]  = ratio_sk2  * celdens[ct][dn][_NA][region];
       }
     }  /* region */
   } /* ratiok == 1 */

 /* possibly block channels */

  for (region=0; region<R_NREGIONS; region ++) {
    celdens[ct][dn][C_NA0] [region] *= (1 - ttxbath);
    celdens[ct][dn][C_NA1] [region] *= (1 - ttxbath);
    celdens[ct][dn][C_NA2] [region] *= (1 - ttxbath);
    celdens[ct][dn][C_NA3] [region] *= (1 - ttxbath);
    celdens[ct][dn][C_NA4] [region] *= (1 - ttxbath);
    celdens[ct][dn][C_NA5][region]  *= (1 - ttxbath);
    celdens[ct][dn][C_NA6][region]  *= (1 - ttxbath);
    if (region==SOMA) {
      celdens[ct][dn][C_NA0] [region]  *= (1 - ttxsoma);
      celdens[ct][dn][C_NA1] [region]  *= (1 - ttxsoma);
      celdens[ct][dn][C_NA2] [region]  *= (1 - ttxsoma);
      celdens[ct][dn][C_NA3] [region]  *= (1 - ttxsoma);
      celdens[ct][dn][C_NA4] [region]  *= (1 - ttxsoma);
      celdens[ct][dn][C_NA5][region]  *= (1 - ttxsoma);
      celdens[ct][dn][C_NA6][region]  *= (1 - ttxsoma);
    }
    if (region==DENDD || region==DEND || region==DENDP) {
      celdens[ct][dn][C_NA0] [region]  *= (1 - ttxdend);
      celdens[ct][dn][C_NA1] [region]  *= (1 - ttxdend);
      celdens[ct][dn][C_NA2] [region]  *= (1 - ttxdend);
      celdens[ct][dn][C_NA3] [region]  *= (1 - ttxdend);
      celdens[ct][dn][C_NA4] [region]  *= (1 - ttxdend);
      celdens[ct][dn][C_NA5][region]  *= (1 - ttxdend);
      celdens[ct][dn][C_NA6][region]  *= (1 - ttxdend);
    }
    celdens[ct][dn][C_K0][region] *= (1 - tea);
    celdens[ct][dn][C_K1][region] *= (1 - tea);
    celdens[ct][dn][C_K2][region] *= (1 - tea);
    celdens[ct][dn][C_K3][region]  *= (1 - fourap);			// KA
    celdens[ct][dn][C_K6][region]  *= (1 - tea);
    celdens[ct][dn][C_K7][region]  *= (1 - tea);

    celdens[ct][dn][C_K4][region]  *= (1 - zd7288);

    // joesterle begin
    celdens[ct][dn][C_HCN1][region]  *= (1 - zd7288);
    celdens[ct][dn][C_HCN2][region]  *= (1 - zd7288);
    celdens[ct][dn][C_HCN3][region]  *= (1 - zd7288);
    celdens[ct][dn][C_HCN4][region]  *= (1 - zd7288);
    // joesterle end
    
    celdens[ct][dn][C_KCA0][region]  *= (1 - ibtox);
    celdens[ct][dn][C_KCA2][region]  *= (1 - ibtox);
    celdens[ct][dn][C_KCA3][region]  *= (1 - ibtox);
  }

// /* When densities in sub-regions are zero, fill them in using
//     densities of neighboring regions.  This allows user
//     to override densities for special cases.
//     Extra sub-regions are not "extra parameters" since by default
//     they are a function of neighboring regions. 
// */
// 
//  
//   int ch;
//  for (ch=0; ch<NCHANS; ch++) {
//    if (ch<=ECHANS) {   			/* ECHANS is last channel in table */
//      if (celdens[ct][ch][R_DEND_PROX]==0)
//	celdens[ct][ch][R_DEND_PROX] = (celdens[ct][dn][ch][R_DEND] +
//				      celdens[ct][ch][dn][SOMA]) * .5;
//      if (celdens[ct][ch][AXON_DIST]==0)
//	celdens[ct][ch][AXON_DIST] =  celdens[ct][ch][dn][AXON];
//    }
//  }

  /* Finally, convert to ndensity (n/um2) from mS/cm2 */

    for (region=0; region<R_NREGIONS; region++) {
      for (ch=0; ch<NCHANS; ch++) 
        if (ch<=ECHANS) {		    /* ECHANS is last channel in table */
	double qc = qcond(chanunit[ch]);
	if (qc > 0) 
	  celdens[ct][dn][ch][region] = celdens[ct][dn][ch][region] / qc;
       }
     }
   } /* for (nd;;)  */
 } /* for (ct;;)  possibly modify channel densities */

} /* modchandens */

/*------------------------------------------------------*/


void printchaninfo() {

  double udens, prden;

  int ct,dn;
  for (ct=0; ct<nceltypes; ct++) {

    prden =  (info_chan) && !disp && getn(ct,BIOPHYS) && getn(ct,MAKE);

    if (prden) {
      printf ("#\n");
      printf ("#   Unitary channel conductances:\n");
      printf ("#\n");
      double qc = exp(log(dqc)*(tempcel-dbasetc)/10)  * 1e12;
      printf ("#   %-5.3g pS   %s\n",dnau  * qc,chname[_NA]);
      printf ("#   %-5.3g pS   %s\n",dnau  * qc,chname[_NA5]);
      printf ("#   %-5.3g pS   %s\n",dnau  * qc,chname[_NA6]);
      printf ("#   %-5.3g pS   %s\n",dnau  * qc,chname[_NA8]);
      printf ("#   %-5.3g pS   %s\n",dku   * qc,chname[_KDR]);
      printf ("#   %-5.3g pS   %s\n",dkau  * qc,chname[_KA]);
      printf ("#   %-5.3g pS   %s\n",dkihu * qc,chname[_KH]);
      printf ("#   %-5.3g pS   %s\n",dkihu * qc,chname[_KHZ]);
      printf ("#   %-5.3g pS   %s\n",dkihu * qc,chname[_KV3]);
      printf ("#   %-5.3g pS   %s\n",dkihu * qc,chname[_K7]);
      printf ("#   %-5.3g pS   %s\n",dkcasu* qc,chname[_SKCA1]);
      printf ("#   %-5.3g pS   %s\n",dkcasu* qc,chname[_SKCA2]);
      printf ("#   %-5.3g pS   %s\n",dkcabu* qc,chname[_BKCA]);
      printf ("#   %-5.3g pS   %s\n",dkcabu* qc,chname[_CLCA1]);
      printf ("#   %-5.3g pS   %s\n",dkcabu* qc,chname[_CLCA2]);
      printf ("#   %-5.3g pS   %s\n",dcalu * qc,chname[_CA]);
      printf ("#   %-5.3g pS   %s\n",dcalu * qc,chname[_CA5]);
      printf ("#   %-5.3g pS   %s\n",dcalu * qc,chname[_CA6]);
      printf ("#   %-5.3g pS   %s\n",dcalu * qc,chname[_CA7]);
      printf ("#   %-5.3g pS   %s\n",dampau * qc,chname[_AMPA1]);
      printf ("#   %-5.3g pS   %s\n",dampau * qc,chname[_AMPA2]);
      printf ("#   %-5.3g pS   %s\n",dampau * qc,chname[_AMPA3]);
      printf ("#   %-5.3g pS   %s\n",dampau * qc,chname[_AMPA4]);
      printf ("#   %-5.3g pS   %s\n",dampau * qc,chname[_AMPA5]);
      printf ("#   %-5.3g pS   %s\n",dgabau * qc,chname[_GABA1]);
      printf ("#   %-5.3g pS   %s\n",dgabau * qc,chname[_GABA2]);
      printf ("#   %-5.3g pS   %s\n",dgabau * qc,chname[_GABA3]);
      printf ("#   %-5.3g pS   %s\n",dgabau * qc,chname[_GABA4]);
      printf ("#\n#\n");
      if (getn(ct,RATIOK)==1) {
        printf ("#   Setting K chans for '%s' by ratio from Na density\n", cname[ct]);
        printf ("# \n");
        printf ("#   ratio_kdr = %g\n", ratio_kdr);
        printf ("#   ratio_ka  = %g\n", ratio_ka);
        printf ("#   ratio_sk1 = %g\n", ratio_sk1);
        printf ("#   ratio_sk2 = %g\n", ratio_sk2);
        printf ("#   ratio_bk  = %g\n", ratio_bk);
      };
      printdens(ct,dn=0,"mS/cm2",udens=0);
      printdens(ct,dn=1,"mS/cm2",udens=0);
      printdens(ct,dn=0,"N/um2",udens=1);
      printdens(ct,dn=1,"N/um2",udens=1);
      print_chan_offset(ct);
    };
  }; /* for (ct;;) */
};   /* proc printchaninfo */

/*------------------------------------------------------*/

void set_chan_offset(int ct, int chvalIndex, int chrate, double val)

{
   chval[ct][chvalIndex][chrate] = val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void set_chan_offsets(chattrib *a, int ct, int chvalIndex)

{
   // fprintf (stderr, "# set_chan_offsets ct = %d, chvalIndex= %d\n", ct, chvalIndex);
   a->voffsm   = chval[ct][chvalIndex][CHOFFM];
   a->voffsh   = chval[ct][chvalIndex][CHOFFH];
   a->taua     = chval[ct][chvalIndex][CHTAUA];
   a->taub     = chval[ct][chvalIndex][CHTAUB];
   a->tauc     = chval[ct][chvalIndex][CHTAUC];
   a->taud     = chval[ct][chvalIndex][CHTAUD];
}

/*------------------------------------------------------*/

void addcicr(cattrib *ca, int ct, int cn, int region) {

/* add CICR and IP3 parameters */

  if ((ca->vm2 = celdens[ct][ndens[ct][cn]][C_VM2][region]) > 0) {		/* if CICR */

        double cas, vm3, c1cicr, kfcicr, kacicr, krcicr, k1cicr, k2cicr, ncicr, mcicr, pcicr;

	ca->cicr = 1;
	if ((c1cicr = celdens[ct][ndens[ct][cn]][C_C1CICR][region]) <=0 || 
	    (c1cicr = celdens[ct][ndens[ct][cn]][C_C1CICR][region]) >=1) c1cicr = c1CICR;
	ca->c1cicr = c1cicr;

	if ((cas    = celdens[ct][ndens[ct][cn]][C_CAS][region]) <=0)	    cas = CASSTART;
	ca->cas = cas;
	if ((vm3    = celdens[ct][ndens[ct][cn]][C_VM3][region]) <=0)	    vm3 = vm3CICR;
	ca->vm3 = vm3;
	if ((kfcicr = celdens[ct][ndens[ct][cn]][C_KFCICR][region]) <=0) kfcicr = kfCICR;
	ca->kfcicr = kfcicr;
	if ((kacicr = celdens[ct][ndens[ct][cn]][C_KACICR][region]) <=0) kacicr = kaCICR;
	ca->kacicr = kacicr;
	if ((krcicr = celdens[ct][ndens[ct][cn]][C_KRCICR][region]) <=0) krcicr = krCICR;
	ca->krcicr = krcicr;
	if ((k1cicr = celdens[ct][ndens[ct][cn]][C_K1CICR][region]) <=0) k1cicr = k1CICR;
	ca->k1cicr = k1cicr;
	if ((k2cicr = celdens[ct][ndens[ct][cn]][C_K2CICR][region]) <=0) k2cicr = k2CICR;
	ca->k2cicr = k2cicr;
		// remove for now
//		if ((ncicr  = celdens[ct][ndens[ct][cn]][C_NCICR][region])  <=0) ncicr  = nCICR;
		ca->ncicr = nCICR;
//		if ((mcicr  = celdens[ct][ndens[ct][cn]][C_MCICR][region])  <=0) mcicr  = mCICR;
		ca->mcicr = mCICR;
//		if ((pcicr  = celdens[ct][ndens[ct][cn]][C_PCICR][region])  <=0) pcicr  = pCICR;
		ca->pcicr = pCICR;

	if ((ca->bip3 = celdens[ct][ndens[ct][cn]][C_BIP3][region]) > 0) {  /* if constant IP3 with CICR */
	double cas2, ip3i, vip3; 

	if ((cas2   = celdens[ct][ndens[ct][cn]][C_CAS2][region]) <=0)	cas2 = CAS2START;
	      ca->cas2 = cas2;
	      if ((ip3i   = celdens[ct][ndens[ct][cn]][C_IP3I][region]) <=0 )	ip3i = IP3ISTART;
	      ca->ip3i = ip3i;
	      if ((vip3   = celdens[ct][ndens[ct][cn]][C_VIP3][region]) <=0)	vip3 = v1IP3;
	      ca->vip3 = vip3;
	}
  }
  if ((ca->v3ip3 = celdens[ct][ndens[ct][cn]][C_V3IP3][region]) > 0) {  /* if dynamic IP3 */

        double cas2, ip3i, v3ip3,v2ip3,v4ip3,d1ip3,d2ip3,d3ip3,d4ip3,a2ip3,a3ip3,b2ip3,k3ip3,c1cicr;

	ca->cicr = 1;
	if ((c1cicr = celdens[ct][ndens[ct][cn]][C_C1CICR][region]) <=0 || 
	    (c1cicr = celdens[ct][ndens[ct][cn]][C_C1CICR][region]) >= 1) c1cicr = c1CICR;
	ca->c1cicr = c1cicr;
       	if ((cas2 = celdens[ct][ndens[ct][cn]][C_CAS2][region]) <=0) cas2 = CAS2START;
       	ca->cas2  = cas2;
       	if ((ip3i = celdens[ct][ndens[ct][cn]][C_IP3I][region]) <=0) ip3i = IP3ISTART;
       	ca->ip3i  = ip3i;
       	if ((v2ip3 = celdens[ct][ndens[ct][cn]][C_V2IP3][region]) <=0) v2ip3 = v2IP3;
       	ca->v2ip3 = v2ip3;
       	if ((v3ip3 = celdens[ct][ndens[ct][cn]][C_V3IP3][region]) <=0) v3ip3 = v3IP3;
       	ca->v3ip3 = v3ip3;
       	if ((v4ip3 = celdens[ct][ndens[ct][cn]][C_V4IP3][region]) <=0) v4ip3 = v4IP3;
       	ca->v4ip3 = v4ip3;

       	//override constant if dynamic rate v3ip3 set
       	ca->bip3 = 0;

       	if ((d1ip3 = celdens[ct][ndens[ct][cn]][C_D1IP3][region]) <=0)  d1ip3 = d1IP3;
       	ca->d1ip3 = d1ip3;
       	if ((d2ip3 = celdens[ct][ndens[ct][cn]][C_D2IP3][region]) <=0)  d2ip3 = d2IP3;
       	ca->d2ip3 = d2ip3;
       	if ((d3ip3 = celdens[ct][ndens[ct][cn]][C_D3IP3][region]) <=0)  d3ip3 = d3IP3;
       	ca->d3ip3 = d3ip3;
       	if ((d4ip3 = celdens[ct][ndens[ct][cn]][C_D4IP3][region]) <=0)  d4ip3 = d4IP3;
       	ca->d4ip3 = d4ip3;
       	if ((a2ip3 = celdens[ct][ndens[ct][cn]][C_A2IP3][region]) <=0)  a2ip3 = a2IP3;
       	ca->a2ip3 = a2ip3;
       	if ((a3ip3 = celdens[ct][ndens[ct][cn]][C_A3IP3][region]) <=0)  a3ip3 = a3IP3;
       	ca->a3ip3 = a3ip3;

       	if ((b2ip3 = celdens[ct][ndens[ct][cn]][C_B2IP3][region]) <=0)  b2ip3 = b2IP3;
       	ca->b2ip3 = b2ip3;
       	if ((k3ip3 = celdens[ct][ndens[ct][cn]][C_K3IP3][region]) <=0)  k3ip3 = k3IP3;
       	ca->k3ip3 = k3ip3;

       	ca->mtypeip3 = mtypeIP3;
       	ca->oip3 = oIP3;
  }
}

/*------------------------------------------------------*/

void addChannelToElem(int ct, int cn, int region, elem* e, int chanType, int chanSubType,
				      int chvalIndex, int chanNoise, double somadist, bool makeEname, 
					bool makeCattribs = false)
{	
       double egrad, lgrad, gradmul, gradadd, ndensity;

	if (EGRAD >= 0) {					    // if EGRAD column is set 
    egrad = celdens[ct][ndens[ct][cn]][chvalIndex][EGRAD];    // get expon grad multiplier for channel
	  if (egrad == 0) gradmul = 1;
	  else gradmul = exp (somadist / egrad);  	            // multiply by distance to soma
	}
	else gradmul = 1;

	if (LGRAD >= 0) {					    // if LGRAD column is set 
    lgrad = celdens[ct][ndens[ct][cn]][chvalIndex][LGRAD];    // get grad multiplier for channel
    gradadd = lgrad * somadist;		  	            // multiply by distance to soma
	}
	else gradadd = 0;

	ndensity = celdens[ct][ndens[ct][cn]][chvalIndex][region] * gradmul + gradadd;
	if (ndensity > 0 || e->attpnt==NULL) {
	   chattrib* a = make_chan(e, chanType, chanSubType);
	   a->ndensity = ndensity;
	   a->voffsm   = chval[ct][chvalIndex][CHOFFM];
	   a->voffsh   = chval[ct][chvalIndex][CHOFFH];
	   a->taua     = chval[ct][chvalIndex][CHTAUA];
	   a->taub     = chval[ct][chvalIndex][CHTAUB];
	   a->tauc     = chval[ct][chvalIndex][CHTAUC];
	   a->taud     = chval[ct][chvalIndex][CHTAUD];
	   if(chanNoise)  nattrib* n = make_chnoise(e);

	   if (makeCattribs && (chanType == CA)) {				// don't use, makeCattribs = F
		/* make capump if Ca chan is defined */
		      double pkm, bmax, bkd, btot, btoti, cashell;
		cattrib* ca = (cattrib*) a;

	   if ((ca->vmax = celdens[ct][ndens[ct][cn]][C_CAPUMP][region]) > 0) {
	      ca->pump = 1;
	       if ((pkm = celdens[ct][ndens[ct][cn]][C_CAPKM] [region]) <= 0) pkm = dcapkm;
	            ca->pkm = pkm;
	   }
 	   if ((ca->btot = celdens[ct][ndens[ct][cn]][C_CABTOT][region]) > 0) { 	// set btot = dcabt in dens file
	  	ca->cabuf = 1;
		if ((bmax = celdens[ct][ndens[ct][cn]][C_CABVMAX][region]) <=0) bmax = dcabf;
		ca->bmax = bmax;
		if ((bkd = celdens[ct][ndens[ct][cn]][C_CABKD][region]) <=0) bkd = dcabr/dcabf;
		ca->bkd = bkd;
		if ((btoti = celdens[ct][ndens[ct][cn]][C_CABTOTI][region]) <=0) btoti = ca->btot; // same as others
		ca->btoti = btoti;
	   }
	   if ((ca->kex = celdens[ct][ndens[ct][cn]][C_CAEXCH][region]) > 0) {
	       ca->exch = 1;
	   }
	   if ((cashell = celdens[ct][ndens[ct][cn]][C_CASHELL][region]) <=0) cashell = nshell;
	   ca->cashell = cashell;

	   addcicr(ca, ct, cn, region);		/* add CICR, IP3 parameters */
        }
     }
	// if (makeEname) ename(e, &chanelem[ct][chvalIndex]);
}


void addChannelsToSphere(int ct, int cn, int nod1, int region, double somadist)
{
	int i;
	double cadens;
	elem *e;
	node *npnt;

	npnt = ndn(ct,cn,nod1);

	e = at(npnt, CHAN);		/* Na v1.2 */
	addChannelToElem(ct, cn, region, e, NA, 0, C_NA0, 0, somadist, true);		
	e = at(npnt, CHAN);		/* Na v1.2 */
	addChannelToElem(ct, cn, region, e, NA, 1, C_NA1, na1nois, somadist, true);		
	e = at(npnt, CHAN);		/* Na v1.2 */
	addChannelToElem(ct, cn, region, e, NA, 2, C_NA2, 0, somadist, true);		
	e = at(npnt, CHAN);		/* Na v1.1 */
	addChannelToElem(ct, cn, region, e, NA, 5, C_NA5, na5nois, somadist, true);	
	e = at(npnt, CHAN);		/* Na v1.6, recurrent */
	addChannelToElem(ct, cn, region, e, NA, 6, C_NA6, na6nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, NA, 8, C_NA8, na8nois, somadist, true);	
	e = at(npnt, CHAN);			/* Kdr */
	addChannelToElem(ct, cn, region, e, K, 0, C_K0, 0, somadist, true);	
	e = at(npnt, CHAN);			/* Kdr */
	addChannelToElem(ct, cn, region, e, K, 1, C_K1, k1nois, somadist, true);	
	e = at(npnt, CHAN);			/* KA */
	addChannelToElem(ct, cn, region, e, K, 2, C_K2, k1nois, somadist, true);	
	e = at(npnt, CHAN);			/* KA */
	addChannelToElem(ct, cn, region, e, K, 3, C_K3, k3nois, somadist, true);	
	e = at(npnt, CHAN);			/* KH */
	addChannelToElem(ct, cn, region, e, K, 4, C_K4, k4nois, somadist, true);	
  
  // begin joesterle
  e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, K, 11, C_HCN1, hcn1nois, somadist, true);
  e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, K, 12, C_HCN2, hcn2nois, somadist, true);	
  e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, K, 13, C_HCN3, hcn3nois, somadist, true);	
  e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, K, 14, C_HCN4, hcn4nois, somadist, true);	
  // end joesterle
  
	e = at(npnt, CHAN);			/* Kir*/
	addChannelToElem(ct, cn, region, e, K, 5, C_K5, k5nois, somadist, true);	
	e = at(npnt, CHAN);			/* Kv3*/
	addChannelToElem(ct, cn, region, e, K, 6, C_K6, k6nois, somadist, true);	
	e = at(npnt, CHAN);			/* Kv3*/
	addChannelToElem(ct, cn, region, e, K, 7, C_K7, k7nois, somadist, true);	
	e = at(npnt, CHAN);			/* KCa */
	addChannelToElem(ct, cn, region, e, KCa, 3, C_KCA3, kca3nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, KCa, 4, C_KCA4, kca4nois, somadist, true);	
	e = at(npnt, CHAN);			
	addChannelToElem(ct, cn, region, e, KCa, 5, C_KCA5, kca5nois, somadist, true);
	e = at(npnt, CHAN);			
	addChannelToElem(ct, cn, region, e, ClCa, 1, C_CLCA1, clcanois, somadist, true);
	e = at(npnt, CHAN);			
	addChannelToElem(ct, cn, region, e, ClCa, 2, C_CLCA2, clcanois, somadist, true);
	e = at(npnt, CHAN);			/* Ca */
	addChannelToElem(ct, cn, region, e, CA, 0, C_CA0, ca1nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 1, C_CA1, ca1nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 2, C_CA2, ca1nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 3, C_CA3, ca3nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 4, C_CA4, ca1nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 5, C_CA5, ca5nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 6, C_CA6, ca6nois, somadist, true);	
	e = at(npnt, CHAN);
	addChannelToElem(ct, cn, region, e, CA, 7, C_CA7, ca7nois, somadist, true);	

	// check for Ca sensitivity in any of the channels in this region
	//
	for (cadens=0,i=0; i<NCHANS; i++) {
	    if (i>=C_KCA0 && i<C_CAPUMP) {
              cadens += celdens[ct][ndens[ct][cn]][i][region];
	    }
	}
	if (cadens > 0) { // If any channels need a Ca comp, make one
			  // (AMPA channels make cacomp separately)
		 double pkm, bmax, bkd, btot, btoti, cashell;

	  e = at(npnt, CACOMP);		
	  cattrib* ca = (cattrib *) make_chan(e, CACOMP, 0);
	  if ((ca->vmax = celdens[ct][ndens[ct][cn]][C_CAPUMP][region]) > 0) {

		ca->pump = 1;
		if ((pkm = celdens[ct][ndens[ct][cn]][C_CAPKM] [region]) <= 0) pkm = dcapkm;
		ca->pkm = pkm;
	  }
	  if ((ca->btot = celdens[ct][ndens[ct][cn]][C_CABTOT][region]) > 0) { 		// set btot = dcabt in dens file
	  	ca->cabuf = 1;
		if ((bmax = celdens[ct][ndens[ct][cn]][C_CABVMAX][region]) <=0) bmax = dcabf;
		ca->bmax = bmax;
		if ((bkd = celdens[ct][ndens[ct][cn]][C_CABKD][region]) <=0) bkd = dcabr/dcabf;
		ca->bkd = bkd;
		if ((btoti = celdens[ct][ndens[ct][cn]][C_CABTOTI][region]) <=0) btoti = ca->btot; // same as others
		ca->btoti = btoti;
	  }
	  if ((ca->kex = celdens[ct][ndens[ct][cn]][C_CAEXCH][region]) > 0) {
	        ca->exch = 1;
	  }
	  if ((cashell = celdens[ct][ndens[ct][cn]][C_CASHELL][region]) <=0) cashell = nshell_soma;
	  ca->cashell = cashell;

	  addcicr(ca, ct, cn,  region);		/* add CICR, IP3 parameters */
     }
}

void addChannelsToCable(int ct, int cn, elem* e, int region, double somadist)
{
	int i;
	double cadens;

	addChannelToElem(ct, cn, region, e, NA, 0, C_NA0, 0, somadist, false);
	addChannelToElem(ct, cn, region, e, NA, 1, C_NA1, na1nois, somadist, false);
	addChannelToElem(ct, cn, region, e, NA, 2, C_NA2, na2nois, somadist, false);
	addChannelToElem(ct, cn, region, e, NA, 5, C_NA5, na5nois, somadist, false);
	addChannelToElem(ct, cn, region, e, NA, 6, C_NA6, na6nois, somadist, false);
	addChannelToElem(ct, cn, region, e, NA, 8, C_NA8, na8nois, somadist, false);
	addChannelToElem(ct, cn, region, e, K,  0, C_K0, 0,       somadist, false);
	addChannelToElem(ct, cn, region, e, K,  1, C_K1, k1nois, somadist, false);
	addChannelToElem(ct, cn, region, e, K,  2, C_K2,  0,      somadist, false);
	addChannelToElem(ct, cn, region, e, K,  3, C_K3,  k3nois, somadist, false);
	addChannelToElem(ct, cn, region, e, K,  4, C_K4,  k4nois, somadist, false);
  
  // joesterle begin
  addChannelToElem(ct, cn, region, e, K,  11, C_HCN1, hcn1nois, somadist, false);
  addChannelToElem(ct, cn, region, e, K,  12, C_HCN2, hcn2nois, somadist, false);
  addChannelToElem(ct, cn, region, e, K,  13, C_HCN3, hcn3nois, somadist, false);
  addChannelToElem(ct, cn, region, e, K,  14, C_HCN4, hcn4nois, somadist, false);
  // joesterle end 
  
	addChannelToElem(ct, cn, region, e, K,  5, C_K5, k5nois, somadist, false);
	addChannelToElem(ct, cn, region, e, K,  6, C_K6, k6nois, somadist, false);
	addChannelToElem(ct, cn, region, e, K,  7, C_K7, k7nois, somadist, false);
	addChannelToElem(ct, cn, region, e, KCa, 3, C_KCA3, kca3nois, somadist, false);
	addChannelToElem(ct, cn, region, e, KCa, 4, C_KCA4, kca4nois, somadist, false);
	addChannelToElem(ct, cn, region, e, KCa, 5, C_KCA5, kca5nois, somadist, false);
	addChannelToElem(ct, cn, region, e, ClCa,1, C_CLCA1, clcanois, somadist, false);
	addChannelToElem(ct, cn, region, e, ClCa,2, C_CLCA2, clcanois, somadist, false);
	addChannelToElem(ct, cn, region, e, CA, 0, C_CA0, ca1nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 1, C_CA1, ca1nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 2, C_CA2, ca1nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 3, C_CA3, ca1nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 4, C_CA4, ca1nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 5, C_CA5, ca5nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 6, C_CA6, ca6nois, false, somadist, true);
	addChannelToElem(ct, cn, region, e, CA, 7, C_CA7, ca7nois, false, somadist, true);

	// check for Ca sensitivity in any of the channels in this region
	//
	for (cadens=0,i=0; i<NCHANS; i++) {
	    if (i>=C_KCA0 && i<C_CAPUMP) {
              cadens += celdens[ct][ndens[ct][cn]][i][region];
	    }
	}
	if (cadens > 0) { // If any channels need a Ca comp, make one
			  // (AMPA channels make cacomp separately)
		 double pkm, bmax, bkd, btot, btoti, cashell;

	  cattrib* ca = (cattrib *) make_chan(e, CACOMP, 0);
	  if ((ca->vmax = celdens[ct][ndens[ct][cn]][C_CAPUMP][region]) > 0) {
		ca->pump = 1;
		if ((pkm = celdens[ct][ndens[ct][cn]][C_CAPKM] [region]) <= 0) pkm = dcapkm;
		ca->pkm = pkm;
	  }
	  if ((ca->btot = celdens[ct][ndens[ct][cn]][C_CABTOT][region]) > 0) { 		// set btot = dcabt in dens file
	  	ca->cabuf = 1;
		if ((bmax = celdens[ct][ndens[ct][cn]][C_CABVMAX][region]) <=0) bmax = dcabf;
		ca->bmax = bmax;
		if ((bkd = celdens[ct][ndens[ct][cn]][C_CABKD][region]) <=0) bkd = dcabr/dcabf;
		ca->bkd = bkd;
		if ((btoti = celdens[ct][ndens[ct][cn]][C_CABTOTI][region]) <=0) btoti = ca->btot; // same as others
		ca->btoti = btoti;
	  }

	  if ((ca->kex = celdens[ct][ndens[ct][cn]][C_CAEXCH][region]) > 0) {
	        ca->exch = 1;
	  }
	  if ((cashell = celdens[ct][ndens[ct][cn]][C_CASHELL][region]) <=0) cashell = nshell_soma;
	  ca->cashell = cashell;

		/* CICR */

	  addcicr(ca, ct, cn, region);		/* add CICR, IP3 parameters */
     }
}


/*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

void make_celseg(int ct, int cn, int nod1, 
		 int nod2, double cd1, double cd2, int region, double cplam)
{
    double x, y, z, dist, local_cplam, local_dia;
    double vstart, cvrev, cri, crm, ccm, somadist;
    elem *e;
    sphere *s;
    cable *c;
    chattrib *a;
    cattrib *ca;
    nattrib *n;

  if (ninfo>=5) fprintf(stderr,"# Make_celseg: nod1 %d nod2 %d reg %d\n", nod1, nod2, region);

  //if (nod1 < nod2) {swap=nod1; nod1=nod2; nod2=swap;}; /* now: nod2 <=nod1 */
  // if (ninfo >=5) fprintf(stderr,"Make_celseg after swap: nod1 %d nod2 %d\n", nod1, nod2);

/*
  if (region==SOMA)      printf("Soma       %d %d\n", nod1,nod2);
  if (region==DEND)      printf("Dend       %d %d\n", nod1,nod2);
  if (region==HILLOCK)   printf("Hillock    %d %d\n", nod1,nod2);
  if (region==AXON_THIN) printf("Axon_thin: %d %d\n", nod1,nod2);
  if (region==AXON)      printf("Axon       %d %d\n", nod1,nod2);
  if (region==AXON_DIST) printf("Axon_long  %d %d\n", nod1,nod2);
 /* */

/*
  x = node [ct][cn][nod1]->xloc - node [gc][cn][nod2]->xloc;
  y = node [ct][cn][nod1]->yloc - node [gc][cn][nod2]->yloc;
  z = node [ct][cn][nod1]->zloc - node [gc][cn][nod2]->zloc;
  dist = sqrt (x*x + y*y);
  printf("r:%g\n",dist /sqrt (x*x + y*y ) );
*/
   
  if (cd2<0) cd2 = cd1;	/* set default no taper */

  if (cplam > 0) local_cplam = cplam;
  else local_cplam = celdens[ct][ndens[ct][cn]][C_CPLAM][region];		/* look in dens_xxx.n for cplam */
  if (local_cplam <= 0) local_cplam = getn(ct,COMPLAM); /* otherwise, use nval_xxx.n */ 
  if (local_cplam <= 0) local_cplam = complam;

  /* dia factor for cell regions */

  local_dia = celdens[ct][ndens[ct][cn]][C_DDIA][region];		/* look in dens_xxx.n for ddia */
  if (local_dia <= 0) local_dia = cell_dia_factor; 			/* otherwise, use dia factor */ 
  else local_dia *= cell_dia_factor;

  /* allow changing diameter in dendrites */

  if (region==DEND || region==DEND_PROX || region==DEND_DIST) {
     local_dia *= dend_dia_factor;
  }

  if (region==DEND)  {
     local_dia *= dendi_dia_factor;
  }

  if (region==DEND_PROX)  {
     local_dia *= dendp_dia_factor;
  }

  if (region==DEND_DIST)  {
     local_dia *= dendd_dia_factor;
  }
 				/* allow changing diameter in axon */

  if (region==AXONT || region==AXON || region==AXOND || region==HILLOCK) {
     local_dia *= ax_dia_factor;
  }

  if (local_dia <= 0) local_dia = 1.0; 
  cd1 *= local_dia;
  cd2 *= local_dia;

  if (cd1 < dia_min_factor) cd1 = dia_min_factor;
  if (cd2 < dia_min_factor) cd2 = dia_min_factor;

/*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

  if (getn(ct,BIOPHYS)==0) {		/* no channels, passive membrane */

    vstart = getn(ct,VSTART);
    cvrev  = getn(ct,VREV);
    crm    = getn(ct,NRM);
    cri    = getn(ct,NRI);
    if (crm==0) crm = drm;
    if (cri==0) cri = dri;

    crm = crm/exp(log(dqrm)*(tempcel-dbasetc)/10.0);	// default dqrm, dqri = 1, no effect
    cri = cri/exp(log(dqri)*(tempcel-dbasetc)/10.0);

    if (ninfo >=5)
      fprintf(stderr, "# Make_celseg %s, biophys=%g;vstart=%g;vrev=%g;rm=%g\n",
                                  cname[ct], getn(ct,BIOPHYS), vstart, cvrev, crm);

    if (nod1==nod2){
      s = make_sphere(nd(ct, cn, nod1), cd1, crm, dcm);
      s->dia = cd1;
      s->vrest = vstart;
      s->vrev = cvrev;
      s->elabl = (char *)regname[region];
      s->region = region;
    }
      //at [ct][cn][nod1] sphere dia=cd1 vrest=vstart vrev=cvrev rm=crm
				//elabl regname[region]
    else {
      if (region!=PSOMA) {		/* if "point soma", then don't make cable */
        cable *c = (cable *)conn(ct, cn, nod1, ct, cn, nod2, CABLE);
        c->cplam = local_cplam;
        c->vrest = vstart;
        c->vrev = cvrev;
        c->Rm = crm;
        c->Ri = cri;
        c->elabl = (char *)regname[region];
        c->region = region;
        c->dia = cd1;
        c->dia2 = cd2;
      }
    }
  }
  else {

    vstart = celdens[ct][ndens[ct][cn]][C_VSTART][region];
    cvrev  = celdens[ct][ndens[ct][cn]][C_VREV][region];
    crm    = celdens[ct][ndens[ct][cn]][C_RM][region];
    cri    = celdens[ct][ndens[ct][cn]][C_RI][region];
    ccm    = celdens[ct][ndens[ct][cn]][C_CM][region];
    if (cri==0) cri = dri;
    if (crm==0) crm = drm;
    if (ccm==0) ccm = dcm;

    crm = crm/exp(log(dqrm)*(tempcel-dbasetc)/10.0);	// default dqrm, dqri = 1, no effect
    cri = cri/exp(log(dqri)*(tempcel-dbasetc)/10.0);

    if (cri < 0.01) {
	fprintf (stderr,"%s: make_celseg: cri  %g too small, stopping...\n",progname,cri);
	return;
    }
    if (crm < 1.0) {
	fprintf (stderr,"%s: make_celseg: crm  %g too small, stopping...\n",progname,crm);
	return;
    }

    if (ninfo >=5) {
      fprintf(stderr, "Make_celseg %s biophys %g region %s node %d %d %d ",
                      cname[ct], getn(ct,BIOPHYS), regname[region], ct, cn, nod1);
      fprintf(stderr,"vstart = %g cvrev %g crm %g cri %g ccm %g \n",vstart,cvrev,crm,cri,ccm);
    }

   /* if radial gradient set, find dist to soma */

     if (EGRAD >= 0 || LGRAD >=0) {			      // if EGRAD or LGRAD column is set 
         somadist = dist3d(ndn(ct,cn,nod1),ndn(ct,cn,soma));  // get distance to soma
     }

    if (nod1==nod2) {
  	// fprintf(stderr,"Sphere nod1 %d region %d\n",nod1,region);
      
      s = make_sphere(nd(ct, cn, nod1), cd1, crm,ccm);
      s->vrest = vstart;
      s->vrev = cvrev;
      s->elabl = (char *)regname[region];
      s->region = region;
      addChannelsToSphere(ct, cn, nod1, region, somadist);

    }	/* if region is soma */

    else { /* cable segment, dendrites or axon */

      if (ninfo >=5) {
        fprintf(stderr, "Make_celseg region is cable; conn nod1 %d %d %d to nod2 %d %d %d\n",
	                 ct, cn, nod1, ct, cn, nod2);
        fprintf(stderr, "ct %s: local_cplam %g complam: %g\n", cname[ct], local_cplam, complam);
        fprintf(stderr, "diameters: cd1=%g cd2=%g\n", cd1, cd2);
      }

      if (region!=PSOMA) {		/* if "point soma", then don't make cable */

        c = make_cable(nd(ct, cn, nod1), nd(ct, cn, nod2));
        c->dia = cd1;
        c->dia2 = cd2;
        c->cplam = local_cplam;
        c->vrest = vstart;
        c->vrev = cvrev;
        c->Rm = crm;
        c->Ri = cri;
        c->Cm = ccm;
        c->elabl = regname[region];
        c->region = region;
        e = (elem *)c;

        addChannelsToCable(ct, cn, e, region, somadist);
      }
    } /* if region is cable  _*/
  } /*if biophys           _*/
} /* proc make_celseg    _*/

/*--------------------------------------------------------*/
