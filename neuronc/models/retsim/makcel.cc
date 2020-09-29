/*  makcel.cc for script retsim.cc */

/* Script to create different cell types.  */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "ncio.h"
#include "retsim.h"
#include "retsim_var.h"

void addChannelsToSphere (int,int,int,int);

extern int *cellnums[NCELTYPES];      /* cell numbers indexed by [ct][cn] */

double *cellmorphs[NCELTYPES] = {NULL}; /* morpology arrays in memory */
int cell_lines[NCELTYPES] = {0};	/* size of morphologies */

void rmcelconn(int ct,int cn);		/* erase connection info */
node *getnpnt (nodeint nodenm1, nodeint nodenm2);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* void printneurvals(void)

 // print out neuron values to check

{
   int i,j;

  printf (" ");
  for (i=0; i<nceltypes; i++ ) {
     printf (" %7s", cname[i]);
  }
  printf ("\n\n");
  for (j=0; j<NPARAMS; j++) {		// print nval array
    printf (" ");
    for (i=0; i<nceltypes; i++ ) {
     printf (" %7.3g", nval[i][j]);
    }
    printf ("\n");
  } 
} /* */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void initneurvals(void)

/* Set the parameters describing the neurons and their circuits */

{
    int c,ct,n,p, nrows, ncols, ncolsx, pindx, pnum;
    double *tnval;
    char *nvfil;

  cname [xrod]   = "rod";		/* set the names */
  cname [xcone]  = "cone";
  cname [hbat]   = "hbat";
  cname [ha]     = "ha";
  cname [hb]     = "hb";
  cname [rbp]    = "rbp";
  cname [dbp1]   = "dbp1";
  cname [dbp2]   = "dbp2";
  cname [dbp3]   = "dbp3";
  cname [dbp4]   = "dbp4";
  cname [hbp1]   = "hbp1";
  cname [hbp2]   = "hbp2";
  cname [aii]    = "aii";
  cname [am]     = "am";
  cname [am2]    = "am2";
  cname [am3]    = "am3";
  cname [am4]    = "am4";
  cname [amh]    = "amh";
  cname [amh2]   = "amh2";
  cname [ams]    = "ams";
  cname [amhs]   = "amhs";
  cname [sbac]   = "sbac";
  cname [a17]    = "a17";
  cname [gca]    = "gca";
  cname [gcb]    = "gcb";
  cname [dsgc]   = "dsgc";
  cname [gcaoff] = "gcaoff";
  cname [gcboff] = "gcboff";

  rname [xglut]	   = "glut";
  rname [xampa]	   = "ampa";
  rname [xampa1]   = "ampa";
  rname [xampa2]   = "ampa2";
  rname [xampa3]   = "ampa3";
  rname [xampa4]   = "ampa4";
  rname [xampa5]   = "ampa5";
  rname [xkainate] = "kainate";
  rname [xmglur6]  = "mGluR6";
  rname [xgaba]	   = "gaba";
  rname [xgaba1]   = "gaba1";
  rname [xgaba2]   = "gaba2";
  rname [xgly]	   = "gly";
  rname [xgapj]	   = "gj";
  rname [xdyad]	   = "dyad";

  if (notinit(nvalfile)) nvalfile = (char *)"nval.n";	/* default neuron params */

  nvfil = emalloc(min(FILNAMSIZ,strlen(confdir)) + min(FILNAMSIZ,strlen(nvalfile)) + 2);
  if (streq(confdir,"")) sprintf (nvfil,"%.50s",nvalfile);
  else                   sprintf (nvfil,"%.50s/%.50s",confdir,nvalfile);

  if (ninfo >=3) fprintf(stderr,"#   makcel: reading neural parameters from '%s'\n", nvalfile);
				  
  if ((tnval=fread (nvfil,&nrows,&ncols))==NULL) {  /* read nval file, empty values are zero */
     disperror();
     fprintf(stderr,"# retsim: makcel: can't read nval file '%s'\n", nvfil);
  }
  else {
    if (nrows > (NPARAMS+1))
     fprintf(stderr,"#   makcel: nval file is too long, %d rows\n", nrows);

    if (ncols > (NCELTYPES+1))
     fprintf(stderr,"#   makcel: nval file is too wide %d \n", ncols);

    ncolsx = ncols - 1;
   for (p=1; p<nrows; p++) { /* transpose and re-index array to make cell the first index */
      pindx = p * ncols;
      pnum  = *(tnval+pindx+ncols-1);     /* get the parameter from the first column */
      // fprintf (stderr,"param %d\n",pnum);
      for (c=0; c<ncolsx; c++) {
        ct = *(tnval + c);        /* get cell type from first row */
        nval[ct][pnum] = *(tnval+pindx+c);
      }
    }
    efree (tnval);
    // printneurvals();
  }
}

/*-----------------------------------------*/

void setmorphfiles(void)
{
  /* Set the filenames for realistic neuron files */

  anatfiles[xrod]  = (char *)"";		/* set the names */
  anatfiles[xcone] = (char *)"";
  anatfiles[hbat]  = (char *)hbat_file;
  anatfiles [ha]   = (char *)ha_file;
  anatfiles [hb]   = (char *)hb_file;
  anatfiles [rbp]  = (char *)rbp_file;
  anatfiles [dbp1] = (char *)dbp1_file;
  anatfiles [dbp2] = (char *)dbp2_file;
  anatfiles [dbp3] = (char *)dbp3_file;
  anatfiles [dbp4] = (char *)dbp4_file;
  anatfiles [hbp1] = (char *)hbp1_file;
  anatfiles [hbp2] = (char *)hbp2_file;
  anatfiles [aii]  = (char *)aii_file;
  anatfiles [sbac] = (char *)sbac_file;
  anatfiles [am]   = (char *)am_file;
  anatfiles [am2]  = (char *)am2_file;
  anatfiles [am3]  = (char *)am3_file;
  anatfiles [am4]  = (char *)am4_file;
  anatfiles [amh]  = (char *)amh_file;
  anatfiles [amh2] = (char *)amh2_file;
  anatfiles [ams]  = (char *)"";
  anatfiles [amhs] = (char *)"";
  anatfiles [a17]  = (char *)a17_file;
  anatfiles [gca]  = (char *)gca_file;
  anatfiles [gcb]  = (char *)gcb_file;
  anatfiles [dsgc] = (char *)dsgc_file;
  anatfiles [gcaoff] = (char *)gcaoff_file;
  anatfiles [gcboff] = (char *)gcboff_file;
}

/*-----------------------------------------*/

void printfilenames(FILE *filout)

{
   int ct, cn, found, morph, maxnum;

  ncfprintf (filout,"#   confdir:           %s\n", confdir);
  ncfprintf (filout,"#   nvalfile:          %s\n", nvalfile);
  for (ct=0; ct<nceltypes; ct++) {      /* get density file names */
    if (getn(ct, MAKE)) {
       if ((morph=int(getn(ct,MORPH)))==0) {
         if (anatfiles[ct]!="")
          ncfprintf (filout,"#   %-5s morph:       %-8s", cname[ct],anatfiles[ct]);
         else {
          setn(ct,MORPH,morph=1);
          ncfprintf (filout,"#   %-5s morph   set:  %-8d",cname[ct],morph);
         };
       }
       else
          ncfprintf (filout,"#   %-5s morph:       artif %d", cname[ct],morph);
       if (getn (ct,BIOPHYS)) {
	   maxnum = getn(ct,MAXNUM);
	   for (cn=found=0; cn<maxnum; cn++) {
	       if (ndens[ct][cn]==0) found=1;
	   }
           if (found) {
		   if (densfil[ct][0]==NULL)
		     ncfprintf (filout,",   densities: %s", "zerodens");
		   else
		     ncfprintf (filout,",   densities: %s", densfil[ct][0]);
	   }
	   for (cn=found=0; cn<maxnum; cn++) {
	       if (ndens[ct][cn]==1) found=1;
	   }
           if (found && densfil[ct][1])
		 ncfprintf (filout,", %s", densfil[ct][1]);
	}
       ncfprintf (filout,"\n");
    };
  };
  ncfprintf    (filout,"#   chanparams file:   %s\n", chanparamsfile);
  if (!notinit(syn_savefile)) 
     ncfprintf (filout,"#   syn_savefile:      %s\n", syn_savefile);
  if (!notinit(syn_restorefile)) 
     ncfprintf (filout,"#   syn_restorefile:   %s\n", syn_restorefile);

}

/*-----------------------------------------*/

int find_ct(const char *ct_name)

{
   int n,found=0;

 if (ct_name) {
   for (n=0; n<nceltypes; n++) {
     if (streq(cname[n],ct_name)) {
       found = 1;
       break;
     }
   }
   if (found) return (n);
 }
 return (-1);
}

/*-----------------------------------------*/

void make_ct(int ct)

{
  switch (ct) {
    case XROD:  make_rods=1;  break;
    case XCONE: make_cones=1; break;
    case HBAT:  make_hbat=1;  break;
    case HA:    make_ha=1;    break;
    case HB:    make_hb=1;    break;
    case RBP:   make_rbp=1;   break;
    case DBP1:  make_dbp1=1;  break;
    case DBP2:  make_dbp2=1;  break;
    case DBP3:  make_dbp3=1;  break;
    case DBP4:  make_dbp4=1;  break;
    case HBP1:  make_hbp1=1;  break;
    case HBP2:  make_hbp2=1;  break;
    case A17:   make_a17=1;   break;
    case AII:   make_aii=1;   break;
    case SBAC:  make_sbac=1;  break;
    case AM:    make_am=1;    break;
    case AM2:   make_am2=1;   break;
    case AM3:   make_am3=1;   break;
    case AM4:   make_am4=1;   break;
    case AMH:   make_amh=1;   break;
    case AMH2:  make_amh2=1;  break;
    case AMS:   make_ams=1;   break;
    case AMHS:  make_amhs=1;  break;
    case GCA:   make_gca=1;   break;
    case GCB:   make_gcb=1;   break;
    case DSGC:  make_dsgc=1;  break;
    case GCAOFF: make_gcaoff=1; break;
    case GCBOFF: make_gcboff=1; break;
  }
}

/*-----------------------------------------*/

void set_ncel(int ct, int n)

{
  switch (ct) {
    case XROD:  make_rods=1; if (n_rods<0)  n_rods=nrods=n;  break;
    case XCONE: make_cones=1;if (n_cones<0) n_cones=ncones=n; break;
    case HBAT:  make_hbat=1; if (n_hbat<0)  n_hbat=nhbat=n;  break;
    case HA:    make_ha=1;   if (n_ha<0)    n_ha=nha=n;    break;
    case HB:    make_hb=1;   if (n_hb<0)    n_hb=nhb=n;    break;
    case RBP:   make_rbp=1;  if (n_rbp<0)   n_rbp=nrbp=n;   break;
    case DBP1:  make_dbp1=1; if (n_dbp1<0)  n_dbp1=ndbp1=n;  break;
    case DBP2:  make_dbp2=1; if (n_dbp2<0)  n_dbp2=ndbp2=n;  break;
    case DBP3:  make_dbp3=1; if (n_dbp3<0)  n_dbp3=ndbp3=n;  break;
    case DBP4:  make_dbp4=1; if (n_dbp4<0)  n_dbp4=ndbp4=n;  break;
    case HBP1:  make_hbp1=1; if (n_hbp1<0)  n_hbp1=nhbp1=n;  break;
    case HBP2:  make_hbp2=1; if (n_hbp2<0)  n_hbp2=nhbp2=n;  break;
    case A17:   make_a17=1;  if (n_a17<0)   n_a17 =na17=n;  break;
    case AII:   make_aii=1;  if (n_aii<0)   n_aii =naii=n;  break;
    case SBAC:  make_sbac=1; if (n_sbac<0)  n_sbac=nsbac=n;  break;
    case AM:    make_am=1;   if (n_am<0)    n_am  =nam=n;  break;
    case AM2:   make_am2=1;  if (n_am2<0)   n_am2 =nam2=n;  break;
    case AM3:   make_am3=1;  if (n_am3<0)   n_am3 =nam3=n;  break;
    case AM4:   make_am4=1;  if (n_am4<0)   n_am2 =nam4=n;  break;
    case AMH:   make_amh=1;  if (n_amh<0)   n_amh =namh=n;  break;
    case AMH2:  make_amh2=1; if (n_amh2<0)  n_amh2=namh2=n;  break;
    case AMS:   make_ams=1;  if (n_ams<0)   n_ams =nams=n;  break;
    case AMHS:  make_amhs=1; if (n_amhs<0)  n_amhs=namhs=n;  break;
    case GCA:   make_gca=1;  if (n_gca<0)   n_gca =ngca=n;  break;
    case GCB:   make_gcb=1;  if (n_gcb<0)   n_gcb =ngcb=n;  break;
    case DSGC:  make_dsgc=1; if (n_dsgc<0)  n_dsgc=ndsgc=n;  break;
    case GCAOFF: make_gcaoff=1; if (n_gcaoff<0) n_gcaoff=ngcaoff=n; break;
    case GCBOFF: make_gcboff=1; if (n_gcboff<0) n_gcboff=ngcboff=n; break;
  }
}

/*-----------------------------------------*/

int paramcheck(int var, int limit)

{
   static int err=0;

  if (var>=limit) {
	if (!err) {
	  ncfprintf(stderr,"retsim: wrong value for param '%d', limit %d\n",
							   var, limit);
	  err = 1;
	}
	return 0;
  }
  else  return 1;
}

/*-----------------------------------------*/

void setn (int cel, int var, double val)

/* Set a value in the neural params table */

{
  if (paramcheck(var,NPARAMS))
    nval[cel][var] = val;
}

/*-----------------------------------------*/

double getn (int cel,int var)

/* Get a value from the neural params table */

{
   static int err=0;

  if (paramcheck(var,NPARAMS))
   return nval[cel][var];
  else return 0;
}


/*-----------------------------------------*/

int anysetn (int var)

/* Check if any cell made has a value set for a param in the neural params table */

{
   int ct, found;

  if (paramcheck(var,NPARAMS)) {
    for (found=ct=0; ct<nceltypes; ct++) {
      if (nval[ct][MAKE]>0 && nval[ct][var]!= 0) {
        found = 1;
        break;
      };
    };
    if (found) return 1;
    else       return 0;
  }
  else return 0;
}

/*-----------------------------------------*/

double getcv(int ctype, int var, int n)

/* Get a value from the synaptic input connection table. */
/*  "var" ranges from CELPRE (0) to CELCONV (NCONNP-1) */

{
  if (paramcheck(var,NCONNP)) {
    if (n>0)
    return nval[ctype][CELPRE1+var+(n-1)*NCONNP];
    else return 0;
  }
  else return 0;
}

/*-----------------------------------------*/

void setcv(int ctype, int var, int n, double val)

/* Set a value from the synaptic input connection table. */
/*  "var" ranges from CELPRE (0) to CELCONV (NCONNP-1) */

{
  if (paramcheck(var,NCONNP)) {
    nval[ctype][CELPRE1+var+(n-1)*NCONNP] = val;
  }
}

/*-----------------------------------------*/

double getsv(int ctype,int var,int n)

/* Given presynaptic celltype and synapse type (n),
   look up a certain parameter "var" in the table*/
/*  "var" ranges from CELPOST (0) to SVREV (NSYNP-1) */

{
  if (paramcheck(var,NSYNP)) {
    if (n>0)
    return nval[ctype][CELPOST1+var+(n-1)*NSYNP];
    else return 0;
  }
  else return 0;
}

/*-----------------------------------------*/

void setsv(int ctype,int var,int n,double val)

/* Set a value from the synaptic parameter table. */
/*  "var" ranges from CELPOST (0) to SVREV (NSYNP-1) */

{
  if (paramcheck(var,NSYNP)) {
    nval[ctype][CELPOST1+var+(n-1)*NSYNP] = val;
  }
}

/*-----------------------------------------*/

void initsynconn(void)

  /* find and index all connections */

{
    int i,j,ct,st;

  for (ct=0; ct<nceltypes; ct++)
    for (j=0; j<nceltypes; j++) {
     for (st=0; st<NSYNTYPES; st++) {
      cellconn[ct][j][st] = 0;
     }
    };

  for (ct=0; ct<nceltypes; ct++) {		/* For all cell types */
    for (j=1; j<=NCONNO; j++) {			/* check all connections. */
      if (getsv(ct,AUTAPSE,j) == 0) { 		/* If not autapse */
        if (getsv(ct,CELPOST,j) >= 0) {		/* If ct cell connects, */
 	  for (st=0;st<NSYNTYPES && cellconn[ct][int(getsv(ct,CELPOST,j))][st] > 0; st++) { } /* check connection. */
	  if (st>=NSYNTYPES) { ncfprintf 
		  (stderr,"initsynconn: %s: %s has too many types of synaptic connections %d %d\n", 
		   					nvalfile,cname[ct],st,j);
		  break;
	  }
 	  cellconn[ct][int(getsv(ct,CELPOST,j))][st] = j;  /* save connection num. */
        };
      };
    };
  };

/*
 printf("#   Connections between cells\n");
 printf("\n");
 for (j=0; j<nceltypes; j++) printf (" %-13s ",cname[j]);
 printf("\n");
 for (j=0; j<nceltypes; j++) {
   for (i=0; i<nceltypes; i++) { 
     for (st=0; st<NSYNTYPES; st++) { 
       printf ("%d ",cellconn[i][j][st]);
     }
     printf(" ");
   }
   printf("\n");
 }
/* */

}

/*-----------------------------------------*/

int getconn(int a,int b)

/* Find the synaptic connection number from cell a to cell b. */

{
  if (a>=0 && b>=0) {
    return cellconn[a][b][0];
  }
  else return 0;
}

/*-----------------------------------------*/

int getconns(int a,int b, int i)

/* Find the synaptic connection number from cell a to cell b. */

{
  if (a>=0 && b>=0) {
    return cellconn[a][b][i];
  }
  else return 0;
}

/*------------------------------------------*/

void enlarge_array (int celltype, int cellnum)

{
  find_maxmin(celltype,cellnum);        /* get size of array, update array size */
  xarrsiz=xmax-xmin;
  yarrsiz=ymax-ymin;
  arrcentx=(xmax+xmin)/2;
  arrcenty=(ymax+ymin)/2;
  arrsiz=max(xarrsiz,yarrsiz);
}

/*--------------------------------------------------------------*/

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, 
		double arrcentx, double arrcenty, 
		double *sxarr, double *syarr, double *sthyarr, double *stharr, int *snarr, int first_center) 

/* create array of cells with specified spacing and regularity */

{
  int i,cn;
  int rsd, ginfo;
  double spacing, density, reg, thetay, thetaz;
  double *xarr,*yarr;

 if (notinit(first_center)) first_center = 0;

 if (sxarr==NULL) {

  /* arrsiz, num of cells not set -> make 1 cell if xarrsiz not init */

   if (notinit(lxarrsiz)) lxarrsiz = xarrsiz;
   if (notinit(lyarrsiz)) lyarrsiz = yarrsiz;
 
   if (ninfo>=2) {
	if (ct!=ams && ct!=amhs) ncfprintf (stderr,"#\n# %ss:\n",cname[ct]);
	else                     ncfprintf (stderr,"#\n# %s's:\n",cname[ct]);
   }
   if (ninfo>=3) {
	   if (ncells < 0) 
	   	ncfprintf (stderr,"# %ss: make %s cells in arrsiz x=%g y=%g\n",
		  cname[ct],"?",lxarrsiz,lyarrsiz);
	   else if (ncells == 1) {
		if (notinit(lxarrsiz)) {
	   	   ncfprintf (stderr,"# %ss: make %d cell in undef array size\n",
		      cname[ct],ncells,lxarrsiz,lyarrsiz);
	        } else {
	   	   ncfprintf (stderr,"# %ss: make %d cell in arrsiz x=%g y=%g\n",
		      cname[ct],ncells,lxarrsiz,lyarrsiz);
		}
	   }
           else
	   	ncfprintf (stderr,"# %ss: make %d cells in arrsiz x=%g y=%g\n",
		  cname[ct],ncells,lxarrsiz,lyarrsiz);
   }
   if (notinit(lxarrsiz)) {		/* ncells, density are set; xarrsiz, lxarrsiz are not */
     if (ncells < 0) ncells = 1;
     ncells = gausnn (&xarr, &yarr, ncells, density=getn(ct,DENS)*1e-6, rsd=ct^rseed,
                                reg=getn(ct,REGU), arrcentx, arrcenty, first_center, ginfo=info);
   }

  /* arrsiz set, num of cells isn't */

   else if (ncells < 0) {

     ncells = gausnn (&xarr, &yarr, density=getn(ct,DENS)*1e-6,
                         rsd=ct^rseed,reg=getn(ct,REGU),lxarrsiz,lyarrsiz,
			arrcentx, arrcenty, first_center, ginfo=info);
   }

   /* arrsiz and nr of cells set by user */

    else {

      ncells = gausnn (&xarr, &yarr, ncells, rsd=ct^rseed, reg=getn(ct,REGU),
                  lxarrsiz,lyarrsiz,arrcentx,arrcenty,first_center, ginfo=info);
    }

 }       /* if (sxarr!=NULL) */

 else {  /* use sxarr, syarr, set by user */
   if (syarr != NULL) {
     xarr = sxarr;
     yarr = syarr;
   }
 }  /* use sxarr, syarr */

   if (ncells > getn(ct, MAXNUM))
        ncells = int(getn(ct, MAXNUM));

   if (ninfo>=2 && ct < a17) {
     spacing = 1 / sqrt (getn(ct,DENS) * 1e-6); /* convert cell/mm^2 to um/cell*/
     ncfprintf (stderr,"# %s spacing: %.3g um\n",cname[ct], spacing);
   };
   if (ninfo>=2 && ct < a17) /* print number of photoreceptors or bipolar cells */
	   ncfprintf (stderr,"# number of %ss: %d\n", cname[ct], ncells);

  if (ncells==1) {
    i = 0;
    if (snarr!=NULL) cn = snarr[i];
    else             cn = i+1;
    if (stharr != NULL) 
         thetaz = stharr[i];
    else thetaz = 0;
    if (sthyarr != NULL)
         thetay = sthyarr[i];
    else thetay = 0;
    if (syarr !=NULL)  makcell (ct, cn, xarr[i], yarr[i],thetay,thetaz,flip);
    else	       makcell (ct, cn, arrcentx, arrcenty,thetay,thetaz,flip);
  }
  else
    for (i=0; i<ncells; i++) {
      if (stharr!=NULL)
	   thetaz = stharr[i];
      else thetaz = rrange(0,360);
      if (sthyarr != NULL)
           thetay = sthyarr[i];
      else thetay = 0;
      if (snarr!=NULL) cn = snarr[i];
      else             cn = i+1;
      if (ninfo>=4)
        ncfprintf (stderr,"# %s #%d, x=%g  y=%g r=%.3g\n",
                cname[ct],cn,xarr[i],yarr[i],thetaz);

      makcell (ct, cn, xarr[i], yarr[i],thetay,thetaz,flip);
    }

  find_maxmin(-1,-1);   /* get size of array, update array size */

  if (setarrsiz==0 && setxarrsiz==0) {
    xarrsiz=xmax-xmin;
    yarrsiz=ymax-ymin;
    zarrsiz=zmax-zmin;
    arrcentx=(xmax+xmin)/2;
    arrcenty=(ymax+ymin)/2;
    arrcentz=(zmax+zmin)/2;
    arrsiz=max(xarrsiz,yarrsiz);
  } 
  setn (ct,NMADE,ncells);
  //if (ninfo>=3) ncfprintf (stderr,"# %ss: made %d cells in arrsiz %g, x=%g y=%g\n",
  //		  cname[ct],ncells,arrsiz,xarrsiz,yarrsiz);
  return ncells;

} /* func setupcells() */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *stharr) 

{
 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, arrcentx, arrcenty, sxarr,syarr,NULL,stharr,NULL,0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *stharr, int *snarr) 

{
 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, arrcentx, arrcenty, sxarr,syarr,NULL,stharr,snarr,0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, double arrcentx, double arrcenty, double *sxarr, double *syarr, double *stharr) 

{
 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, arrcentx, arrcenty, sxarr,syarr,NULL,stharr,NULL,0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz) 

{
    double larrcentx, larrcenty;

   if (notinit(arrcentx)) larrcentx = lxarrsiz * 0.5; 
   else                   larrcentx = arrcentx;
   if (notinit(arrcenty)) larrcenty = lyarrsiz * 0.5;
   else                   larrcenty = arrcenty;

 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, larrcentx, larrcenty, NULL, NULL, NULL, NULL, NULL, 0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, double arrcentx, double arrcenty) 

{
 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, arrcentx, arrcenty, NULL, NULL, NULL, NULL, NULL, 0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, 
				double arrcentx, double arrcenty, int first_cent) 

{
 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, arrcentx, arrcenty, NULL, NULL, NULL, NULL, NULL, first_cent);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double larrsiz) 

{
 return setupcells(ct, ncells, larrsiz, LARGENUM);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells) 

{
 return setupcells(ct, ncells, LARGENUM, LARGENUM);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* first_cent=1 => sets first cell at (0,0) */

int setupcells(int ct,int ncells, int first_cent) 

{
 return setupcells(ct, ncells, LARGENUM, LARGENUM, arrcentx, arrcenty, NULL, NULL, NULL, NULL, NULL, first_cent);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double *sxarr, double *syarr) 

{
 return setupcells(ct, ncells, LARGENUM, LARGENUM, sxarr, syarr, NULL);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double *sxarr, double *syarr, double *stharr) 

{
 return setupcells(ct, ncells, LARGENUM, LARGENUM, sxarr, syarr, stharr);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double *sxarr, double *syarr, double *stharr, int *snarr) 

{
 return setupcells(ct, ncells, LARGENUM, LARGENUM, sxarr, syarr, stharr, snarr);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* sthyarr allows different y rotations for stereo pairs */

int setupcells(int ct,int ncells, double lxarrsiz, double lyarrsiz, double *sxarr, double *syarr, double *sthyarr, double *stharr, int *snarr) 

{
 return setupcells(ct, ncells, lxarrsiz, lyarrsiz, arrcentx, arrcenty, sxarr,syarr,sthyarr,stharr,snarr,0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setupcells(int ct,int ncells, double *sxarr, double *syarr, double *sthyarr, double *stharr, int *snarr) 

{
 return setupcells(ct, ncells, LARGENUM, LARGENUM, arrcentx, arrcenty, sxarr,syarr,sthyarr,stharr,snarr,0);
}

/*------------ procs to prune connections ------------------------*/



void erase_cell(int ct, int cn)

/* Erase a cell and remove its connections */

{
     register node *np,*npn;

  ncell_erased[ct]++;
  for (np=getnpnt(ct,cn); np=foreachn(np, ct, cn); np=npn) {
       npn = np->hcnext;
       erasenode (np);               /* erase nodes */
  }
  rmcelconn(ct,cn);                   /* erase connection info */
  if (ninfo>=2) ncfprintf(stderr,".");
}

/*----------------------------------------------------------------*/

void checkcellin (int ct)

/* Erase cells that don't have any inputs */
/*  Useful for bipolar and amacrine cells */
/*  Ignores gap junctions to same type */

{
    int ncell_in,cn,nn;
    node *npnt,*np,*tnpnt;

 for (npnt=nodepnt; npnt=foreach (npnt, ct,-1,soma, NULL,&cn,NULL); npnt=tnpnt) {
    tnpnt = npnt->next;
    ncell_in  = tot_ncel_ind(ct,cn);       /* get number of inputs */
    if (ninfo >= 3) ncfprintf (stderr,"checkcellin: %s %d: %d\n",cname[ct],cn,ncell_in);
    if (ncell_in == 0) {
       erase_cell(ct,cn);			/* erase unconnected cells */
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void checkcellin (int ct, int cti)

/* Erase cells that don't have inputs from another type */
/*  Useful for bipolar and amacrine cells */
/*  Ignores gap junctions to same type */

{
    int ncell_in,cn,nn;
    node *npnt,*np,*tnpnt;

 for (npnt=nodepnt; npnt=foreach (npnt, ct,-1,soma, NULL,&cn,NULL); npnt=tnpnt) {
    tnpnt = npnt->next;
    ncell_in  = ncel_in(ct,cn,cti);       /* get number of inputs */
    if (ninfo >= 3) ncfprintf (stderr,"checkcellin: %s %d: %d\n",cname[ct],cn,ncell_in);
    if (ncell_in == 0) {
       erase_cell(ct,cn);			/* erase unconnected cells */
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void checkcellout (int ct)

/* Erase cells that don't have any outputs */
/*  Useful for bipolar and amacrine cells */
/*  Ignores gap junctions to same type */

{
    int ncell_out,cn,nn;
    node *npnt,*np, *tnpnt;

 for (npnt=nodepnt; npnt=foreach (npnt, ct,-1,soma, NULL,&cn,NULL); npnt=tnpnt) {
    tnpnt = npnt->next;
    ncell_out  = tot_ncel_outd(ct,cn);    /* get number of outputs */
    if (ninfo >= 3) ncfprintf (stderr,"checkcellout: %s %d: %d\n",cname[ct],cn,ncell_out);
    if (ncell_out == 0) {
       erase_cell(ct,cn);			/* erase unconnected cells */
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void checkcellout (int ct, int cto)

/* Erase cells that don't have outputs to another type */
/*  Useful for bipolar and amacrine cells */
/*  Ignores gap junctions to same type */

{
    int ncell_out,cn,nn;
    node *npnt,*np, *tnpnt;

 for (npnt=nodepnt; npnt=foreach (npnt, ct,-1,soma, NULL,&cn,NULL); npnt=tnpnt) {
    tnpnt = npnt->next;
    ncell_out  = ncel_out(ct,cn,cto);    /* get number of outputs */
    if (ninfo >= 3) ncfprintf (stderr,"checkcellout: %s %d: %d\n",cname[ct],cn,ncell_out);
    if (ncell_out == 0) {
       erase_cell(ct,cn);			/* erase unconnected cells */
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void checkcellout (int ct, int ct2, int ct3)

/* Erase cells that don't converge to a third type, */
/*  either directly or through a second type */
/*  Useful for bipolar and amacrine cells that converge to a gc */
/*  Ignores gap junctions to same type */

{
    int	cn,nn;
    int connected,converging;
    node *npnt,*np, *tnpnt;

 for (npnt=nodepnt; npnt=foreach (npnt, ct,-1,soma, NULL,&cn,NULL); npnt=tnpnt) {
    tnpnt = npnt->next;
    connected  = ncel_out(ct,cn,ct3);     /* determine if directly connected */
    converging = connected2 (ct,cn,ct2,ct3); /* determine if converging */
    if (ninfo >= 3) ncfprintf (stderr,"checkcellout2: %s %d: %d %d\n",
		       cname[ct],cn,connected,converging);
    if (!connected && !converging) {
       erase_cell(ct,cn);			/* erase unconnected cells */
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void checkcellout (int ct, int ct2, int ct3, int ct4)

/* Erase cells that don't converge to a fourth type, */
/*  either directly or indirectly through a second or third type */
/*  Useful for bipolar and amacrine cells that converge to a gc */
/*  Ignores gap junctions to same type */

{
    int	cn,nn;
    int connected,converging2, converging3;
    node *npnt,*np, *tnpnt;

 for (npnt=nodepnt; npnt=foreach (npnt, ct,-1,soma, NULL,&cn,NULL); npnt=tnpnt) {
    tnpnt = npnt->next;
    connected  = ncel_out(ct,cn,ct4);     /* determine if directly connected */
    converging2 = connected2 (ct,cn,ct2,ct4); /* determine if converging */
    converging3 = connected3 (ct,cn,ct2,ct3,ct4); /* determine if converging */
    if (ninfo >= 3) ncfprintf (stderr,"checkcellout3: %s %d: %d %d %d\n",
		       cname[ct],cn,connected,converging2, converging3);
    if (!connected && !converging2 && !converging3) {
       erase_cell(ct,cn);			/* erase unconnected cells */
    }
  }
}

/*--------------------------------------------------------------*/

void set_rcolors(int ct, int cn)

/* retrieve and set region colors, in "disprcolor()" for display(). */

{
  rcolors = celdens[ct][ndens[ct][cn]][C_COLOR];
}

/*--------------------------------------------------------------*/

void unset_rcolors(void)

/* retrieve and set region colors, in "disprcolor()" for display(). */

{
  rcolors = NULL;
}

/*--------------------------------------------------------------*/

/*
void dispcelltype(int ct,int n,double scal,double nscal) 

{
    int colval, nod, i, q;
    int dsgcnodq, sbcelq, sbnodq;
    int sbcolor[nsb+1] = {1};

  for (i=1;i<=nsb;i++)
    sbcolor[i]=getn(sb,NCOLOR)-1+i;      // give each sb a different color

  for (i=1; i<=n; i++) {
    if ( !(colval=getn(ct,NCOLOR))) colval = i;
    if (ct==dsgc) {
      for (q=0;q<sb_dsconns;q++) {              // show synapses
        dsgcnodq =  dsgc_sbconn_nodenrs[i][q];  // get nod with sb input
        if (dsgcnodq>=0 && dsgc_in_syns[i][dsgcnodq][0]==sb) {
          sbcelq=dsgc_in_syns[i][dsgcnodq][1];  // sb cell providing input
          sbnodq=dsgc_in_syns[i][dsgcnodq][2];  // sb node providing input
          display synapse matching [dsgc][i][sbnodq] color sbcolor[sbcelq];
        }
      }
      display matching [dsgc][i][-1] only color colval dscale scal;
    }
    else if (ct==sb)
      display matching [ct][i][-1] only color sbcolor[i] dscale scal
    else
      display matching [ct][i][-1] only color colval dscale scal;

    display comps matching [ct][1][-1];
    display node matching [ct][1][-1] dscale nscal;
  }
}
/* */

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

void displaycell(int ct,int cn,int colval, double scal,double nscal)
{
    int color,ns,d;
    double dscale;

    ns = int(nscal);
    display (SYNAPSE,MATCHING,ndt(ct,cn,-1), color=5, 	  dscale=synap_scale);    /* show synapses magenta*/
    display (SPHERE,MATCHING, ndt(ct,cn,-1), color=colval,dscale=scal);
    display (CABLE,MATCHING,  ndt(ct,cn,-1), color=colval,dscale=scal);
    display (GJ,MATCHING,     ndt(ct,cn,-1), color=6,     dscale=gj_scale);
    display (ROD,MATCHING,    ndt(ct,cn,-1), color=colval,dscale=scal);
    display (CONE,MATCHING,   ndt(ct,cn,-1), color=colval,dscale=scal);
    display (COMPS,MATCHING,  ndt(ct,cn,-1), color=1,     dscale=scal);
    if (ns<0) {
      switch (ns) {
       case -1: display (NODE, MATCHING,  ndt(ct,cn,soma), color=white, dscale=nscal); break;
       case -2: display (NODE, MATCHING,  ndt(ct,cn,soma), color=white, dscale=nscal); break;
       case -3: display (NODE, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -4: display (NODE, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -5: display (NODE, MATCHING,  ndt(-1,-1,soma), color=white, dscale=nscal); break;
       case -6: display (NODE, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -7: display (NODE, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -9: break;
       default:
                display (NODE, MATCHING, ndt(ct,cn,soma), color=white, dscale=nscal); /* display soma node */
      }
    } 
    else display (NODE, MATCHING, ndt(ct,cn,-1), color=white, dscale=nscal);  /* positive nscal */
    // dscale = abs(nscal);
    // d = int(dscale);
    // dscale = -(6.0 + dscale - d);   // -6 -> display nodenm2, nodenm3 
    // display (LABEL, MATCHING, ndt(ct,cn,-1), color=white, dscale);  /* always display labels */

      switch (ns) {
       case -1: display (LABEL, MATCHING,  ndt(ct,cn,soma), color=white, dscale=nscal); break;
       case -2: display (LABEL, MATCHING,  ndt(ct,cn,soma), color=white, dscale=nscal); break;
       case -3: display (LABEL, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -4: display (LABEL, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -5: display (LABEL, MATCHING,  ndt(-1,-1,soma), color=white, dscale=nscal); break;
       case -6: display (LABEL, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -7: display (LABEL, MATCHING,  ndt(ct,cn,-1), color=white, dscale=nscal); break;
       case -9: break;
       default:
                display (LABEL, MATCHING, ndt(ct,cn,soma), color=white, dscale=nscal); /* display soma node */
      }
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

void dispcelltype(int ct, int start, int stop, double scal, double nscal) 

{
    int colval, i;

 if (disp>0 && disp<DMOVIE)
  for (i=start; i<=stop; i++) {
    if ((colval=int(getn(ct,NCOLOR)))<0) colval = i % brtwht + 0;
    if (colval == RCOLOR) set_rcolors(ct,i);		    /* set region colors */
    if (!notinit(disp_ct)) {		/* disp_ct and disp_cn override node displays */
      if (!notinit(disp_cn)) {
	 if (ct==disp_ct && i==disp_cn) 
	    displaycell(disp_ct,i,colval,scal,nscal);
      }
      else if (ct==disp_ct) 
	    displaycell(disp_ct,i,colval,scal,nscal);
    }
    else {  displaycell(ct,i,colval,scal,nscal);
    }
  }
  unset_rcolors();					/* unset region colors */
}

/*--------------------------------------------------------------*/

void dispcelltype(int ct, int start, int stop, double scal, double nscal, double zmax, double zmin) 

{
   if (!notinit(zmax) || !notinit(zmin)) display_z(zmax, zmin);
   dispcelltype(ct, start, stop, scal, nscal);
   display_z(LARGENODE, -LARGENODE);
}

/*--------------------------------------------------------------*/

void makhden(int ctype, int cn, int nden,double xpos,double ypos)

/* make Hz cell with dendrites according to parameters */

/* diabs  = starting size of dendrite for taper */
/* diaspc = space constant of dendrite taper */

{
   int i,region,sn;
   double angl,sumangl;
   double xden, yden;
   double leng, thisangl;
   double crm,cvrest;
   double rad,dia,diabs,diaspc,cplam;
   sphere *s;
   node *np;

/* if (ctype==ha) printf ("ha # %d nden %d diabs %g diaspc %g\n",
		cn,nden,diabs,diaspc);
if (ctype==hb) printf ("hb # %d nden %d diabs %g diaspc %g\n",
		cn,nden,diabs,diaspc);
/* */

   rad    = getn(ctype,DTREEDIA)/2;
   diabs  = getn(ctype,TAPERABS);	/* starting taper diameter */
   diaspc = getn(ctype,TAPERSPC);	/* taper space constant */
   celnode[ctype][cn] = 0;
   cplam = -1;

   np = loc(nd(ctype,cn,soma),xpos,ypos,getn(ctype,SOMAZ));
   dia=getn(ctype,SOMADIA);
   make_celseg(ctype,cn,soma,soma,dia,0,region=SOMA,cplam);

   angl = 360 / nden;
   sumangl = 0;
   for (i=1; i<=nden; i++) {
     thisangl = angl * rrange (0.8,1.2);
     sumangl += thisangl;
     leng = rad * rrange (0.8, 1.2) * 0.2;
//printf ("angle %g length %g\n", angl, leng);

     xden = cos(sumangl/DEG) * leng + xpos;
     yden = sin(sumangl/DEG) * leng + ypos;
//printf ("loc %g %gi\n", xden, yden);
     np = nd(ctype,cn,++celnode[ctype][cn]);
     np = loc(np,xden,yden,getn(ctype,DENDARBZ)+2);

     sn = celnode[ctype][cn];		/* sequence node */   
     celnode[ctype][cn] = taperdistden (ctype,cn,0,celnode[ctype][cn],
				diabs,diaspc,DEND_PROX, celnode[ctype][cn]);

     if (leng < rad * 0.5) {                  /* make first bifurcation */
           double shortang, shortlen, shsumang;
           double xden1, yden1;

        shortlen  = rrange(0.5,1) * (rad-leng) * 0.5;
	if (shortlen+leng > rad) shortlen = rad - leng;
        shortang = angl * rrange(0.4, 1.0) * 0.5;
        shsumang = sumangl - shortang;		/* less angle */

        xden1 = cos(shsumang/DEG) * shortlen + xden;
        yden1 = sin(shsumang/DEG) * shortlen + yden;
        np = nd(ctype,cn,++celnode[ctype][cn]);
        np = loc(np,xden1,yden1,getn(ctype,DENDARBZ)+2);
        celnode[ctype][cn] = taperdistden(ctype,cn,sn,celnode[ctype][cn],
				diabs,diaspc,DEND_DIST,celnode[ctype][cn]);

        shortlen  = rrange(0.5,1) * (rad-leng) * 0.5;
	if (shortlen+leng > rad) shortlen = rad - leng;
        shortang = angl * rrange(0.4, 1.0) * 0.5;
        shsumang = sumangl + shortang;		/* greater angle */

        xden1 = cos(shsumang/DEG) * shortlen + xden;
        yden1 = sin(shsumang/DEG) * shortlen + yden;
        np = nd(ctype,cn,++celnode[ctype][cn]);
        np = loc(np,xden1,yden1,getn(ctype,DENDARBZ)+2);
        celnode[ctype][cn] = taperdistden(ctype,cn,sn,celnode[ctype][cn],
				diabs,diaspc,DEND_DIST,celnode[ctype][cn]);
     }
   }     /* for (i=0;i<nden;) */
}

/*-----------------------------------------*/

void extend_gc_branch(int brnum, int obrnum, int ctype, int n, int sn,
				double cumleng, double angl, int ord)

  /* Extend an existing branch on a ganglion cell dendrite. */
  /*  Includes competition between branches of same cell type. */

{
	int i, do_branch, nconn,nelem;
	int cnum, cnod, cbrnum;
        int nsn, brdir;
	int n1a, n1b, n1c, n1d;
	int n2a, n2b, n2c, n2d, csn;
	double leng, cmindist, dist, mdist, brprob_dia, brprob_dist;
	double sdia, shortang, shortlen, shsumang, sangl;
	double xdiff, ydiff;
	double xm, ym, zm;
	double xden, yden,xden1,yden1;
        double brdia,edia,mdia,growthresh,cplam;
	double crm, cvrest, rad, diabs, diaspc, stratd;
	double cd1, cd2;
	elem *e,*epnt;
	node *np;

   crm    = getn(ctype,NRM);
   cvrest = getn(ctype,VSTART);
   rad    = getn(ctype,DTREEDIA)/2;
   diabs  = getn(ctype,TAPERABS);
   diaspc = getn(ctype,TAPERSPC);
   stratd  = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
   cplam = -1;

   brdir = ffs();		/* branch direction */
   mdia   = rrange(0.6, 0.9);	/* dia of "main" side of branch */
   brdia  = rrange(0.3, 0.6);	/* dia of side branch */

	/* find branches at distal end of segment */

   e = get_elempnt(brnum);

   n1a = e->node1a;
   n1b = e->node1b;
   n1c = e->node1c;
   n1d = e->node1d;

   n2a = e->node2a;
   n2b = e->node2b;
   n2c = e->node2c;
   n2d = e->node2d;
   sdia = ((cable *)e)->dia;

   nconn = int(get_nfield(ndt(n2a,n2b,n2c,n2d),NUMCONN));
   xden  = get_nfield(ndt(n2a,n2b,n2c,n2d),XLOC);
   yden  = get_nfield(ndt(n2a,n2b,n2c,n2d),YLOC);
   xden1 = get_nfield(ndt(n1a,n1b,n1c,n1d),XLOC);
   yden1 = get_nfield(ndt(n1a,n1b,n1c,n1d),YLOC);

   xdiff = xden - xden1;
   ydiff = yden - yden1;
   if (xdiff==0) xdiff=1e-6;
   sangl = atan(ydiff/xdiff) * DEG;
   if (xdiff < 0) sangl += 180;

   /* if (ninfo>=3) printf ("xden %d %g %g %g %g %g\n",
			n1c,xden1,yden1,n2c,xden,yden); /* */
   csn = n2c;

   if (ninfo>=3) printf ("tracing from node %d %d %d\n",
			n2a,n2b,n2c);
   if (ninfo>=3) printf ("nconn %d\n",nconn);
   if (ninfo>=3) printf ("sangl %g\n",sangl);

	/* follow each branch to its end */

   if (nconn>1) {
     for (i=1; i<=nconn; i++) {
        nelem = int(get_nfield(ndt(n2a,n2b,n2c),ELNUM,i));
        if (ninfo>=3) printf ("conn %d: elem %d\n", i,nelem);
	if (nelem==brnum) continue;
        if (get_efield(nelem,NTYPE) != CABLE) continue;
        if (ninfo>=3) printf ("following conn %d: elem %d\n", i,nelem);
	leng = get_efield(nelem,LENGTH);
        extend_gc_branch(nelem,brnum,ctype,n,n2c,cumleng+leng,angl,ord);
     };
   }
   else {	/* If this is the end, add to it */


     if (ninfo>=3) printf ("found the end at node %d %d %d, %d\n", n2a,n2b,n2c,sn);

	/* compute location of new tip of main dendrite */

     shortlen  = rrange (0.8,1.2) * rad * 0.2; //(rad*1.5-cumleng) *  0.5;
     if (shortlen+cumleng > rad) shortlen = rad - cumleng;
     shortang = angl * rrange(0.6,1) * 0.6;
     shsumang = sangl + shortang * brdir;
     if (shsumang < 0) shsumang += 360;
  if (ninfo>=3) printf ("shsumang = %g %g %g\n",shsumang,sangl,shortang);
        if (ninfo>=3) printf ("xden yden %g %g\n", xden,yden);
     xm = cos(shsumang/DEG) * shortlen + xden;
     ym = sin(shsumang/DEG) * shortlen + yden;
     zm = getn(ctype,DENDARBZ);
     nsn = ++celnode[ctype][n];
     np = loc(nd(ctype,n,nsn),xm,ym,zm);

		/* Check proximity of other branches of same cell type. */

     cmindist = mdist = 1e10; /* now find closest dendrite */
     cbrnum = -1;
     for (epnt=elempnt; epnt=foreach(epnt,CABLE,ctype,-1,-1,-1,NULL,&cnum,&cnod,NULL);
			epnt = epnt->next) {
       cbrnum = epnt->elnum;
       if (cbrnum==brnum) continue;		/* ignore this cable */
       if (cbrnum==obrnum) continue;		/* ignore this cable */
//       if ((dist=e3dist([n2a][n2b][n2c],cbrnum)) < cmindist) {
//           cmindist = dist;			/* remember this dist */
//       };
       if ((dist=endist3d(cbrnum, ndn(ctype,n,nsn))) < mdist) {
           cmindist = dist;			/* remember this dist */
       };

 //if (ninfo>=3) printf ("cnod %d n2c %d dist %g\n", cnod,n2c,dist);
     };   /* foreach cable, remember dist. */

		/* Extend this branch if other branches are far enough away. */

     if (ninfo>=3) printf ("closest branch %g\n", cmindist);

     if (sdia > getn(ctype,DTIPDIA)*1.5) {	/* only grow if fat enough */

       growthresh = getn (ctype,GROWTHR);
       if (cmindist > growthresh) { 		/* grow, maybe make branch */

 	 do_branch = 0;

         if (cmindist > growthresh*2 &&
            sdia > getn(ctype,DTIPDIA)*3) {  /* only make branch if fat enough */
	    brprob_dist =  cmindist / growthresh / 30;
	    brprob_dia  =  sdia / getn(ctype,DTIPDIA)/20;
     if (ninfo>=3) printf ("branch probability %g %g\n", brprob_dia, brprob_dist);
	    //if (brprob_dist < brprob_dia) brprob_dia = brprob_dist;
	    do_branch = (drand() < (brprob_dia));
//	    if (mdist < growthresh*2) do_branch = 0;

	 };
				/* extend main dendrite */
	 if (! do_branch) {
	    brdir = 0;
	 }


/*  If end of new dendritic tip is too close to cells of same type, */
/*   move tip to make dendrite shorter. */

	 if (cmindist < growthresh*2) {
            shsumang = sangl + shortang * brdir;
     	    xm = cos(shsumang/DEG) * shortlen * 0.5 + xden;
     	    ym = sin(shsumang/DEG) * shortlen * 0.5 + yden;
	    np = loc(nd(ctype,n,nsn),xm,ym,zm);
	 }

      if (ninfo>=3) printf ("extending branch from %d to %d\n", csn, nsn);

	 edia = taperdia(ctype,n,csn,ctype,n,celnode[ctype][n],sdia,diaspc);
         make_celseg (ctype, n, csn, celnode[ctype][n],cd1=sdia*mdia,cd2=edia,DEND,cplam);

         if (do_branch) {		/* make side branch */

	/* Make side branch if other branches are further than threshold. */

	   shortlen  = rrange (0.6,1.0) * (rad*1.2-cumleng)*0.3;
	   if (shortlen+cumleng > rad) shortlen = rad - cumleng;
	   shortang = angl * rrange(0.6,1) * 0.5;
	   shsumang = sangl - shortang * brdir;

	   xden1 = cos(shsumang/DEG) * shortlen + xden;
	   yden1 = sin(shsumang/DEG) * shortlen + yden;
	   np = loc(nd(ctype,n,++celnode[ctype][n]),xden1,yden1,getn(ctype,DENDARBZ));
	   nsn = celnode[ctype][n];
     if (ninfo>=3) printf ("making side branch to %d\n", nsn);
	   edia = taperdia(ctype,n,csn,ctype,n,celnode[ctype][n],sdia,diaspc);
           make_celseg (ctype, n, csn, celnode[ctype][n],cd1=sdia*brdia,
							cd2=edia*brdia,DEND,cplam);
        }

     }		/* cmindist > growthresh */

    }	/* if (sdia > ) */
    else  erase (ndn(ctype,n,nsn));

   }   /* this is end */

   if (ninfo>=3) printf ("extend_gc_branch end\n");
}

/*-----------------------------------------*/

void mak_gc_branch (int ctype, int n, int sn, 
			double cumleng, double angl,double sumangl,
			double cxden,double cyden,int ord,double sdia)

  /* Make a branch on a ganglion cell dendrite */

{
        double shortang, shortlen, shsumang;
        int brdir,nsn;
	double cd1, cd2;
	double xden1, yden1;
        double brdia,edia,mdia;
	double crm, cvrest, rad;
	double diabs, diaspc, stratd,cplam;
	node *np;

 if (ord < 5) {

   crm    = getn(ctype,NRM);
   cvrest = getn(ctype,VSTART);
   rad    = getn(ctype,DTREEDIA)/2;
   diabs  = getn(ctype,TAPERABS);
   diaspc = getn(ctype,TAPERSPC);
   stratd  = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
   cplam = -1;

   brdir = ffs();		/* branch direction */
   mdia   = rrange(0.6, 0.9);	/* dia of "main" side of branch */
   brdia  = rrange(0.3, 0.6);	/* dia of side branch */

	/* make first order branch */

   shortlen  = rrange (0.6,1) * (rad-cumleng) * (1-stratd) * 0.6;
   shortlen  = rrange (0.6,1) * (rad*1.5-cumleng) *  0.3;
   if (shortlen+cumleng > rad) shortlen = rad - cumleng;
   shortang = angl * rrange(0.3,1) * 0.5; // * (rad*0.4 + rad)/shortlen;
   shsumang = sumangl + shortang * brdir;

   xden1 = cos(shsumang/DEG) * shortlen + cxden;
   yden1 = sin(shsumang/DEG) * shortlen + cyden;
   np = loc(nd(ctype,n,++celnode[ctype][n]),xden1,yden1,getn(ctype,DENDARBZ));
   nsn = celnode[ctype][n];
   edia = taperdia(ctype,n,sn,ctype,n,celnode[ctype][n],sdia,diaspc);
           
   make_celseg (ctype, n, sn, celnode[ctype][n],	cd1=sdia*mdia,
						cd2=edia*mdia,DEND,cplam);

     /* make higher order dendrite */

    if (cumleng+shortlen < rad * 0.9)
      mak_gc_branch (ctype,n,nsn,cumleng+shortlen,angl,
			shsumang,xden1,yden1,ord+1,edia);

	/* make side branch */

   shortlen  = rrange (0.6,1.0) * (rad*1.2-cumleng)*1.0;
   if (shortlen+cumleng > rad) shortlen = rad - cumleng;
   shortang = angl * rrange(0.5,1) * 0.5; // * (rad*0.4)/shortlen;
   shsumang = sumangl - shortang * brdir;

   xden1 = cos(shsumang/DEG) * shortlen + cxden;
   yden1 = sin(shsumang/DEG) * shortlen + cyden;
   np = loc(nd(ctype,n,++celnode[ctype][n]),xden1,yden1,getn(ctype,DENDARBZ));
   nsn = celnode[ctype][n];
   edia = taperdia(ctype,n,sn,ctype,n,celnode[ctype][n],sdia,diaspc);

   make_celseg (ctype, n, sn, celnode[ctype][n],	cd1=sdia*brdia,
						cd2=edia*brdia,DEND,cplam);

     /* make higher order dendrite */

 //    if (cumleng+shortlen < rad * 1.0)
 //      mak_gc_branch (ctype,n,nsn,cumleng+shortlen,angl,
 //    			shsumang,xden1,yden1,ord+3,edia);

 }
}

/*-----------------------------------------*/

void mak_am_branch (int ctype, int n, int sn, 
	double cumleng, double angl,double sumangl,double cxden,double cyden,
	int ord,double sdia)

  /* Make a branch on an amacrine cell dendrite */

{
        double shortang, shortlen, shsumang;
        int brdir,nsn;
	double cd1, cd2;
	double xden1, yden1;
        double brdia,edia,mdia;
	double crm, cvrest, rad;
	double diabs, diaspc, stratd,cplam;
	node *np;

 if (ord < 5) {

   crm    = getn(ctype,NRM);
   cvrest = getn(ctype,VSTART);
   rad    = getn(ctype,DTREEDIA)/2;
   diabs  = getn(ctype,TAPERABS);
   diaspc = getn(ctype,TAPERSPC);
   stratd  = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
   cplam = -1;

   brdir = ffs();		/* branch direction */
   mdia   = rrange(0.6, 0.9);	/* dia of "main" side of branch */
   brdia  = rrange(0.3, 0.6);	/* dia of side branch */

	/* make first order branch */

   shortlen  = rrange (0.6,1) * (rad-cumleng) * (1-stratd) * 0.6;
   shortlen  = rrange (0.6,1) * (rad*1.5-cumleng) *  0.3;
   if (shortlen+cumleng > rad) shortlen = rad - cumleng;
   shortang = angl * rrange(0.3,1) * 0.5; // * (rad*0.4 + rad)/shortlen;
   shsumang = sumangl + shortang * brdir;

   xden1 = cos(shsumang/DEG) * shortlen + cxden;
   yden1 = sin(shsumang/DEG) * shortlen + cyden;
   np = loc(nd(ctype,n,++celnode[ctype][n]),xden1,yden1,getn(ctype,DENDARBZ));
   nsn = celnode[ctype][n];
   edia = taperdia(ctype,n,sn,ctype,n,celnode[ctype][n],sdia*mdia,diaspc);

   make_celseg (ctype, n, sn, celnode[ctype][n],cd1=sdia*mdia,
						cd2=edia,DEND,cplam);
     /* make higher order dendrite */

    if (cumleng+shortlen < rad * 0.8)
      mak_gc_branch (ctype,n,nsn,cumleng+shortlen,angl,
			shsumang,xden1,yden1,ord+1,edia);

	/* make side branch */

   shortlen  = rrange (0.6,1.0) * (rad-cumleng)*.8;
   if (shortlen+cumleng > rad) shortlen = rad - cumleng;
   shortang = angl * rrange(0.5,1) * 0.5; // * (rad*0.4)/shortlen;
   shsumang = sumangl - shortang * brdir;

   xden1 = cos(shsumang/DEG) * shortlen + cxden;
   yden1 = sin(shsumang/DEG) * shortlen + cyden;
   np = loc(nd(ctype,n,++celnode[ctype][n]),xden1,yden1,getn(ctype,DENDARBZ));
   nsn = celnode[ctype][n];
   edia = taperdia(ctype,n,sn,ctype,n,celnode[ctype][n],sdia*brdia,diaspc);
   make_celseg (ctype, n, sn, celnode[ctype][n],	cd1=sdia*brdia,cd2=edia,DEND,cplam);

     /* make higher order dendrite */

     if (cumleng+shortlen < rad * 1.0)
       mak_gc_branch (ctype,n,nsn,cumleng+shortlen,angl,
     			shsumang,xden1,yden1,ord+3,edia);

 }
}

/*-------------------------------------------------------------*/

#define cabldata(row,col) (*(ctcabldata+row*ccols+col))

void make_gc_comps(int ct, int cn, int region, int ct2)

/* For labeled presyn terminals in a cell region, 
   make single GC compartments with cellnum = presyn nodenum, below each cbp presynaptic terminal */

{    
    int n,rnode;
    int i, found, clines, ccols;
    double *ctcabldata;
    double dia = 1.0;
    double  dz = -2.0;
    node *npnt, *np;
    sphere *s;

   ccols = C_DENDN+1; 			// setup for "cabldata()"
   clines = cell_lines[ct];
   ctcabldata = cellmorphs[ct];

  for(n=1,npnt=nodepnt; npnt=foreach(npnt,ct,cn,-1,NULL,NULL,&rnode); npnt=npnt->next) {
      if (npnt->region==region) {
	    for (found=i=0; i<clines; i++) { 			// check if it has label
	           if ((cabldata(i,C_NODE)==rnode) && 
	               (cabldata(i,C_DENDN)>0)) {
			    // fprintf (stderr,"dendn %g\n",cabldata(i,C_DENDN));
	                  found = 1; break;
	           }
	    }
	    if (!found) continue;
	    np = loc(nd(ct2,n,soma),npnt->xloc,npnt->yloc,npnt->zloc+dz); /*locate the nodes */
	    s = make_sphere(nd(ct2, n, soma), dia, drm);
	    s->elabl = (char *)regname[region];
	    s->region = region;
            celnode[ct2][n] = 0;
	    n++;
      }
   }
  n--;
  ncfprintf (stderr,"# n%s %d\n",cname[ct2],n);
  set_ncel(ct2,n);
  setn(ct2,MAKE,1);
  setn(ct2,NMADE,n);
  made_gc_comps = 1;
}

/*-------------------------------------------------------------*/

/* Diameter of AXON_THIN region is set to "ath_dia" in anatfile.
   Then when read in, the value of ath_dia determines the
   diameter used.  The advantage of this is that we can test
   the effect of different diameters easily. */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void mak_real_cel(int ct, int cn, double xpos, double ypos, double thetax, double thetay, double thetaz, int flip, int dend1, int dend2, int dend3, int dend4, int dend5, double dend_cplam, int taper)

/* make real cell from anatomy file */

{
    int i, dend, clines=0,ccols=0, region;
    int nod, maxnod, aflag=0, somafound=0;
    int p, par, par_indx;
    double somaoffsetx=0, somaoffsety=0, somaoffsetz=0;
    double cd1, cd2, arbor_scale, dia_factor,d_factor;
    double *ctcabldata, par_dia, cplam;
    double ax, ay, az, dx, dy, dz, dx2, dy2, dz2;
    node *np;
    const char *anatfname;
    char *tanatfname;
    FILE *fp;


  if (notinit(ath_dia)) ath_dia = 0.6667;    /* default dia of thin segment */
  arbor_scale = getn(ct,ARBSCALE);	     /* xy scale for dend tree */
  if (arbor_scale == 0) arbor_scale = 1.0;
  dia_factor  = getn(ct,DENDDIA);	     /* dend dia (thickness) factor */
  if (dia_factor == 0) dia_factor = 1.0;

  if (fp = fopen(anatfiles[ct], "r")) {   /* Check for user-specified anatomy file */
    anatfname = anatfiles[ct];
    aflag = 0;
    fclose(fp);
  } else {

    //prepend confdir, check again
    tanatfname = emalloc(min(FILNAMSIZ,strlen(confdir)) + min(FILNAMSIZ,strlen(anatfiles[ct])) + 2);
    aflag = 1;
    sprintf(tanatfname, "%.50s/%.50s", confdir, anatfiles[ct]);
    anatfname = tanatfname;
  }
  if (ninfo>2) fprintf(stderr, "# %s %-3d Using anatomy file '%s'.\n", cname[ct], cn, anatfname);

  ctcabldata = fread (anatfname, &clines, &ccols);	/* read anatomy file */
  cellmorphs[ct] = ctcabldata;

  if (clines<=0) {
     fprintf (stderr,"Missing anatomy file '%s'.\n",anatfname);
     if (aflag) efree(tanatfname);
     clines = 0;
     return;
  };
  if (aflag) efree(tanatfname);
  cell_lines[ct] = clines;

  thetax *= PI / 180.0;
  thetay *= PI / 180.0;
  thetaz *= PI / 180.0;

  somafound = 0;
  somaoffsetx = somaoffsety = somaoffsetz = 0;
  for (i=0; i<clines; i++) {	/*find soma, and set offset (used below)*/
    if(cabldata(i,C_REG)==SOMA && !somafound) {
      somaoffsetx = cabldata(i,C_AX);
      somaoffsety = cabldata(i,C_AY);
      somaoffsetz = cabldata(i,C_AZ);
      somafound = 1;
    }
  }

  maxnod = 0;
  for (i=0; i<clines; i++) {				/* locate nodes first */

    dend   = int(cabldata(i,C_DENDN));
    region = int(cabldata(i,C_REG));

    if ((region==DEND) || (region==DEND_PROX)) 
    	if (!getn(ct,MAKE_DEND)) continue;		/* if skip dendrites */

    if ((region==AXON) || (region==AXON_THIN) ||
        (region==AXON_DIST) || (region==HILLOCK)) 
       if (!getn(ct,MAKE_AXON)) continue;		/* if skip axon */

    if (region==AXON_DIST)
        if (!getn(ct,MAKE_DIST)) continue;		/* if skip long axon */

    /* If any of "dend1-5" are set non-negative, and "dend_cplam" is negative */
    /* skip all dendrites except the one specified by dend1-5 */

    /* If "dend1-5" is set non-negative, and "dend_cplam" is positive */
    /* set complam for dendrites other than dend equal to dend_cplam */

    cplam = -1;
    if (ccols>C_DENDN) {  // if dendrite number column exists in morphology file
       if (dend_cplam <= 0) {
          if ((dend > soma) && 
              (!notinit(dend1) || !notinit(dend2) || !notinit(dend3) || !notinit(dend4) || !notinit(dend5)) && 
	      (dend1!=dend && dend2 != dend && dend3 != dend && dend4 != dend && dend5 != dend)) continue;	
        } else {
          if ((dend > soma) && 
              (!notinit(dend1) || !notinit(dend2) || !notinit(dend3) || !notinit(dend4) || !notinit(dend5)) && 
	      (dend1!=dend && dend2 != dend && dend3 != dend && dend4 != dend && dend5 != dend)) cplam = dend_cplam;	
        }
    }

    
    nod= int(cabldata(i,C_NODE));		/* node number of segment */
    dx = ax = (cabldata(i,C_AX)-somaoffsetx)*arbor_scale;	/* make soma zero offset */
    dy = ay = (cabldata(i,C_AY)-somaoffsety)*arbor_scale;
    dz = az = (cabldata(i,C_AZ)-somaoffsetz);			/* make soma zero relative z coord */

    // dx = (ax * cos(thetaz)  + ay * sin(thetaz));
    // dy = (ax * -sin(thetaz) + ay * cos(thetaz));
 
    // dz = (az * cos(thetay)  + dx * sin(thetay));
    // dx = (az * -sin(thetay) + dx * cos(thetay));

    // dz = (dz * cos(thetax)  + dy * sin(thetax));
    // dy = (dz * -sin(thetax) + dy * cos(thetax));

    dy  = (ay * cos(thetax) - az * sin(thetax));	   /* x-axis rotation */
    dz  = (ay * sin(thetax) + az * cos(thetax));

    dz2 = (dz * cos(thetay) - ax * sin(thetay));    	   /* y-axis rotation */
    dx  = (dz * sin(thetay) + ax * cos(thetay));
    dz = dz2;

    if (flip>0) dx = -dx;

    dx2 = (dx * cos(thetaz)  - dy * sin(thetaz));  /* do z rotation last after cell is tilted correctly */
    dy  = (dx * sin(thetaz)  + dy * cos(thetaz));
    dx = dx2;
 
    // fprintf (stderr,"nod %-6d ax %-6.4g ay %-6.4g az %-6.4g sx %-6.4g sy %-6.4g sz %-6.4g\n",
    // 		    nod, ax,ay,az,somaoffsetx,somaoffsety,somaoffsetz);

    np = loc(nd(ct,cn,nod),xpos+dx,ypos+dy,getn(ct,SOMAZ)+dz,region,dend); /*locate the nodes */

   /* For starburst amacrine, set additional dia factor */
   /*   but only for sb1 for now. */

    if (ct==sbac) {
        double r;

      if (notinit(dia_prox_rad)) dia_prox_rad = 25;
      if (notinit(dia_prox_factor)) dia_prox_factor = 2;

      r = sqrt(ax*ax + ay*ay);
      if (i>0 && r < dia_prox_rad && sbac_file=="morph_sb1")
        dia_factor *= dia_prox_factor;
      else
        dia_factor *= 1.0;
    } else
        dia_factor *= 1.0;

    //fprintf(stderr, "dsgc nod %d: x=%g  y=%g z=%g\n",
    //                      nod, xpos+dx, ypos+dy, getn(ct,SOMAZ)-az);

    /* somaoffset is added to C_AZ, making it zero for soma;  */
    /*   adding getn then makes z-coords. */


    if (taper>0) {
       par_indx = -1;
       par = cabldata(i,C_PAR);		 /* find parent node */
       for (p=0; p<clines; p++) {	
 	  if (par==cabldata(p,C_NODE)) {	 /* get index to parent node */
	     if (par!=cabldata(p,C_PAR))	 /* if not a sphere */
	         par_indx = p;
          }
       }
       if (par_indx == -1) par_dia = cabldata(i,C_DIA);  	     /* no taper */
       else 		   par_dia = cabldata(par_indx,C_DIA);       /* taper */

    /* don't taper dendrite going to soma */

       if (par_indx==-1) par_indx = 0;
       if (int(cabldata(i,C_REG)!=SOMA) && int(cabldata(par_indx,C_REG)==SOMA)) 
		    par_dia = cabldata(i,C_DIA); 	     /* no taper */ 
    }
    else (par_dia = -1);	/* no taper */

    /* make the dendritic tree and axon with taper */

    if ((region==DEND_DIST) || (region==DEND) ||
        (region==DEND_PROX))   d_factor = dia_factor;  /* scale dendrites only */
    else           d_factor = 1.0;

    make_celseg (ct, cn, nod, int(cabldata(i,C_PAR)),
                 cd1=cabldata(i,C_DIA)*d_factor,
                 cd2=par_dia*d_factor,
                 int(cabldata(i,C_REG)),cplam);
    if (maxnod < nod) maxnod = nod;

  } /* for (i;;) */

  celnode[ct][cn] = maxnod;
}

/*-----------------------------------------*/

int dendn_node (int ct, int dendn)

/* Look up first instance of a number in dendn column of morph file,
   and return its node number. This allows a more general index 
   useful to find a cell node in cells that have different node numbers. */
{
    int i, clines, ccols, dend_node, found;
    double *ctcabldata;

   ccols = C_DENDN+1; 
   clines = cell_lines[ct];
   ctcabldata = cellmorphs[ct];
   for (found=i=0; i<clines; i++) {
	   if (cabldata(i,C_DENDN)==dendn) {
		  dend_node = cabldata(i,C_NODE);
		  found = 1;
		  break;
	   }
   }
   if (!found) {
	 dend_node = NULLVAL;
	//ncfprintf(stderr,"retsim: can't find label %d\n",dendn);
   }
   return (dend_node);
}

/*-----------------------------------------*/

int make_axon_branch (int ct, int cn, int startnode, int newnode, double zoffs, double zlen, double angle, double zangle,
			double brdia, double len, double vdia, int nsegs, int reg1, int reg2) 

{

   int i,parent;
   double dx,dy,dz;
   double xpos,ypos,zpos,seglen;
   double cplam;
   node *np;


  xpos = nd(ct,cn,startnode)->xloc;
  ypos = nd(ct,cn,startnode)->yloc;
  zpos = nd(ct,cn,startnode)->zloc + zoffs;

  parent = startnode;
  for (i=1; i<=nsegs; i++) {
      if (i==1) seglen = len + zlen;
      else      seglen = len;
      dx = cos(angle/DEG) * cos(zangle/DEG)* seglen;
      dy = sin(angle/DEG) * cos(zangle/DEG)* seglen;
      dz = sin(zangle/DEG) * seglen;
      // fprintf (stderr,"axbr %g x %g y %g\n",axbr,x,y);

      xpos += dx; 
      ypos += dy;
      zpos += dz;

      np = loc(nd(ct,cn,newnode), xpos,ypos,zpos);
      make_celseg (ct, cn, parent,  newnode, brdia, brdia, reg1, cplam=-1); /* make 1 seg of axon branch */
      if (vdia>0)
	make_celseg (ct, cn, newnode, newnode, vdia,  0,   reg2, cplam=-1); /* make varicosity */
      parent = newnode++;
  } 
  return parent;
}

/*-----------------------------------------*/

void mak_bp(int ct, int cn, int nden, double xpos, double ypos, double thetax, 
		double thetay, double thetaz, int flip, int morph)

{
   int i,sn,pd,dendr,nsegs,naxbr,axnode,oldnode,naxseg;
   int next,next2;
   int lob_append;
   double cd1, cd2, somadia, axbrdia, seglen, vdia;
   double dtheta, start_angle, sumangle, brangle, zangle;
   double xden, yden, somaz, axarbz, axbrzoffs, axnodedist, axloc;
   double crm,cvrest,rad;
   double dia, diabs,diaspc,edia;
   double stratd, cplam;
   sphere *s;
   node *np;
 
   crm    = getn(ct,NRM);
   cvrest = getn(ct,VSTART);
   rad    = getn(ct,DTREEDIA)/2;
   diabs  = getn(ct,TAPERABS);
   diaspc = getn(ct,TAPERSPC);
   stratd  = getn(ct,STRATDIA); /* inside radius of stratified dend. annulus */
   cplam = -1;

  switch (morph) {
  
    case MORPH_REAL:  	 	  /* realistic morphology */

      if (notinit(den_dia)) den_dia = 0.3;  /* set dendrite thin dia in file */
      if (notinit(taper)) taper = 1;	  /* make cables with taper */
      mak_real_cel(ct,cn,xpos,ypos,thetax,thetay,thetaz,flip,
		    bp_dend1,bp_dend1,bp_dend1,bp_dend1,bp_dend1,cplam,taper);
      break;

    case MORPH_A1: 		/* bp-like artif morphology */

      np = loc(nd(ct,cn,soma), xpos,ypos,getn(ct,SOMAZ));
      make_celseg (ct, cn, soma, soma, cd1=getn(ct,SOMADIA),cd2=0,SOMA,cplam);
      np = loc(nd(ct,cn,axtrm), xpos,ypos,getn(ct,AXARBZ));
      make_celseg (ct, cn, soma, axtrm, cd1=0.5 ,cd2=0.4, AXON,cplam); /* connect to soma */
      make_celseg (ct, cn, axtrm, axtrm, cd1=1,cd2=0, AXON,cplam);   /* sphere at tip */
      np = loc(nd(ct,cn,dendr=2), xpos,ypos,getn(ct,DENDARBZ));
      make_celseg (ct, cn, soma, dendr, cd1=0.5 ,cd2=0.3, DEND,cplam); /* connect to soma */
      celnode[ct][cn] = 2;
      break;
  
    case MORPH_A2: 		/* just soma, artif morphology */

      if (getn(ct,BIOPHYS)==0) cplam = 0.5;		/* if no biophys, just one comp */
      np = loc(nd(ct,cn,soma), xpos,ypos,getn(ct,SOMAZ));
      make_celseg (ct, cn, soma, soma, cd1=getn(ct,SOMADIA),cd2=0,SOMA,cplam);
      np = loc(nd(ct,cn,axtrm), xpos,ypos,getn(ct,AXARBZ));
      make_celseg (ct, cn, soma, axtrm, cd1=0.5 ,cd2=0.4, AXON,cplam); /* connect to soma */
      make_celseg (ct, cn, axtrm, axtrm, cd1=1,cd2=0, AXOND,cplam);   /* sphere at tip */
      celnode[ct][cn] = 1;
      break;

    case MORPH_T4: 		/* artif morphology with branched axon, sparse, stratified, t4 */

      // somadia = getn(ct,SOMADIA);
      // axarbz = getn(ct,AXARBZ);
      somadia = 8;
      axbrzoffs = 0;		// Z offset for start of branches
      somaz = getn(ct,SOMAZ);	// zloc of soma
      axarbz = somaz - 30;	// bottom of proximal axon

      naxbr = 4; 		// number of axon branches
      nsegs = 2;		// number of segments per branch
      seglen = 5;		// length of each segment
      start_angle = thetaz;	// angle of first branch
      axbrdia = 0.2;		// diameter of axon branches
      brangle = 0;		// angle at branchpoint
      zangle = -2;

      np = loc(nd(ct,cn,soma), xpos,ypos,somaz);
      make_celseg (ct, cn, soma, soma, cd1=somadia,cd2=0,SOMA,cplam);
      np = loc(nd(ct,cn,dendr=2), xpos,ypos,somaz+getn(ct,DENDARBZ));
      make_celseg (ct, cn, soma, dendr, cd1=0.5 ,cd2=0.3, DEND,cplam); /* make dends, connect to soma */
      celnode[ct][cn] = 2;

      oldnode = soma;
      axnode = 10;
      naxseg = 20;
      axnodedist = (somaz - axarbz - axbrzoffs) / 5;  // make 5 nodes along proximal axon to just above terminal
      axloc = somaz-axnodedist;
      for (i=0; i<naxseg; i++) {
        np = loc(nd(ct,cn,axnode), xpos,ypos,axloc);
        make_celseg (ct, cn, oldnode, axnode, cd1=0.5, cd2=0.5, R4+i, cplam); /* make prox axon, connect to soma */
        axloc -= axnodedist;
        oldnode = axnode++;
        if (axloc < axarbz) break;
      }

      dtheta = 360 / naxbr;
      sumangle = start_angle;
      for (i=1; i<=naxbr; i++) {
         next  = make_axon_branch (ct,cn,oldnode, i*100, -axbrzoffs, 0, sumangle, zangle= -20, axbrdia, seglen=5, vdia=2.0, nsegs=1, R9, R10);
         next2 = make_axon_branch (ct,cn,next, next+1,  0, 0, sumangle, zangle= -87, axbrdia, seglen=5, vdia=2.0, nsegs=2, R9, R10);
         next  = make_axon_branch (ct,cn,next, next2+1, 0, 0, sumangle, zangle= -20, axbrdia, seglen=5, vdia=2.0, nsegs=1, R9, R10);
         next2 = make_axon_branch (ct,cn,next, next+1,  0, -1, sumangle, zangle= -87, axbrdia, seglen=5, vdia=2.0, nsegs=2, R9, R10);
         sumangle += dtheta;
         celnode[ct][cn] += nsegs*2;
      }
      break;

    case MORPH_T8: 		/* artif morphology with branched axon, sparse, stratified, t8 */

      // somadia = getn(ct,SOMADIA);
      // axarbz = getn(ct,AXARBZ);
      somadia = 6;
      axbrzoffs = 5;		// Z offset for start of branches
      somaz = getn(ct,SOMAZ);
      axarbz = somaz - 80;		// bottom of proximal axon

      naxbr = 3; 		// number of axon branches
      nsegs = 2;		// number of segments per branch
      seglen = 8;		// length of each segment
      start_angle = thetaz;	// angle of first branch
      axbrdia = 0.2;		// diameter of axon branches
      brangle = 30;		// angle at branchpoint
      zangle = -3;

      np = loc(nd(ct,cn,soma), xpos,ypos,somaz);
      make_celseg (ct, cn, soma, soma, cd1=somadia,cd2=0,SOMA,cplam);
      np = loc(nd(ct,cn,dendr=2), xpos,ypos,getn(ct,DENDARBZ));
      make_celseg (ct, cn, soma, dendr, cd1=0.5 ,cd2=0.3, R1,cplam); /* make dends, connect to soma */
      celnode[ct][cn] = 2;

      oldnode = soma;
      axnode = 10;
      naxseg = 20;
      axnodedist = (somaz - axarbz - axbrzoffs) / 5;  // make 5 nodes along proximal axon to just above terminal
      axloc = somaz-axnodedist;
      for (i=0; i<naxseg; i++) {
        np = loc(nd(ct,cn,axnode), xpos,ypos,axloc);
        make_celseg (ct, cn, oldnode, axnode, cd1=0.5, cd2=0.5, R4+i, cplam); /* make prox axon, connect to soma */
        axloc -= axnodedist;
        oldnode = axnode++;
        if (axloc < axarbz) break;
      }

      dtheta = 360 / naxbr;
      sumangle = start_angle;
      for (i=1; i<=naxbr; i++) {
         next  = make_axon_branch (ct,cn,oldnode, i*100, -axbrzoffs, 0, sumangle, zangle, axbrdia, seglen, vdia=2.0, nsegs, R9, R10);
         next2 = make_axon_branch (ct,cn,next, next+1, 0,  0, sumangle - brangle, zangle, axbrdia, seglen, vdia=2.0, nsegs+1, R9, R10);
         next2 = make_axon_branch (ct,cn,next, next2+1, 0, 0, sumangle + brangle, zangle, axbrdia, seglen, vdia=2.0, nsegs+1, R9, R10);
         sumangle += dtheta;
         celnode[ct][cn] += nsegs*3+2;
      }
      break;

  } /* switch morph */
}

/*-----------------------------------------*/

void mak_hz(int ctype, int cn, int nden, double xpos, double ypos, 
		double thetax, double thetay, double thetaz, int flip, int morph)

{  
  int i,region,sn;
  double angl,sumangl;
  double xden, yden;
  double leng, thisangl;
  double crm,cvrest;
  double rad,dia,diabs,diaspc,cplam;
  sphere *s;
  node *np;

  
  if (morph==MORPH_REAL) { 		/* realistic morphology */
    if (notinit(den_dia)) den_dia = 0.3;   /* set dendrite thin dia in file */
    if (notinit(taper)) taper = 1;	  /* make cables with taper */
    mak_real_cel(ctype,cn,xpos,ypos,thetax,thetay,thetaz,flip,
		    hz_dend1,hz_dend2,hz_dend3,hz_dend4,hz_dend5,cplam,taper);
  }
  else if (morph==MORPH_A1) {
    rad    = getn(ctype,DTREEDIA)/2;
    diabs  = getn(ctype,TAPERABS);	/* starting taper diameter */
    diaspc = getn(ctype,TAPERSPC);	/* taper space constant */
    celnode[ctype][cn] = 0;
    cplam = -1;

    np = loc(nd(ctype,cn,soma),xpos,ypos,getn(ctype,SOMAZ));
    dia=getn(ctype,SOMADIA);
    make_celseg(ctype,cn,soma,soma,dia,0,region=SOMA,cplam);

    angl = 360 / nden;
    sumangl = 0;
    for (i=1; i<=nden; i++) {
      thisangl = angl * rrange (0.8,1.2);
      sumangl += thisangl;
      leng = rad * rrange (0.8, 1.2) * 0.2;
//printf ("angle %g length %g\n", angl, leng);

      xden = cos(sumangl/DEG) * leng + xpos;
      yden = sin(sumangl/DEG) * leng + ypos;
//printf ("loc %g %gi\n", xden, yden);
      np = nd(ctype,cn,++celnode[ctype][cn]);
      np = loc(np,xden,yden,getn(ctype,DENDARBZ)+2);

      sn = celnode[ctype][cn];		/* sequence node */   
      celnode[ctype][cn] = taperdistden (ctype,cn,0,celnode[ctype][cn],
				diabs,diaspc,DEND_PROX, celnode[ctype][cn]);

      if (leng < rad * 0.5) {                  /* make first bifurcation */
        double shortang, shortlen, shsumang;
        double xden1, yden1;

        shortlen  = rrange(0.5,1) * (rad-leng) * 0.5;
	if (shortlen+leng > rad) shortlen = rad - leng;
        shortang = angl * rrange(0.4, 1.0) * 0.5;
        shsumang = sumangl - shortang;		/* less angle */

        xden1 = cos(shsumang/DEG) * shortlen + xden;
        yden1 = sin(shsumang/DEG) * shortlen + yden;
        np = nd(ctype,cn,++celnode[ctype][cn]);
        np = loc(np,xden1,yden1,getn(ctype,DENDARBZ)+2);
        celnode[ctype][cn] = taperdistden(ctype,cn,sn,celnode[ctype][cn],
				diabs,diaspc,DEND_DIST,celnode[ctype][cn]);

        shortlen  = rrange(0.5,1) * (rad-leng) * 0.5;
	if (shortlen+leng > rad) shortlen = rad - leng;
        shortang = angl * rrange(0.4, 1.0) * 0.5;
        shsumang = sumangl + shortang;		/* greater angle */

        xden1 = cos(shsumang/DEG) * shortlen + xden;
        yden1 = sin(shsumang/DEG) * shortlen + yden;
        np = nd(ctype,cn,++celnode[ctype][cn]);
        np = loc(np,xden1,yden1,getn(ctype,DENDARBZ)+2);
        celnode[ctype][cn] = taperdistden(ctype,cn,sn,celnode[ctype][cn],
				diabs,diaspc,DEND_DIST,celnode[ctype][cn]);
      }
    }     /* for (i=0;i<nden;) */
  }
}


/*-----------------------------------------*/

void mak_aii(int ctype, int cn, int nden, double xpos, double ypos, 
		double thetax, double thetay, double thetaz, int flip, int morph)

{
   int i,sn,pd;
   int lob_append;
   double cd1, cd2;
   double angl,sumangl;
   double xden, yden;
   double leng, thisangl;
   double crm,cvrest,rad;
   double dia, diabs,diaspc,edia;
   double stratd, cplam;
   sphere *s;
   node *np;
 
   crm    = getn(ctype,NRM);
   cvrest = getn(ctype,VSTART);
   rad    = getn(ctype,DTREEDIA)/2;
   diabs  = getn(ctype,TAPERABS);
   diaspc = getn(ctype,TAPERSPC);
   stratd  = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
   cplam = -1;

  if (morph==MORPH_REAL) { 		/* realistic morphology */
    if (notinit(taper)) taper = 1;	  /* make cables with taper */
    den_dia = 0.3;			/* set dendrite thin dia in file */
    mak_real_cel(ctype,cn,xpos,ypos,thetax,thetay,thetaz,flip,
		    aii_dend1,aii_dend2,aii_dend3,aii_dend4,aii_dend5,cplam,taper);
  }
  else if (morph==MORPH_A1) {		/* aii-like artif morphology */
   celnode[ctype][cn] = pd = 0;
   np = loc(nd(ctype,cn,soma),xpos,ypos,getn(ctype,SOMAZ));
   s = make_sphere(np,dia=10);
   make_celseg (ctype, cn, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);

   // Add lobular appendage 

   np = loc(nd(ctype,cn,lob_append= ++celnode[ctype][cn]),xpos+10,ypos,getn(ctype,SOMAZ-5));
   make_celseg (ctype, cn, lob_append, lob_append, cd1=6.0,cd2=0,DEND_DIST,cplam);
   make_celseg (ctype, cn, soma, lob_append, cd1=0.3,cd2=0.3,DEND_PROX,cplam);

   angl = 360 / nden;
   sumangl = 0;
   for (i=1; i<=nden; i++) {
     thisangl = angl * rrange (0.9, 1.1);
     sumangl += thisangl;
     leng = rrange(0.5,1.5) * rad * stratd;
//printf ("angle %g length %g\n", angl, leng);

     xden = cos(sumangl/DEG) * leng + xpos;
     yden = sin(sumangl/DEG) * leng + ypos;
//printf ("loc %g %gi\n", xden, yden);
     np = loc(nd(ctype,cn,++celnode[ctype][cn]),xden,yden,getn(ctype,DENDARBZ));
     sn = celnode[ctype][cn];
     edia = taperdia(ctype,cn,pd,ctype,cn,celnode[ctype][cn],diabs,diaspc);
     make_celseg (ctype, cn, pd, celnode[ctype][cn],cd1=diabs,cd2=edia,DEND_PROX,cplam);

     if (leng < rad * 0.5) {                  /* make secondary dendrite */
	   mak_am_branch (ctype,cn,sn, leng,angl,sumangl,xden,yden,2,edia);
      }

   }
  } /* make artificial cell morphology */
}

/*-----------------------------------------*/

void mak_amac(int ctype, int n, int nden, double xpos, double ypos, 
		double thetax, double thetay, double thetaz, int flip, int morph)

{
   int i,j, pd,sn, nvaricos;
   double ax, ay, az, dx, dy, dz, dx2, dy2;
   double cd1, cd2, dvaricos,seglen;
   double angl,sumangl;
   double xden, yden, xden2, yden2;
   double leng, leng2, thisangl;
   double crm,cri,cvrest,rad;
   double diabs,diaspc,edia,cplam;
   double stratd,dendarbor;
   node *np;

   if (notinit(dia_prox_rad))    dia_prox_rad    = 25; /* radius for prox dia */
   if (notinit(dia_prox_factor)) dia_prox_factor = 2; /* dia factor */
   if (notinit(taper)) taper = 1; /* make cables with taper */

   if      (morph==MORPH_REAL) {
              mak_real_cel (ctype,n,xpos,ypos,thetax,thetay,thetaz,flip,
                           am_dend1,am_dend2,am_dend3,am_dend4,am_dend5,
			   am_dend_cplam,taper);  /* use morph file */
						   }
  else if (morph==MORPH_A1) {

    crm    = getn(ctype,NRM);
    cri    = getn(ctype,NRI);
    cvrest = getn(ctype,VSTART);
    rad    = getn(ctype,DTREEDIA)/2;
    diabs  = getn(ctype,TAPERABS);
    diaspc = getn(ctype,TAPERSPC);
    stratd = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
    cplam = -1;			   // use am_dia_cplam if wanted here, see sbac_dend_cplam

    thetax /= DEG;
    thetay /= DEG;
    thetaz /= DEG;

    celnode[ctype][n] = pd = 0;
    np = loc(nd(ctype,n,soma),xpos,ypos,getn(ctype,SOMAZ));
    make_celseg (ctype, n, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
 
    angl = 360 / nden;
    sumangl = 0;
    for (i=1; i<=nden; i++) {
//      thisangl = angl * rrange (0.9, 1.1);
      thisangl = angl;
      // leng = rrange(0.9,1.1) * rad * stratd;
      leng = rad * stratd;
//fprintf (stderr,"angle %g length %g nden %d\n", angl, leng, nden);

      xden = cos(sumangl/DEG) * leng;
      yden = sin(sumangl/DEG) * leng;

      dx = ax = xden;
      dy = ay = yden;
      dz = az = 0;	/* make soma zero relative z coord */
  
      dz  = (az * cos(thetay)  + ax * sin(thetay));	/* do xy rotation first to get cell tilted correctly */
      dx2 = (az * -sin(thetay) + ax * cos(thetay));

      dz  = (dz * cos(thetax)  + ay * sin(thetax));
      dy2 = (dz * -sin(thetax) + ay * cos(thetax));

      dx = (dx2 *  cos(thetaz)  + dy2 * sin(thetaz));    /* do z rotation last after cell is tilted correctly */
      dy = (dx2 * -sin(thetaz)  + dy2 * cos(thetaz));
 
    // fprintf (stderr,"nod %-6.4d ax %-6.4g ay %-6.4g az %-6.4g sx %-6.4g sy %-6.4g sz %-6.4g\n",
    // 		    nod, ax,ay,az,somaoffsetx,somaoffsety,somaoffsetz);

      np = loc(nd(ctype,n,++celnode[ctype][n]),xpos+dx,ypos+dy,getn(ctype,DENDARBZ)+dz);
      edia = taperdia(ctype,n,pd,ctype,n,celnode[ctype][n],diabs,diaspc);
      make_celseg (ctype, n, pd, celnode[ctype][n],cd1=diabs,cd2=diabs,DEND,cplam);

      if ((seglen=getn(ctype,SEGLEN))<=0) seglen = 10;  /* dendritic segment length */
      dendarbor = rad-leng;
      nvaricos = dendarbor/seglen;               /* varicosities on a dendrite spaced seglen apart */
      dvaricos=1.5;                              /* dia of varicosity */
      // leng2 = rrange(0.9,1.1) * rad/nvaricos; /* inter-varicosity distance */
      leng2 = int(dendarbor/nvaricos);           /* inter-varicosity distance */
      xden2 = xden;
      yden2 = yden;
     
      for (j=1; j<=nvaricos; j++) {             /* make varicosities on each dendrite */

        xden2 += cos(sumangl/DEG) * leng2;
        yden2 += sin(sumangl/DEG) * leng2;

        dx = ax = xden2;			/* make soma zero offset */
        dy = ay = yden2;
        dz = az = 0;	/* make soma zero relative z coord */
  
        dz  = (az * cos(thetay)  + ax * sin(thetay));	/* do xy rotation first to get cell tilted correctly */
        dx2 = (az * -sin(thetay) + ax * cos(thetay));

        dz  = (dz * cos(thetax)  + ay * sin(thetax));
        dy2 = (dz * -sin(thetax) + ay * cos(thetax));

        dx = (dx2 *  cos(thetaz)  + dy2 * sin(thetaz));    /* do z rotation last after cell is tilted correctly */
        dy = (dx2 * -sin(thetaz)  + dy2 * cos(thetaz));
 
        sn = celnode[ctype][n];
        np = loc(nd(ctype,n,++celnode[ctype][n]),xpos+dx,ypos+dy,getn(ctype,DENDARBZ)+dz);
        make_celseg (ctype, n, celnode[ctype][n], celnode[ctype][n],cd1=dvaricos,cd2=dvaricos,VARICOS,cplam);
        make_celseg (ctype, n, sn, celnode[ctype][n],cd1=diabs,cd2=diabs,DEND,cplam);
      }
      sumangl += thisangl;
    }     /* for (i=0;i<nden;) */
  }
  else if (morph==MORPH_A2) {

    crm    = getn(ctype,NRM);
    cvrest = getn(ctype,VSTART);
    rad    = getn(ctype,DTREEDIA)/2;
    diabs  = getn(ctype,TAPERABS);
    diaspc = getn(ctype,TAPERSPC);
    stratd  = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
    cplam = -1;			   // use am_dia_cplam if wanted here, see sbac_dend_cplam

    celnode[ctype][n] = pd = 0;
    np = loc(nd(ctype,n,soma),xpos,ypos,getn(ctype,SOMAZ));
    make_celseg (ctype, n, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
 
    angl = 360 / nden;
    sumangl = 0;
    for (i=1; i<=nden; i++) {
      thisangl = angl * rrange (0.9, 1.1);
      sumangl += thisangl;
      leng = rrange(0.9,1.1) * rad * stratd;
//printf ("angle %g length %g\n", angl, leng);

      xden = cos(sumangl/DEG) * leng + xpos;
      yden = sin(sumangl/DEG) * leng + ypos;
//printf ("loc %g %gi\n", xden, yden);
      np = loc(nd(ctype,n,++celnode[ctype][n]),xden,yden,getn(ctype,DENDARBZ));
      sn = celnode[ctype][n];
      edia = taperdia(ctype,n,pd,ctype,n,celnode[ctype][n],diabs,diaspc);
      make_celseg (ctype, n, pd, celnode[ctype][n],cd1=diabs,cd2=edia,DEND,cplam);
 
      if (leng < rad * 0.5) {                  /* make secondary dendrite */
	   mak_am_branch (ctype,n,sn, leng,angl,sumangl,xden,yden,2,edia);
      }
    }     /* for (i=0;i<nden;) */
  }
  else fprintf (stderr,"retsim: mak_amac, unknown morph, %d\n",morph);
}

/*-----------------------------------------*/

void mak_a17(int ctype, int n, int nden, double xpos, double ypos, 
		double thetax, double thetay, double thetaz, int flip, int morph)

{
   int i,j, pd,sn, nvaricos;
   double cd1, cd2, dvaricos;
   double angl,sumangl;
   double xden, yden, xden2, yden2;
   double leng, leng2, thisangl;
   double crm,cvrest,rad;
   double diabs,diaspc,edia,edia2,cplam;
   double stratd;
   node *np;

   if (notinit(dia_prox_rad))    dia_prox_rad    = 25; /* radius for prox dia */
   if (notinit(dia_prox_factor)) dia_prox_factor = 2; /* dia factor */
   if (notinit(taper)) taper = 1; /* make cables with taper */

   if      (morph==MORPH_REAL) {
              mak_real_cel (ctype,n,xpos,ypos,thetax,thetay,thetaz,flip,
                           am_dend1,am_dend2,am_dend3,am_dend4,am_dend5,
			   am_dend_cplam,taper);  /* use morph file */
						   }
  else if (morph==MORPH_A1) {

    crm    = getn(ctype,NRM);
    cvrest = getn(ctype,VSTART);
    rad    = getn(ctype,DTREEDIA)/2;
    diabs  = getn(ctype,TAPERABS);
    diaspc = getn(ctype,TAPERSPC);
    stratd  = getn(ctype,STRATDIA); /* inside radius of stratified dend. annulus */
    cplam = -1;			   // use am_dia_cplam if wanted here, see sbac_dend_cplam

    celnode[ctype][n] = pd = 0;
    np = loc(nd(ctype,n,soma),xpos,ypos,getn(ctype,SOMAZ));
    make_celseg (ctype, n, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
 
    angl = 360 / nden;
    sumangl = 0;
    for (i=1; i<=nden; i++) {
//      thisangl = angl * rrange (0.9, 1.1);
      thisangl = angl;
      leng = rrange(0.9,1.1) * rad * stratd;
//fprintf (stderr,"angle %g length %g nden %d\n", angl, leng, nden);

      xden = cos(sumangl/DEG) * leng + xpos;
      yden = sin(sumangl/DEG) * leng + ypos;
//printf ("loc %g %gi\n", xden, yden);
      np = loc(nd(ctype,n,++celnode[ctype][n]),xden,yden,getn(ctype,DENDARBZ));
      edia = taperdia(ctype,n,pd,ctype,n,celnode[ctype][n],diabs,diaspc);
      make_celseg (ctype, n, pd, celnode[ctype][n],cd1=diabs,cd2=edia,DEND,cplam);

       
      nvaricos=20;                            /* number of varicosities on a dendrite */
      dvaricos=1.5;                           /* dia of varicosity */
      leng2 = rrange(0.9,1.1) * rad/nvaricos; /* inter-varicosity distance */
      xden2 = xden;
      yden2 = yden;
     
      for (j=1; j<=nvaricos; j++) {             /* make varicosities on each dendrite */
        xden2 += cos(sumangl/DEG) * leng2;
        yden2 += sin(sumangl/DEG) * leng2;

      sn = celnode[ctype][n];
      np = loc(nd(ctype,n,++celnode[ctype][n]),xden2,yden2,getn(ctype,DENDARBZ));
      make_celseg (ctype, n, celnode[ctype][n], celnode[ctype][n],cd1=dvaricos,cd2=dvaricos,VARICOS,cplam);
      make_celseg (ctype, n, sn, celnode[ctype][n],cd1=edia,cd2=edia,DEND,cplam);
      }
      sumangl += thisangl;
      
      //      leng2 = rrange(0.9,1.1) * rad;
      //      xden2 = cos(sumangl/DEG) * leng2 + xpos;
      //      yden2 = sin(sumangl/DEG) * leng2 + ypos;
      //      sn = celnode[ctype][n];
      //      np = loc(nd(ctype,n,++celnode[ctype][n]),xden2,yden2,getn(ctype,DENDARBZ));
      //      edia2 = taperdia(ctype,n,sn,ctype,n,celnode[ctype][n],diabs,diaspc);
      //      make_celseg (ctype, n, sn, celnode[ctype][n],cd1=edia,cd2=edia2,DEND,cplam);
      
      // if (leng < rad * 0.5) {                  /* make secondary dendrite */
      //       mak_am_branch (ctype,n,sn, leng,angl,sumangl,xden,yden,2,edia);
      // }
    }     /* for (i=0;i<nden;) */
  }
  else fprintf (stderr,"retsim: mak_a17, unknown morph, %d\n",morph);
}

/*-------------------------------------------------------------*/

void mak_gc(int ctype, int n, int nden, double xpos, double ypos, 
		double thetax, double thetay, double thetaz, int flip, int morph)

/* make gc dendrites according to parameters */

/* diabs  = starting size of dendrite for taper */
/* diaspc = space constant of dendrite taper */

{
   int i, brnum;
   int pd,sn,csteps,nconn;
   double d1,d2, cd1,cd2;
   double angl,sumangl;
   double xden,yden;
   double leng,thisangl;
   double diabs,diaspc,edia,dend_dia;
   double stratd,cplam;
   double crm, cvrest;
   double rad,pdz;
   node *np;

  if (notinit(gc_dend_cplam)) gc_dend_cplam = -1; 
  if (notinit(thetaz)) thetaz = 0; 
  if (notinit(taper)) taper = 1;	  /* make cables with taper */
  if (morph==MORPH_REAL) { 		/* realistic gc morphology */
    mak_real_cel(ctype,n,xpos,ypos,thetax,thetay,thetaz,flip,gc_dend1,
			gc_dend2,gc_dend3,gc_dend4,gc_dend5,gc_dend_cplam,taper);
  }
  else if (morph==MORPH_A1) {		/* alpha-like gc artif morphology */
    crm    = getn(ctype,NRM);
    cvrest = getn(ctype,VSTART);
    rad    = getn(ctype,DTREEDIA)/2;
    diabs  = getn(ctype,TAPERABS);
    diaspc = getn(ctype,TAPERSPC);
    stratd  = getn(ctype,STRATDIA); /* inside rad of stratified dend. annulus */
    cplam = -1;

    celnode[ctype][n] = 0;
    np = loc(nd(ctype,n,soma),xpos,ypos,getn(ctype,SOMAZ));

    /* make the dendritic tree and axon */
    make_celseg (ctype, n, soma, soma, d1=getn(ctype,SOMADIA),d2=0,SOMA,cplam);

    /* make primary dendrite */
    pd = celnode[ctype][n] = 1;		/* soma is start of primary dend */
    pdz = (getn(ctype,DENDARBZ) + getn(ctype,SOMAZ)) * 0.5;

    np = loc(nd(ctype,n,pd),xpos,ypos,pdz);
    celnode[ctype][n] = taperdistden (ctype,n,0,pd,diabs,diaspc,DEND_PROX,pd);

    angl = 360 / nden;
    sumangl = 0;
    for (i=1; i<=nden; i++) {
      thisangl = angl * rrange (0.9, 1.1);
      sumangl += thisangl;
      leng = rrange(0.3,1) * rad * stratd;
//printf ("angle %g length %g\n", angl, leng);

      xden = cos(sumangl/DEG) * leng + xpos;
      yden = sin(sumangl/DEG) * leng + ypos;
      np = loc(nd(ctype,n,++celnode[ctype][n]),xden,yden,getn(ctype,DENDARBZ));
      sn = celnode[ctype][n];
//printf ("sn %d loc %g %gi\n", sn, xden, yden);
      edia = taperdia(ctype,n,pd,ctype,n,celnode[ctype][n],diabs,diaspc);
      make_celseg (ctype, n, pd, celnode[ctype][n],cd1=diabs,cd2=edia,DEND,cplam);

      if (leng < rad * 0.5)   		/* make secondary dendrite */
 mak_gc_branch (ctype,n,sn, leng,angl,sumangl,xden,yden,2,edia);
    }     /* for (i=0;i<nden;) */

    /* This is the new part for "extending" gc branches */
/*
    csteps = 4;
    for (j=0; j<csteps; j++) {
      printf("cell step %d\n",j);
      nconn = node [ctype][n][0] -> numconn;
      for (i=1; i<=nconn; i++) {
        brnum = node [ctype][n][0] -> i;
        if (element brnum->ntype != ntype(cable))
	  continue;
        printf("den %d of %d, elem %d\n",i,nconn,brnum);
        xden = node [ctype][n][soma] -> xloc;
        yden = node [ctype][n][soma] -> yloc;
        extend_gc_branch (brnum,brnum,ctype,n,0,leng,angl,2);
      }
    }
*/

  } 	/* morph==MORPH_A1) */

  else if (morph==MORPH_A2) {			/* ds-like gc artif morphology */
  }

  else if (morph==MORPH_SIMP) {			/* simple straight cable */
      double somadia;

    np = loc(nd(ctype,n,soma),xpos+40,ypos,getn(ctype,SOMAZ));
    make_celseg(ctype, n, soma, soma, somadia=8, somadia=8, SOMA, cplam);
    np = loc(nd(ctype,n,1),xpos+15,ypos,getn(ctype,DENDARBZ));
    make_celseg(ctype, n, soma, 1, 0.5, 0.5, DEND, cplam);
  }
  else if (morph==MORPH_SOMA) {			/* simple sphere */
      double somadia;
    np = loc(nd(ctype,n,soma),xpos,ypos,getn(ctype,SOMAZ));
    make_celseg(ctype, n, soma, soma, somadia=8, somadia=8, SOMA, cplam);
  }
}

/*-------------------------------------------------------------*/

void makcell (int ctype, int cn, double xpos, double ypos, double thetay, double thetaz, int flip)

/* Procedure to make the beginnings of a neuron. */

{
    int ct,ccount=0;
    int nden,aden;
    double crm,cvrest,cvrev;
    double rad,dia,cd1,cd2;
    double bginten;
    double thetax, cplam;
    node *np;
    photorec *p;
    sphere *s;

  ct = ctype;
  rad = getn(ctype,DTREEDIA) / 2;
  if (rad <= 0) rad=1;
  crm    = getn(ctype,NRM);
  cvrev  = getn(ctype,VREV);
  cvrest = getn(ctype,VSTART);
  cplam = -1;
  thetax = 0;
  if (notinit(flip)) flip = 0;
  if (ctype==xrod) { 		/* procedure to make rod */

    if (notinit(rod_maxcond)) rod_maxcond = 10e-12;

    /* simple rod: one compartment: */
    np = loc(nd(xrod,cn,soma), xpos,ypos,0);
    make_celseg (xrod, cn, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
    p = make_rod(np,dia=2); p->xpos=xpos; p->ypos=ypos; p->timec1=0.3; 
	p->maxcond=10e-12; p->photnoise=pnoise; p->darknoise=dnoise; p->linit=1;
	if (!notinit(rod_timec)) p->timec1=rod_timec; 
	if (!notinit(rod_loopg)) p->loopg=rod_loopg;
	//p->dkrseed = cn^737591; p->phrseed = cn^451234; /* set from command line */
    np = loc(nd(xrod,cn,axtrm),xpos,ypos,getn(xrod,AXARBZ));
    make_celseg (xrod, cn, axtrm, axtrm, cd1=2,cd2=0,AXON_DIST,cplam);
    make_celseg (xrod, cn, soma, axtrm, cd1=0.2,cd2=0.2,AXON,cplam);
    celnode[ctype][cn] = axtrm;
  }

  else if (ctype==xcone) {	/* procedure to make cone */

    /* simple cone: one compartment: */

    if (notinit(bg_inten)) bginten = 1e4;
    else                   bginten = bg_inten;

    if (notinit(cone_maxcond)) cone_maxcond = 480e-12;
    if (notinit(cone_pigm))    cone_pigm    = 14;

    if (getn(ct,MORPH) == 1) {
      np = loc(nd(xcone,cn,soma), xpos,ypos,getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
      p = make_cone(np,dia=2.5);
      p->xpos=xpos;
      p->ypos=ypos;
	    p->pigm = cone_pigm;
      p->pathl = 7.5;
      p->attf = 0.9;
	    p->maxcond=cone_maxcond;
      p->linit=bginten; 
	    if (!notinit(cone_timec)) p->timec1=cone_timec; //p->timec1=0.2;
	    if (!notinit(cone_loopg)) p->loopg=cone_loopg;
	    p->photnoise=pnoise; p->darknoise=dnoise;
      np = loc(nd(xcone,cn,axtrm),xpos,ypos,getn(xcone,AXARBZ)+getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, axtrm, cd1=0.2,cd2=0.2,AXON,cplam);
      make_celseg (xcone, cn, axtrm, axtrm, cd1=3,cd2=0,AXON_DIST,cplam);
      celnode[ctype][cn] = soma;
    }

    else if (getn(ct,MORPH) == 2) {		/* single compartment cone, cplam = 1 */
      np = loc(nd(xcone,cn,soma), xpos,ypos,getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
      p = (photorec*)make_transducer(np, xpos, ypos);
      np = loc(nd(xcone,cn,axtrm),xpos,ypos,getn(xcone,AXARBZ)+getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, axtrm, cd1=0.2,cd2=0.2,AXON,cplam=1);
      make_celseg (xcone, cn, axtrm, axtrm, cd1=3,cd2=0,AXON_DIST,cplam);
      celnode[ctype][cn] = soma;
    }

    else if (getn(ct,MORPH) == 3) {		/* single compartment cone, cplam = 1 */
      np = loc(nd(xcone,cn,soma), xpos,ypos,getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, soma, cd1=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
      p = (photorec*)make_itransducer(np, xpos, ypos);
      np = loc(nd(xcone,cn,axtrm),xpos,ypos,getn(xcone,AXARBZ)+getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, axtrm, cd1=0.2,cd2=0.2,AXON,cplam=1);
      make_celseg (xcone, cn, axtrm, axtrm, cd1=3,cd2=0,AXON_DIST,cplam);
      celnode[ctype][cn] = soma;
    }

    /* more realistic cone: OS, axon, terminal: 2 compartments */
   
    else if (getn(ct,MORPH) == 4) {
      // Define location of soma.
      np = loc(nd(xcone,cn,soma), xpos,ypos,getn(xcone,SOMAZ));
      // Make soma compartment.
      make_celseg (xcone, cn, soma, soma, dia=getn(ctype,SOMADIA),cd2=0,SOMA,cplam);
      p = make_cone(np,dia=3); p->xpos = xpos; p->ypos=ypos;
      // Link photoreceptor properties to global variables.
	    p->maxcond=cone_maxcond; p->linit=bginten;
	    p->photnoise=pnoise; p->darknoise=dnoise;
      
	    if (!notinit(cone_timec)) p->timec1=cone_timec; //p->timec1=0.2;
	    if (!notinit(cone_loopg)) p->loopg=cone_loopg;
      // Define distance of axon terminal with respect to soma. Can be set in nval file.
      np = loc(nd(xcone,cn,axtrm), xpos,ypos,getn(xcone,AXARBZ)+getn(xcone,SOMAZ));
      // Create axon.
      make_celseg (xcone, cn, soma, axtrm, cd1=0.2,cd2=0.2,AXON,cplam);
      // Create axon terminal.
      make_celseg (xcone, cn, axtrm, axtrm, cd1=3,cd2=0,AXON_DIST,cplam);
      celnode[ctype][cn] = axtrm;
    }
    
    else if (getn(ct,MORPH) == 5) {
      np = loc(nd(xcone,cn,soma), xpos,ypos,getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, soma, cd1=5.13,cd2=0,SOMA,cplam);
      p = make_cone(np); 
      p->xpos=xpos;
      p->ypos=ypos;
	    p->dia=0.505; // --> A = pi * r^2 = 0.2 um² <-> d = 2*r = 0.505 [Nikonov 2006] // joesterle
      p->pathl = 13.4; // https://onlinelibrary.wiley.com/doi/epdf/10.1002/cne.901880204
      p->pigm = cone_pigm;
      p->attf = 1;
	    p->maxcond=cone_maxcond;
      p->linit=bginten; 
      
	    if (!notinit(cone_timec)) p->timec1=cone_timec; //p->timec1=0.2;
	    if (!notinit(cone_loopg)) p->loopg=cone_loopg;
	    p->photnoise=pnoise;
      p->darknoise=dnoise;
      
      np = loc(nd(xcone,cn,axtrm),xpos,ypos,getn(xcone,AXARBZ)+getn(xcone,SOMAZ));
      make_celseg (xcone, cn, soma, axtrm, cd1=1.3,cd2=1.3,AXON,cplam);
      make_celseg (xcone, cn, axtrm, axtrm, cd1=6,cd2=0,AXON_DIST,cplam);
      celnode[ctype][cn] = soma;
    }
  }

  else if (ctype==rbp) {    		/* ctype == rbp ------------ */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(rbp_thetax)) thetax += rbp_thetax; 
    if (!notinit(rbp_thetay)) thetay += rbp_thetay;
    if (!notinit(rbp_thetaz)) thetaz +=  rbp_thetaz; 
    if (!notinit(rbp_flip))     flip +=  rbp_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==dbp1) {        /* ctype == dbp1 ----------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(dbp1_thetax)) thetax += dbp1_thetax;
    if (!notinit(dbp1_thetay)) thetay += dbp1_thetay;
    if (!notinit(dbp1_thetaz)) thetaz += dbp1_thetaz; 
    if (!notinit(dbp1_flip))     flip +=  dbp1_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==dbp2) {        /* ctype == dbp2 ----------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(dbp2_thetax)) thetax += dbp2_thetax;
    if (!notinit(dbp2_thetay)) thetay += dbp2_thetay;
    if (!notinit(dbp2_thetaz)) thetaz += dbp2_thetaz; 
    if (!notinit(dbp2_flip))     flip += dbp2_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==dbp3) {        /* ctype == dbp3 ----------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(dbp3_thetax)) thetax += dbp3_thetax;
    if (!notinit(dbp3_thetay)) thetay += dbp3_thetay;
    if (!notinit(dbp3_thetaz)) thetaz += dbp3_thetaz; 
    if (!notinit(dbp3_flip))     flip += dbp3_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==dbp4) {        /* ctype == dbp4 ----------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(dbp4_thetax)) thetax += dbp4_thetax;
    if (!notinit(dbp4_thetay)) thetay += dbp4_thetay;
    if (!notinit(dbp4_thetaz)) thetaz += dbp4_thetaz; 
    if (!notinit(dbp4_flip))     flip += dbp4_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==hbp1) {        /* ctype == hbp1 ----------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(hbp1_thetax)) thetax += hbp1_thetax;
    if (!notinit(hbp1_thetay)) thetay += hbp1_thetay;
    if (!notinit(hbp1_thetaz)) thetaz += hbp1_thetaz; 
    if (!notinit(hbp1_flip))     flip += hbp1_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==hbp2) {        /* ctype == hbp2 ----------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(hbp2_thetax)) thetax += hbp2_thetax;
    if (!notinit(hbp2_thetay)) thetay += hbp2_thetay; 
    if (!notinit(hbp2_thetaz)) thetaz += hbp2_thetaz; 
    if (!notinit(hbp2_flip))     flip += hbp2_flip; 
    mak_bp(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==ha) {		/* ha --------------------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3)+0.5);
    if (!notinit(ha_thetax)) thetax += ha_thetax;
    if (!notinit(ha_thetay)) thetay += ha_thetay; 
    if (!notinit(ha_thetaz)) thetaz += ha_thetaz; 
    if (!notinit(ha_flip))     flip += ha_flip; 
    mak_hz(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==hb) {		/* hb --------------------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.6,nden*1.3)+0.5);
    if (!notinit(hb_thetax)) thetax += hb_thetax;
    if (!notinit(hb_thetay)) thetay += hb_thetay; 
    if (!notinit(hb_thetaz)) thetaz += hb_thetaz; 
    if (!notinit(hb_flip))     flip += hb_flip; 
    mak_hz(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==hbat) {	/* hbat --------------------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3)+0.5);
    if (!notinit(hbat_thetax)) thetax += hbat_thetax;
    if (!notinit(hbat_thetay)) thetay += hbat_thetay; 
    if (!notinit(hbat_thetaz)) thetaz += hbat_thetaz; 
    if (!notinit(hbat_flip))     flip += hbat_flip; 
    mak_hz(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==aii) { 		 /* aii  -------------------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(aii_thetax)) thetax += aii_thetax;
    if (!notinit(aii_thetay)) thetay += aii_thetay;
    if (!notinit(aii_thetaz)) thetaz += aii_thetaz; 
    if (!notinit(aii_flip))   flip += aii_flip; 
    mak_aii(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==a17) { 		 /* amac  -------------------- */
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(a17_thetax)) thetax += a17_thetax;
    if (!notinit(a17_thetay)) thetay += a17_thetay;
    if (!notinit(a17_thetaz)) thetaz += a17_thetaz;
    if (!notinit(a17_flip))   flip += a17_flip; 
    //    mak_a17(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,int(getn(ctype,MORPH)));
    mak_a17(ctype,cn,nden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==sbac) { 		    /* sbac  -------------------- */
    if (!notinit(sbtheta)) thetaz = sbtheta;
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.8,nden*1.3));
    if (!notinit(sbac_thetax)) thetax += sbac_thetax; 
    if (!notinit(sbac_thetay)) thetay += sbac_thetay;
    if (!notinit(sbac_thetaz)) thetaz += sbac_thetaz; 
    if (!notinit(sbac_flip))     flip += sbac_flip; 
    mak_sbac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==am) { 		   /* am  -------------------- */
    nden = int(getn(ctype,NDENDR));
    if (nden > 2) aden = int(rrange(nden*.8,nden*1.3));
    else aden = nden;
    if (!notinit(am_thetax)) thetax += am_thetax;
    if (!notinit(am_thetay)) thetay += am_thetay;
    if (!notinit(am_thetaz)) thetaz += am_thetaz; 
    if (!notinit(am_flip))     flip += am_flip; 
    mak_amac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==am2) { 		   /* am2  ------------------- */
    nden = int(getn(ctype,NDENDR));
    if (nden > 2) aden = int(rrange(nden*.8,nden*1.3));
    else aden = nden;
    if (!notinit(am2_thetax)) thetax += am2_thetax;
    if (!notinit(am2_thetay)) thetay += am2_thetay;
    if (!notinit(am2_thetaz)) thetaz += am2_thetaz; 
    if (!notinit(am2_flip))     flip += am2_flip; 
    mak_amac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==am3) { 		   /* am3  ------------------- */
    nden = int(getn(ctype,NDENDR));
    if (nden > 2) aden = int(rrange(nden*.8,nden*1.3));
    else aden = nden;
    if (!notinit(am3_thetax)) thetax += am3_thetax;
    if (!notinit(am3_thetay)) thetay += am3_thetay;
    if (!notinit(am3_thetaz)) thetaz += am3_thetaz; 
    if (!notinit(am3_flip))     flip += am3_flip; 
    mak_amac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==am4) { 		   /* am4  ------------------- */
    nden = int(getn(ctype,NDENDR));
    if (nden > 2) aden = int(rrange(nden*.8,nden*1.3));
    else aden = nden;
    if (!notinit(am4_thetax)) thetax += am4_thetax;
    if (!notinit(am4_thetay)) thetay += am4_thetay;
    if (!notinit(am4_thetaz)) thetaz += am4_thetaz; 
    if (!notinit(am4_flip))     flip += am4_flip; 
    mak_amac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==amh) { 		   /* amh  -------------------- */
    nden = int(getn(ctype,NDENDR));
    if (nden > 2) aden = int(rrange(nden*.8,nden*1.3));
    else aden = nden;
    if (!notinit(amh_thetax)) thetax += amh_thetax;
    if (!notinit(amh_thetay)) thetay += amh_thetay;
    if (!notinit(amh_thetaz)) thetaz += amh_thetaz; 
    if (!notinit(amh_flip))     flip += amh_flip; 
    mak_amac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==amh2) { 		   /* amh2  ------------------- */
    nden = int(getn(ctype,NDENDR));
    if (nden > 2) aden = int(rrange(nden*.8,nden*1.3));
    else aden = nden;
    if (!notinit(amh2_thetax)) thetax += amh2_thetax;
    if (!notinit(amh2_thetay)) thetay += amh2_thetay;
    if (!notinit(amh2_thetaz)) thetaz += amh2_thetaz; 
    if (!notinit(amh2_flip))     flip += amh2_flip; 
    mak_amac(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }

  else if (ctype==ams) {        	   /* ctype == ams, small amac ---- */
    					   /* small-field amacrine, like cbp */
    if (getn(ct,MORPH) == MORPH_A1) {
      np = loc(nd(ct,cn,soma), xpos,ypos,getn(ctype,SOMAZ));
      make_celseg (ct, cn, soma, soma, cd1=getn(ct,SOMADIA),cd2=0,SOMA,cplam);
      np = loc(nd(ct,cn,axtrm), xpos,ypos,getn(ctype,AXARBZ));/* axon tip */
      make_celseg (ct, cn, soma, axtrm, cd1=0.8 ,cd2=0.8, AXON,cplam); /* connect to soma */
      make_celseg (ct, cn, axtrm, axtrm, cd1=1,cd2=0, AXON,cplam);     /* sphere at tip */
      celnode[ctype][cn] = axtrm;
    }
    else if (getn(ct,MORPH) == MORPH_A2) {	/* small-field morph, only soma */
      np = loc(nd(ct,cn,soma), xpos,ypos,getn(ctype,SOMAZ));
      make_celseg (ct, cn, soma, soma, cd1=getn(ct,SOMADIA),cd2=0,SOMA,cplam);
      celnode[ctype][cn] = soma;
    }
  }
  else if (ctype==amhs) {        	   /* ctype == amhs, small hyp amac */
    					   /* small-field amacrine, like cbp */
    if (getn(ct,MORPH) == MORPH_A1) {
      np = loc(nd(ct,cn,soma), xpos,ypos,getn(ctype,SOMAZ));
      make_celseg (ct, cn, soma, soma, cd1=getn(ct,SOMADIA),cd2=0,SOMA,cplam);
      np = loc(nd(ct,cn,axtrm), xpos,ypos,getn(ctype,AXARBZ));/* axon tip */
      make_celseg (ct, cn, soma, axtrm, cd1=0.8 ,cd2=0.8, AXON,cplam); /* connect to soma */
      make_celseg (ct, cn, axtrm, axtrm, cd1=1,cd2=0, AXON,cplam);     /* sphere at tip */
      celnode[ctype][cn] = axtrm;
    }
    else if (getn(ct,MORPH) == MORPH_A2) {	/* small-field morph, only soma */
      np = loc(nd(ct,cn,soma), xpos,ypos,getn(ctype,SOMAZ));
      make_celseg (ct, cn, soma, soma, cd1=getn(ct,SOMADIA),cd2=0,SOMA,cplam);
      celnode[ctype][cn] = soma;
    }
  }

  else if (ctype >=gca && ctype <= gcboff) {    /* gca  --------------------- */
       if (ctype==dsgc) {
          if (!notinit(dsgc_thetax)) thetax += dsgc_thetax;
          if (!notinit(dsgc_thetay)) thetay += dsgc_thetay;
          if (!notinit(dsgc_thetaz)) thetaz += dsgc_thetaz;
          if (!notinit(dsgc_flip))     flip += dsgc_flip; 
       } else {
          if (!notinit(gctheta)) thetaz = gctheta;
          if (!notinit(gc_thetax)) thetax += gc_thetax;
          if (!notinit(gc_thetay)) thetay += gc_thetay;
          if (!notinit(gc_thetaz)) thetaz += gc_thetaz;
          if (!notinit(gc_flip))     flip += gc_flip; 
       }
    nden = int(getn(ctype,NDENDR));
    aden = int(rrange(nden*.85,nden*1.2));
    mak_gc(ctype,cn,aden,xpos,ypos,thetax,thetay,thetaz,flip,int(getn(ctype,MORPH)));
  }
  else
    fprintf(stderr, "makcel.n: proc makcell: cell type %d not available.", ctype);

  ccount = cellnums[ct][NCELLS]+1;		// allows multiple calls to setupcells()
  cellnums[ct][ccount] = cn;
  cellnums[ct][NCELLS] = ccount;
}

/*-------------------------------------------------------------*/

int trace_node (int ct, int cn, int sn, int tn)

/* Given starting node "sn", trace back towards soma and stop 
   if trace node "tn" is found. 
   Based on the first element at a node being the parent.
*/ 

{
   int nn, nelem, onn;
   int nod1, nod2;
   elem *epnt;

   for (onn=nn=sn; nn!=tn; sn ) {  		/* trace path back to soma */
     epnt = get_elempnt(int(get_nfield(ndn(ct,cn,nn),ELNUM,1)));
     nod1=  int(get_efield(epnt,NODE1C));     /* first node */
     nod2=  int(get_efield(epnt,NODE2C));     /* second node */
     //fprintf (stderr,"nod1 %d nod2 %d\n",nod1,nod2);
     if (nod2 == nn)
          nn = nod1;
     else nn = nod2;
     //printf ("tracing from node %d to %d\n",onn,nn);
     onn = nn;
     if (nn==soma) break;
   }  /* for (nn=sn;;) */
   if (nn==tn) return 1;
   else        return 0;
}

/*-------------------------------------------------------------*/

int trace_node_alt (int ct, int cn, int sn, int tn)

/* Based on the first node on a cable connecting to the parent. */

{
   int i,n,nelem,nconn;
   int nod1, nod2;
   elem *epnt;

 for (n=sn; n!=tn; n) {
   nconn = int(get_nfield(ndn(ct,cn,n),NUMCONN));
   for (i=1; i<=nconn; i++) {  /* Find path back to soma */
     epnt = get_elempnt(int(get_nfield(ndn(ct,cn,n),ELNUM,i)));

     if (epnt->ctype != CABLE) continue;

     if (ninfo>2)
      printf ("tracing from node %d conn %d\n",n,i);

     nod1= int(get_efield(epnt,NODE1C));     /* first node */
     nod2= int(get_efield(epnt,NODE2C));     /* second node */
     fprintf (stderr,"nod1 %d nod2 %d n %d\n",nod1,nod2,n);
     if (nod2 == n)
       continue;             /* descendent node */
     if (nod1 == n) {
       n = nod2;             /* follow parent node */
       break;
     }
   }  /* for (i;;) */
   if (n==soma) break;
 }
 if (n==tn) return 1;
 else       return 0;
}

/*-------------------------------------------------------------*/

void mak_sharp_electrode(int ct, int cn)

/* Make circuit analog to sharp electrode. */

/*  With sharp electrode recordings, spike height is attenuated,
    spikes have little or no after-hyperpolarization, and spikes
    are observed during rising edge of PSP.  

    With patch recordings, spikes always have after-hyperpolarization,
    and each spike leaves membrane voltage at the same potential
    so spiking prevents the membrane voltage from rising.

    The sharp electrode analog reduces spike height and after-hyperpolarization
    but does not affect the low frequency response.
*/

{
   int recpnt2, recpnt3;
   double vstart;
   resistor *rpnt;
   capac *cpnt;
   diode *dpnt;

#define RECPNT 10000

  recpnt2 = recpnt+2;
  recpnt3 = recpnt+3;
 
  if (notinit(elec_resist1)) elec_resist1 = 100e6;
  if (notinit(elec_resist2)) elec_resist2 = 100e6;
  if (notinit(elec_capac1))  elec_capac1 =  1e-12;
  if (notinit(elec_capac2))  elec_capac2 =  100e-12;
  vstart = getn(ct,VREST);
  //conn (nd(ct,soma),nd(ct,recpnt3),VBUF); 
  //cpnt = (capac*)at(nd(ct,recpnt3),GNDCAP); cpnt->c = elec_capac1, cpnt->vrest=vstart;
  rpnt = (resistor*)conn (nd(ct,soma),nd(ct,recpnt),RESISTOR); rpnt->r = elec_resist1;
  cpnt = (capac*)at(nd(ct,recpnt),GNDCAP); cpnt->c = elec_capac1, cpnt->vrest=vstart;
  //dpnt = conn (nd(gc,recpnt3),nd(gc,recpnt),DIODE); rpnt->r = elec_resist1;

  rpnt = (resistor*)conn(nd(ct,recpnt),nd(ct,recpnt2),RESISTOR); rpnt->r = elec_resist2;
  cpnt = (capac*)at(nd(ct,recpnt2),GNDCAP); cpnt->c = elec_capac2, cpnt->vrest=vstart;
}

/*-------------------------------------------------------------*/
