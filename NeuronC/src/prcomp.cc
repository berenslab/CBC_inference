/* Module prcomp in Program nc */

/* displays compartment list */

/*        Feb 91                  R.G. Smith */

extern "C" {

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

}

#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "ncelem.h"
#include "control.h"
#include "ncio.h"
#include <vector>

#define DEBUG 1
#ifdef DEBUG
#include "ncdebug.h"
#endif


extern comp *compnt;
extern conn *connpnt;

extern double xrot;
extern double yrot;
extern double zrot;
extern double rxcent,rycent,rzcent;
extern double dxcent,dycent,dzcent;
extern double dsize;
extern int cumchan;

extern FILE *compout;
extern int prcomps;

int testlin(int nprint);
void setrot(double xrot, double yrot, double zrot, 
		double xcent, double ycent, double zcent, 
		double rxcent, double rycent, double rzcent, double scal);
void drcomp(double xloc, double yloc, double zloc, double dia, int color, 
			double vmax, double vmin);
void drcconn (double x1, double y1, double z1, 
		double x2, double y2, double z2);
conn *printconns(comp *pnt, int nprint);
char *findsym(int);
comp *othercomp (comp *pnt,conn *cpnt);
char *prnode (int n1, int n2, int n3, int n4); 
ntcomp *findntcomp (comp *cpnt, int type);

/*------------------------------------*/

void prcomp(void)

/* print a list of all compartments, with node number,
    capacitances, and resistances. Label photrecs "rec" and
    the node number of the photorec.
*/
  
{
   comp *pnt;
   conn *cpnt;
   conlst *lpnt;
   node *npnt;
   int n,nprint,missing;
   double dia,rm;


if (debug & NCPRCOMP) ncfprintf (stderr,"prcomp start prmap %d disp %d\n",prmap,disp);

 if (cumchan) {
  ncfprintf (stdout,"#c Temperature %g deg Celsius\n",tempcel);
ncfprintf (stdout,"#c T Q10 coefficients dqc %g dqm %g dqh %g dqn %g dqkca %g dqca %g\n",
		dqc,dqm,dqh,dqn,dqkca,dqca);
 }

if (prmap&PCOND) ncfprintf
 (stdout,
   "# comp node     nodrm       nodcap        connections (siemens)\n");
else if (prmap&PCOMP) ncfprintf
 (stdout,
   "# comp node     equiv comp dia, rm        connections (ohms)\n");

if (prmap&(PCOND|PCOMP)) 
 for (pnt=compnt; pnt; pnt=pnt->next) {
   ncfprintf(stdout,"#c %-3d ",pnt->num);

						/* print the nodes */
  for (n=0,lpnt=pnt->nodlst; lpnt; lpnt=lpnt->next,n++) {
        if (n) ncfprintf (stdout,
		"\n#c %-3d ",pnt->num);    /* second node and up */
        if (lpnt->conpnt) npnt = (node *)lpnt->conpnt;
        if (npnt->nodenm4 > NULLVAL)
		 ncfprintf (stdout,"[%-d %-d %-d %-d]",
		   npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4);
        else if (npnt->nodenm3 > NULLVAL)
		 ncfprintf (stdout,"[%-d %-d %-d]",
			npnt->nodenm1,npnt->nodenm2,npnt->nodenm3);
        else if (npnt->nodenm2 > NULLVAL)
		 ncfprintf (stdout,"[%-d %-d]",npnt->nodenm1,npnt->nodenm2);
	else
                 ncfprintf (stdout,"[%-4d]   ",npnt->nodenm1);
  }

  if (!n) ncfprintf (stdout,"        ");
                                                /*   print capac, resist */
  if (dcm==0) dcm = 1e-6;			/* assume 1 uf / cm2 */
  dia = 2. * sqrt(pnt->cap/(dcm*4*MPI)) * 1e4; /* find sphere dia from cap */
  if ((rm=pnt->rm)==0) rm = 1e-20;
  rm  = MPI * dia * dia / (rm * 1e8);
  missing = (pnt->rm > (LARGERES * .9) || pnt->cap < SMALLCAP * 1.1);

   if (missing) {
     if (prmap&PCOMP) 
       ncfprintf (stdout,
		"MISSING dia %-5.3g Rm %-6g ",dia,rm); /* sphere printout*/
     else
      ncfprintf(stdout,
	" MISSING gm %-8.3g cap %-8.3g  ",pnt->rm,pnt->cap);/* standard printout */
   }
    else {
     if (pnt->jnoise) {
      if (prmap&PCOMP) 
        ncfprintf (stdout,
	  " sphere dia %-.3g, Rm %-g J",dia,rm); /* sphere printout*/
      else
       ncfprintf(stdout,
	  " Gm %-8.3g J,     cap %-8.3g",pnt->rm,pnt->cap); /* std. print*/
      }
     else {
      if (prmap&PCOMP) 
       ncfprintf (stdout,
		" sphere dia %-.3g, Rm %-g",dia,rm); /* sphere printout*/
      else
       ncfprintf(stdout,
	  " Gm %-8.3g,     cap %-8.3g",pnt->rm,pnt->cap); /* std. printout*/
      }
     }
   nprint = 3;                                  /* things printed on line */
   cpnt = printconns(pnt,nprint);
   ncfprintf (stdout,"\n");
   }
  fflush (stdout);
  if (debug & NCPRCOMP) ncfprintf (stderr,"prcomp end\n");
}

/*------------------------------------*/

conn *printconns(comp *pnt, int nprint)

/* Print the connections to a compartment.  */

{
   const char *sc;
   conlst *lpnt;
   comp *pnt2;
   chan *cpt;
   conn *cpnt;
   ntcomp *npt;
   double cond,gjcond,ri,area;
   int nchan,n;

   if (dcm==0) dcm = 1e-6;			/* assume 1 uf / cm2 */
   area = pnt->cap/dcm;				/* find sphere area from cap */
   if (area == 0) area = 1e-20;
   for (n=0,lpnt=pnt->clst; lpnt; lpnt=lpnt->next) {  /* find connections */
      n++;
      if (!(cpnt=lpnt->conpnt)) continue;
      pnt2 = cpnt->comp1;

      switch (cpnt->ctype) {
        case SYN2:
        case AMPA:
        case KAINATE:
        case NMDA:
        case GABA:
        case GLY:
        case CGMP:
        case NA:
        case K:
        case ClCa:
        case KCa:
        case CA:
          if (prmap&PCOMP)
	      cond = ((chan *)cpnt)->maxcond / area;
          else
	      cond = ((chan *)cpnt)->maxcond;
	  break;
	 default:
	   break;
      }
      switch (cpnt->ctype) {
        case ROD:
                nprint = testlin(nprint);
                if (((photrec *)cpnt)->recnm3 != NULLVAL)
		  ncfprintf (stdout,
			"rod %-4d %-2d %-2d",((photrec *)cpnt)->recnm1,
				((photrec *)cpnt)->recnm2,
				((photrec *)cpnt)->recnm3);
                else if (((photrec *)cpnt)->recnm2 != NULLVAL)
		  ncfprintf (stdout,"rod %-4d %-4d",((photrec *)cpnt)->recnm1,
						 ((photrec *)cpnt)->recnm2);
		else
                  ncfprintf (stdout,"rod %-4d",((photrec *)cpnt)->recnm1);
                break;
        case CONE:
                nprint = testlin(nprint);
                if (((photrec *)cpnt)->recnm3 != NULLVAL)
		  ncfprintf (stdout,
			"cone %-4d %-2d %-2d",((photrec *)cpnt)->recnm1,
				 ((photrec *)cpnt)->recnm2,
				 ((photrec *)cpnt)->recnm3);
                else if (((photrec *)cpnt)->recnm2 != NULLVAL)
		  ncfprintf (stdout,
			"cone %-4d %-4d",((photrec *)cpnt)->recnm1,
						 ((photrec *)cpnt)->recnm2);
		else
                  ncfprintf (stdout,
			"cone %-4d",((photrec *)cpnt)->recnm1);
                break;
        case CHR:
                nprint = testlin(nprint);
                if (((photrec *)cpnt)->recnm3 != NULLVAL)
		  ncfprintf (stdout,
			"chr %-4d %-2d %-2d",((photrec *)cpnt)->recnm1,
				 ((photrec *)cpnt)->recnm2,
				 ((photrec *)cpnt)->recnm3);
                else if (((photrec *)cpnt)->recnm2 != NULLVAL)
		  ncfprintf (stdout,
			"chr %-4d %-4d",((photrec *)cpnt)->recnm1,
						 ((photrec *)cpnt)->recnm2);
		else
                  ncfprintf (stdout,
			"chr %-4d",((photrec *)cpnt)->recnm1);
                break;
        case VTRANSDUCER:
                nprint = testlin(nprint);
                if (((photrec *)cpnt)->recnm3 != NULLVAL)
		  ncfprintf (stdout,
			"tranduc %-4d %-2d %-2d",((photrec *)cpnt)->recnm1,
				 ((photrec *)cpnt)->recnm2,
				 ((photrec *)cpnt)->recnm3);
                else if (((photrec *)cpnt)->recnm2 != NULLVAL)
		  ncfprintf (stdout,
			"tranduc %-4d %-4d",((photrec *)cpnt)->recnm1,
						 ((photrec *)cpnt)->recnm2);
		else
                  ncfprintf (stdout,
			"tranduc %-4d",((photrec *)cpnt)->recnm1);
                break;
        case ITRANSDUCER:
                nprint = testlin(nprint);
                if (((photrec *)cpnt)->recnm3 != NULLVAL)
		  ncfprintf (stdout,
			"itranduc %-4d %-2d %-2d",((photrec *)cpnt)->recnm1,
				 ((photrec *)cpnt)->recnm2,
				 ((photrec *)cpnt)->recnm3);
                else if (((photrec *)cpnt)->recnm2 != NULLVAL)
		  ncfprintf (stdout,
			"itranduc %-4d %-4d",((photrec *)cpnt)->recnm1,
						 ((photrec *)cpnt)->recnm2);
		else
                  ncfprintf (stdout,
			"itranduc %-4d",((photrec *)cpnt)->recnm1);
                break;
        case LOAD:
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"load %-8.3g",((load *)cpnt)->conduct);
                break;
	case CGMP:
	case SYN2:
	case AMPA:
	case KAINATE:
	case NMDA:
	case GABA:
	case GLY:
        case NA:
        case CA:
        case K:
        case KCa:
        case ClCa:
                nprint = testlin(nprint);
		if (nchan=int(((chan *)cpnt)->nchan))
                     ncfprintf (stdout,
			"%s(%d) %-.3g n=%d",
		  findsym(cpnt->ctype),((chan *)(long)cpnt->stype),cond,nchan);
                else ncfprintf (stdout,
			"%s(%d) %-.3g",
			findsym(cpnt->ctype),((chan *)(long)cpnt->stype),cond);
                if (((chan *)cpnt)->compe!=NULL) 
			ncfprintf (stdout,"pvext %-d",((chan *)cpnt)->compe->num);
                if (((chan *)cpnt)->compp!=NULL) 
			ncfprintf (stdout,"pvext %-d",((chan *)cpnt)->compp->num);
                break;
        case SYNAPSE:
      	      pnt2 = ((synap *)cpnt)->comp1;
              if (pnt2 == pnt) {		/* presynaptic comp */
                pnt2 = ((synap *)cpnt)->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"syn(%.3g) to -> %-4d",((synap *)cpnt)->maxcond,pnt2->num);
              }
              else {				/* postsynaptic comp */
                nprint = testlin(nprint);
		if (cpt=((synap*)cpnt)->resp1) {
                   ncfprintf (stdout,"syn from <- %-2d",pnt2->num);
		}
		else if (cpt=((synap*)cpnt)->resp2) {
                   ncfprintf (stdout,"syn from <- %-2d",pnt2->num);
		}
		else {
		   ncfprintf (stdout,
			"syn(%.3g) from <- %-2d",((synap*)cpnt)->maxcond,pnt2->num);
		}
		if (cpt=((synap*)cpnt)->resp1) {
                  ncfprintf (stdout,
			"->%s", findsym(cpt->ctype));
		}
		else if (cpt=((synap*)cpnt)->resp2) {
                  ncfprintf (stdout,
			"->%s", findsym(cpt->ctype));
		}
		else if (((synap*)cpnt)->mesg1) {
                    ncfprintf (stdout,
			"->%s",findsym(((synap*)cpnt)->mesg1));
  		}
		else if (((synap*)cpnt)->mesg2) {
                    ncfprintf (stdout,
			"->%s",findsym(((synap*)cpnt)->mesg2));
  		}
              }
                break;
        case GJ:
                if (pnt2 == pnt) pnt2 = cpnt->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
		if (prmap & PCOMP) {
	          if (cpnt->conduct==0) gjcond = 1e-30;
		  else gjcond = ((gj *)cpnt)->maxcond;
		  ri = 1/gjcond;		
                  ncfprintf (stdout,
			"%8.3g gj %-d",ri,pnt2->num);
	        }
		else {
                  ncfprintf (stdout,
			"gj (%.3g) <-> %-d",((gj *)cpnt)->maxcond,pnt2->num);
		}
                break;
        case BUF:
              if (pnt2 == pnt) {
                pnt2 = ((dbuf *)cpnt)->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"vbuf to   -> %-4d",pnt2->num);
              }
              else {
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"vbuf from -> %-4d",pnt2->num);
              }
              break;
        case NBUF:
              if (pnt2 == pnt) {
                pnt2 = ((ndbuf *)cpnt)->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"nbuf to   -> %-4d",pnt2->num);
              }
              else {
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"nbuf from -> %-4d",pnt2->num);
              }
              break;
        case DIODE:
              if (pnt2 == pnt) {
                pnt2 = ((dbuf *)cpnt)->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"diode to   -> %-4d",pnt2->num);
              }
              else {
                nprint = testlin(nprint);
                ncfprintf (stdout,
			"diode from -> %-4d",pnt2->num);
              }
              break;
        case AXIALRES:
                if (pnt2 == pnt) pnt2 = cpnt->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
		if (prmap&PCOMP) { 
		  if (cpnt->conduct==0) gjcond = 1e-30;
		  else gjcond = cpnt->conduct;
		  ri = 1/gjcond;		
                  ncfprintf (stdout,
			"res ax %8.3g -> %-d",ri,pnt2->num);
		}
		else
                  ncfprintf (stdout,
			"res ax %8.3g -> %-d",cpnt->conduct,pnt2->num);
                break;
        default:
                if (pnt2 == pnt) pnt2 = cpnt->comp2;  
                if (!pnt2) break;
                nprint = testlin(nprint);
		if (prmap&PCOMP) { 
		  if (cpnt->conduct==0) gjcond = 1e-30;
		  else gjcond = cpnt->conduct;
		  ri = 1/gjcond;		
                  ncfprintf (stdout,
			"%8.3g -> %-d",ri,pnt2->num);
		}
		else
                  ncfprintf (stdout,
			"%8.3g -> %-d",cpnt->conduct,pnt2->num);
                break;
      }
   }
  if (pnt->capnt) {
		 int nshell;

                nprint = testlin(nprint);
		nshell = pnt->capnt->cashell - 1;
                if (nshell > 1) sc="s";
                else            sc="";
		ncfprintf (stdout,
			"Cai comp %-d shell%s dr %g", nshell,sc,pnt->capnt->dr);
                nprint = testlin(nprint);
		nshell = pnt->capnt->caoshell - 1;
                if (nshell > 1) sc="s";
                else            sc="";
		ncfprintf (stdout,
			"Cao comp %-d shell%s", nshell,sc);
		if (pnt->capnt->vmax > 0) {
                  nprint = testlin(nprint);
		  ncfprintf (stdout,"Ca pump %-g", pnt->capnt->vmax);
		}
		if (pnt->capnt->kexi > 0) {
                  nprint = testlin(nprint);
		  ncfprintf (stdout,"Ca exch %-g", pnt->capnt->kexi);
		}
		if (pnt->capnt->cabf > 0) {
                  nprint = testlin(nprint);
		  ncfprintf (stdout,"Ca buf %-g %-g", pnt->capnt->cabf, pnt->capnt->cabr);
		}
  }
  if (npt=findntcomp(pnt,CAMP)) {
		  int n;

                nprint = testlin(nprint);
		
		if ((n=npt->n) > 1)
                  ncfprintf (stdout,
			"cAMP comp %-d inputs",n);
		else
                  ncfprintf (stdout,
			"cAMP comp",n);
  }
  if (npt=findntcomp(pnt,CGMP)) {
		  int n;

                nprint = testlin(nprint);
		if ((n=npt->n) > 1)
                  ncfprintf (stdout,
			"cGMP comp %-d inputs",n);
		else
                  ncfprintf (stdout,
			"cGMP comp",n);
  }

  return (cpnt);	/* return pointer to the last connection */
}

/*------------------------------------*/

int nodematch(node *npnt, int n1a, int n1b, int n1c, int n1d, int except, 
			  int na, int nb, int nc, int nd)
{
   int match,nomatch;

 match=1;
 if (n1a != NULLVAL) 			/* check node 1 first */	
    if (npnt->nodenm1!=n1a) match=0;
 if (n1b != NULLVAL)
    if (npnt->nodenm2!=n1b) match=0;
 if (n1c != NULLVAL)
    if (npnt->nodenm3!=n1c) match=0;
 if (n1d != NULLVAL)
    if (npnt->nodenm4!=n1d) match=0;

 if (match) {
   nomatch = 0;
   if (na!=NULLVAL || nb!=NULLVAL || nc!=NULLVAL) {
     nomatch = 1;
     if (na!=NULLVAL && npnt->nodenm1!=na) nomatch=0;     
     if (nb!=NULLVAL && npnt->nodenm2!=nb) nomatch=0;     
     if (nc!=NULLVAL && npnt->nodenm3!=nc) nomatch=0;     
     if (nd!=NULLVAL && npnt->nodenm4!=nd) nomatch=0;     
   }
   if (nomatch) match = 0;
 }
 return match;
}

/*------------------------------------*/

void dcomp(int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype, 
	int na, int nb, int nc, int nd, 
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int,int),
	double vmax, double vmin, 
	double dscale, int hide, int excl)
{
   comp *pnt,*pnt1,*pnt2;
   conn *cpnt,*cpnt1,*cpnt2,*cnpnt;
   conlst *lpnt,*nlpnt,*tlpnt;
   node *npnt;
   int j,n,n1,n2;
   int numcomp,num1,num2;
   double dia;
   double avgx,avgy,avgz;
   double xloc,yloc,zloc;
   double x1,y1,z1;
   double x2,y2,z2;
   double oldx,oldy,oldz;
   int match;

#ifdef DEBUG
  if ((debug & NCPRCOMP) && (debugz & 2))
	 ncfprintf (stderr,
		"dcomp n1a %d n1b %d n1c %d\n",n1a,n1b,n1c);
#endif

if (disp & DCOMP) {			/* display compartments */

/* to find the location of a compartment, first see if it 
    represents one or more nodes.  If so, find node locations
    and average their location for comp's location.  If comp
    is not at a node, then it must be in a cable.  On the first
    such compartment, find both ends of the cable (comps that
    are at nodes) and calculate the positions of comp by its
    relative location in the cable.  Since compartments inside
    cables are always created sequentially, all the remaining
    compartments can be displayed immediately (in a local "for"
    loop) after calculating their positions.  Note that the
    beginning and end compartments (the ones at nodes) may have
    been created previously, so they are not necessarily sequential
    with the compartments inside the cable.  They must be displayed
    separately in the main loop. 
*/

 if (dcm==0) dcm = 1e-6;			/* assume 1 uf / cm2 */
 //setrot(xrot,yrot,zrot,dxcent,dycent,dzcent,rxcent,rycent,rzcent,dsize); 

  std::vector<int> comp_nodes_ct;
  std::vector<int> comp_nodes_cn;
  std::vector<int> comp_nodes_nd;
 
 for (pnt=compnt; pnt; pnt=pnt->next) {	/* for all comps */

    comp_nodes_ct.clear();
    comp_nodes_cn.clear();
    comp_nodes_nd.clear();
 
    /* ncfprintf (stderr,"comp # %d\n",pnt->num);  /* */

    avgx=avgy=avgz=0.0;
	/* find nodes that point to comp */
    for (match=n=0,lpnt=pnt->nodlst; lpnt; lpnt=lpnt->next,n++) {
      if (lpnt->conpnt) {
        npnt = (node *)(lpnt->conpnt);
	      match |= nodematch(npnt,n1a,n1b,n1c,n1d,exceptype,na,nb,nc,nd);
        avgx += npnt->xloc;			/* find average location */
        avgy += npnt->yloc;
        avgz += npnt->zloc;
		
		    // Save one node to compartment for external voltage.
		    comp_nodes_ct.push_back(npnt->nodenm1);
		    comp_nodes_cn.push_back(npnt->nodenm2);
		    comp_nodes_nd.push_back(npnt->nodenm3);
		    //node_to_comp[3] = npnt->nodenm4;
      } 
    }
	/* if there are one or more nodes, draw the comp */
	/*   at the nodes' average location: */

    if (n) {
      if (match) {
        xloc = avgx/n;
        yloc = avgy/n;
        zloc = avgz/n;
        dia = 2. * sqrt(pnt->cap/(dcm*4*MPI)) * 1e4; /* sphere dia from cap */
        drcomp(xloc,yloc,zloc,dia*dscale,color,vmax,vmin);     /* */
          /* ncfprintf (stderr,"comp %d nodes %d xloc %-5.3g yloc %-5.3g dia %-5.3g\n",
                    pnt->num,n,xloc,yloc,dia);   /* */
	    if (prcomps) {
        ncfprintf (compout,"comp %d nodes %d nodes_to_comp [", pnt->num, n);
        for (int nd_idx=0; nd_idx<(int)comp_nodes_ct.size(); nd_idx++) {
          ncfprintf (compout,"(%d,%d,%d);", comp_nodes_ct[nd_idx],comp_nodes_cn[nd_idx],comp_nodes_nd[nd_idx]);
        }
        ncfprintf (compout,"] xloc %-5.3g yloc %-5.3g zloc %-5.3g dia %-5.3g Cm %-5.3g Rm %-5.3g\n",xloc,yloc,zloc,dia,pnt->cap,pnt->rm);
          
        //ncfprintf (compout,"comp %d nodes %d node_to_comp %i %i %i xloc %-5.3g yloc %-5.3g zloc %-5.3g dia %-5.3g Cm %-5.3g Rm %-5.3g\n",
				//	pnt->num,n,node_to_comp[0],node_to_comp[1],node_to_comp[2],xloc,yloc,zloc,dia,pnt->cap,pnt->rm); /* new to write morphology to file */
	      fflush(compout);
	    }
      }
    } else {

      /* otherwise, comp must be inside cable. */
	  /*   traverse list of comp -> conn -> comp to find */
	  /*   end nodes of cable and then compute location  */
	  /*   of comp from end nodes' locs                  */
	  /*   Ignore conns that are not of type AXIALRES    */
      
      
      
      match = 0;
      avgx=avgy=avgz=0.0;
      num1=num2=0;
      for (j=0,lpnt=pnt->clst; lpnt; lpnt=lpnt->next,j++)  {/* find first node */
      
        if (!(cpnt1=lpnt->conpnt)) break;
        if (cpnt1->ctype!=AXIALRES) {j--;continue;}
        if ((pnt1=cpnt1->comp1)==pnt) pnt1 = cpnt1->comp2;  /* start it off */
         
        for (n=0; pnt1->nodlst==NULL; n++) {		/* traverse cable */
      
          tlpnt = pnt1->clst;
          if (!(cpnt2=tlpnt->conpnt)) break;		  /* find conn */
          if (cpnt1==cpnt2) {
	        tlpnt = tlpnt->next;
	        cpnt2 = tlpnt->conpnt; 			/* wrong conn? */
          }
          if (!cpnt2) break;
          while (cpnt2->ctype!=AXIALRES) {
            tlpnt = tlpnt->next;
	        cpnt2 = tlpnt->conpnt;
            if (!cpnt2) break;
          }
      
          if ((pnt2=cpnt2->comp1)==pnt1) pnt2 = cpnt2->comp2;
          pnt1=pnt2;
          cpnt1=cpnt2;
        }
        if (j==0) {		 /* get avg loc of first end compartment */
          num1 = n;
          avgx=avgy=avgz=0.0;
          for (n=0,nlpnt=pnt1->nodlst; nlpnt; nlpnt=nlpnt->next,n++) {
      
            if (nlpnt->conpnt) {
              npnt = (node *)(nlpnt->conpnt);
	          match |= nodematch(npnt,n1a,n1b,n1c,n1d,exceptype,na,nb,nc,nd);
              avgx += npnt->xloc;			/* find average location */
              avgy += npnt->yloc;
              avgz += npnt->zloc;
            } 
          }
          if (n) {
            x1 = avgx / n;			/* average loc of first compartment */
            y1 = avgy / n;
            z1 = avgz / n;
          }
        }  /* if (j==0) */
      
        else if (j==1) {		 /* add loc of second end compartment */
      
          num2 = n;
          avgx=avgy=avgz=0.0;
          for (n=0,nlpnt=pnt1->nodlst; nlpnt; nlpnt=nlpnt->next,n++) {
      
            if (nlpnt->conpnt) {
              if (!(npnt=(node *)(nlpnt->conpnt))) {n--; continue;}
	          match |= nodematch(npnt,n1a,n1b,n1c,n1d,exceptype,na,nb,nc,nd);
              avgx = avgx + npnt->xloc;		/* find average location */
              avgy += npnt->yloc;
              avgz += npnt->zloc;
            } 
          }
          if (n) {
            x2 = y2 = z2 = 0.0;
            if (avgx) x2 = avgx / n;	/* average loc of second end comp */
            if (avgy) y2 = avgy / n;
            if (avgz) z2 = avgz / n;
          }
        }   /* if (j==1) */
      }    /* for (j=0;;) */
      
      if (!match) continue;		/* don't draw comp; check next one */
      
      numcomp = num1 + num2 + 2;		/* number of segments (conns) */
      
	  		/* Now display all the compartments in the cable. */
 	  		/* Interpolate location from end compartment locs. */
      oldx=x1;
      oldy=y1;
      oldz=z1;
      for (j=1; j<numcomp; j++) {
        xloc = (x1 * (numcomp-j) + x2 * j) / numcomp;
        yloc = (y1 * (numcomp-j) + y2 * j) / numcomp;
        zloc = (z1 * (numcomp-j) + z2 * j) / numcomp;
        if (match) {
          dia = 2. * sqrt(pnt->cap/(dcm*4*MPI)) * 1e4; /* sphere dia from cap */
          drcomp(xloc,yloc,zloc,dia*dscale,color,vmax,vmin); /* draw the comp */
          if (disp&DCONN) drcconn(oldx,oldy,oldz,xloc,yloc,zloc);
        }						 /* draw connection */
        oldx=xloc;
        oldy=yloc;
        oldz=zloc;
      
/*     ncfprintf (stderr,"comp %d seg %d xloc %-5.3g yloc %-5.3g dia %-5.3g\n",
	  			pnt->num,j,xloc,yloc,dia); /* */
        if (j>1) pnt = pnt->next;
        if (prcomps) {
          ncfprintf (compout,"comp %d nodes %d nodes_to_comp [", pnt->num, n);
          for (int nd_idx=0; nd_idx<(int)comp_nodes_ct.size(); nd_idx++) {
            ncfprintf (compout,"(%d,%d,%d);", comp_nodes_ct[nd_idx],comp_nodes_cn[nd_idx],comp_nodes_nd[nd_idx]);
          }
          ncfprintf (compout,"] xloc %-5.3g yloc %-5.3g zloc %-5.3g dia %-5.3g Cm %-5.3g Rm %-5.3g\n",xloc,yloc,zloc,dia,pnt->cap,pnt->rm);
          fflush(compout);
        }
      } 
      if (match && disp&DCONN) drcconn(oldx,oldy,oldz,x2,y2,z2);
			        /* draw the last conn in cable */

    } 	/* else if (!n) */
 }     /* for (pnt=compnt;;)  */
}    /* if (disp&DCOMP) */

				/* Now display all the connections */
				/*  associated with nodes. */
				/* Problem is that we don't want to draw */
				/* them more than once, so we must use */
				/* a loop based on the connection list, */
				/* not the compartment list. */
if (disp & DCONN) {
  for (cnpnt=connpnt; cnpnt; cnpnt=cnpnt->next) {

    if (!(pnt1=cnpnt->comp1)) continue;
    if (!(pnt2=cnpnt->comp2)) continue;

/* ncfprintf (stderr,"conn c1 %d c2 %d nod %d\n",
			pnt1->num,pnt2->num,pnt1->nodlst);   /* */

    if (prcomps) {
      ncfprintf (compout,"conn c1 %d c2 %d nod %d Ri %-5.3g\n",
			pnt1->num,pnt2->num,pnt1->nodlst,1/cnpnt->conduct);   /* new to write morphology to file */      
      fflush(compout);
    }
     
    avgx=avgy=avgz=0.0;
    for (match=n1=0,nlpnt=pnt1->nodlst; nlpnt; nlpnt=nlpnt->next,n1++) {

      if (nlpnt->conpnt) {
        npnt = (node *)(nlpnt->conpnt);
	match |= nodematch(npnt,n1a,n1b,n1c,n1d,exceptype,na,nb,nc,nd);
        avgx += npnt->xloc;			/* find average location */
        avgy += npnt->yloc;
        avgz += npnt->zloc;
      } 
    }
    if (n1) {
      x1 = avgx / n1;			/* average loc of first end comp */
      y1 = avgy / n1;
      z1 = avgz / n1;
    }

    avgx=avgy=avgz=0.0;
    for (n2=0,nlpnt=pnt2->nodlst; nlpnt; nlpnt=nlpnt->next,n2++) {

      if (nlpnt->conpnt) {
        npnt=(node *)(nlpnt->conpnt);
	match |= nodematch(npnt,n1a,n1b,n1c,n1d,exceptype,na,nb,nc,nd);
        avgx += npnt->xloc;			/* find average location */
        avgy += npnt->yloc;
        avgz += npnt->zloc;
      } 
    }
    if (n2) {
      x2 = avgx / n2;			/* average loc of second end comp */
      y2 = avgy / n2;
      z2 = avgz / n2;
    }
    if (n1 && n2) 
       if (match) drcconn(x1,y1,z1,x2,y2,z2);        /* draw the connection */

  }     /* for (cnpnt=connpnt;;) */
 }    /* if (disp&DCONN) */

#ifdef DEBUG
  if ((debug & NCPRCOMP) && (debugz & 2))
		ncfprintf (stderr,"dcomp end\n");
#endif

}
/*------------------------------------*/

int testlin(int nprint)
{
  ncfprintf (stdout, ", ");
  if (++nprint > 4) {
    nprint=2;
    ncfprintf (stdout,"\n#c        ");
  }
  else ncfprintf (stdout,"");
  return (nprint);
}

/*------------------------------------*/

void prelem (elem *epnt)

/* print parameter values for an element. */

{
    int i;
    cable *cepnt;
    sphere *sepnt;
    synapse *sypnt;
    FILE *str;

  str = stderr;
  switch (epnt->ctype) {

    case CABLE:
	cepnt = (cable *)epnt;
	ncfprintf (str,"'cable'   %d",cepnt->elnum);
	if (cepnt->length != NULLVAL) ncfprintf (str," length %-6.3g",cepnt->length);
	if (cepnt->dia    != NULLVAL) ncfprintf (str," dia %-6.3g",cepnt->dia);
	if (cepnt->dia2   != NULLVAL) ncfprintf (str," dia2 %-6.3g",cepnt->dia2);
	if (cepnt->Rm     != NULLVAL) ncfprintf (str," Rm %-6.4g",cepnt->Rm);
	if (cepnt->Ri     != NULLVAL) ncfprintf (str," Ri %-6.4g",cepnt->Ri);
	if (cepnt->cplam  != NULLVAL) ncfprintf (str," cplam %-6.3g",cepnt->cplam);
        break;

    case SPHERE:
	sepnt = (sphere *)epnt;
	ncfprintf (str,"'sphere'  %d",sepnt->elnum);
	if (sepnt->dia    != NULLVAL) ncfprintf (str," dia %-6.3g",sepnt->dia);
	if (sepnt->Rm     != NULLVAL) ncfprintf (str," Rm %-6.4g",sepnt->Rm);
        break;

    case SYNAPSE:
	sypnt = (synapse *)epnt;
	ncfprintf (str,"'synapse' %d",sypnt->elnum);
	if (sypnt->ntact  != NULLVAL) {
	   if (sypnt->ntact  == OPEN) 
		 ncfprintf (str," open");
	   else if (sypnt->ntact  == OPEN) 
		 ncfprintf (str," close");
	}
	if (sypnt->vrev   != NULLVAL) ncfprintf (str," vrev %-6.3g",sypnt->vrev);
	if (sypnt->thresh != NULLVAL) ncfprintf (str," thresh %-6.3g",sypnt->thresh);
	if (sypnt->nfilt1 != NULLVAL) {
	        ncfprintf(str," timec1");
		for (i=0; i<sypnt->nfilt1; i++) {	  
	         if (sypnt->timec1[i]!=0) ncfprintf(str," %-6.3g",sypnt->timec1[i]);
		}
	}
	if (sypnt->nfilt2 != NULLVAL) {
	        ncfprintf(str," timec2");
		for (i=0; i<sypnt->nfilt2; i++) {	  
	         if (sypnt->timec2[i]!=0)ncfprintf(str," %-6.3g",sypnt->timec2[i]);
		}
	}
	if (sypnt->tfall2 != NULLVAL) ncfprintf (str," tfall2 %-6.3g",sypnt->tfall2);
	if (sypnt->nfilt3 != NULLVAL) {
	        ncfprintf(str," timec3");
		for (i=0; i<sypnt->nfilt3; i++) {	  
	         if (sypnt->timec3[i]!=0)ncfprintf(str," %-6.3g",sypnt->timec3[i]);
		}
	}
	if (sypnt->tfall3 != NULLVAL) ncfprintf (str," tfall3 %-6.3g",sypnt->tfall3);
	if (sypnt->vgain  != NULLVAL) ncfprintf (str," vgain  %-6.3g",sypnt->vgain);
	if (sypnt->ngain  != NULLVAL) ncfprintf (str," ngain  %-6.3g",sypnt->ngain);
	if (sypnt->cgain  != NULLVAL) ncfprintf (str," cgain  %-6.3g",sypnt->cgain);
	if (sypnt->coff   != NULLVAL) ncfprintf (str," coff   %-6.3g",sypnt->coff);
	if (sypnt->maxcond!= NULLVAL) ncfprintf(str," maxcond %-6.3g",sypnt->maxcond);
	if (sypnt->nkd     != NULLVAL) ncfprintf(str," nkd %-6.3g",sypnt->nkd);
	if (sypnt->dyadelem > 0)
		 ncfprintf(str,"dyad  %-6.3g",sypnt->dyadelem);
        else
	if (sypnt->curve  != NULLVAL) {
	    if (sypnt->curve == LINEAR)
		 ncfprintf(str,"linear  %-6.3g",sypnt->ngain);
	    else if (sypnt->curve == EXPON)
		 ncfprintf(str,"expon  %-6.3g",sypnt->ngain);
	}
	if (sypnt->mesg1 != NULLVAL) {
		 ncfprintf(str," %s out     ",findsym(sypnt->mesg1));
	}
	if (sypnt->mesg2 != NULLVAL) {
		 ncfprintf(str," %s out     ",findsym(sypnt->mesg2));
	}
        break;

    case BUF: {
                vbuf *bepnt;

	bepnt = (vbuf *)epnt;
	ncfprintf (str,"'vbuf'    %d",bepnt->elnum);
	if (bepnt->delay  != NULLVAL) ncfprintf (str," delay %-6.3g",bepnt->delay);
	if (bepnt->offset != NULLVAL) ncfprintf (str," offset %-6.3g",bepnt->offset);
	if (bepnt->gain   != NULLVAL) ncfprintf (str," gain %-6.3g",bepnt->gain);
	if (bepnt->tau    != NULLVAL) ncfprintf (str," tau %-6.3g",bepnt->tau);
	if (bepnt->lphp   != NULLVAL) ncfprintf (str," lphp %-6d",bepnt->lphp);
        }
        break;

    case NBUF: {
                nbuf *bepnt;

	bepnt = (nbuf *)epnt;
	ncfprintf (str,"'nbuf'    %d",bepnt->elnum);
	if (bepnt->gain    != NULLVAL) ncfprintf (str," gain %-6.3g",bepnt->gain);
	if (bepnt->offset  != NULLVAL) ncfprintf (str," offset %-6.3g",bepnt->offset);
	if (bepnt->ntoffset!= NULLVAL) ncfprintf (str," ntoffset %-6.3g",bepnt->ntoffset);
	if (bepnt->ntrans  != NULLVAL) ncfprintf (str," ntrans %-6d",bepnt->ntrans);
        }
        break;

    case CACOMP: 
    case CHAN: {
		attrib *apnt;

        if (!(apnt=epnt->attpnt)) {
		ncfprintf (str,"No channel parameters available.\n");
		break;
        }
	else {
         chattrib *chapnt = (chattrib *)apnt;
         
	 ncfprintf (str,"'chan'    %d",epnt->elnum);
	 if (apnt->ctype  != NULLVAL) ncfprintf (str," %s type %d",
			findsym(apnt->ctype), apnt->stype);
	 if (chapnt->vrev   != NULLVAL) ncfprintf (str," vrev %g",chapnt->vrev);
	 if (chapnt->maxcond!= NULLVAL) ncfprintf (str," maxcond %g",chapnt->maxcond);
	}
       }
       break;
    case ROD: 
	 ncfprintf (str,"'rod'    %d",epnt->elnum);
         break;
    case CONE: 
	 ncfprintf (str,"'cone'    %d",epnt->elnum);
    case CHR: 
	 ncfprintf (str,"'chr'    %d",epnt->elnum);
       break;


    default: break;
  }
  if (epnt->nocondens) ncfprintf (str," saved %-6d",epnt->nocondens);
  if (epnt->modif) ncfprintf (str," modif %-6d",epnt->modif);
  ncfprintf (str,"\n");
}

/*------------------------------------*/
