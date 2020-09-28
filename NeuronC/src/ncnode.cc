/* ncnode.cc */

/* Routines to make and access nodes */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncelem.h"

// #define NHASHSIZ 199                /* prime number to make better hash */
// #define NHASHSIZ 4919               /* prime number to make better hash */
// #define NHASHSIZ 19937              /* prime number to make better hash */
// #define NHASHSIZ 99929              /* prime number to make better hash */
 #define NHASHSIZ 199373                /* prime number to make better hash */

// #define NHASHSIZC 29927             /* prime number to make better hash */
#define NHASHSIZC 4919                  /* prime number to make better hash */

#define NHASHSIZMIN 199                /* prime number to make better hash */

node **nhashtab = NULL;                   /* node hash table */
node **nhashtabe = NULL;                  /* node hash end table */
node **nhashtabc = NULL;                  /* node hash table for ct, cn only */
node **nhashtabce = NULL;                 /* node hash end table for ct, cn only */
     
int nhashsiz  = NHASHSIZ;
int nhashsizc = NHASHSIZC;

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef __cplusplus
}
#endif

#include "ncio.h"

extern node *nodepnt;
extern node *nodend;
extern int cumnode;
extern char *progname;
extern int ehashsiz;
extern elem *ehashtab[];
extern elem *enhashtab[];
extern int cumelem;

char *emalloc(unsigned int n);
void execerror(const char *s, const char *t);
void efree(void *ptr);

/*---------------------------------------------------*/

/*
void xdebugf(void)

{
   ncfprintf (stderr,"nhashtab %d\n",nhashtab[5]);
}
*/

/*---------------------------------------------------*/

void ninithash(void)
{
   int i;

   if (nhashtab !=NULL) efree (nhashtab);
   if (nhashtabe!=NULL) efree (nhashtabe);

   if ((nhashtab=(node **)emalloc(nhashsiz*sizeof(node*)))==NULL)
       ncfprintf (stderr,"ninithash: can't allocate hash table %d\n",nhashsiz);

   if ((nhashtabe=(node **)emalloc(nhashsiz*sizeof(node*)))==NULL)
       ncfprintf (stderr,"ninithash: can't allocate hash table %d\n",nhashsiz);

   for (i=0; i<nhashsiz; i++) {
     nhashtab[i] = NULL;
     nhashtabe[i] = NULL;
   }
}

/*---------------------------------------------------*/

void ninithashc(void)
{
   int i;

   if (nhashtabc !=NULL) efree (nhashtabc);
   if (nhashtabce!=NULL) efree (nhashtabce);
   if ((nhashtabc=(node **)emalloc(nhashsizc*sizeof(node*)))==NULL)
       ncfprintf (stderr,"ninithashc: can't allocate hashc table %d\n",nhashsizc);

   if ((nhashtabce=(node **)emalloc(nhashsizc*sizeof(node*)))==NULL)
       ncfprintf (stderr,"ninithashc: can't allocate hashc table %d\n",nhashsizc);

   for (i=0; i<nhashsizc; i++) {
     nhashtabc[i] = NULL;
     nhashtabce[i] = NULL;
   }
}

/*---------------------------------------------------*/

int nodehash(nodeint nodea, nodeint nodeb, 
	     nodeint nodec, nodeint noded)

/* Make disordered index of node number */

{
   int i;

   i = 1;
//   if (nodea > 1) i *= nodea; 
//   if (nodeb > 1) i *= nodeb; 
//   if (nodec > 1) i *= nodec; 
//   if (noded > 1) i *= noded; 
   if     (noded > 0) {
	              i += noded; 
       if (nodec > 0) i += nodec*100; 
       if (nodeb > 0) i += nodeb*10000; 
       if (nodea > 0) i += nodea*1000000; 
   } else {
       if (nodec > 0) i += nodec; 
       if (nodeb > 0) i += nodeb*10000; 
       if (nodea > 0) i += nodea*1000000; 
   }
   if (i<0) i = -i;
   return (i % nhashsiz); 
}

/*---------------------------------------------------*/

int nodehashc(nodeint nodea, nodeint nodeb)

/* Make disordered index of (ct,cn) node number */

{
   int i;

   i = 1;
//   if (nodea > 0) i *= nodea; 
//   if (nodeb > 0) i *= nodeb; 
   if (nodea > 0) i += nodea * 10000; 
   if (nodeb > 0) i += nodeb; 
   return (i % nhashsizc); 
}

/*---------------------------------------------------*/

node *getnpnt (nodeint nodenm1, nodeint nodenm2)

/* get starting node pointer for (ct,cn) hash */

{
  return (nhashtabc[nodehashc(nodenm1,nodenm2)]);
}

/*---------------------------------------------------*/

void ninstall (node *newnod)
              
/* Install a node in hash table using all 4 node numbers */
/*  Uses "direct chaining" with field "hnext". */

{
    int i;
    register node *npnt,*nlast;
    static int ninitfl=0;

   if (!ninitfl) {			/* initialize table once at start */
        ninitfl = 1;
        ninithash();
   }
   newnod->hnext = NULL;
   newnod->hlast = NULL;
   					
   i=nodehash(newnod->nodenm1,newnod->nodenm2,	/* initial index into nhashtab*/ 
	      newnod->nodenm3,newnod->nodenm4);

   if ((npnt=nhashtab[i])==NULL) {      /* install directly in table */
      nhashtab[i]  = newnod;
      nhashtabe[i] = newnod;
   }
   else {                               /* otherwise, go to end of list */
//      for (; npnt; npnt=npnt->hnext)    /* find last symbol */
//          nlast = npnt;    
//      nlast->hnext = newnod;            /* put new symbol at end of list */
//      newnod->hlast = nlast;

     nlast = nhashtabe[i];
     nlast->hnext = newnod;
     newnod->hlast = nlast;
     nhashtabe[i] = newnod;
   }
}

/*---------------------------------------------------*/

void ninstallc (node *newnod)
              
/* Install a node in hash table using ct, cn node numbers */
/*  Uses "direct chaining" with field "hcnext". */

{
    int i;
    register node *npnt,*nlast;
    static int ninitfl=0;

   if (!ninitfl) {			/* initialize table once at start */
        ninitfl = 1;
        ninithashc();
   }
   newnod->hcnext = NULL;
   newnod->hclast = NULL;
   i=nodehashc(newnod->nodenm1,newnod->nodenm2); /* initial index into nhashtab*/ 
   if ((npnt=nhashtabc[i])==NULL) {      /* install directly in table */
      nhashtabc[i]  = newnod;
      nhashtabce[i] = newnod;
   }
   else {                               /* otherwise, go to end of list */
//      for (; npnt; npnt=npnt->hcnext)    /* find last symbol */
//          nlast = npnt;    
//      nlast->hcnext = newnod;            /* put new symbol at end of list */
//      newnod->hclast = nlast;

      nlast = nhashtabce[i];
      nlast->hcnext = newnod;
      newnod->hclast = nlast;
      nhashtabce[i] = newnod;
   }
}

/*---------------------------------------------------*/

void printhash(void)

{
#define HBINSIZ 50
#define HNBINSIZ 50
     int i,j,n;
     int ncable;
     node *npnt;
     elem *epnt;
     int hbin[HBINSIZ]  = {0};
     int hcbin[HBINSIZ] = {0};
     int ebin[HBINSIZ]  = {0};
     int enbin[HNBINSIZ]  = {0};

   for (n=i=0; i<nhashsiz; i++) {
      for (j=0,npnt=nhashtab[i]; npnt; j++, npnt=npnt->hnext) ;
      n += j;
      if (j>=HBINSIZ) j = HBINSIZ-1;
      hbin[j]++;
   }
   fprintf(stderr,"# node hash table size %d\n", nhashsiz);
   fprintf(stderr,"# total nodes %d\n", n);
   for (i=0; i<HBINSIZ; i++) {
      fprintf(stderr,"# %d %d\n", i, hbin[i]);
   }

   fprintf(stderr,"# \n");
   for (n=i=0; i<nhashsizc; i++) {
      for (j=0,npnt=nhashtabc[i]; npnt; j++, npnt=npnt->hcnext) ;
      n += j;
      if (j>=HBINSIZ) j = HBINSIZ-1;
      hcbin[j]++;
   }
   fprintf(stderr,"# nodec hash table size %d\n", nhashsizc);
   fprintf(stderr,"# total nodesc %d\n", n);
   for (i=0; i<HBINSIZ; i++) {
      fprintf(stderr,"# %d %d\n", i, hcbin[i]);
   }

   fprintf(stderr,"# \n");
   for (n=i=0; i<ehashsiz; i++) {
      for (j=0,epnt=ehashtab[i]; epnt; j++, epnt=epnt->hnext) ;
      n += j;
      if (j>=HBINSIZ) j = HBINSIZ-1;
      ebin[j]++;
   }
   fprintf(stderr,"# elem hash table size %d\n", ehashsiz);
   fprintf(stderr,"# total elemsn %d cumelem %d\n", n, cumelem);
   for (i=0; i<HBINSIZ; i++) {
      fprintf(stderr,"# %d %d\n", i, ebin[i]);
   }

   fprintf(stderr,"# \n");
   ncable = 0;
   for (n=i=0; i<ehashsiz; i++) {
      for (j=0,epnt=enhashtab[i]; epnt; j++, epnt=epnt->hnnext) 
      	if (epnt->ctype==CABLE) ncable++;
      n += j;
      // fprintf(stderr,"# elem %d %d\n", i,j);
      if (j>=HNBINSIZ) j = HNBINSIZ-1;
      enbin[j]++;
   }
   fprintf(stderr,"# elem hash table size %d\n", ehashsiz);
   fprintf(stderr,"# total elemsn %d cumelem %d ncables %d\n", n, cumelem, ncable);
   for (i=0; i<HNBINSIZ; i++) {
      fprintf(stderr,"# %d %d\n", i, enbin[i]);
   }
}
/*---------------------------------------------------*/

node *maknod(nodeint nodea, nodeint nodeb, 
	     nodeint nodec, nodeint noded)
                          
/* make a new node, and link it to the
   beginning of the node list. */

{
    node *npnt;
 
  if ((npnt=(node *)emalloc(sizeof(node))) == NULL) {
     ncfprintf (stderr,"no space left for node %d\n",cumnode+1);
     return ((node *)NULL);  
  }
  npnt->next = nodepnt;		/* Set up "next" and "last" to connect all nodes. */
  npnt->last = NULL;
  if (nodepnt!=NULL) nodepnt->last = npnt;
  nodepnt = npnt;  		/*  We use "hnext" for hash table. */
  if (nodend==NULL) nodend = npnt;
  nodepnt->nodenm1 = nodea;
  nodepnt->nodenm2 = nodeb;
  nodepnt->nodenm3 = nodec;
  nodepnt->nodenm4 = noded;

  ninstall(npnt);	/* install in hash table */
  ninstallc(npnt);	/* install in ct, cn hash table */

  nodepnt->elemlst = (conlst *)NULL;
  nodepnt->comptr = (comp *)NULL;
  nodepnt->compnum = 0;
  nodepnt->ctype = NODE;
/*
  if (nodea > MAXSHORT || nodeb > MAXSHORT || 
      nodec > MAXSHORT || noded > MAXSHORT)   {
    ncfprintf (stderr,"Node number %d %d %d %d too large\n",
			nodea,nodeb,nodec,noded);
    execerror ("Bad node value: ","stopping...");
  }
*/
  if (nodea < 0 || nodeb < 0 || nodec < 0 || nodec < 0) {
    if (nodea > NULLVAL && nodeb > NULLVAL && 
	nodec > NULLVAL && noded > NULLVAL) { 
    ncfprintf (stderr,"Node number %d %d %d %d is negative\n",
			nodea,nodeb,nodec,noded);
    execerror ("Bad node value: ","stopping...");
   }
  }
  nodepnt->xloc = LARGENODE;
  nodepnt->yloc = LARGENODE;
  nodepnt->zloc = LARGENODE;
  nodepnt->region = NULLVAL;
  nodepnt->dendn = NULLVAL;
  nodepnt->label = NULLVAL;
  nodepnt->labeltext = NULL;
  cumnode++;                    /* increment total */
  return (nodepnt); 
}

/*---------------------------------------------------*/

void rehashnodes(void)

/* Erase initial hash tables and rehash all nodes */
/*  Need some way to expand hash tables during the time that */
/*  nodes are being created, i.e. dynamic allocation. */

{
    int i;
    register node* npnt;

  //nhashsiz = cumnode / 10.0;	/* calculate new hash table size */
  nhashsiz = cumnode / 5.0;	/* calculate new hash table size */
  if (nhashsiz < NHASHSIZMIN) nhashsiz = NHASHSIZMIN;

  ninithash();		/* erase and reinitialize hash tables */
  ninithashc();

  for (npnt=nodepnt; npnt; npnt=npnt->next) {

       ninstall(npnt);	/* install in hash table */
       ninstallc(npnt);	/* install in ct, cn hash table */
  } 
}

/*---------------------------------------------------*/

void delnode (node *npnt)

/* Delete a node, given a pointer to it. */
/* Don't delete node's element list because this is done by "erasenode()" */

{
   int i,found;

  if (!npnt) return;
  if (npnt->comptr) {
     ncfprintf (stderr,"delnode: can't delete node: it has a compartment.\n");
	  return;			/* stop if compartment */ 
  }

/* ncfprintf (stderr," %d %d %d %d\n", 
	npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4); */

	/* First, delete from hash table: */

  i = nodehash (npnt->nodenm1,npnt->nodenm2,npnt->nodenm3,npnt->nodenm4);
  
  if (npnt->hnext!=NULL) npnt->hnext->hlast = npnt->hlast;     /* delete hash next pointer */
  if (npnt->hlast!=NULL) npnt->hlast->hnext = npnt->hnext;     /* delete hash last pointer */
  else nhashtab[i] = npnt->hnext;
  if (nhashtabe[i]==npnt) nhashtabe[i] = npnt->hlast;

	/* Delete from ct, cn hash table: */

  i = nodehashc (npnt->nodenm1,npnt->nodenm2);
  
  if (npnt->hcnext!=NULL) npnt->hcnext->hclast = npnt->hclast;     /* delete hash next pointer */
  if (npnt->hclast!=NULL) npnt->hclast->hcnext = npnt->hcnext;     /* delete hash last pointer */
  else nhashtabc[i] = npnt->hcnext;
  if (nhashtabce[i]==npnt) nhashtabce[i] = npnt->hclast;

	/* Patch next, last pointers */

  if (npnt->last!=NULL) npnt->last->next = npnt->next;     /* patch last pointer */
  if (npnt->next!=NULL) npnt->next->last = npnt->last;     /* patch next pointer */
  if (nodepnt==npnt) nodepnt = npnt->next; 
  if (nodend==npnt) nodend = npnt->last; 
  efree (npnt);
  cumnode--;
}

/*---------------------------------------------------*/

// long int total_lookups=0;
// long int total_finds=0;

node *findnode(nodeint num1, nodeint num2, 
		nodeint num3, nodeint num4, const char *str)
                      
/* Find node among list of all nodes.
   Nodes are placed in hash table, which provides
   faster access than possible with one sequential list.
   The node itself is not moved, therefore only the 
   pointers of the node list need rearranging.
   Pointers in other lists that point to nodes remain correct.
*/

{
   register node *npnt;
   int i,j,found;
   char sbuf[100];
   static int err=0;

  i = nodehash(num1,num2,num3,num4);
  for (j=0,found=0,npnt = nhashtab[i]; npnt; npnt = npnt->hnext, j++) {
    if ((npnt->nodenm2==num2) && (npnt->nodenm3==num3) &&
        (npnt->nodenm1==num1) && (npnt->nodenm4==num4)) {
       found = 1;
       break;
    }
  }

//   total_lookups += ++j;
//   total_finds ++;

  if (found) return npnt;
  else {
    if (str) {
      if (num4 != NULLVAL)
         sprintf (sbuf,"%.40s: can't find node %d %d %d %d\n",str,num1,num2,num3,num4);
      else if (num3 != NULLVAL)
         sprintf (sbuf,"%.40s: can't find node %d %d %d\n",str,num1,num2,num3);
      else if (num2 != NULLVAL)
         sprintf (sbuf,"%.40s: can't find node %d %d\n",str,num1,num2);
      else 
         sprintf (sbuf,"%.40s: can't find node %d\n",str,num1);

      if (err<5) {
         execerror("nc: ",sbuf);
	 // ncfprintf (stderr,sbuf);
	 err++;
      }
    }
   return (node *)NULL; 
  }
}

/*---------------------------------------------------*/

node *findnode(nodeint num1, nodeint num2, nodeint num3, const char *str)
                      
/* Find node among list of all nodes.
   Nodes are placed in hash table, which provides
   faster access than possible with one sequential list.
   The node itself is not moved, therefore only the 
   pointers of the node list need rearranging.
   Pointers in other lists that point to nodes remain correct.
*/


{
   register node *npnt;
   int i,j,found;
   char sbuf[100];
   static int err=0;

  i = nodehash(num1,num2,num3,-1);
  for (npnt = nhashtab[i]; npnt; npnt = npnt->hnext) {
    if (npnt->nodenm3!=num3) continue;
    if (npnt->nodenm2!=num2) continue;
    if (npnt->nodenm1!=num1) continue;
    return (npnt);
  }

//   total_lookups += ++j;
//   total_finds ++;

 if (str) {
     if (num3 != NULLVAL)
         sprintf (sbuf,"findnode: can't find node %d %d %d\n",num1,num2,num3);
     else if (num2 != NULLVAL)
         sprintf (sbuf,"findnode: can't find node %d %d\n",num1,num2);
     else  sprintf (sbuf,"findnode: can't find node %d\n",num1);
     if (err<5) {
        execerror(sbuf,NULL);
        // ncfprintf (stderr,sbuf);
        err++;
     }
   }
   return (node *)NULL; 
}

/*---------------------------------------------------*/

node *findnode(nodeint num1, nodeint num2, nodeint num3)

/* Find node among list of all nodes.
   Nodes are placed in hash table, which provides
   faster access than possible with one sequential list.
   The node itself is not moved, therefore only the 
   pointers of the node list need rearranging.
   Pointers in other lists that point to nodes remain correct.
*/


{
   register node *npnt;
   int i,j,found;
   char sbuf[100];
   static int err=0;

  i = nodehash(num1,num2,num3,-1);
  for (npnt = nhashtab[i]; npnt; npnt = npnt->hnext) {
    if (npnt->nodenm3!=num3) continue;
    if (npnt->nodenm2!=num2) continue;
    if (npnt->nodenm1!=num1) continue;
    return (npnt);
  }
  return (node *)NULL; 
}

/*---------------------------------------------------*/
/*
node *findnode(nodeint num1, nodeint num2, nodeint num3)

{
   findnode(num1, num2, num3, NULLVAL, NULL);
}
/* */

/*---------------------------------------------------*/

/*
node *findnode(nodeint num1, nodeint num2, nodeint num3, const char *str)

{
   findnode(num1, num2, num3, NULLVAL, str);
}
/* */

/*---------------------------------------------------*/

 
node *findnode(nodeint num1, nodeint num2)

{
   register node *npnt;
   int i,j,found;
   char sbuf[100];
   static int err=0;

  i = nodehashc(num1,num2);
  for (j=0,found=0,npnt = nhashtabc[i]; npnt; npnt = npnt->hcnext, j++) {
      if (npnt->nodenm2!=num2) continue;
      if (npnt->nodenm1!=num1) continue; 
      return (npnt);
  }

   if (num2 != NULLVAL)
         sprintf (sbuf,"findnode: can't find node %d %d\n",num1,num2);
   else  sprintf (sbuf,"findnode: can't find node %d\n",num1);
   if (err<5) {
        execerror(sbuf,NULL);
        // ncfprintf (stderr,sbuf);
        err++;
   }
   return (node *)NULL; 
}


/*---------------------------------------------------*/
/*
node *findnode(nodeint num1, nodeint num2, const char *str)

{
   findnode(num1, num2, NULLVAL, NULLVAL, str);
}

*/

/*---------------------------------------------------*/

node *findnode(nodeint num1)

{
   findnode(num1, NULLVAL);
}

/*---------------------------------------------------*/

