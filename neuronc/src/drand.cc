/* module drand in program nc */

/* random number function */

#include <stdint.h>

#include "ncsub.h"
#include "ncomp.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
int getpid(void);

#ifdef __cplusplus
}
#endif

#include "ncio.h"

extern int rseed;

/* This set of rand() functions can be set to run either      */
/*  the glibc2 version (in rnd_glibc.cc) or one of the GSL    */
/*  set of rng functions, "taus" (in "rnd_taus.cc").          */
/*  Also to "mt" (mersenne twister) (in "rnd_mt.cc")          */

/* A good random number generator is hard to find. The glibc2 */
/* generator has correlations for state sizes less than 128,  */ 
/*  as tested by the "diehard" program:                       */
/*                                                            */
/*    http://stat.fsu.edu/~geo/diehard.html                   */

/* The GSL "taus" generator is smaller, just as fast, and is  */
/* supposed to be one of the best in terms of its randomness. */

/* To set which version to use, uncomment one of the sets     */
/*  of definitions below, and set RNDSIZ in "ncomp.h" to      */
/*  the same size as RANDSTATSIZ (or bigger)                  */
/*  Then, remember to "make clean" and "make" in nc/src       */

/* For "rnd_glibc.cc" (from glibc2 library, use these definitions: */

/*
#define rinitstate    r_initstate
#define rsetstate     r_setstate
#define random        r_random
#define random_double r_random_double
#define RANDSTATSIZ 256

char *r_initstate(unsigned long int seed, void *state, int n);
char *r_setstate(void *state);
int r_random(void);
double r_random_double(void);
/* */

/* for"rnd_taus.cc" (from Gnu Sci Library), use these definitions: */
/*
#define rinitstate    taus_init
#define rsetstate     taus_setstate
#define random        taus_get
#define random_double taus_get_double

#define RANDSTATSIZ (3*sizeof(long int))

void taus_init(unsigned long int seed, void *state, int n);
void taus_setstate(void *state);
unsigned long taus_get(void);
double taus_get_double(void);

/* */

/* for"rnd113_taus.cc" (from Gnu Sci Library), use these definitions: */

#define rinitstate    taus113_init
#define rsetstate     taus113_setstate
#define random        taus113_get
#define random_double taus113_get_double

#define RANDSTATSIZ (4*sizeof(long int))

void taus113_init(unsigned long int seed, void *state, int n);
void taus113_setstate(void *state);
unsigned long taus113_get(void);
double taus113_get_double(void);

/* */

/* for"rnd_tinymt64.cc", use these definitions: */
/*
#define rinitstate    tinymt_init
#define rsetstate     tinymt_setstate
#define random        tinymt_get
#define random_double tinymt_get_double

#define RANDSTATSIZ (8*sizeof(uint32_t))

void tinymt_init(unsigned long long int seed, void *state, int n);
void tinymt_setstate(void *state);
unsigned long tinymt_get(void);
double tinymt_get_double(void);

/* */

/* for"rnd_mt.cc" (from Gnu Sci Library), use these definitions: */
/*
#define rinitstate    mt_init
#define rsetstate     mt_setstate
#define random        mt_get
#define random_double mt_get_double

#define RANDSTATSIZ (625*sizeof(long int))

void mt_init(unsigned long int seed, void *state, int n);
void mt_setstate(void *state);
unsigned long mt_get(void);
double mt_get_double(void);

/* */

void execerror(char *s, char *t);

/*-----------------------------------------------------------*/

/* Call setrand(seed) first to initialize.  */

static char _rstate[RANDSTATSIZ] = {0};

double drand(void)
{
  return (random_double());
}

/*-----------------------------------------------------------*/

double drand_pos(void)

{
   double posval;

 for (posval=0; (posval=random_double()) <= 0; ); 
  return (posval);
}


/*-----------------------------------------------------------*/

unsigned long irand(void)
{
  return (random());
}

/*-----------------------------------------------------------*/

/* Call this function to invoke a numbered RNG */
 
double rrand(int ngen)
{
  randstate *rpnt;
  randstate *findrand(int ngen);
  void restorstate();
  double retval;

  rpnt = findrand(ngen);
  if (rpnt) {
    rsetstate (rpnt->rstate);
    retval = random_double();
    restorstate();
  }
  else { 
   retval = -10;
  }
  return retval;
}


/*-----------------------------------------------------------*/

void initstate(unsigned long int val, char *rstate, int size)

/* Call this procedure to initialize random state. */

{
  rinitstate (val,(void*)rstate,size);
}

/*-----------------------------------------------------------*/

void drand_setstate(char *rstate)

/* Call this procedure to set random state. */

{
  rsetstate ((void*)rstate);
}

/*-----------------------------------------------------------*/

void restorstate(void)

/* Call this procedure to restore default random state. */

{
  rsetstate (_rstate);
}

/*-----------------------------------------------------------*/

void setrand(int val)

/* Call this prodecure with seed to initialize entire RNG. */

{
  if (val<=0) val = getpid();
  rinitstate (val,_rstate,RANDSTATSIZ); 
  rseed = val;
}

/*-----------------------------------------------------------*/

#define RHASHSIZ 199             /* prime number to make better hash */

randstate *rhashtab[RHASHSIZ] = {0};  /* noise generator hash table */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef __cplusplus
}
#endif

extern int cumrand;

char *emalloc(unsigned int n);

/*---------------------------------------------------*/

void rinithash(void)
{
   int i;

   for (i=0; i<RHASHSIZ; i++) 
     rhashtab[i] = 0;
}

/*---------------------------------------------------*/

int randhash(int n)

/* Make disordered index of rand noise generator number */

{
   if (n<0) n = -n;
   return (n % RHASHSIZ); 
}

/*---------------------------------------------------*/

void initrand (int ngen, int nrseed)
              
/* Initialize and install a noise generator in hash table */
/*  Uses "direct chaining" with field "next". */

{
    int i,found;
    randstate *rpnt,*rlast,*stpnt;
    static int rinitfl=0;

  /* look first for existing random noise generator */

  i = randhash(ngen);
  for (found=0,rpnt = rhashtab[i]; rpnt; rpnt = rpnt->next) {
    if (rpnt->ngen==ngen) {
       found = 1;
       break;
    }
  }

  if (!found) {
     if ((stpnt=(randstate *)emalloc(sizeof(randstate))) == NULL) {
       ncfprintf (stderr,"No space left for rand noise gen %d\n",ngen);
     }
   }
   else stpnt = rpnt;

   if (stpnt) {
	if (nrseed <= 0) nrseed = 45833^ngen^rseed;
	rinitstate (nrseed,stpnt->rstate,RNDSIZ);
   }
   ++cumrand;

   if (!rinitfl) {			/* initialize table once at start */
        rinitfl = 1;
        rinithash();
   }
   stpnt->next = 0;
   stpnt->ngen = ngen;
   i=randhash(ngen); 			/* initial index into nhashtab*/ 
   if (!(rpnt=rhashtab[i])) {           /* install directly in table */
      rhashtab[i] = stpnt;
   }
   else {                               /* otherwise, go to end of list */
      for (; rpnt; rpnt=rpnt->next)     /* find last one */
          rlast = rpnt;    
      rlast->next = stpnt;              /* put new one at end of list */
   }
}

/*---------------------------------------------------*/

void delrand (randstate *rpnt)

/* Delete a random noise generator, given a pointer to it. */

{
   int i,found;
   randstate *rlast,*rpt;
 
  if (!rpnt) return;

	/* First, delete from hash table: */

  i = randhash(rpnt->ngen);
  rlast = (randstate*)NULL;
  for (found=0, rpt=rhashtab[i]; rpt; rlast=rpt, rpt=rpt->next) {
    if (rpt==rpnt) {
        found = 1;
        break;
    }
  }
  if (found) {
     if (rlast) rlast->next = rpt->next;     /* delete hash pointer */
     else rhashtab[i] = rpt->next;
  }
}

/*---------------------------------------------------*/

randstate *findrand(int ngen)
                      
/* Find random noise generator among list of them. 
   They are placed in hash table, which provides
   faster access than possible with one sequential list.
*/

{
   randstate *rpnt;
   int i,found;
        
  i = randhash(ngen);
  for (found=0,rpnt = rhashtab[i]; rpnt; rpnt = rpnt->next) {
    if (rpnt->ngen==ngen) {
       found = 1;
       break;
    }
  }

  if (found) return rpnt;
  else {
     ncfprintf (stderr,"findrand: can't find random noise gen %d\n",ngen);
  }
  return (randstate *)NULL; 
}

/*---------------------------------------------------*/

