
/* rndtest, a simple program to test the random number generator. */

/* Use this to generate a 10 MB binary file for use with  */
/*  "diehard" (die.c in nc/tests, a tester for randomness). */

#include <stdio.h>

/* Use one of three random number generators: */
/*

 glibc  (several sizes, small and fast or large and hi-quality 

 taus   (only 3 4-byte integers, 12 bytes, fastest in gsl)

 mt     (mersenne twister, longest random period, fast but not fastest 

 See http://www.gnu.org/software/gsl for more on these random number generators.
 Uncomment below the one that you wish to test. 

*/

/* Use these definitions for "random.cc" (from glibc2 library) */

/*
#define rinitstate    r_initstate
#define rsetstate     r_setstate
#define random        r_random
#define random_double r_random_double
#define RANDSTATSIZ 128

char *r_initstate(unsigned long int seed, void *state, int n);
char *r_setstate(void *state);
int r_random(void);
double r_random_double(void);
/* */

/* Use these definitions for "taus.cc" (from gnu sci library) */

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

/* for"rnd113_taus.cc" (from Gnu Sci Library), use these
 * definitions: */
/*
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

main()

{
   long int i,n,x,x1;
   static char buf[RANDSTATSIZ];

  rinitstate(347321,buf,RANDSTATSIZ);
  n = 10000000;
  for (i=0; i<n; i++) {
     x=random();
     fwrite (&x,1,1,stdout);
  }
}

