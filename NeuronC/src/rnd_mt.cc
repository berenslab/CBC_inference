/* rnd_mt.cc in "nc" */

/* From GSL rng/mt.c */

/* This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this program;
   if not, write to the Free Foundation, Inc., 59 Temple Place, Suite
   330, Boston, MA 02111-1307 USA

   Original implementation was copyright (C) 1997 Makoto Matsumoto and
   Takuji Nishimura. Coded by Takuji Nishimura, considering the
   suggestions by Topher Cooper and Marc Rieffel in July-Aug. 1997, "A
   C-program for MT19937: Integer version (1998/4/6)"

   This implementation copyright (C) 1998 Brian Gough. I reorganized
   the code to use the module framework of GSL.  The license on this
   implementation was changed from LGPL to GPL, following paragraph 3
   of the LGPL, version 2.

   Update:

   The seeding procedure has been updated to match the 10/99 release
   of MT19937.

   Update:

   The seeding procedure has been updated again to match the 2002
   release of MT19937

   The original code included the comment: "When you use this, send an
   email to: matumoto@math.keio.ac.jp with an appropriate reference to
   your work".

   Makoto Matsumoto has a web page with more information about the
   generator, http://www.math.keio.ac.jp/~matumoto/emt.html. 

   The paper below has details of the algorithm.

   From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
   623-dimensionally equidistributerd uniform pseudorandom number
   generator". ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1 (Jan. 1998), Pages 3-30

   You can obtain the paper directly from Makoto Matsumoto's web page.

   The period of this generator is 2^{19937} - 1.

*/

/* Modified by RGS 2007/02/16

   changes

   Removed mt_1998_set()
    and    mt_1999_set()
   Renamed mt_set() to mt_init()
     changed order of args and added extra dummy arg.
   Also note that the size of the state set by mt_setstate
   must be correctly set in ncomp.h (RNDSIZ).

*/
 
// #include <config.h>
#include <stdlib.h>
// #include <gsl/gsl_rng.h>

#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;   

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;   

typedef struct
  {
    unsigned long mt[N];
    int mti;
  }
mt_state_t;

/* make global state variable pointer so rng remembers */

static mt_state_t mt_state_buf = {0};
static mt_state_t *mt_state = &mt_state_buf;

/*------------------------------------------------------*/

void mt_setstate (void *vstate)
{
  mt_state = (mt_state_t *) vstate;
}

/*------------------------------------------------------*/

unsigned long mt_get (void)
{
  static unsigned long int k,y;
  static int kk;
  static mt_state_t *state;
  static unsigned long int *mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

  state = (mt_state_t *) mt_state;
  mt = state->mt;
  if (state->mti >= N)
    {   /* generate N words at one time */

      for (kk = 0; kk < N - M; kk++)
        {
          y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
      for (; kk < N - 1; kk++)
        {
          y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

      {
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
      }

      state->mti = 0;
    }

  /* Tempering */
  
  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);

  state->mti++;

  return k;
}

/*------------------------------------------------------*/

double mt_get_double (void)
{
  return mt_get () / 4294967296.0 ;
}

/*------------------------------------------------------*/

void mt_init (unsigned long int s, void *vstate, int siz)
{
  static int i;
  static mt_state_t *state;

  mt_setstate(vstate);
  state = (mt_state_t *) mt_state;

  if (s == 0)
    s = 4357;   /* the default seed is 4357 */

  state->mt[0]= s & 0xffffffffUL;

  for (i = 1; i < N; i++)
    {
      /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
         Ed. p.106 for multiplier. */

      state->mt[i] =
        (1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
      
      state->mt[i] &= 0xffffffffUL;
    }

  state->mti = i;
}

