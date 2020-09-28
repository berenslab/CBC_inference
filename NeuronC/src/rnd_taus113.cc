/* rnd_taus113.cc in "nc" */


/* rng/taus113.c
 *
 * Copyright (C) 2002 Atakan Gurkan
 * Based on the file taus.c which has the notice
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdlib.h>

/* This is a maximally equidistributed combined, collision free 
   Tausworthe generator, with a period ~2^{113}. The sequence is,

   x_n = (z1_n ^ z2_n ^ z3_n ^ z4_n)  

   b = (((z1_n <<  6) ^ z1_n) >> 13)
   z1_{n+1} = (((z1_n & 4294967294) << 18) ^ b)
   b = (((z2_n <<  2) ^ z2_n) >> 27)
   z2_{n+1} = (((z2_n & 4294967288) <<  2) ^ b)
   b = (((z3_n << 13) ^ z3_n) >> 21)
   z3_{n+1} = (((z3_n & 4294967280) <<  7) ^ b)
   b = (((z4_n <<  3)  ^ z4_n) >> 12)
   z4_{n+1} = (((z4_n & 4294967168) << 13) ^ b)

   computed modulo 2^32. In the formulas above '^' means exclusive-or 
   (C-notation), not exponentiation. 
   The algorithm is for 32-bit integers, hence a bitmask is used to clear 
   all but least significant 32 bits, after left shifts, to make the code 
   work on architectures where integers are 64-bit.

   The generator is initialized with 
   zi = (69069 * z{i+1}) MOD 2^32 where z0 is the seed provided
   During initialization a check is done to make sure that the initial seeds 
   have a required number of their most significant bits set.
   After this, the state is passed through the RNG 10 times to ensure the
   state satisfies a recurrence relation.

   References:
   P. L'Ecuyer, "Tables of Maximally-Equidistributed Combined LFSR Generators",
   Mathematics of Computation, 68, 225 (1999), 261--269.
     http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
   P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe Generators", 
   Mathematics of Computation, 65, 213 (1996), 203--213.
     http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps
   the online version of the latter contains corrections to the print version.
*/

/* Modified by R.G. Smith, Sept, 2016 */

/*  This taus113 RNG method has a longer period than, but is only
 *  15% slower than the taus2 method. */

/* Changes:

  Renamed taus113_set() to taus113_init(),
    changed order of args and added extra dummy arg to taus_init()
  Added taus113_setstate
  Added static global definition of the state variable pointer
    static taus113_state_t *taus113_state = NULL;
    #define state taus113_state

  Removed gsl types at end:

     static const gsl_rng_type taus113_type =
     const gsl_rng_type *gsl_rng_taus = &taus113_type;

*/

typedef struct
{
  unsigned long int z1, z2, z3, z4;
}
taus113_state_t;

/* make global state variable pointer so rng remembers */

static taus113_state_t taus113_state_buf = {0UL,0UL,0UL,0UL};
static taus113_state_t *taus113_state = &taus113_state_buf;

/*-----------------------------------------------------------*/

void taus113_setstate (void *vstate)
{
  taus113_state = (taus113_state_t *) vstate;
}

/*-----------------------------------------------------------*/

unsigned long taus113_get (void)
{
  static taus113_state_t *state;
  static unsigned long b1, b2, b3, b4;

#define MASK 0xffffffffUL

  state = (taus113_state_t *) taus113_state;
  b1 = ((((state->z1 << 6UL) & MASK) ^ state->z1) >> 13UL);
  state->z1 = ((((state->z1 & 4294967294UL) << 18UL) & MASK) ^ b1);

  b2 = ((((state->z2 << 2UL) & MASK) ^ state->z2) >> 27UL);
  state->z2 = ((((state->z2 & 4294967288UL) << 2UL) & MASK) ^ b2);

  b3 = ((((state->z3 << 13UL) & MASK) ^ state->z3) >> 21UL);
  state->z3 = ((((state->z3 & 4294967280UL) << 7UL) & MASK) ^ b3);

  b4 = ((((state->z4 << 3UL) & MASK) ^ state->z4) >> 12UL);
  state->z4 = ((((state->z4 & 4294967168UL) << 13UL) & MASK) ^ b4);

  return (state->z1 ^ state->z2 ^ state->z3 ^ state->z4);

}

/*-----------------------------------------------------------*/

double taus113_get_double (void)
{
  return taus113_get() / 4294967296.0;
}

/*-----------------------------------------------------------*/

void taus113_init (unsigned long int s, void *vstate, int siz)
{
  static taus113_state_t *state;

  taus113_setstate(vstate);
  state = taus113_state;
 
  if (!s)
    s = 1UL;                    /* default seed is 1 */

#define LCG(n) ((69069UL * n) & 0xffffffffUL)

  state->z1 = LCG (s);
  if (state->z1 < 2UL)
    state->z1 += 2UL;
  state->z2 = LCG (state->z1);
  if (state->z2 < 8UL)
    state->z2 += 8UL;
  state->z3 = LCG (state->z2);
  if (state->z3 < 16UL)
    state->z3 += 16UL;
  state->z4 = LCG (state->z3);
  if (state->z4 < 128UL)
    state->z4 += 128UL;

  /* Calling RNG ten times to satify recurrence condition */
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();
  taus113_get ();

  return;
}

