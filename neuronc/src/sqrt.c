/* sqrt() for the DEC Alpha */

/* Copyright (C) 1996 Free Software Foundation, Inc.
   Contributed by David Mosberger (davidm@cs.arizona.edu).
   Based on C source that is Copyright 1995 Linus Torvalds.

This file is part of the GNU C Library.

The GNU C Library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The GNU C Library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.

You should have received a copy of the GNU Library General Public
License along with the GNU C Library; see the file COPYING.LIB.  If
not, write to the Free Software Foundation, Inc., 675 Mass Ave,
Cambridge, MA 02139, USA.  */

/* This version is much faster than generic sqrt implementation, but
it doesn't handle exceptional values or the inexact flag.  Don't use
this if _IEEE_FP or _IEEE_FP_INEXACT is in effect. */

/*
 * This code is in the public domain, copyright 1995 Linus Torvalds.
 */
#include <errno.h>

static struct sqrt_data_struct {
	const int T2[64];
	const unsigned long half, one_and_a_half, two_to_minus_30, one;
	const unsigned long dn, up, Nan;
} sqrt_data = { {
	0x1500, 0x2ef8, 0x4d67, 0x6b02, 0x87be, 0xa395, 0xbe7a, 0xd866,
	0xf14a, 0x1091b,0x11fcd,0x13552,0x14999,0x15c98,0x16e34,0x17e5f,
	0x18d03,0x19a01,0x1a545,0x1ae8a,0x1b5c4,0x1bb01,0x1bfde,0x1c28d,
	0x1c2de,0x1c0db,0x1ba73,0x1b11c,0x1a4b5,0x1953d,0x18266,0x16be0,
	0x1683e,0x179d8,0x18a4d,0x19992,0x1a789,0x1b445,0x1bf61,0x1c989,
	0x1d16d,0x1d77b,0x1dddf,0x1e2ad,0x1e5bf,0x1e6e8,0x1e654,0x1e3cd,
	0x1df2a,0x1d635,0x1cb16,0x1be2c,0x1ae4e,0x19bde,0x1868e,0x16e2e,
	0x1527f,0x1334a,0x11051,0xe951, 0xbe01, 0x8e0d, 0x5924, 0x1edd },
	0x3fe0000000000000,	/* half */
	0x3ff8000000000000,	/* one_and_a_half */
	0x3e10000000000000,	/* two_to_minus_30 */
	0x3ff0000000000000,	/* one */
	0x3fefffffffffffff,	/* __dn = nextafter(1,-Inf) */
	0x3ff0000000000001,	/* __up = nextafter(1,+Inf) */
	0xffffffffffffffff	/* __Nan */	
};

#define lobits(x) (((unsigned int *)&x)[0])
#define hibits(x) (((unsigned int *)&x)[1])

static inline double initial_guess(double x, unsigned int k,
	struct sqrt_data_struct * ptr)
{
	double ret = 0.0;

	k = 0x5fe80000 - (k >> 1);
	k = k - ptr->T2[63&(k>>14)];
	hibits(ret) = k;
	return ret;
}

/* up = nextafter(1,+Inf), dn = nextafter(1,-Inf) */

#define __half			(ptr->half)
#define __one_and_a_half	(ptr->one_and_a_half)
#define __two_to_minus_30	(ptr->two_to_minus_30)
#define __one			(ptr->one)
#define __up			(ptr->up)
#define __dn			(ptr->dn)
#define __Nan			(ptr->Nan)

#define Double(x) (*(double *)&x)

/* Multiply with chopping rounding.. */
#define choppedmul(a,b,c) \
__asm__("multc %1,%2,%0":"=f" (c):"f" (a), "f" (b))

double
sqrt(double x)
{
  struct sqrt_data_struct * ptr;
  unsigned long k, bits;
  double y, z, zp, zn;
  double dn, up, low, high;
  double half, one_and_a_half, one, two_to_minus_30;

  ptr = &sqrt_data;
  *(double *)&bits = x;
  k = bits;

  /* Negative or NaN or Inf */
  if ((k >> 52) >= 0x7ff)
    goto special;
  y = initial_guess(x, k >> 32, ptr);
  half = Double(__half);
  one_and_a_half = Double(__one_and_a_half);
  y = y*(one_and_a_half - half*x*y*y);
  dn = Double(__dn);
  two_to_minus_30 = Double(__two_to_minus_30);
  y = y*((one_and_a_half - two_to_minus_30) - half*x*y*y);
  up = Double(__up);
  z = x*y;
  one = Double(__one);
  z = z + half*z*(one-z*y);

  choppedmul(z,dn,zp);
  choppedmul(z,up,zn);

  choppedmul(z,zp,low);
  low = low - x;
  choppedmul(z,zn,high);
  high = high - x;

  /* I can't get gcc to use fcmov's.. */
  __asm__("fcmovge %2,%3,%0"
	  :"=f" (z)
	  :"0" (z), "f" (low), "f" (zp));
  __asm__("fcmovlt %2,%3,%0"
	  :"=f" (z)
	  :"0" (z), "f" (high), "f" (zn));
  return z;	/* Argh! gcc jumps to end here */
special:
  /* throw away sign bit */
  k <<= 1;
  /* -0 */
  if (!k)
    return x;
  /* special? */
  if ((k >> 53) == 0x7ff) {
    /* NaN? */
    if (k << 11)
      return x;
    /* sqrt(+Inf) = +Inf */
    if (x > 0)
      return x;
  }
  /* -ve or -Inf */
  errno = EDOM;
  x = Double(__Nan);
  return x;
}

