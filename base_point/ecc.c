

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/times.h>
//#include <iostream>
//#include "tommath.h"
#include <time.h>
#include <unistd.h>


#define BIT_LEN 800 
#define KEY_LONG 128  //私钥比特长
#define P_LONG 200    //有限域P比特长
#define EN_LONG 40    //一次取明文字节数(x,20)(y,20)

#define  Max 0xfff
//-------------------------------------------
//-------------------------------------------

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "ecc.h"
//#include <time.h>

int KARATSUBA_MUL_CUTOFF = 74;
int KARATSUBA_SQR_CUTOFF = 124;
int TOOM_MUL_CUTOFF = 350;
int TOOM_SQR_CUTOFF = 400;
 
const char *mp_s_rmap = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/";  //??? ?????


const mp_digit ltm_prime_tab[] = {
	0x0002, 0x0003, 0x0005, 0x0007, 0x000B, 0x000D, 0x0011, 0x0013,
	0x0017, 0x001D, 0x001F, 0x0025, 0x0029, 0x002B, 0x002F, 0x0035,
	0x003B, 0x003D, 0x0043, 0x0047, 0x0049, 0x004F, 0x0053, 0x0059,
	0x0061, 0x0065, 0x0067, 0x006B, 0x006D, 0x0071, 0x007F,
#ifndef MP_8BIT
	0x0083,
	0x0089, 0x008B, 0x0095, 0x0097, 0x009D, 0x00A3, 0x00A7, 0x00AD,
	0x00B3, 0x00B5, 0x00BF, 0x00C1, 0x00C5, 0x00C7, 0x00D3, 0x00DF,
	0x00E3, 0x00E5, 0x00E9, 0x00EF, 0x00F1, 0x00FB, 0x0101, 0x0107,
	0x010D, 0x010F, 0x0115, 0x0119, 0x011B, 0x0125, 0x0133, 0x0137,

	0x0139, 0x013D, 0x014B, 0x0151, 0x015B, 0x015D, 0x0161, 0x0167,
	0x016F, 0x0175, 0x017B, 0x017F, 0x0185, 0x018D, 0x0191, 0x0199,
	0x01A3, 0x01A5, 0x01AF, 0x01B1, 0x01B7, 0x01BB, 0x01C1, 0x01C9,
	0x01CD, 0x01CF, 0x01D3, 0x01DF, 0x01E7, 0x01EB, 0x01F3, 0x01F7,
	0x01FD, 0x0209, 0x020B, 0x021D, 0x0223, 0x022D, 0x0233, 0x0239,
	0x023B, 0x0241, 0x024B, 0x0251, 0x0257, 0x0259, 0x025F, 0x0265,
	0x0269, 0x026B, 0x0277, 0x0281, 0x0283, 0x0287, 0x028D, 0x0293,
	0x0295, 0x02A1, 0x02A5, 0x02AB, 0x02B3, 0x02BD, 0x02C5, 0x02CF,

	0x02D7, 0x02DD, 0x02E3, 0x02E7, 0x02EF, 0x02F5, 0x02F9, 0x0301,
	0x0305, 0x0313, 0x031D, 0x0329, 0x032B, 0x0335, 0x0337, 0x033B,
	0x033D, 0x0347, 0x0355, 0x0359, 0x035B, 0x035F, 0x036D, 0x0371,
	0x0373, 0x0377, 0x038B, 0x038F, 0x0397, 0x03A1, 0x03A9, 0x03AD,
	0x03B3, 0x03B9, 0x03C7, 0x03CB, 0x03D1, 0x03D7, 0x03DF, 0x03E5,
	0x03F1, 0x03F5, 0x03FB, 0x03FD, 0x0407, 0x0409, 0x040F, 0x0419,
	0x041B, 0x0425, 0x0427, 0x042D, 0x043F, 0x0443, 0x0445, 0x0449,
	0x044F, 0x0455, 0x045D, 0x0463, 0x0469, 0x047F, 0x0481, 0x048B,

	0x0493, 0x049D, 0x04A3, 0x04A9, 0x04B1, 0x04BD, 0x04C1, 0x04C7,
	0x04CD, 0x04CF, 0x04D5, 0x04E1, 0x04EB, 0x04FD, 0x04FF, 0x0503,
	0x0509, 0x050B, 0x0511, 0x0515, 0x0517, 0x051B, 0x0527, 0x0529,
	0x052F, 0x0551, 0x0557, 0x055D, 0x0565, 0x0577, 0x0581, 0x058F,
	0x0593, 0x0595, 0x0599, 0x059F, 0x05A7, 0x05AB, 0x05AD, 0x05B3,
	0x05BF, 0x05C9, 0x05CB, 0x05CF, 0x05D1, 0x05D5, 0x05DB, 0x05E7,
	0x05F3, 0x05FB, 0x0607, 0x060D, 0x0611, 0x0617, 0x061F, 0x0623,
	0x062B, 0x062F, 0x063D, 0x0641, 0x0647, 0x0649, 0x064D, 0x0653
#endif
};

int mp_init(mp_int * a);


/* pre-calculate the value required for Barrett reduction
 * For a given modulus "b" it calulates the value required in "a"
 */
int mp_reduce_setup (mp_int * a, mp_int * b)
{
  int     res;
  
  if ((res = mp_2expt (a, b->used * 2 * DIGIT_BIT)) != MP_OKAY) {
    return res;
  }
  return mp_div (a, b, a, NULL);
}


/* computes a = 2**b 
 *
 * Simple algorithm which zeroes the int, grows it then just sets one bit
 * as required.
 */
int
mp_2expt (mp_int * a, int b)
{
  int     res;

  /* zero a as per default */
  mp_zero (a);

  /* grow a to accomodate the single bit */
  if ((res = mp_grow (a, b / DIGIT_BIT + 1)) != MP_OKAY) {
    return res;
  }

  /* set the used count of where the bit will go */
  a->used = b / DIGIT_BIT + 1;

  /* put the single bit in its place */
  a->dp[b / DIGIT_BIT] = ((mp_digit)1) << (b % DIGIT_BIT);

  return MP_OKAY;
}

 int mp_reduce_2k(mp_int *a, mp_int *n, mp_digit d)
 {
    mp_int q;
    int    p, res;
 
    if ((res = mp_init(&q)) != MP_OKAY) {
       return res;
    }
 
    p = mp_count_bits(n);
 top:
    /* q = a/2**p, a = a mod 2**p */
    if ((res = mp_div_2d(a, p, &q, a)) != MP_OKAY) {
       goto ERR;
    }
 
    if (d != 1) {
       /* q = q * d */
       if ((res = mp_mul_d(&q, d, &q)) != MP_OKAY) {
          goto ERR;
       }
    }
 
    /* a = a + q */
    if ((res = s_mp_add(a, &q, a)) != MP_OKAY) {
       goto ERR;
    }
 
    if (mp_cmp_mag(a, n) != MP_LT) {
       s_mp_sub(a, n, a);
       goto top;
    }
 
 ERR:
    mp_clear(&q);
    return res;
}

 int
 mp_dr_reduce (mp_int * x, mp_int * n, mp_digit k)
 {
   int      err, i, m;
   mp_word  r;
   mp_digit mu, *tmpx1, *tmpx2;
 
   /* m = digits in modulus */
   m = n->used;
 
   /* ensure that "x" has at least 2m digits */
   if (x->alloc < m + m) {
     if ((err = mp_grow (x, m + m)) != MP_OKAY) {
       return err;
     }
   }
 
 /* top of loop, this is where the code resumes if
  * another reduction pass is required.
  */
 top:
   /* aliases for digits */
   /* alias for lower half of x */
   tmpx1 = x->dp;
 
   /* alias for upper half of x, or x/B**m */
   tmpx2 = x->dp + m;
 
   /* set carry to zero */
   mu = 0;
 
   /* compute (x mod B**m) + k * [x/B**m] inline and inplace */
   for (i = 0; i < m; i++) {
       r         = ((mp_word)*tmpx2++) * ((mp_word)k) + *tmpx1 + mu;
       *tmpx1++  = (mp_digit)(r & MP_MASK);
       mu        = (mp_digit)(r >> ((mp_word)DIGIT_BIT));
   }
 
   /* set final carry */
   *tmpx1++ = mu;
 
   /* zero words above m */
   for (i = m + 1; i < x->used; i++) {
       *tmpx1++ = 0;
   }
 
   /* clamp, sub and return */
   mp_clamp (x);
 
   /* if x >= n then subtract and reduce again
    * Each successive "recursion" makes the input smaller and smaller.
    */
   if (mp_cmp_mag (x, n) != MP_LT) {
     s_mp_sub(x, n, x);
     goto top;
   }
   return MP_OKAY;
 }

 void mp_dr_setup(mp_int *a, mp_digit *d)
 {
    /* the casts are required if DIGIT_BIT is one less than
     * the number of bits in a mp_digit [e.g. DIGIT_BIT==31]
     */
    *d = (mp_digit)((((mp_word)1) << ((mp_word)DIGIT_BIT)) -
         ((mp_word)a->dp[0]));
 }

 int
 mp_montgomery_setup (mp_int * n, mp_digit * rho) 
 {     
   mp_digit x, b;
    
 /* fast inversion mod 2**k
  * 
  * Based on the fact that
  *    
  * XA = 1 (mod 2**n)  =>  (X(2-XA)) A = 1 (mod 2**2n)
  *                    =>  2*X*A - X*X*A*A = 1
  *                    =>  2*(1) - (1)     = 1
  */
   b = n->dp[0]; 
       
   if ((b & 1) == 0) {
     return MP_VAL;
   }          
          
   x = (((b + 2) & 4) << 1) + b; /* here x*a==1 mod 2**4 */
   x *= 2 - b * x;               /* here x*a==1 mod 2**8 */
 #if !defined(MP_8BIT)
   x *= 2 - b * x;               /* here x*a==1 mod 2**16 */
 #endif
 #if defined(MP_64BIT) || !(defined(MP_8BIT) || defined(MP_16BIT))
   x *= 2 - b * x;               /* here x*a==1 mod 2**32 */
 #endif
 #ifdef MP_64BIT
   x *= 2 - b * x;               /* here x*a==1 mod 2**64 */
 #endif
 
   /* rho = -1/m mod b */
   *rho = (((mp_word)1 << ((mp_word) DIGIT_BIT)) - x) & MP_MASK;
 
   return MP_OKAY;
 }

/* determines the setup value */
int mp_reduce_2k_setup(mp_int *a, mp_digit *d)
{
   int res, p;
   mp_int tmp;
   
   if ((res = mp_init(&tmp)) != MP_OKAY) {
      return res;
   }
   
   p = mp_count_bits(a);
   if ((res = mp_2expt(&tmp, p)) != MP_OKAY) {
      mp_clear(&tmp);
      return res;
   }
   
   if ((res = s_mp_sub(&tmp, a, &tmp)) != MP_OKAY) {
      mp_clear(&tmp);
      return res;
   }
   
   *d = tmp.dp[0];
   mp_clear(&tmp);
   return MP_OKAY;
}
 int
 mp_montgomery_reduce (mp_int * x, mp_int * n, mp_digit rho)
 {
   int     ix, res, digs;
   mp_digit mu;
 
   /* can the fast reduction [comba] method be used?
    *
    * Note that unlike in mul you're safely allowed *less*
    * than the available columns [255 per default] since carries
    * are fixed up in the inner loop.
    */
   digs = n->used * 2 + 1;
   if ((digs < MP_WARRAY) &&
       n->used <
       (1 << ((CHAR_BIT * sizeof (mp_word)) - (2 * DIGIT_BIT)))) {
     return fast_mp_montgomery_reduce (x, n, rho);
   }
 
   /* grow the input as required */
   if (x->alloc < digs) {
     if ((res = mp_grow (x, digs)) != MP_OKAY) {
       return res;
     }
   }
   x->used = digs;
 
   for (ix = 0; ix < n->used; ix++) {
     /* mu = ai * rho mod b
      *
      * The value of rho must be precalculated via
      * montgomery_setup() such that
      * it equals -1/n0 mod b this allows the
      * following inner loop to reduce the
      * input one digit at a time
      */
     mu = (mp_digit) (((mp_word)x->dp[ix]) * ((mp_word)rho) & MP_MASK);
 
     /* a = a + mu * m * b**i */
     {
       register int iy;
       register mp_digit *tmpn, *tmpx, u;
       register mp_word r;
 
       /* alias for digits of the modulus */
       tmpn = n->dp;
 
       /* alias for the digits of x [the input] */
       tmpx = x->dp + ix;
 
       /* set the carry to zero */
       u = 0;
 
       /* Multiply and add in place */
       for (iy = 0; iy < n->used; iy++) {
         /* compute product and sum */
         r       = ((mp_word)mu) * ((mp_word)*tmpn++) +
                   ((mp_word) u) + ((mp_word) * tmpx);
 
         /* get carry */
         u       = (mp_digit)(r >> ((mp_word) DIGIT_BIT));
 
         /* fix digit */
         *tmpx++ = (mp_digit)(r & ((mp_word) MP_MASK));
       }
       /* At this point the ix'th digit of x should be zero */
 
 
       /* propagate carries upwards as required*/
       while (u) {
         *tmpx   += u;
         u        = *tmpx >> DIGIT_BIT;
         *tmpx++ &= MP_MASK;
       }
     }
   }
 
   /* at this point the n.used'th least
    * significant digits of x are all zero
    * which means we can shift x to the
    * right by n.used digits and the
    * residue is unchanged.
    */
 
   /* x = x/b**n.used */
   mp_clamp(x);
   mp_rshd (x, n->used);
 
   /* if x >= n then x = x - n */
   if (mp_cmp_mag (x, n) != MP_LT) {
     return s_mp_sub (x, n, x);
   }
 
   return MP_OKAY;
 }

int fast_mp_montgomery_reduce (mp_int * x, mp_int * n, mp_digit rho)
{
  int     ix, res, olduse;
  mp_word W[MP_WARRAY];

  /* get old used count */
  olduse = x->used;

  /* grow a as required */
  if (x->alloc < n->used + 1) {
    if ((res = mp_grow (x, n->used + 1)) != MP_OKAY) {
      return res;
    }
  }

  /* first we have to get the digits of the input into
 *    * an array of double precision words W[...]
 *       */
  {
    register mp_word *_W;
    register mp_digit *tmpx;

    /* alias for the W[] array */
    _W   = W;

    /* alias for the digits of  x*/
    tmpx = x->dp;

    /* copy the digits of a into W[0..a->used-1] */
    for (ix = 0; ix < x->used; ix++) {
      *_W++ = *tmpx++;
    }

    /* zero the high words of W[a->used..m->used*2] */
    for (; ix < n->used * 2 + 1; ix++) {
      *_W++ = 0;
    }
  }

  /* now we proceed to zero successive digits
 *    * from the least significant upwards
 *       */
  for (ix = 0; ix < n->used; ix++) {
    /* mu = ai * m' mod b
 *      *
 *           * We avoid a double precision multiplication (which isn't required)
 *                * by casting the value down to a mp_digit.  Note this requires
 *                     * that W[ix-1] have  the carry cleared (see after the inner loop)
 *                          */
    register mp_digit mu;
    mu = (mp_digit) (((W[ix] & MP_MASK) * rho) & MP_MASK);

    /* a = a + mu * m * b**i
 *      *
 *           * This is computed in place and on the fly.  The multiplication
 *                * by b**i is handled by offseting which columns the results
 *                     * are added to.
 *                          *
 *                               * Note the comba method normally doesn't handle carries in the
 *                                    * inner loop In this case we fix the carry from the previous
 *                                         * column since the Montgomery reduction requires digits of the
 *                                              * result (so far) [see above] to work.  This is
 *                                                   * handled by fixing up one carry after the inner loop.  The
 *                                                        * carry fixups are done in order so after these loops the
 *                                                             * first m->used words of W[] have the carries fixed
 *                                                                  */
    {
      register int iy;
      register mp_digit *tmpn;
      register mp_word *_W;

      /* alias for the digits of the modulus */
      tmpn = n->dp;

      /* Alias for the columns set by an offset of ix */
      _W = W + ix;

      /* inner loop */
      for (iy = 0; iy < n->used; iy++) {
          *_W++ += ((mp_word)mu) * ((mp_word)*tmpn++);
      }
    }

    /* now fix carry for next digit, W[ix+1] */
    W[ix + 1] += W[ix] >> ((mp_word) DIGIT_BIT);
  }
  /* now we have to propagate the carries and
 *    * shift the words downward [all those least
 *       * significant digits we zeroed].
 *          */
  {
    register mp_digit *tmpx;
    register mp_word *_W, *_W1;

    /* nox fix rest of carries */

    /* alias for current word */
    _W1 = W + ix;

    /* alias for next word, where the carry goes */
    _W = W + ++ix;

    for (; ix <= n->used * 2 + 1; ix++) {
      *_W++ += *_W1++ >> ((mp_word) DIGIT_BIT);
    }

    /* copy out, A = A/b**n
 *      *
 *           * The result is A/b**n but instead of converting from an
 *                * array of mp_word to mp_digit than calling mp_rshd
 *                     * we just copy them in the right order
 *                          */

    /* alias for destination word */
    tmpx = x->dp;

    /* alias for shifted double precision result */
    _W = W + n->used;

    for (ix = 0; ix < n->used + 1; ix++) {
      *tmpx++ = (mp_digit)(*_W++ & ((mp_word) MP_MASK));
    }

    /* zero oldused digits, if the input a was larger than
 *      * m->used+1 we'll have to clear the digits
 *           */
    for (; ix < olduse; ix++) {
      *tmpx++ = 0;
    }
  }

  /* set the max used and clamp */
  x->used = n->used + 1;
  mp_clamp (x);

  /* if A >= m then A = A - m */
  if (mp_cmp_mag (x, n) != MP_LT) {
    return s_mp_sub (x, n, x);
  }
  return MP_OKAY;
}

/* determines the setup value */
int mp_reduce_2k_setup_l(mp_int *a, mp_int *d)
{
   int    res;
   mp_int tmp;
   
   if ((res = mp_init(&tmp)) != MP_OKAY) {
      return res;
   }
   
   if ((res = mp_2expt(&tmp, mp_count_bits(a))) != MP_OKAY) {
      goto ERR;
   }
   
   if ((res = s_mp_sub(&tmp, a, d)) != MP_OKAY) {
      goto ERR;
   }
   
ERR:
   mp_clear(&tmp);
   return res;
}


/* low level squaring, b = a*a, HAC pp.596-597, Algorithm 14.16 */
int s_mp_sqr (mp_int * a, mp_int * b)
{
  mp_int  t;
  int     res, ix, iy, pa;
  mp_word r;
  mp_digit u, tmpx, *tmpt;

  pa = a->used;
  if ((res = mp_init_size (&t, 2*pa + 1)) != MP_OKAY) {
    return res;
  }

  /* default used is maximum possible size */
  t.used = 2*pa + 1;

  for (ix = 0; ix < pa; ix++) {
    /* first calculate the digit at 2*ix */
    /* calculate double precision result */
    r = ((mp_word) t.dp[2*ix]) +
        ((mp_word)a->dp[ix])*((mp_word)a->dp[ix]);

    /* store lower part in result */
    t.dp[ix+ix] = (mp_digit) (r & ((mp_word) MP_MASK));

    /* get the carry */
    u           = (mp_digit)(r >> ((mp_word) DIGIT_BIT));

    /* left hand side of A[ix] * A[iy] */
    tmpx        = a->dp[ix];

    /* alias for where to store the results */
    tmpt        = t.dp + (2*ix + 1);
    
    for (iy = ix + 1; iy < pa; iy++) {
      /* first calculate the product */
      r       = ((mp_word)tmpx) * ((mp_word)a->dp[iy]);

      /* now calculate the double precision result, note we use
       * addition instead of *2 since it's easier to optimize
       */
      r       = ((mp_word) *tmpt) + r + r + ((mp_word) u);

      /* store lower part */
      *tmpt++ = (mp_digit) (r & ((mp_word) MP_MASK));

      /* get carry */
      u       = (mp_digit)(r >> ((mp_word) DIGIT_BIT));
    }
    /* propagate upwards */
    while (u != ((mp_digit) 0)) {
      r       = ((mp_word) *tmpt) + ((mp_word) u);
      *tmpt++ = (mp_digit) (r & ((mp_word) MP_MASK));
      u       = (mp_digit)(r >> ((mp_word) DIGIT_BIT));
    }
  }

  mp_clamp (&t);
  mp_exch (&t, b);
  mp_clear (&t);
  return MP_OKAY;
}

/* reduces a modulo n where n is of the form 2**p - d 
   This differs from reduce_2k since "d" can be larger
   than a single digit.
*/
int mp_reduce_2k_l(mp_int *a, mp_int *n, mp_int *d)
{
   mp_int q;
   int    p, res;
   
   if ((res = mp_init(&q)) != MP_OKAY) {
      return res;
   }
   
   p = mp_count_bits(n);    
top:
   /* q = a/2**p, a = a mod 2**p */
   if ((res = mp_div_2d(a, p, &q, a)) != MP_OKAY) {
      goto ERR;
   }
   
   /* q = q * d */
   if ((res = mp_mul(&q, d, &q)) != MP_OKAY) { 
      goto ERR;
   }
   
   /* a = a + q */
   if ((res = s_mp_add(a, &q, a)) != MP_OKAY) {
      goto ERR;
   }
   
   if (mp_cmp_mag(a, n) != MP_LT) {
      s_mp_sub(a, n, a);
      goto top;
   }
   
ERR:
   mp_clear(&q);
   return res;
}



/* reduces x mod m, assumes 0 < x < m**2, mu is 
 * precomputed via mp_reduce_setup.
 * From HAC pp.604 Algorithm 14.42
 */
int mp_reduce (mp_int * x, mp_int * m, mp_int * mu)
{
  mp_int  q;
  int     res, um = m->used;

  /* q = x */
  if ((res = mp_init_copy (&q, x)) != MP_OKAY) {
    return res;
  }

  /* q1 = x / b**(k-1)  */
  mp_rshd (&q, um - 1);         

  /* according to HAC this optimization is ok */
  if (((unsigned long) um) > (((mp_digit)1) << (DIGIT_BIT - 1))) {
    if ((res = mp_mul (&q, mu, &q)) != MP_OKAY) {
      goto CLEANUP;
    }
  } else {
#ifdef BN_S_MP_MUL_HIGH_DIGS_C
    if ((res = s_mp_mul_high_digs (&q, mu, &q, um)) != MP_OKAY) {
      goto CLEANUP;
    }
#elif defined(BN_FAST_S_MP_MUL_HIGH_DIGS_C)
    if ((res = fast_s_mp_mul_high_digs (&q, mu, &q, um)) != MP_OKAY) {
      goto CLEANUP;
    }
#else 
    { 
      res = MP_VAL;
      goto CLEANUP;
    }
#endif
  }

  /* q3 = q2 / b**(k+1) */
  mp_rshd (&q, um + 1);         

  /* x = x mod b**(k+1), quick (no division) */
  if ((res = mp_mod_2d (x, DIGIT_BIT * (um + 1), x)) != MP_OKAY) {
    goto CLEANUP;
  }

  /* q = q * m mod b**(k+1), quick (no division) */
  if ((res = s_mp_mul_digs (&q, m, &q, um + 1)) != MP_OKAY) {
    goto CLEANUP;
  }

  /* x = x - q */
  if ((res = mp_sub (x, &q, x)) != MP_OKAY) {
    goto CLEANUP;
  }

  /* If x < 0, add b**(k+1) to it */
  if (mp_cmp_d (x, 0) == MP_LT) {
    mp_set (&q, 1);
    if ((res = mp_lshd (&q, um + 1)) != MP_OKAY)
      goto CLEANUP;
    if ((res = mp_add (x, &q, x)) != MP_OKAY)
      goto CLEANUP;
  }

  /* Back off if it's too big */
  while (mp_cmp (x, m) != MP_LT) {
    if ((res = s_mp_sub (x, m, x)) != MP_OKAY) {
      goto CLEANUP;
    }
  }
  
CLEANUP:
  mp_clear (&q);

  return res;
}

int
mp_toom_sqr(mp_int *a, mp_int *b)
{
    mp_int w0, w1, w2, w3, w4, tmp1, a0, a1, a2;
    int res, B;

    /* init temps */
    if ((res = mp_init_multi(&w0, &w1, &w2, &w3, &w4, &a0, &a1, &a2, &tmp1, NULL)) != MP_OKAY) {
       return res;
    }

    /* B */
    B = a->used / 3;

    /* a = a2 * B**2 + a1 * B + a0 */
    if ((res = mp_mod_2d(a, DIGIT_BIT * B, &a0)) != MP_OKAY) {
       goto ERR;
    }

    if ((res = mp_copy(a, &a1)) != MP_OKAY) {
       goto ERR;
    }
    mp_rshd(&a1, B);
    mp_mod_2d(&a1, DIGIT_BIT * B, &a1);

    if ((res = mp_copy(a, &a2)) != MP_OKAY) {
       goto ERR;
    }
    mp_rshd(&a2, B*2);

    /* w0 = a0*a0 */
    if ((res = mp_sqr(&a0, &w0)) != MP_OKAY) {
       goto ERR;
    }

    /* w4 = a2 * a2 */
    if ((res = mp_sqr(&a2, &w4)) != MP_OKAY) {
       goto ERR;
    }

    /* w1 = (a2 + 2(a1 + 2a0))**2 */
    if ((res = mp_mul_2(&a0, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul_2(&tmp1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a2, &tmp1)) != MP_OKAY) {
       goto ERR;
    }

    if ((res = mp_sqr(&tmp1, &w1)) != MP_OKAY) {
       goto ERR;
    }

    /* w3 = (a0 + 2(a1 + 2a2))**2 */
    if ((res = mp_mul_2(&a2, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul_2(&tmp1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a0, &tmp1)) != MP_OKAY) {
       goto ERR;
    }

    if ((res = mp_sqr(&tmp1, &w3)) != MP_OKAY) {
       goto ERR;
    }


    /* w2 = (a2 + a1 + a0)**2 */
    if ((res = mp_add(&a2, &a1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a0, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_sqr(&tmp1, &w2)) != MP_OKAY) {
       goto ERR;
    }

    /* now solve the matrix

       0  0  0  0  1
       1  2  4  8  16
       1  1  1  1  1
       16 8  4  2  1
       1  0  0  0  0

       using 12 subtractions, 4 shifts, 2 small divisions and 1 small multiplication.
     */

     /* r1 - r4 */
     if ((res = mp_sub(&w1, &w4, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - r0 */
     if ((res = mp_sub(&w3, &w0, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r1/2 */
     if ((res = mp_div_2(&w1, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3/2 */
     if ((res = mp_div_2(&w3, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r2 - r0 - r4 */
     if ((res = mp_sub(&w2, &w0, &w2)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w2, &w4, &w2)) != MP_OKAY) {
        goto ERR;
     }
     /* r1 - r2 */
     if ((res = mp_sub(&w1, &w2, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - r2 */
     if ((res = mp_sub(&w3, &w2, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r1 - 8r0 */
     if ((res = mp_mul_2d(&w0, 3, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w1, &tmp1, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - 8r4 */
     if ((res = mp_mul_2d(&w4, 3, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w3, &tmp1, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* 3r2 - r1 - r3 */
     if ((res = mp_mul_d(&w2, 3, &w2)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w2, &w1, &w2)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w2, &w3, &w2)) != MP_OKAY) {
        goto ERR;
     }
     /* r1 - r2 */
     if ((res = mp_sub(&w1, &w2, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - r2 */
     if ((res = mp_sub(&w3, &w2, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r1/3 */
     if ((res = mp_div_3(&w1, &w1, NULL)) != MP_OKAY) {
        goto ERR;
     }
     /* r3/3 */
     if ((res = mp_div_3(&w3, &w3, NULL)) != MP_OKAY) {
        goto ERR;
     }

     /* at this point shift W[n] by B*n */
     if ((res = mp_lshd(&w1, 1*B)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_lshd(&w2, 2*B)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_lshd(&w3, 3*B)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_lshd(&w4, 4*B)) != MP_OKAY) {
        goto ERR;
     }

     if ((res = mp_add(&w0, &w1, b)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_add(&w2, &w3, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_add(&w4, &tmp1, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_add(&tmp1, b, b)) != MP_OKAY) {
        goto ERR;
     }

ERR:
     mp_clear_multi(&w0, &w1, &w2, &w3, &w4, &a0, &a1, &a2, &tmp1, NULL);
     return res;
}




int mp_karatsuba_sqr (mp_int * a, mp_int * b)
{
  mp_int  x0, x1, t1, t2, x0x0, x1x1;
  int     B, err;

  err = MP_MEM;

  /* min # of digits */
  B = a->used;

  /* now divide in two */
  B = B >> 1;

  /* init copy all the temps */
  if (mp_init_size (&x0, B) != MP_OKAY)
    goto ERR;
  if (mp_init_size (&x1, a->used - B) != MP_OKAY)
    goto X0;

  /* init temps */
  if (mp_init_size (&t1, a->used * 2) != MP_OKAY)
    goto X1;
  if (mp_init_size (&t2, a->used * 2) != MP_OKAY)
    goto T1;
  if (mp_init_size (&x0x0, B * 2) != MP_OKAY)
    goto T2;
  if (mp_init_size (&x1x1, (a->used - B) * 2) != MP_OKAY)
    goto X0X0;

  {
    register int x;
    register mp_digit *dst, *src;

    src = a->dp;

    /* now shift the digits */
    dst = x0.dp;
    for (x = 0; x < B; x++) {
      *dst++ = *src++;
    }

    dst = x1.dp;
    for (x = B; x < a->used; x++) {
      *dst++ = *src++;
    }
  }

  x0.used = B;
  x1.used = a->used - B;

  mp_clamp (&x0);

  /* now calc the products x0*x0 and x1*x1 */
  if (mp_sqr (&x0, &x0x0) != MP_OKAY)
    goto X1X1;           /* x0x0 = x0*x0 */
  if (mp_sqr (&x1, &x1x1) != MP_OKAY)
    goto X1X1;           /* x1x1 = x1*x1 */

  /* now calc (x1-x0)**2 */
  if (mp_sub (&x1, &x0, &t1) != MP_OKAY)
    goto X1X1;           /* t1 = x1 - x0 */
  if (mp_sqr (&t1, &t1) != MP_OKAY)
    goto X1X1;           /* t1 = (x1 - x0) * (x1 - x0) */

  /* add x0y0 */
  if (s_mp_add (&x0x0, &x1x1, &t2) != MP_OKAY)
    goto X1X1;           /* t2 = x0x0 + x1x1 */
  if (mp_sub (&t2, &t1, &t1) != MP_OKAY)
    goto X1X1;           /* t1 = x0x0 + x1x1 - (x1-x0)*(x1-x0) */

  /* shift by B */
  if (mp_lshd (&t1, B) != MP_OKAY)
    goto X1X1;           /* t1 = (x0x0 + x1x1 - (x1-x0)*(x1-x0))<<B */
  if (mp_lshd (&x1x1, B * 2) != MP_OKAY)
    goto X1X1;           /* x1x1 = x1x1 << 2*B */

  if (mp_add (&x0x0, &t1, &t1) != MP_OKAY)
    goto X1X1;           /* t1 = x0x0 + t1 */
  if (mp_add (&t1, &x1x1, b) != MP_OKAY)
    goto X1X1;           /* t1 = x0x0 + t1 + x1x1 */

  err = MP_OKAY;

X1X1:mp_clear (&x1x1);
X0X0:mp_clear (&x0x0);
T2:mp_clear (&t2);
T1:mp_clear (&t1);
X1:mp_clear (&x1);
X0:mp_clear (&x0);
ERR:
  return err;
}

int fast_s_mp_sqr (mp_int * a, mp_int * b)
{
  int       olduse, res, pa, ix, iz;
  mp_digit   W[MP_WARRAY], *tmpx;
  mp_word   W1;

  /* grow the destination as required */
  pa = a->used + a->used;
  if (b->alloc < pa) {
    if ((res = mp_grow (b, pa)) != MP_OKAY) {
      return res;
    }
  }

  /* number of output digits to produce */
  W1 = 0;
  for (ix = 0; ix < pa; ix++) { 
      int      tx, ty, iy;
      mp_word  _W;
      mp_digit *tmpy;

      /* clear counter */
      _W = 0;

      /* get offsets into the two bignums */
      ty = MIN(a->used-1, ix);
      tx = ix - ty;

      /* setup temp aliases */
      tmpx = a->dp + tx;
      tmpy = a->dp + ty;

      /* this is the number of times the loop will iterrate, essentially
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MIN(a->used-tx, ty+1);

      /* now for squaring tx can never equal ty 
       * we halve the distance since they approach at a rate of 2x
       * and we have to round because odd cases need to be executed
       */
      iy = MIN(iy, (ty-tx+1)>>1);

      /* execute loop */
      for (iz = 0; iz < iy; iz++) {
         _W += ((mp_word)*tmpx++)*((mp_word)*tmpy--);
      }

      /* double the inner product and add carry */
      _W = _W + _W + W1;

      /* even columns have the square term in them */
      if ((ix&1) == 0) {
         _W += ((mp_word)a->dp[ix>>1])*((mp_word)a->dp[ix>>1]);
      }

      /* store it */
      W[ix] = (mp_digit)(_W & MP_MASK);

      /* make next carry */
      W1 = _W >> ((mp_word)DIGIT_BIT);
  }

  /* setup dest */
  olduse  = b->used;
  b->used = a->used+a->used;

  {
    mp_digit *tmpb;
    tmpb = b->dp;
    for (ix = 0; ix < pa; ix++) {
      *tmpb++ = W[ix] & MP_MASK;
    }

    /* clear unused digits [that existed in the old copy of c] */
    for (; ix < olduse; ix++) {
      *tmpb++ = 0;
    }
  }
  mp_clamp (b);
  return MP_OKAY;
}

int mp_toom_mul(mp_int *a, mp_int *b, mp_int *c)
{
    mp_int w0, w1, w2, w3, w4, tmp1, tmp2, a0, a1, a2, b0, b1, b2;
    int res, B;
        
    /* init temps */
    if ((res = mp_init_multi(&w0, &w1, &w2, &w3, &w4, 
                             &a0, &a1, &a2, &b0, &b1, 
                             &b2, &tmp1, &tmp2, NULL)) != MP_OKAY) {
       return res;
    }
    
    /* B */
    B = MIN(a->used, b->used) / 3;
    
    /* a = a2 * B**2 + a1 * B + a0 */
    if ((res = mp_mod_2d(a, DIGIT_BIT * B, &a0)) != MP_OKAY) {
       goto ERR;
    }

    if ((res = mp_copy(a, &a1)) != MP_OKAY) {
       goto ERR;
    }
    mp_rshd(&a1, B);
    mp_mod_2d(&a1, DIGIT_BIT * B, &a1);

    if ((res = mp_copy(a, &a2)) != MP_OKAY) {
       goto ERR;
    }
    mp_rshd(&a2, B*2);
    
    /* b = b2 * B**2 + b1 * B + b0 */
    if ((res = mp_mod_2d(b, DIGIT_BIT * B, &b0)) != MP_OKAY) {
       goto ERR;
    }

    if ((res = mp_copy(b, &b1)) != MP_OKAY) {
       goto ERR;
    }
    mp_rshd(&b1, B);
    mp_mod_2d(&b1, DIGIT_BIT * B, &b1);

    if ((res = mp_copy(b, &b2)) != MP_OKAY) {
       goto ERR;
    }
    mp_rshd(&b2, B*2);
    
    /* w0 = a0*b0 */
    if ((res = mp_mul(&a0, &b0, &w0)) != MP_OKAY) {
       goto ERR;
    }
    
    /* w4 = a2 * b2 */
    if ((res = mp_mul(&a2, &b2, &w4)) != MP_OKAY) {
       goto ERR;
    }
    
    /* w1 = (a2 + 2(a1 + 2a0))(b2 + 2(b1 + 2b0)) */
    if ((res = mp_mul_2(&a0, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul_2(&tmp1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a2, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    
    if ((res = mp_mul_2(&b0, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp2, &b1, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul_2(&tmp2, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp2, &b2, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    
    if ((res = mp_mul(&tmp1, &tmp2, &w1)) != MP_OKAY) {
       goto ERR;
    }
    
    /* w3 = (a0 + 2(a1 + 2a2))(b0 + 2(b1 + 2b2)) */
    if ((res = mp_mul_2(&a2, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul_2(&tmp1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a0, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    
    if ((res = mp_mul_2(&b2, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp2, &b1, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul_2(&tmp2, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp2, &b0, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    
    if ((res = mp_mul(&tmp1, &tmp2, &w3)) != MP_OKAY) {
       goto ERR;
    }
    

    /* w2 = (a2 + a1 + a0)(b2 + b1 + b0) */
    if ((res = mp_add(&a2, &a1, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp1, &a0, &tmp1)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&b2, &b1, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_add(&tmp2, &b0, &tmp2)) != MP_OKAY) {
       goto ERR;
    }
    if ((res = mp_mul(&tmp1, &tmp2, &w2)) != MP_OKAY) {
       goto ERR;
    }
    
    /* now solve the matrix 
    
       0  0  0  0  1
       1  2  4  8  16
       1  1  1  1  1
       16 8  4  2  1
       1  0  0  0  0
       
       using 12 subtractions, 4 shifts, 
              2 small divisions and 1 small multiplication 
     */
     
     /* r1 - r4 */
     if ((res = mp_sub(&w1, &w4, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - r0 */
     if ((res = mp_sub(&w3, &w0, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r1/2 */
     if ((res = mp_div_2(&w1, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3/2 */
     if ((res = mp_div_2(&w3, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r2 - r0 - r4 */
     if ((res = mp_sub(&w2, &w0, &w2)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w2, &w4, &w2)) != MP_OKAY) {
        goto ERR;
     }
     /* r1 - r2 */
     if ((res = mp_sub(&w1, &w2, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - r2 */
     if ((res = mp_sub(&w3, &w2, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r1 - 8r0 */
     if ((res = mp_mul_2d(&w0, 3, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w1, &tmp1, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - 8r4 */
     if ((res = mp_mul_2d(&w4, 3, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w3, &tmp1, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* 3r2 - r1 - r3 */
     if ((res = mp_mul_d(&w2, 3, &w2)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w2, &w1, &w2)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_sub(&w2, &w3, &w2)) != MP_OKAY) {
        goto ERR;
     }
     /* r1 - r2 */
     if ((res = mp_sub(&w1, &w2, &w1)) != MP_OKAY) {
        goto ERR;
     }
     /* r3 - r2 */
     if ((res = mp_sub(&w3, &w2, &w3)) != MP_OKAY) {
        goto ERR;
     }
     /* r1/3 */
     if ((res = mp_div_3(&w1, &w1, NULL)) != MP_OKAY) {
        goto ERR;
     }
     /* r3/3 */
     if ((res = mp_div_3(&w3, &w3, NULL)) != MP_OKAY) {
        goto ERR;
     }
     
     /* at this point shift W[n] by B*n */
     if ((res = mp_lshd(&w1, 1*B)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_lshd(&w2, 2*B)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_lshd(&w3, 3*B)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_lshd(&w4, 4*B)) != MP_OKAY) {
        goto ERR;
     }     
     
     if ((res = mp_add(&w0, &w1, c)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_add(&w2, &w3, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_add(&w4, &tmp1, &tmp1)) != MP_OKAY) {
        goto ERR;
     }
     if ((res = mp_add(&tmp1, c, c)) != MP_OKAY) {
        goto ERR;
     }     
     
ERR:
     mp_clear_multi(&w0, &w1, &w2, &w3, &w4, 
                    &a0, &a1, &a2, &b0, &b1, 
                    &b2, &tmp1, &tmp2, NULL);
     return res;
}

int fast_s_mp_mul_digs (mp_int * a, mp_int * b, mp_int * c, int digs)
{
  int     olduse, res, pa, ix, iz;
  mp_digit W[MP_WARRAY];
  register mp_word  _W;

  /* grow the destination as required */
  if (c->alloc < digs) {
    if ((res = mp_grow (c, digs)) != MP_OKAY) {
      return res;
    }
  }

  /* number of output digits to produce */
  pa = MIN(digs, a->used + b->used);

  /* clear the carry */
  _W = 0;
  for (ix = 0; ix < pa; ix++) { 
      int      tx, ty;
      int      iy;
      mp_digit *tmpx, *tmpy;

      /* get offsets into the two bignums */
      ty = MIN(b->used-1, ix);
      tx = ix - ty;

      /* setup temp aliases */
      tmpx = a->dp + tx;
      tmpy = b->dp + ty;

      /* this is the number of times the loop will iterrate, essentially 
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MIN(a->used-tx, ty+1);

      /* execute loop */
      for (iz = 0; iz < iy; ++iz) {
         _W += ((mp_word)*tmpx++)*((mp_word)*tmpy--);
      }

      /* store term */
      W[ix] = ((mp_digit)_W) & MP_MASK;

      /* make next carry */
      _W = _W >> ((mp_word)DIGIT_BIT);
  }

  /* store final carry */
  W[ix] = (mp_digit)(_W & MP_MASK);

  /* setup dest */
  olduse  = c->used;
  c->used = pa;

  {
    register mp_digit *tmpc;
    tmpc = c->dp;
    for (ix = 0; ix < pa+1; ix++) {
      /* now extract the previous digit [below the carry] */
      *tmpc++ = W[ix];
    }

    /* clear unused digits [that existed in the old copy of c] */
    for (; ix < olduse; ix++) {
      *tmpc++ = 0;
    }
  }
  mp_clamp (c);
  return MP_OKAY;
}

#define s_mp_mul(a, b, c) s_mp_mul_digs(a, b, c, (a)->used + (b)->used + 1)


int s_mp_mul_digs (mp_int * a, mp_int * b, mp_int * c, int digs)
{
  mp_int  t;
  int     res, pa, pb, ix, iy;
  mp_digit u;
  mp_word r;
  mp_digit tmpx, *tmpt, *tmpy;

  /* can we use the fast multiplier? */
  if (((digs) < MP_WARRAY) &&
      MIN (a->used, b->used) < 
          (1 << ((CHAR_BIT * sizeof (mp_word)) - (2 * DIGIT_BIT)))) {
    return fast_s_mp_mul_digs (a, b, c, digs);
  }

  if ((res = mp_init_size (&t, digs)) != MP_OKAY) {
    return res;
  }
  t.used = digs;

  /* compute the digits of the product directly */
  pa = a->used;
  for (ix = 0; ix < pa; ix++) {
    /* set the carry to zero */
    u = 0;

    /* limit ourselves to making digs digits of output */
    pb = MIN (b->used, digs - ix);

    /* setup some aliases */
    /* copy of the digit from a used within the nested loop */
    tmpx = a->dp[ix];
    
    /* an alias for the destination shifted ix places */
    tmpt = t.dp + ix;
    
    /* an alias for the digits of b */
    tmpy = b->dp;

    /* compute the columns of the output and propagate the carry */
    for (iy = 0; iy < pb; iy++) {
      /* compute the column as a mp_word */
      r       = ((mp_word)*tmpt) +
                ((mp_word)tmpx) * ((mp_word)*tmpy++) +
                ((mp_word) u);

      /* the new column is the lower part of the result */
      *tmpt++ = (mp_digit) (r & ((mp_word) MP_MASK));

      /* get the carry word from the result */
      u       = (mp_digit) (r >> ((mp_word) DIGIT_BIT));
    }
    /* set carry if it is placed below digs */
    if (ix + iy < digs) {
      *tmpt = u;
    }
  }

  mp_clamp (&t);
  mp_exch (&t, c);

  mp_clear (&t);
  return MP_OKAY;
}



int fast_mp_invmod (mp_int * a, mp_int * b, mp_int * c)
{
  mp_int  x, y, u, v, B, D;
  int     res, neg;

  /* 2. [modified] b must be odd   */
  if (mp_iseven (b) == 1) {
    return MP_VAL;
  }

  /* init all our temps */
  if ((res = mp_init_multi(&x, &y, &u, &v, &B, &D, NULL)) != MP_OKAY) {
     return res;
  }

  /* x == modulus, y == value to invert */
  if ((res = mp_copy (b, &x)) != MP_OKAY) {
    goto LBL_ERR;
  }

  /* we need y = |a| */
  if ((res = mp_mod (a, b, &y)) != MP_OKAY) {
    goto LBL_ERR;
  }

  /* 3. u=x, v=y, A=1, B=0, C=0,D=1 */
  if ((res = mp_copy (&x, &u)) != MP_OKAY) {
    goto LBL_ERR;
  }
  if ((res = mp_copy (&y, &v)) != MP_OKAY) {
    goto LBL_ERR;
  }
  mp_set (&D, 1);

top:
  /* 4.  while u is even do */
  while (mp_iseven (&u) == 1) {
    /* 4.1 u = u/2 */
    if ((res = mp_div_2 (&u, &u)) != MP_OKAY) {
      goto LBL_ERR;
    }
    /* 4.2 if B is odd then */
    if (mp_isodd (&B) == 1) {
      if ((res = mp_sub (&B, &x, &B)) != MP_OKAY) {
        goto LBL_ERR;
      }
    }
    /* B = B/2 */
    if ((res = mp_div_2 (&B, &B)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }

  /* 5.  while v is even do */
  while (mp_iseven (&v) == 1) {
    /* 5.1 v = v/2 */
    if ((res = mp_div_2 (&v, &v)) != MP_OKAY) {
      goto LBL_ERR;
    }
    /* 5.2 if D is odd then */
    if (mp_isodd (&D) == 1) {
      /* D = (D-x)/2 */
      if ((res = mp_sub (&D, &x, &D)) != MP_OKAY) {
        goto LBL_ERR;
      }
    }
    /* D = D/2 */
    if ((res = mp_div_2 (&D, &D)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }

  /* 6.  if u >= v then */
  if (mp_cmp (&u, &v) != MP_LT) {
    /* u = u - v, B = B - D */
    if ((res = mp_sub (&u, &v, &u)) != MP_OKAY) {
      goto LBL_ERR;
    }

    if ((res = mp_sub (&B, &D, &B)) != MP_OKAY) {
      goto LBL_ERR;
    }
  } else {
    /* v - v - u, D = D - B */
    if ((res = mp_sub (&v, &u, &v)) != MP_OKAY) {
      goto LBL_ERR;
    }

    if ((res = mp_sub (&D, &B, &D)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }

  /* if not zero goto step 4 */
  if (mp_iszero (&u) == 0) {
    goto top;
  }

  /* now a = C, b = D, gcd == g*v */

  /* if v != 1 then there is no inverse */
  if (mp_cmp_d (&v, 1) != MP_EQ) {
    res = MP_VAL;
    goto LBL_ERR;
  }

  /* b is now the inverse */
  neg = a->sign;
  while (D.sign == MP_NEG) {
    if ((res = mp_add (&D, b, &D)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }
  mp_exch (&D, c);
  c->sign = neg;
  res = MP_OKAY;

LBL_ERR:mp_clear_multi (&x, &y, &u, &v, &B, &D, NULL);
  return res;
}


int mp_invmod_slow (mp_int * a, mp_int * b, mp_int * c)
{
  mp_int  x, y, u, v, A, B, C, D;
  int     res;

  /* b cannot be negative */
  if (b->sign == MP_NEG || mp_iszero(b) == 1) {
    return MP_VAL;
  }

  /* init temps */
  if ((res = mp_init_multi(&x, &y, &u, &v, 
                           &A, &B, &C, &D, NULL)) != MP_OKAY) {
     return res;
  }

  /* x = a, y = b */
  if ((res = mp_mod(a, b, &x)) != MP_OKAY) {
      goto LBL_ERR;
  }
  if ((res = mp_copy (b, &y)) != MP_OKAY) {
    goto LBL_ERR;
  }

  /* 2. [modified] if x,y are both even then return an error! */
  if (mp_iseven (&x) == 1 && mp_iseven (&y) == 1) {
    res = MP_VAL;
    goto LBL_ERR;
  }

  /* 3. u=x, v=y, A=1, B=0, C=0,D=1 */
  if ((res = mp_copy (&x, &u)) != MP_OKAY) {
    goto LBL_ERR;
  }
  if ((res = mp_copy (&y, &v)) != MP_OKAY) {
    goto LBL_ERR;
  }
  mp_set (&A, 1);
  mp_set (&D, 1);

top:
  /* 4.  while u is even do */
  while (mp_iseven (&u) == 1) {
    /* 4.1 u = u/2 */
    if ((res = mp_div_2 (&u, &u)) != MP_OKAY) {
      goto LBL_ERR;
    }
    /* 4.2 if A or B is odd then */
    if (mp_isodd (&A) == 1 || mp_isodd (&B) == 1) {
      /* A = (A+y)/2, B = (B-x)/2 */
      if ((res = mp_add (&A, &y, &A)) != MP_OKAY) {
         goto LBL_ERR;
      }
      if ((res = mp_sub (&B, &x, &B)) != MP_OKAY) {
         goto LBL_ERR;
      }
    }
    /* A = A/2, B = B/2 */
    if ((res = mp_div_2 (&A, &A)) != MP_OKAY) {
      goto LBL_ERR;
    }
    if ((res = mp_div_2 (&B, &B)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }

  /* 5.  while v is even do */
  while (mp_iseven (&v) == 1) {
    /* 5.1 v = v/2 */
    if ((res = mp_div_2 (&v, &v)) != MP_OKAY) {
      goto LBL_ERR;
    }
    /* 5.2 if C or D is odd then */
    if (mp_isodd (&C) == 1 || mp_isodd (&D) == 1) {
      /* C = (C+y)/2, D = (D-x)/2 */
      if ((res = mp_add (&C, &y, &C)) != MP_OKAY) {
         goto LBL_ERR;
      }
      if ((res = mp_sub (&D, &x, &D)) != MP_OKAY) {
         goto LBL_ERR;
      }
    }
    /* C = C/2, D = D/2 */
    if ((res = mp_div_2 (&C, &C)) != MP_OKAY) {
      goto LBL_ERR;
    }
    if ((res = mp_div_2 (&D, &D)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }

  /* 6.  if u >= v then */
  if (mp_cmp (&u, &v) != MP_LT) {
    /* u = u - v, A = A - C, B = B - D */
    if ((res = mp_sub (&u, &v, &u)) != MP_OKAY) {
      goto LBL_ERR;
    }

    if ((res = mp_sub (&A, &C, &A)) != MP_OKAY) {
      goto LBL_ERR;
    }

    if ((res = mp_sub (&B, &D, &B)) != MP_OKAY) {
      goto LBL_ERR;
    }
  } else {
    /* v - v - u, C = C - A, D = D - B */
    if ((res = mp_sub (&v, &u, &v)) != MP_OKAY) {
      goto LBL_ERR;
    }

    if ((res = mp_sub (&C, &A, &C)) != MP_OKAY) {
      goto LBL_ERR;
    }

    if ((res = mp_sub (&D, &B, &D)) != MP_OKAY) {
      goto LBL_ERR;
    }
  }

  /* if not zero goto step 4 */
  if (mp_iszero (&u) == 0)
    goto top;

  /* now a = C, b = D, gcd == g*v */

  /* if v != 1 then there is no inverse */
  if (mp_cmp_d (&v, 1) != MP_EQ) {
    res = MP_VAL;
    goto LBL_ERR;
  }

  /* if its too low */
  while (mp_cmp_d(&C, 0) == MP_LT) {
      if ((res = mp_add(&C, b, &C)) != MP_OKAY) {
         goto LBL_ERR;
      }
  }
  
  /* too big */
  while (mp_cmp_mag(&C, b) != MP_LT) {
      if ((res = mp_sub(&C, b, &C)) != MP_OKAY) {
         goto LBL_ERR;
      }
  }
  
  /* C is now the inverse */
  mp_exch (&C, c);
  res = MP_OKAY;
LBL_ERR:mp_clear_multi (&x, &y, &u, &v, &A, &B, &C, &D, NULL);
  return res;
}

/* divide by three (based on routine from MPI and the GMP manual) */
int
mp_div_3 (mp_int * a, mp_int *c, mp_digit * d)
{
  mp_int   q;
  mp_word  w, t;
  mp_digit b;
  int      res, ix;
  
  /* b = 2**DIGIT_BIT / 3 */
  b = (((mp_word)1) << ((mp_word)DIGIT_BIT)) / ((mp_word)3);

  if ((res = mp_init_size(&q, a->used)) != MP_OKAY) {
     return res;
  }
  
  q.used = a->used;
  q.sign = a->sign;
  w = 0;
  for (ix = a->used - 1; ix >= 0; ix--) {
     w = (w << ((mp_word)DIGIT_BIT)) | ((mp_word)a->dp[ix]);

     if (w >= 3) {
        /* multiply w by [1/3] */
        t = (w * ((mp_word)b)) >> ((mp_word)DIGIT_BIT);

        /* now subtract 3 * [w/3] from w, to get the remainder */
        w -= t+t+t;

        /* fixup the remainder as required since
         * the optimization is not exact.
         */
        while (w >= 3) {
           t += 1;
           w -= 3;
        }
      } else {
        t = 0;
      }
      q.dp[ix] = (mp_digit)t;
  }

  /* [optional] store the remainder */
  if (d != NULL) {
     *d = (mp_digit)w;
  }

  /* [optional] store the quotient */
  if (c != NULL) {
     mp_clamp(&q);
     mp_exch(&q, c);
  }
  mp_clear(&q);
  
  return res;
}



void mp_clear_multi(mp_int *mp, ...) 
{
    mp_int* next_mp = mp;
    va_list args;
    va_start(args, mp);
    while (next_mp != NULL) {
        mp_clear(next_mp);
        next_mp = va_arg(args, mp_int*);
    }
    va_end(args);
}


/* determines if a number is a valid DR modulus */
int mp_dr_is_modulus(mp_int *a)
{
   int ix;

   /* must be at least two digits */
   if (a->used < 2) {
      return 0;
   }

   /* must be of the form b**k - a [a <= b] so all
    * but the first digit must be equal to -1 (mod b).
    */
   for (ix = 1; ix < a->used; ix++) {
       if (a->dp[ix] != MP_MASK) {
          return 0;
       }
   }
   return 1;
}


/* determines if mp_reduce_2k can be used */
int mp_reduce_is_2k(mp_int *a)
{
   int ix, iy, iw;
   mp_digit iz;
   
   if (a->used == 0) {
      return MP_NO;
   } else if (a->used == 1) {
      return MP_YES;
   } else if (a->used > 1) {
      iy = mp_count_bits(a);
      iz = 1;
      iw = 1;
    
      /* Test every bit from the second digit up, must be 1 */
      for (ix = DIGIT_BIT; ix < iy; ix++) {
          if ((a->dp[iw] & iz) == 0) {
             return MP_NO;
          }
          iz <<= 1;
          if (iz > (mp_digit)MP_MASK) {
             ++iw;
             iz = 1;
          }
      }
   }
   return MP_YES;
}



/* computes Y == G**X mod P, HAC pp.616, Algorithm 14.85
 *
 * Uses a left-to-right k-ary sliding window to compute the modular exponentiation.
 * The value of k changes based on the size of the exponent.
 *
 * Uses Montgomery or Diminished Radix reduction [whichever appropriate]
 */

#ifdef MP_LOW_MEM
   #define TAB_SIZE 32
#else
   #define TAB_SIZE 256
#endif

int mp_exptmod_fast (mp_int * G, mp_int * X, mp_int * P, mp_int * Y, int redmode)
{
  mp_int  M[TAB_SIZE], res;
  mp_digit buf, mp;
  int     err, bitbuf, bitcpy, bitcnt, mode, digidx, x, y, winsize;

  /* use a pointer to the reduction algorithm.  This allows us to use
   * one of many reduction algorithms without modding the guts of
   * the code with if statements everywhere.
   */
  int     (*redux)(mp_int*,mp_int*,mp_digit);

  /* find window size */
  x = mp_count_bits (X);
  if (x <= 7) {
    winsize = 2;
  } else if (x <= 36) {
    winsize = 3;
  } else if (x <= 140) {
    winsize = 4;
  } else if (x <= 450) {
    winsize = 5;
  } else if (x <= 1303) {
    winsize = 6;
  } else if (x <= 3529) {
    winsize = 7;
  } else {
    winsize = 8;
  }

#ifdef MP_LOW_MEM
  if (winsize > 5) {
     winsize = 5;
  }
#endif

  /* init M array */
  /* init first cell */
  if ((err = mp_init(&M[1])) != MP_OKAY) {
     return err;
  }

  /* now init the second half of the array */
  for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
    if ((err = mp_init(&M[x])) != MP_OKAY) {
      for (y = 1<<(winsize-1); y < x; y++) {
        mp_clear (&M[y]);
      }
      mp_clear(&M[1]);
      return err;
    }
  }

  /* determine and setup reduction code */
  if (redmode == 0) {
#ifdef BN_MP_MONTGOMERY_SETUP_C   
     /* now setup montgomery  */
     if ((err = mp_montgomery_setup (P, &mp)) != MP_OKAY) {
        goto LBL_M;
     }
#else
     err = MP_VAL;
     goto LBL_M;
#endif

     /* automatically pick the comba one if available (saves quite a few calls/ifs) */
#ifdef BN_FAST_MP_MONTGOMERY_REDUCE_C
     if (((P->used * 2 + 1) < MP_WARRAY) &&
          P->used < (1 << ((CHAR_BIT * sizeof (mp_word)) - (2 * DIGIT_BIT)))) {
        redux = fast_mp_montgomery_reduce;
     } else 
#endif
     {
#ifdef BN_MP_MONTGOMERY_REDUCE_C
        /* use slower baseline Montgomery method */
        redux = mp_montgomery_reduce;
#else
        err = MP_VAL;
        goto LBL_M;
#endif
     }
  } else if (redmode == 1) {
#if defined(BN_MP_DR_SETUP_C) && defined(BN_MP_DR_REDUCE_C)
     /* setup DR reduction for moduli of the form B**k - b */
     mp_dr_setup(P, &mp);
     redux = mp_dr_reduce;
#else
     err = MP_VAL;
     goto LBL_M;
#endif
  } else {
#if defined(BN_MP_REDUCE_2K_SETUP_C) && defined(BN_MP_REDUCE_2K_C)
     /* setup DR reduction for moduli of the form 2**k - b */
     if ((err = mp_reduce_2k_setup(P, &mp)) != MP_OKAY) {
        goto LBL_M;
     }
     redux = mp_reduce_2k;
#else
     err = MP_VAL;
     goto LBL_M;
#endif
  }

  /* setup result */
  if ((err = mp_init (&res)) != MP_OKAY) {
    goto LBL_M;
  }

  /* create M table
   *

   *
   * The first half of the table is not computed though accept for M[0] and M[1]
   */

  if (redmode == 0) {
#ifdef BN_MP_MONTGOMERY_CALC_NORMALIZATION_C
     /* now we need R mod m */
     if ((err = mp_montgomery_calc_normalization (&res, P)) != MP_OKAY) {
       goto LBL_RES;
     }
#else 
     err = MP_VAL;
     goto LBL_RES;
#endif

     /* now set M[1] to G * R mod m */
     if ((err = mp_mulmod (G, &res, P, &M[1])) != MP_OKAY) {
       goto LBL_RES;
     }
  } else {
     mp_set(&res, 1);
     if ((err = mp_mod(G, P, &M[1])) != MP_OKAY) {
        goto LBL_RES;
     }
  }

  /* compute the value at M[1<<(winsize-1)] by squaring M[1] (winsize-1) times */
  if ((err = mp_copy (&M[1], &M[1 << (winsize - 1)])) != MP_OKAY) {
    goto LBL_RES;
  }

  for (x = 0; x < (winsize - 1); x++) {
    if ((err = mp_sqr (&M[1 << (winsize - 1)], &M[1 << (winsize - 1)])) != MP_OKAY) {
      goto LBL_RES;
    }
    if ((err = redux (&M[1 << (winsize - 1)], P, mp)) != MP_OKAY) {
      goto LBL_RES;
    }
  }

  /* create upper table */
  for (x = (1 << (winsize - 1)) + 1; x < (1 << winsize); x++) {
    if ((err = mp_mul (&M[x - 1], &M[1], &M[x])) != MP_OKAY) {
      goto LBL_RES;
    }
    if ((err = redux (&M[x], P, mp)) != MP_OKAY) {
      goto LBL_RES;
    }
  }

  /* set initial mode and bit cnt */
  mode   = 0;
  bitcnt = 1;
  buf    = 0;
  digidx = X->used - 1;
  bitcpy = 0;
  bitbuf = 0;

  for (;;) {
    /* grab next digit as required */
    if (--bitcnt == 0) {
      /* if digidx == -1 we are out of digits so break */
      if (digidx == -1) {
        break;
      }
      /* read next digit and reset bitcnt */
      buf    = X->dp[digidx--];
      bitcnt = (int)DIGIT_BIT;
    }

    /* grab the next msb from the exponent */
    y     = (mp_digit)(buf >> (DIGIT_BIT - 1)) & 1;
    buf <<= (mp_digit)1;

    /* if the bit is zero and mode == 0 then we ignore it
     * These represent the leading zero bits before the first 1 bit
     * in the exponent.  Technically this opt is not required but it
     * does lower the # of trivial squaring/reductions used
     */
    if (mode == 0 && y == 0) {
      continue;
    }

    /* if the bit is zero and mode == 1 then we square */
    if (mode == 1 && y == 0) {
      if ((err = mp_sqr (&res, &res)) != MP_OKAY) {
        goto LBL_RES;
      }
      if ((err = redux (&res, P, mp)) != MP_OKAY) {
        goto LBL_RES;
      }
      continue;
    }

    /* else we add it to the window */
    bitbuf |= (y << (winsize - ++bitcpy));
    mode    = 2;

    if (bitcpy == winsize) {
      /* ok window is filled so square as required and multiply  */
      /* square first */
      for (x = 0; x < winsize; x++) {
        if ((err = mp_sqr (&res, &res)) != MP_OKAY) {
          goto LBL_RES;
        }
        if ((err = redux (&res, P, mp)) != MP_OKAY) {
          goto LBL_RES;
        }
      }

      /* then multiply */
      if ((err = mp_mul (&res, &M[bitbuf], &res)) != MP_OKAY) {
        goto LBL_RES;
      }
      if ((err = redux (&res, P, mp)) != MP_OKAY) {
        goto LBL_RES;
      }

      /* empty window and reset */
      bitcpy = 0;
      bitbuf = 0;
      mode   = 1;
    }
  }

  /* if bits remain then square/multiply */
  if (mode == 2 && bitcpy > 0) {
    /* square then multiply if the bit is set */
    for (x = 0; x < bitcpy; x++) {
      if ((err = mp_sqr (&res, &res)) != MP_OKAY) {
        goto LBL_RES;
      }
      if ((err = redux (&res, P, mp)) != MP_OKAY) {
        goto LBL_RES;
      }

      /* get next bit of the window */
      bitbuf <<= 1;
      if ((bitbuf & (1 << winsize)) != 0) {
        /* then multiply */
        if ((err = mp_mul (&res, &M[1], &res)) != MP_OKAY) {
          goto LBL_RES;
        }
        if ((err = redux (&res, P, mp)) != MP_OKAY) {
          goto LBL_RES;
        }
      }
    }
  }

  if (redmode == 0) {
     /* fixup result if Montgomery reduction is used
      * recall that any value in a Montgomery system is
      * actually multiplied by R mod n.  So we have
      * to reduce one more time to cancel out the factor
      * of R.
      */
     if ((err = redux(&res, P, mp)) != MP_OKAY) {
       goto LBL_RES;
     }
  }

  /* swap res with Y */
  mp_exch (&res, Y);
  err = MP_OKAY;
LBL_RES:mp_clear (&res);
LBL_M:
  mp_clear(&M[1]);
  for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
    mp_clear (&M[x]);
  }
  return err;
}



#ifdef MP_LOW_MEM
   #define TAB_SIZE 32
#else
   #define TAB_SIZE 256
#endif

int s_mp_exptmod (mp_int * G, mp_int * X, mp_int * P, mp_int * Y, int redmode)
{
  mp_int  M[TAB_SIZE], res, mu;
  mp_digit buf;
  int     err, bitbuf, bitcpy, bitcnt, mode, digidx, x, y, winsize;
  int (*redux)(mp_int*,mp_int*,mp_int*);

  /* find window size */
  x = mp_count_bits (X);
  if (x <= 7) {
    winsize = 2;
  } else if (x <= 36) {
    winsize = 3;
  } else if (x <= 140) {
    winsize = 4;
  } else if (x <= 450) {
    winsize = 5;
  } else if (x <= 1303) {
    winsize = 6;
  } else if (x <= 3529) {
    winsize = 7;
  } else {
    winsize = 8;
  }

#ifdef MP_LOW_MEM
    if (winsize > 5) {
       winsize = 5;
    }
#endif

  /* init M array */
  /* init first cell */
  if ((err = mp_init(&M[1])) != MP_OKAY) {
     return err; 
  }

  /* now init the second half of the array */
  for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
    if ((err = mp_init(&M[x])) != MP_OKAY) {
      for (y = 1<<(winsize-1); y < x; y++) {
        mp_clear (&M[y]);
      }
      mp_clear(&M[1]);
      return err;
    }
  }

  /* create mu, used for Barrett reduction */
  if ((err = mp_init (&mu)) != MP_OKAY) {
    goto LBL_M;
  }
  
  if (redmode == 0) {
     if ((err = mp_reduce_setup (&mu, P)) != MP_OKAY) {
        goto LBL_MU;
     }
     redux = mp_reduce;
  } else {
     if ((err = mp_reduce_2k_setup_l (P, &mu)) != MP_OKAY) {
        goto LBL_MU;
     }
     redux = mp_reduce_2k_l;
  }    

  /* create M table
   *
   * The M table contains powers of the base, 
   * e.g. M[x] = G**x mod P
   *
   * The first half of the table is not 
   * computed though accept for M[0] and M[1]
   */
  if ((err = mp_mod (G, P, &M[1])) != MP_OKAY) {
    goto LBL_MU;
  }

  /* compute the value at M[1<<(winsize-1)] by squaring 
   * M[1] (winsize-1) times 
   */
  if ((err = mp_copy (&M[1], &M[1 << (winsize - 1)])) != MP_OKAY) {
    goto LBL_MU;
  }

  for (x = 0; x < (winsize - 1); x++) {
    /* square it */
    if ((err = mp_sqr (&M[1 << (winsize - 1)], 
                       &M[1 << (winsize - 1)])) != MP_OKAY) {
      goto LBL_MU;
    }

    /* reduce modulo P */
    if ((err = redux (&M[1 << (winsize - 1)], P, &mu)) != MP_OKAY) {
      goto LBL_MU;
    }
  }

  /* create upper table, that is M[x] = M[x-1] * M[1] (mod P)
   * for x = (2**(winsize - 1) + 1) to (2**winsize - 1)
   */
  for (x = (1 << (winsize - 1)) + 1; x < (1 << winsize); x++) {
    if ((err = mp_mul (&M[x - 1], &M[1], &M[x])) != MP_OKAY) {
      goto LBL_MU;
    }
    if ((err = redux (&M[x], P, &mu)) != MP_OKAY) {
      goto LBL_MU;
    }
  }

  /* setup result */
  if ((err = mp_init (&res)) != MP_OKAY) {
    goto LBL_MU;
  }
  mp_set (&res, 1);

  /* set initial mode and bit cnt */
  mode   = 0;
  bitcnt = 1;
  buf    = 0;
  digidx = X->used - 1;
  bitcpy = 0;
  bitbuf = 0;

  for (;;) {
    /* grab next digit as required */
    if (--bitcnt == 0) {
      /* if digidx == -1 we are out of digits */
      if (digidx == -1) {
        break;
      }
      /* read next digit and reset the bitcnt */
      buf    = X->dp[digidx--];
      bitcnt = (int) DIGIT_BIT;
    }

    /* grab the next msb from the exponent */
    y     = (buf >> (mp_digit)(DIGIT_BIT - 1)) & 1;
    buf <<= (mp_digit)1;

    /* if the bit is zero and mode == 0 then we ignore it
     * These represent the leading zero bits before the first 1 bit
     * in the exponent.  Technically this opt is not required but it
     * does lower the # of trivial squaring/reductions used
     */
    if (mode == 0 && y == 0) {
      continue;
    }

    /* if the bit is zero and mode == 1 then we square */
    if (mode == 1 && y == 0) {
      if ((err = mp_sqr (&res, &res)) != MP_OKAY) {
        goto LBL_RES;
      }
      if ((err = redux (&res, P, &mu)) != MP_OKAY) {
        goto LBL_RES;
      }
      continue;
    }

    /* else we add it to the window */
    bitbuf |= (y << (winsize - ++bitcpy));
    mode    = 2;

    if (bitcpy == winsize) {
      /* ok window is filled so square as required and multiply  */
      /* square first */
      for (x = 0; x < winsize; x++) {
        if ((err = mp_sqr (&res, &res)) != MP_OKAY) {
          goto LBL_RES;
        }
        if ((err = redux (&res, P, &mu)) != MP_OKAY) {
          goto LBL_RES;
        }
      }

      /* then multiply */
      if ((err = mp_mul (&res, &M[bitbuf], &res)) != MP_OKAY) {
        goto LBL_RES;
      }
      if ((err = redux (&res, P, &mu)) != MP_OKAY) {
        goto LBL_RES;
      }

      /* empty window and reset */
      bitcpy = 0;
      bitbuf = 0;
      mode   = 1;
    }
  }

  /* if bits remain then square/multiply */
  if (mode == 2 && bitcpy > 0) {
    /* square then multiply if the bit is set */
    for (x = 0; x < bitcpy; x++) {
      if ((err = mp_sqr (&res, &res)) != MP_OKAY) {
        goto LBL_RES;
      }
      if ((err = redux (&res, P, &mu)) != MP_OKAY) {
        goto LBL_RES;
      }

      bitbuf <<= 1;
      if ((bitbuf & (1 << winsize)) != 0) {
        /* then multiply */
        if ((err = mp_mul (&res, &M[1], &res)) != MP_OKAY) {
          goto LBL_RES;
        }
        if ((err = redux (&res, P, &mu)) != MP_OKAY) {
          goto LBL_RES;
        }
      }
    }
  }

  mp_exch (&res, Y);
  err = MP_OKAY;
LBL_RES:mp_clear (&res);
LBL_MU:mp_clear (&mu);
LBL_M:
  mp_clear(&M[1]);
  for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
    mp_clear (&M[x]);
  }
  return err;
}
void mp_set(mp_int * a, mp_digit b)
{
	mp_zero(a);
	a->dp[0] = b & MP_MASK;
	a->used = (a->dp[0] != 0) ? 1 : 0;
}


/* stores a bignum as a ASCII string in a given radix (2..64)
?????mp_int???a??r?????,????????temp(char?)???????????????
*/
int mp_toradix(mp_int * a, char *str, int radix)
{
	int     res, digs;
	mp_int  t;
	mp_digit d;
	char   *_s = str;

	/* check range of the radix */
	if (radix < 2 || radix > 64) {
		return MP_VAL;
	}

	/* quick out if its zero */
	if (mp_iszero(a) == 1) {
		*str++ = '0';
		*str = '\0';
		return MP_OKAY;
	}

	if ((res = mp_init_copy(&t, a)) != MP_OKAY) {
		return res;
	}

	/* if it is negative output a - */
	if (t.sign == MP_NEG) {
		++_s;
		*str++ = '-';
		t.sign = MP_ZPOS;
	}

	digs = 0;
	while (mp_iszero(&t) == 0) {
		if ((res = mp_div_d(&t, (mp_digit)radix, &t, &d)) != MP_OKAY) {
			mp_clear(&t);
			return res;
		}
		*str++ = mp_s_rmap[d];
		++digs;
	}

	/* reverse the digits of the string.  In this case _s points
	* to the first digit [exluding the sign] of the number]
	*/
	bn_reverse((unsigned char *)_s, digs);

	/* append a NULL so the string is properly terminated */
	*str = '\0';

	mp_clear(&t);
	return MP_OKAY;
}


int mp_init_set_int(mp_int * a, unsigned long b)
{
	int err;
	if ((err = mp_init(a)) != MP_OKAY) 
	{
		return err;
	}
	return mp_set_int(a, b);
}

/* initialize and set a digit 
??b???????
???mp_int????????????mp_int?????????????(1?2?)
*/
int mp_init_set(mp_int * a, mp_digit b)
{
	int err;
	if ((err = mp_init(a)) != MP_OKAY)
	{
		return err;
	}
	mp_set(a, b);
	return err;
}

/* makes a truly random prime of a given size (bits),
*
* Flags are as follows:
*
*   LTM_PRIME_BBS      - make prime congruent to 3 mod 4
*   LTM_PRIME_SAFE     - make sure (p-1)/2 is prime as well (implies LTM_PRIME_BBS)
*   LTM_PRIME_2MSB_OFF - make the 2nd highest bit zero
*   LTM_PRIME_2MSB_ON  - make the 2nd highest bit one
*
* You have to supply a callback which fills in a buffer with random bytes.  "dat" is a parameter you can
* have passed to the callback (e.g. a state or something).  This function doesn't use "dat" itself
* so it can be NULL
*
*/

/* This is possibly the mother of all prime generation functions, muahahahahaha!
*dat ??????NULL
????????,a????,t?????,????10,size?????????
*/
int mp_prime_random_ex(mp_int *a, int t, int size, int flags, ltm_prime_callback cb, void *dat)
{
	unsigned char *tmp, maskAND, maskOR_msb, maskOR_lsb;
	int res, err, bsize, maskOR_msb_offset;

	/* sanity check the input */
	if (size <= 1 || t <= 0) {
		return MP_VAL;
	}

	/* LTM_PRIME_SAFE implies LTM_PRIME_BBS */
	if (flags & LTM_PRIME_SAFE) 
	{
		flags |= LTM_PRIME_BBS;
	}

	/* calc the byte size */
	bsize = (size >> 3) + ((size & 7) ? 1 : 0);

	/* we need a buffer of bsize bytes */
	tmp = OPT_CAST(unsigned char) XMALLOC(bsize);
	if (tmp == NULL) 
	{
		return MP_MEM;
	}

	/* calc the maskAND value for the MSbyte*/
	maskAND = ((size & 7) == 0) ? 0xFF : (0xFF >> (8 - (size & 7)));

	/* calc the maskOR_msb */
	maskOR_msb = 0;
	maskOR_msb_offset = ((size & 7) == 1) ? 1 : 0;
	if (flags & LTM_PRIME_2MSB_ON) 
	{
		maskOR_msb |= 0x80 >> ((9 - size) & 7);
	}

	/* get the maskOR_lsb */
	maskOR_lsb = 1;
	if (flags & LTM_PRIME_BBS) 
	{
		maskOR_lsb |= 3;
	}

	do 
		{
		/* read the bytes */
		if (cb(tmp, bsize, dat) != bsize) {
			err = MP_VAL;
			goto error;
		}

		/* work over the MSbyte */
		tmp[0] &= maskAND;
		tmp[0] |= 1 << ((size - 1) & 7);

		/* mix in the maskORs */
		tmp[maskOR_msb_offset] |= maskOR_msb;
		tmp[bsize - 1] |= maskOR_lsb;

		/* read it in */
		if ((err = mp_read_unsigned_bin(a, tmp, bsize)) != MP_OKAY)     { goto error; }

		/* is it prime? */
		if ((err = mp_prime_is_prime(a, t, &res)) != MP_OKAY)           { goto error; }
		if (res == MP_NO) {
			continue;
		}

		if (flags & LTM_PRIME_SAFE) {
			/* see if (a-1)/2 is prime */
			if ((err = mp_sub_d(a, 1, a)) != MP_OKAY)                    { goto error; }
			if ((err = mp_div_2(a, a)) != MP_OKAY)                       { goto error; }

			/* is it prime? */
			if ((err = mp_prime_is_prime(a, t, &res)) != MP_OKAY)        { goto error; }
		}
	} while (res == MP_NO);

	if (flags & LTM_PRIME_SAFE) {
		/* restore a to the original value */
		if ((err = mp_mul_2(a, a)) != MP_OKAY)                          { goto error; }
		if ((err = mp_add_d(a, 1, a)) != MP_OKAY)                       { goto error; }
	}

	err = MP_OKAY;
error:
	XFREE(tmp);
	return err;
}

/*
?????a?b??,??b???????
*/
int mp_expt_d(mp_int * a, mp_digit b, mp_int * c)
{
	int     res, x;
	mp_int  g;
	
	if ((res = mp_init_copy(&g, a)) != MP_OKAY) {
		return res;
	}

	/* set initial result */
	mp_set(c, 1);

	for (x = 0; x < (int)DIGIT_BIT; x++) {
		/* square */
		if ((res = mp_sqr(c, c)) != MP_OKAY) {
			mp_clear(&g);
			return res;
		}

		/* if the bit is set multiply */
		if ((b & (mp_digit)(((mp_digit)1) << (DIGIT_BIT - 1))) != 0) {
			if ((res = mp_mul(c, &g, c)) != MP_OKAY) {
				mp_clear(&g);
				return res;
			}
		}

		/* shift to next bit */
		b <<= 1;
	}

	mp_clear(&g);
	return MP_OKAY;
}

/*
??mp_int???a,??b=a^2,???????mp_int b???????
*/
int mp_sqr(mp_int * a, mp_int * b)
{
	int     res;

#ifdef BN_MP_TOOM_SQR_C
	/* use Toom-Cook? */
	if (a->used >= TOOM_SQR_CUTOFF) {
		res = mp_toom_sqr(a, b);
		/* Karatsuba? */
	}
	else
#endif
#ifdef BN_MP_KARATSUBA_SQR_C
		if (a->used >= KARATSUBA_SQR_CUTOFF) {
			res = mp_karatsuba_sqr(a, b);
		}
		else
#endif
		{
#ifdef BN_FAST_S_MP_SQR_C
			/* can we use the fast comba multiplier? */
			if ((a->used * 2 + 1) < MP_WARRAY &&
				a->used <
				(1 << (sizeof(mp_word) * CHAR_BIT - 2 * DIGIT_BIT - 1))) {
				res = fast_s_mp_sqr(a, b);
			}
			else
#endif
#ifdef BN_S_MP_SQR_C
				res = s_mp_sqr(a, b);
#else
				res = MP_VAL;
#endif
		}
	b->sign = MP_ZPOS;
	return res;
}

/*
??????????
*/
int mp_sqrt(mp_int *arg, mp_int *ret)
{
	int res;
	mp_int t1, t2;

	/* must be positive */
	if (arg->sign == MP_NEG) {
		return MP_VAL;
	}

	/* easy out */
	if (mp_iszero(arg) == MP_YES) {
		mp_zero(ret);
		return MP_OKAY;
	}

	if ((res = mp_init_copy(&t1, arg)) != MP_OKAY) {
		return res;
	}

	if ((res = mp_init(&t2)) != MP_OKAY) {
		goto E2;
	}

	/* First approx. (not very bad for large arg) */
	mp_rshd(&t1, t1.used / 2);

	/* t1 > 0  */
	if ((res = mp_div(arg, &t1, &t2, NULL)) != MP_OKAY) {
		goto E1;
	}
	if ((res = mp_add(&t1, &t2, &t1)) != MP_OKAY) {
		goto E1;
	}
	if ((res = mp_div_2(&t1, &t1)) != MP_OKAY) {
		goto E1;
	}
	/* And now t1 > sqrt(arg) */
	do {
		if ((res = mp_div(arg, &t1, &t2, NULL)) != MP_OKAY) {
			goto E1;
		}
		if ((res = mp_add(&t1, &t2, &t1)) != MP_OKAY) {
			goto E1;
		}
		if ((res = mp_div_2(&t1, &t1)) != MP_OKAY) {
			goto E1;
		}
		/* t1 >= sqrt(arg) >= t2 at this point */
	} while (mp_cmp_mag(&t1, &t2) == MP_GT);

	mp_exch(&t1, ret);

E1: mp_clear(&t2);
E2: mp_clear(&t1);
	return res;
}

/* multiply by a digit
?????,???????mp_int?????????
??mp_int a?mp_digit b ??c=ab
*/
int mp_mul_d(mp_int * a, mp_digit b, mp_int * c)
{
	mp_digit u, *tmpa, *tmpc;
	mp_word  r;
	int      ix, res, olduse;

	/* make sure c is big enough to hold a*b */
	if (c->alloc < a->used + 1) {
		if ((res = mp_grow(c, a->used + 1)) != MP_OKAY) 
		{
			return res;
		}
	}

	/* get the original destinations used count */
	olduse = c->used;

	/* set the sign */
	c->sign = a->sign;

	/* alias for a->dp [source] */
	tmpa = a->dp;

	/* alias for c->dp [dest] */
	tmpc = c->dp;

	/* zero carry */
	u = 0;

	/* compute columns */
	for (ix = 0; ix < a->used; ix++) {
		/* compute product and carry sum for this term */
		r = ((mp_word)u) + ((mp_word)*tmpa++) * ((mp_word)b);

		/* mask off higher bits to get a single digit */
		*tmpc++ = (mp_digit)(r & ((mp_word)MP_MASK));

		/* send carry into next iteration */
		u = (mp_digit)(r >> ((mp_word)DIGIT_BIT));
	}

	/* store final carry [if any] and increment ix offset  */
	*tmpc++ = u;
	++ix;

	/* now zero digits above the top */
	while (ix++ < olduse) {
		*tmpc++ = 0;
	}

	/* set used count */
	c->used = a->used + 1;
	mp_clamp(c);

	return MP_OKAY;
}

/* high level multiplication (handles sign)
???????????,??????????c=ab
*/
int mp_mul(mp_int * a, mp_int * b, mp_int * c)
{
	int     res, neg;
	neg = (a->sign == b->sign) ? MP_ZPOS : MP_NEG;

	/* use Toom-Cook? */
#ifdef BN_MP_TOOM_MUL_C
	if (MIN(a->used, b->used) >= TOOM_MUL_CUTOFF) {
		res = mp_toom_mul(a, b, c);
	}
	else
#endif
#ifdef BN_MP_KARATSUBA_MUL_C
		/* use Karatsuba? */
		if (MIN(a->used, b->used) >= KARATSUBA_MUL_CUTOFF) {
			res = mp_karatsuba_mul(a, b, c);
		}
		else
#endif
		{
			/* can we use the fast multiplier?
			*
			* The fast multiplier can be used if the output will
			* have less than MP_WARRAY digits and the number of
			* digits won't affect carry propagation
			*/
			int     digs = a->used + b->used + 1;

#ifdef BN_FAST_S_MP_MUL_DIGS_C
			if ((digs < MP_WARRAY) &&
				MIN(a->used, b->used) <=
				(1 << ((CHAR_BIT * sizeof(mp_word)) - (2 * DIGIT_BIT)))) {
				res = fast_s_mp_mul_digs(a, b, c, digs);
			}
			else
#endif
#ifdef BN_S_MP_MUL_DIGS_C
				res = s_mp_mul(a, b, c); /* uses s_mp_mul_digs */
#else
				res = MP_VAL;
#endif

		}
	c->sign = (c->used > 0) ? neg : MP_ZPOS;
	return res;
}

/*
????mp_int???a?b,???????c=a+b
*/
int mp_add(mp_int * a, mp_int * b, mp_int * c)
{
	int     sa, sb, res;

	/* get sign of both inputs */
	sa = a->sign;
	sb = b->sign;

	/* handle two cases, not four */
	if (sa == sb) {
		/* both positive or both negative */
		/* add their magnitudes, copy the sign */
		c->sign = sa;
		res = s_mp_add(a, b, c);
	}
	else {
		/* one positive, the other negative */
		/* subtract the one with the greater magnitude from */
		/* the one of the lesser magnitude.  The result gets */
		/* the sign of the one with the greater magnitude. */
		if (mp_cmp_mag(a, b) == MP_LT) {
			c->sign = sb;
			res = s_mp_sub(b, a, c);
		}
		else {
			c->sign = sa;
			res = s_mp_sub(a, b, c);
		}
	}
	return res;
}

/*
????mp_int???a?b,???????c=a-b
*/
int mp_sub(mp_int * a, mp_int * b, mp_int * c)
{
	int     sa, sb, res;

	sa = a->sign;
	sb = b->sign;

	if (sa != sb) {
		/* subtract a negative from a positive, OR */
		/* subtract a positive from a negative.    */
		/* In either case, ADD their magnitudes,   */
		/* and use the sign of the first number.   */
		c->sign = sa;
		res = s_mp_add(a, b, c);
	}
	else {
		/* subtract a positive from a positive, OR  */
		/* subtract a negative from a negative.     */
		/* First, take the difference between their */
		/* magnitudes, then... */
		if (mp_cmp_mag(a, b) != MP_LT) {
			/* Copy the sign from the first */
			c->sign = sa;
			/* The first has a larger or equal magnitude */
			res = s_mp_sub(a, b, c);
		}
		else {
			/* The result has the *opposite* sign from */
			/* the first number. */
			c->sign = (sa == MP_ZPOS) ? MP_NEG : MP_ZPOS;
			/* The second has a larger magnitude */
			res = s_mp_sub(b, a, c);
		}
	}
	return res;
}

/*
????mp_int???a?b,??c= a mod b,0 <= c <= b
*/
int mp_mod(mp_int * a, mp_int * b, mp_int * c)
{
	mp_int  t;
	int     res;

	if ((res = mp_init(&t)) != MP_OKAY) {
		return res;
	}

	if ((res = mp_div(a, b, NULL, &t)) != MP_OKAY) {
		mp_clear(&t);
		return res;
	}

	if (t.sign != b->sign) {
		res = mp_add(b, &t, c);
	}
	else {
		res = MP_OKAY;
		mp_exch(&t, c);
	}

	mp_clear(&t);
	return res;
}

/*
??????,????mp_int???a?b,??????????(a?b???)
??????MP_GT=1????,MP_EQ=0????,MP_LT=-1????
*/
int mp_cmp(mp_int * a, mp_int * b)
{
	/* compare based on sign */
	if (a->sign != b->sign) {
		if (a->sign == MP_NEG) {
			return MP_LT;
		}
		else {
			return MP_GT;
		}
	}

	/* compare digits */
	if (a->sign == MP_NEG) {
		/* if negative compare opposite direction */
		return mp_cmp_mag(b, a);
	}
	else {
		return mp_cmp_mag(a, b);
	}
}

/*
mp_int???a,??:???a?????
*/
void mp_clear(mp_int * a)
{
	int i;

	/* only do anything if a hasn't been freed previously */
	if (a->dp != NULL) { // if (a->dp != NULL) 
		/* first zero the digits */
		for (i = 0; i < a->used; i++) {
			a->dp[i] = 0;
		}

		/* free ram */
		XFREE(a->dp);

		/* reset members to make debugging easier */
		a->dp = NULL;
		a->alloc = a->used = 0;
		a->sign = MP_ZPOS;
	}
}

/*
??????mp_int??a,???????,mp_int??b??a?????????
mp_int??b???a?????????,mp_int??a??????????mp_int??b??
*/
int mp_copy(mp_int * a, mp_int * b)
{
	int   res, n;

	/* if dst == src do nothing */
	if (a == b) {
		return MP_OKAY;
	}

	/* grow dest */
	if (b->alloc < a->used) {
		if ((res = mp_grow(b, a->used)) != MP_OKAY) {
			return res;
		}
	}

	/* zero b and copy the parameters over */
  {
	  register mp_digit *tmpa, *tmpb;

	  /* pointer aliases */

	  /* source */
	  tmpa = a->dp;

	  /* destination */
	  tmpb = b->dp;

	  /* copy all the digits */
	  for (n = 0; n < a->used; n++) {
		  *tmpb++ = *tmpa++;
	  }

	  /* clear high digits */
	  for (; n < b->used; n++) {
		  *tmpb++ = 0;
	  }
  }

  /* copy used count and sign */
  b->used = a->used;
  b->sign = a->sign;
  return MP_OKAY;
}

/*
???mp_int??a,??b???a
*/
int mp_init_copy(mp_int * a, mp_int * b)
{
	int     res;

	if ((res = mp_init(a)) != MP_OKAY) {
		return res;
	}
	return mp_copy(b, a);
}

/*
??mp_int???a,??a?????
*/
void mp_zero(mp_int * a)
{
	int       n;
	mp_digit *tmp;

	a->sign = MP_ZPOS;
	a->used = 0;

	tmp = a->dp;
	for (n = 0; n < a->alloc; n++) {
		*tmp++ = 0;
	}
}

/*
??mp_int???a?b,(a,b)=1, p >= 2,0<a<p ,
????c=a^-1(mod b).?ab=1(mod p)
*/
int mp_invmod(mp_int * a, mp_int * b, mp_int * c)
{   /* hac 14.61, pp608 */
	/* b cannot be negative */
	if (b->sign == MP_NEG || mp_iszero(b) == 1) {
		return MP_VAL;
	}
#ifdef BN_FAST_MP_INVMOD_C  //????
	/* if the modulus is odd we can use a faster routine instead */
	if (mp_isodd(b) == 1) {
		return fast_mp_invmod(a, b, c);
	}
#endif

#ifdef BN_MP_INVMOD_SLOW_C  //????
	return mp_invmod_slow(a, b, c);
#endif

	return MP_VAL;
}

/*
????? d = a*b(mod c)
*/
int mp_mulmod(mp_int * a, mp_int * b, mp_int * c, mp_int * d)
{  /* d = a * b (mod c) */
	int     res;
	mp_int  t;

	if ((res = mp_init(&t)) != MP_OKAY) {
		return res;
	}

	if ((res = mp_mul(a, b, &t)) != MP_OKAY) {
		mp_clear(&t);
		return res;
	}
	res = mp_mod(&t, c, d);
	mp_clear(&t);
	return res;
}

/* b = -a */
int mp_neg(mp_int * a, mp_int * b)
{
	int     res;
	if (a != b) {
		if ((res = mp_copy(a, b)) != MP_OKAY) {
			return res;
		}
	}

	if (mp_iszero(b) != MP_YES) {
		b->sign = (a->sign == MP_ZPOS) ? MP_NEG : MP_ZPOS;
	}
	else {
		b->sign = MP_ZPOS;
	}

	return MP_OKAY;
}


int mp_submod(mp_int * a, mp_int * b, mp_int * c, mp_int * d)
{
	int     res;
	mp_int  t;


	if ((res = mp_init(&t)) != MP_OKAY) {
		return res;
	}

	if ((res = mp_sub(a, b, &t)) != MP_OKAY) {
		mp_clear(&t);
		return res;
	}
	res = mp_mod(&t, c, d);
	mp_clear(&t);
	return res;
}




int mp_lshd(mp_int * a, int b)
{
	int     x, res;

	/* if its less than zero return */
	if (b <= 0) {
		return MP_OKAY;
	}

	/* grow to fit the new digits */
	if (a->alloc < a->used + b) {
		if ((res = mp_grow(a, a->used + b)) != MP_OKAY) {
			return res;
		}
	}

  {
	  register mp_digit *top, *bottom;

	  /* increment the used by the shift amount then copy upwards */
	  a->used += b;

	  /* top */
	  top = a->dp + a->used - 1;

	  /* base */
	  bottom = a->dp + a->used - 1 - b;

	  /* much like mp_rshd this is implemented using a sliding window
	  * except the window goes the otherway around.  Copying from
	  * the bottom to the top.  see bn_mp_rshd.c for more info.
	  */
	  for (x = a->used - 1; x >= b; x--) {
		  *top-- = *bottom--;
	  }

	  /* zero the lower digits */
	  top = a->dp;
	  for (x = 0; x < b; x++) {
		  *top++ = 0;
	  }
  }
  return MP_OKAY;
}



/* reads a unsigned char array, assumes the msb is stored first [big endian] */
int mp_read_unsigned_bin(mp_int * a, const unsigned char *b, int c)
{
	int     res;

	/* make sure there are at least two digits */
	if (a->alloc < 2) {
		if ((res = mp_grow(a, 2)) != MP_OKAY) {
			return res;
		}
	}

	/* zero the int */
	mp_zero(a);

	/* read the bytes in */
	while (c-- > 0) {
		if ((res = mp_mul_2d(a, 8, a)) != MP_OKAY) {
			return res;
		}

#ifndef MP_8BIT
		a->dp[0] |= *b++;
		a->used += 1;
#else
		a->dp[0] = (*b & MP_MASK);
		a->dp[1] |= ((*b++ >> 7U) & 1);
		a->used += 2;
#endif
	}
	mp_clamp(a);
	return MP_OKAY;
}

void bn_reverse(unsigned char *s, int len)
{
	int     ix, iy;
	unsigned char t;

	ix = 0;
	iy = len - 1;
	while (ix < iy) {
		t = s[ix];
		s[ix] = s[iy];
		s[iy] = t;
		++ix;
		--iy;
	}
}


int mp_cmp_d(mp_int * a, mp_digit b)
{
	/* compare based on sign */
	if (a->sign == MP_NEG) {
		return MP_LT;
	}

	/* compare based on magnitude */
	if (a->used > 1) {
		return MP_GT;
	}

	/* compare the only digit of a to b */
	if (a->dp[0] > b) {
		return MP_GT;
	}
	else if (a->dp[0] < b) {
		return MP_LT;
	}
	else {
		return MP_EQ;
	}
}

int mp_mod_2d(mp_int * a, int b, mp_int * c)
{
	int     x, res;

	/* if b is <= 0 then zero the int */
	if (b <= 0) {
		mp_zero(c);
		return MP_OKAY;
	}

	/* if the modulus is larger than the value than return */
	if (b >= (int)(a->used * DIGIT_BIT)) {
		res = mp_copy(a, c);
		return res;
	}

	/* copy */
	if ((res = mp_copy(a, c)) != MP_OKAY) {
		return res;
	}

	/* zero digits above the last digit of the modulus */
	for (x = (b / DIGIT_BIT) + ((b % DIGIT_BIT) == 0 ? 0 : 1); x < c->used; x++) {
		c->dp[x] = 0;
	}
	/* clear the digit that is not completely outside/inside the modulus */
	c->dp[b / DIGIT_BIT] &=
		(mp_digit)((((mp_digit)1) << (((mp_digit)b) % DIGIT_BIT)) - ((mp_digit)1));
	mp_clamp(c);
	return MP_OKAY;
}


/* shift left by a certain bit count */
int mp_mul_2d(mp_int * a, int b, mp_int * c)
{
	mp_digit d;
	int      res;

	/* copy */
	if (a != c) {
		if ((res = mp_copy(a, c)) != MP_OKAY) {
			return res;
		}
	}

	if (c->alloc < (int)(c->used + b / DIGIT_BIT + 1)) {
		if ((res = mp_grow(c, c->used + b / DIGIT_BIT + 1)) != MP_OKAY) {
			return res;
		}
	}

	/* shift by as many digits in the bit count */
	if (b >= (int)DIGIT_BIT) {
		if ((res = mp_lshd(c, b / DIGIT_BIT)) != MP_OKAY) {
			return res;
		}
	}

	/* shift any bit count < DIGIT_BIT */
	d = (mp_digit)(b % DIGIT_BIT);
	if (d != 0) {
		register mp_digit *tmpc, shift, mask, r, rr;
		register int x;

		/* bitmask for carries */
		mask = (((mp_digit)1) << d) - 1;

		/* shift for msbs */
		shift = DIGIT_BIT - d;

		/* alias */
		tmpc = c->dp;

		/* carry */
		r = 0;
		for (x = 0; x < c->used; x++) {
			/* get the higher bits of the current word */
			rr = (*tmpc >> shift) & mask;

			/* shift the current word and OR in the carry */
			*tmpc = ((*tmpc << d) | r) & MP_MASK;
			++tmpc;

			/* set the carry to the carry bits of the current word */
			r = rr;
		}

		/* set final carry */
		if (r != 0) {
			c->dp[(c->used)++] = r;
		}
	}
	mp_clamp(c);
	return MP_OKAY;
}


static int s_is_power_of_two(mp_digit b, int *p)
{
	int x;

	for (x = 1; x < DIGIT_BIT; x++) {
		if (b == (((mp_digit)1) << x)) {
			*p = x;
			return 1;
		}
	}
	return 0;
}

void mp_exch(mp_int * a, mp_int * b)
{
	mp_int  t;

	t = *a;
	*a = *b;
	*b = t;
}



int mp_grow(mp_int * a, int size)
{
	int     i;
	mp_digit *tmp;

	/* if the alloc size is smaller alloc more ram */
	if (a->alloc < size) {
		/* ensure there are always at least MP_PREC digits extra on top */
		size += (MP_PREC * 2) - (size % MP_PREC);

		/* reallocate the array a->dp
		*
		* We store the return in a temporary variable
		* in case the operation failed we don't want
		* to overwrite the dp member of a.
		*/
		tmp = OPT_CAST(mp_digit) XREALLOC(a->dp, sizeof(mp_digit) * size);
		if (!tmp) {  //if (tmp == NULL)
			/* reallocation failed but "a" is still valid [can be freed] */
			return MP_MEM;
		}

		/* reallocation succeeded so set a->dp */
		a->dp = tmp;

		/* zero excess digits */
		i = a->alloc;
		a->alloc = size;
		for (; i < a->alloc; i++) {
			a->dp[i] = 0;
		}
	}
	return MP_OKAY;
}

int mp_init_size(mp_int * a, int size)
{
	int x;

	/* pad size so there are always extra digits */
	size += (MP_PREC * 2) - (size % MP_PREC);

	/* alloc mem */
	a->dp = OPT_CAST(mp_digit) XMALLOC(sizeof(mp_digit) * size);
	if (a->dp == NULL) {
		return MP_MEM;
	}

	/* set the members */
	a->used = 0;
	a->alloc = size;
	a->sign = MP_ZPOS;

	/* zero the digits */
	for (x = 0; x < size; x++) {
		a->dp[x] = 0;
	}
	return MP_OKAY;
}

void mp_clamp(mp_int * a)
{
	/* decrease used while the most significant digit is
	* zero.
	*/
	while (a->used > 0 && a->dp[a->used - 1] == 0) {
		--(a->used);
	}

	/* reset the sign flag if used == 0 */
	if (a->used == 0) {
		a->sign = MP_ZPOS;
	}
}

int mp_mod_d(mp_int * a, mp_digit b, mp_digit * c)                                  
{
	return mp_div_d(a, b, NULL, c);
}
int mp_div_2d(mp_int * a, int b, mp_int * c, mp_int * d)        
{
	mp_digit D, r, rr;
	int     x, res;
	mp_int  t;


	/* if the shift count is <= 0 then we do no work */
	if (b <= 0) {
		res = mp_copy(a, c);
		if (d != NULL) {  // if (d != NULL)
			mp_zero(d);
		}
		return res;
	}

	if ((res = mp_init(&t)) != MP_OKAY) {
		return res;
	}

	/* get the remainder */
	if (d != NULL) {  // if (d != NULL)
		if ((res = mp_mod_2d(a, b, &t)) != MP_OKAY) {
			mp_clear(&t);
			return res;
		}
	}

	/* copy */
	if ((res = mp_copy(a, c)) != MP_OKAY) {
		mp_clear(&t);
		return res;
	}

	/* shift by as many digits in the bit count */
	if (b >= (int)DIGIT_BIT) {
		mp_rshd(c, b / DIGIT_BIT);
	}

	/* shift any bit count < DIGIT_BIT */
	D = (mp_digit)(b % DIGIT_BIT);
	if (D != 0) {
		register mp_digit *tmpc, mask, shift;

		/* mask */
		mask = (((mp_digit)1) << D) - 1;

		/* shift for lsb */
		shift = DIGIT_BIT - D;

		/* alias */
		tmpc = c->dp + (c->used - 1);

		/* carry */
		r = 0;
		for (x = c->used - 1; x >= 0; x--) {
			/* get the lower  bits of this word in a temp */
			rr = *tmpc & mask;

			/* shift the current word and mix in the carry bits from the previous word */
			*tmpc = (*tmpc >> D) | (r << shift);
			--tmpc;

			/* set the carry to the carry bits of the current word found above */
			r = rr;
		}
	}
	mp_clamp(c);
	if (d != NULL) { // if (d != NULL)
		mp_exch(&t, d);
	}
	mp_clear(&t);
	return MP_OKAY;
}

void mp_rshd(mp_int * a, int b)
{
	int     x;

	/* if b <= 0 then ignore it */
	if (b <= 0) {
		return;
	}

	/* if b > used then simply zero it and return */
	if (a->used <= b) {
		mp_zero(a);
		return;
	}

  {
	  register mp_digit *bottom, *top;

	  /* shift the digits down */

	  /* bottom */
	  bottom = a->dp;

	  /* top [offset into digits] */
	  top = a->dp + b;

	  /* this is implemented as a sliding window where
	  * the window is b-digits long and digits from
	  * the top of the window are copied to the bottom
	  *
	  * e.g.

	  b-2 | b-1 | b0 | b1 | b2 | ... | bb |   ---->
	  /\                   |      ---->
	  \-------------------/      ---->
	  */
	  for (x = 0; x < (a->used - b); x++) {
		  *bottom++ = *top++;
	  }

	  /* zero the top digits */
	  for (; x < a->used; x++) {
		  *bottom++ = 0;
	  }
  }

  /* remove excess digits */
  a->used -= b;
}

int mp_div_d(mp_int * a, mp_digit b, mp_int * c, mp_digit * d)
{
	mp_int  q;
	mp_word w;
	mp_digit t;
	int     res, ix;

	/* cannot divide by zero */
	if (b == 0) {
		return MP_VAL;
	}

	/* quick outs */
	if (b == 1 || mp_iszero(a) == 1) {
		if (d != NULL) { // if (d != NULL)
			*d = 0;
		}
		if (c != NULL) { // if (c != NULL)
			return mp_copy(a, c);
		}
		return MP_OKAY;
	}

	/* power of two ? */
	if (s_is_power_of_two(b, &ix) == 1) {
		if(d != NULL) { // if (d != NULL)
			*d = a->dp[0] & ((((mp_digit)1) << ix) - 1);
		}
		if (c != NULL) { // if (c != NULL)
			return mp_div_2d(a, ix, c, NULL);
		}
		return MP_OKAY;
	}

#ifdef BN_MP_DIV_3_C
	/* three? */
	if (b == 3) {
		return mp_div_3(a, c, d);
	}
#endif

	/* no easy answer [c'est la vie].  Just division */
	if ((res = mp_init_size(&q, a->used)) != MP_OKAY) {
		return res;
	}

	q.used = a->used;
	q.sign = a->sign;
	w = 0;
	for (ix = a->used - 1; ix >= 0; ix--) {
		w = (w << ((mp_word)DIGIT_BIT)) | ((mp_word)a->dp[ix]);

		if (w >= b) {
			t = (mp_digit)(w / b);
			w -= ((mp_word)t) * ((mp_word)b);
		}
		else {
			t = 0;
		}
		q.dp[ix] = (mp_digit)t;
	}

	if (d != NULL) {
		*d = (mp_digit)w;
	}

	if (c != NULL) {
		mp_clamp(&q);
		mp_exch(&q, c);
	}
	mp_clear(&q);

	return res;
}

int mp_cmp_mag(mp_int * a, mp_int * b)
{
	int     n;
	mp_digit *tmpa, *tmpb;

	/* compare based on # of non-zero digits */
	if (a->used > b->used) {
		return MP_GT;
	}

	if (a->used < b->used) {
		return MP_LT;
	}

	/* alias for a */
	tmpa = a->dp + (a->used - 1);

	/* alias for b */
	tmpb = b->dp + (a->used - 1);

	/* compare based on digits  */
	for (n = 0; n < a->used; ++n, --tmpa, --tmpb) {
		if (*tmpa > *tmpb) {
			return MP_GT;
		}

		if (*tmpa < *tmpb) {
			return MP_LT;
		}
	}
	return MP_EQ;
}



int mp_prime_is_divisible(mp_int * a, int *result)
{
	int     err, ix;
	mp_digit res;

	/* default to not */
	*result = MP_NO;

	for (ix = 0; ix < PRIME_SIZE; ix++) {
		/* what is a mod LBL_prime_tab[ix] */
		if ((err = mp_mod_d(a, ltm_prime_tab[ix], &res)) != MP_OKAY) {
			return err;
		}

		/* is the residue zero? */
		if (res == 0) {
			*result = MP_YES;
			return MP_OKAY;
		}
	}

	return MP_OKAY;
}


int mp_prime_is_prime(mp_int * a, int t, int *result)
{
	mp_int  b;
	int     ix, err, res;

	/* default to no */
	*result = MP_NO;

	/* valid value of t? */
	if (t <= 0 || t > PRIME_SIZE) {
		return MP_VAL;
	}

	/* is the input equal to one of the primes in the table? */
	for (ix = 0; ix < PRIME_SIZE; ix++) {
		if (mp_cmp_d(a, ltm_prime_tab[ix]) == MP_EQ) {
			*result = 1;
			return MP_OKAY;
		}
	}

	/* first perform trial division */
	if ((err = mp_prime_is_divisible(a, &res)) != MP_OKAY) {
		return err;
	}

	/* return if it was trivially divisible */
	if (res == MP_YES) {
		return MP_OKAY;
	}

	/* now perform the miller-rabin rounds */
	if ((err = mp_init(&b)) != MP_OKAY) {
		return err;
	}

	for (ix = 0; ix < t; ix++) {
		/* set the prime */
		mp_set(&b, ltm_prime_tab[ix]);

		if ((err = mp_prime_miller_rabin(a, &b, &res)) != MP_OKAY) {
			goto LBL_B;
		}

		if (res == MP_NO) {
			goto LBL_B;
		}
	}

	/* passed the test */
	*result = MP_YES;
LBL_B:mp_clear(&b);
	return err;
}

 int mp_reduce_is_2k_l(mp_int *a)
 {
    int ix, iy;
 
    if (a->used == 0) {
       return MP_NO;
    } else if (a->used == 1) {
       return MP_YES;
    } else if (a->used > 1) {
       /* if more than half of the digits are -1 we're sold */
       for (iy = ix = 0; ix < a->used; ix++) {
           if (a->dp[ix] == MP_MASK) {
               ++iy;
           }
       }
       return (iy >= (a->used/2)) ? MP_YES : MP_NO;
 
    }
    return MP_NO;
 }

/* this is a shell function that calls either the normal or Montgomery
* exptmod functions.  Originally the call to the montgomery code was
* embedded in the normal function but that wasted alot of stack space
* for nothing (since 99% of the time the Montgomery code would be called)
*/
int mp_exptmod(mp_int * G, mp_int * X, mp_int * P, mp_int * Y)
{
	int dr;

	/* modulus P must be positive */
	if (P->sign == MP_NEG) {
		return MP_VAL;
	}

	/* if exponent X is negative we have to recurse */
	if (X->sign == MP_NEG) {
#ifdef BN_MP_INVMOD_C  //????????
		mp_int tmpG, tmpX;
		int err;

		/* first compute 1/G mod P */
		if ((err = mp_init(&tmpG)) != MP_OKAY) {
			return err;
		}
		if ((err = mp_invmod(G, P, &tmpG)) != MP_OKAY) {
			mp_clear(&tmpG);
			return err;
		}

		/* now get |X| */
		if ((err = mp_init(&tmpX)) != MP_OKAY) {
			mp_clear(&tmpG);
			return err;
		}
		if ((err = mp_abs(X, &tmpX)) != MP_OKAY) {
			mp_clear_multi(&tmpG, &tmpX, NULL);
			return err;
		}

		/* and now compute (1/G)**|X| instead of G**X [X < 0] */
		err = mp_exptmod(&tmpG, &tmpX, P, Y);
		mp_clear_multi(&tmpG, &tmpX, NULL);
		return err;
#else 
		/* no invmod */
		return MP_VAL;
#endif
	}

	/* modified diminished radix reduction */
//#if defined(BN_MP_REDUCE_IS_2K_L_C) && defined(BN_MP_REDUCE_2K_L_C)
#if 1
	if (mp_reduce_is_2k_l(P) == MP_YES) {
		return s_mp_exptmod(G, X, P, Y, 1);
	}
#endif

#ifdef BN_MP_DR_IS_MODULUS_C
	/* is it a DR modulus? */
	dr = mp_dr_is_modulus(P);
#else
	/* default to no */
	dr = 0;
#endif

#ifdef BN_MP_REDUCE_IS_2K_C
	/* if not, is it a unrestricted DR modulus? */
	if (dr == 0) {
		dr = mp_reduce_is_2k(P) << 1;
	}
#endif

	/* if the modulus is odd or dr != 0 use the montgomery method */
#ifdef BN_MP_EXPTMOD_FAST_C
	if (mp_isodd(P) == 1 || dr != 0) {
		return mp_exptmod_fast(G, X, P, Y, dr);
	}
	else {
#endif
#ifdef BN_S_MP_EXPTMOD_C
		/* otherwise use the generic Barrett reduction technique */
		return s_mp_exptmod(G, X, P, Y, 0);
#else
		/* no exptmod for evens */
		return MP_VAL;
#endif
#ifdef BN_MP_EXPTMOD_FAST_C
	}
#endif
}

/* returns the number of bits in an int */
int mp_count_bits(mp_int * a)
{
	int     r;
	mp_digit q;

	/* shortcut */
	if (a->used == 0) {
		return 0;
	}

	/* get number of digits and add that */
	r = (a->used - 1) * DIGIT_BIT;

	/* take the last digit and count the bits in it */
	q = a->dp[a->used - 1];
	while (q > ((mp_digit)0)) {
		++r;
		q >>= ((mp_digit)1);
	}
	return r;
}

/* b = a*2 */
int mp_mul_2(mp_int * a, mp_int * b)
{
	int     x, res, oldused;

	/* grow to accomodate result */
	if (b->alloc < a->used + 1) {
		if ((res = mp_grow(b, a->used + 1)) != MP_OKAY) {
			return res;
		}
	}

	oldused = b->used;
	b->used = a->used;

	{
		register mp_digit r, rr, *tmpa, *tmpb;

		/* alias for source */
		tmpa = a->dp;

		/* alias for dest */
		tmpb = b->dp;

		/* carry */
		r = 0;
		for (x = 0; x < a->used; x++) {

			/* get what will be the *next* carry bit from the
			* MSB of the current digit
			*/
			rr = *tmpa >> ((mp_digit)(DIGIT_BIT - 1));

			/* now shift up this digit, add in the carry [from the previous] */
			*tmpb++ = ((*tmpa++ << ((mp_digit)1)) | r) & MP_MASK;

			/* copy the carry that would be from the source
			* digit into the next iteration
			*/
			r = rr;
		}

		/* new leading digit? */
		if (r != 0) {
			/* add a MSB which is always 1 at this point */
			*tmpb = 1;
			++(b->used);
		}

		/* now zero any excess digits on the destination
		* that we didn't write to
		*/
		tmpb = b->dp + b->used;
		for (x = b->used; x < oldused; x++) {
			*tmpb++ = 0;
		}
	}
	b->sign = a->sign;
	return MP_OKAY;
}

/* b = a/2 */
int mp_div_2(mp_int * a, mp_int * b)
{
	int     x, res, oldused;

	/* copy */
	if (b->alloc < a->used) {
		if ((res = mp_grow(b, a->used)) != MP_OKAY) {
			return res;
		}
	}

	oldused = b->used;
	b->used = a->used;
	{
		register mp_digit r, rr, *tmpa, *tmpb;

		/* source alias */
		tmpa = a->dp + b->used - 1;

		/* dest alias */
		tmpb = b->dp + b->used - 1;

		/* carry */
		r = 0;
		for (x = b->used - 1; x >= 0; x--) {
			/* get the carry for the next iteration */
			rr = *tmpa & 1;

			/* shift the current digit, add in carry and store */
			*tmpb-- = (*tmpa-- >> 1) | (r << (DIGIT_BIT - 1));

			/* forward carry to next iteration */
			r = rr;
		}

		/* zero excess digits */
		tmpb = b->dp + b->used;
		for (x = b->used; x < oldused; x++) {
			*tmpb++ = 0;
		}
	}
	b->sign = a->sign;
	mp_clamp(b);
	return MP_OKAY;
}



/* low level addition, based on HAC pp.594, Algorithm 14.7 */
int s_mp_add(mp_int * a, mp_int * b, mp_int * c)
{
	mp_int *x;
	int     olduse, res, min, max;

	/* find sizes, we let |a| <= |b| which means we have to sort
	* them.  "x" will point to the input with the most digits
	*/
	if (a->used > b->used) {
		min = b->used;
		max = a->used;
		x = a;
	}
	else {
		min = a->used;
		max = b->used;
		x = b;
	}

	/* init result */
	if (c->alloc < max + 1) {
		if ((res = mp_grow(c, max + 1)) != MP_OKAY) {
			return res;
		}
	}

	/* get old used digit count and set new one */
	olduse = c->used;
	c->used = max + 1;

	{
		register mp_digit u, *tmpa, *tmpb, *tmpc;
		register int i;

		/* alias for digit pointers */

		/* first input */
		tmpa = a->dp;

		/* second input */
		tmpb = b->dp;

		/* destination */
		tmpc = c->dp;

		/* zero the carry */
		u = 0;
		for (i = 0; i < min; i++) {
			/* Compute the sum at one digit, T[i] = A[i] + B[i] + U */
			*tmpc = *tmpa++ + *tmpb++ + u;

			/* U = carry bit of T[i] */
			u = *tmpc >> ((mp_digit)DIGIT_BIT);

			/* take away carry bit from T[i] */
			*tmpc++ &= MP_MASK;
		}

		/* now copy higher words if any, that is in A+B
		* if A or B has more digits add those in
		*/
		if (min != max) {
			for (; i < max; i++) {
				/* T[i] = X[i] + U */
				*tmpc = x->dp[i] + u;

				/* U = carry bit of T[i] */
				u = *tmpc >> ((mp_digit)DIGIT_BIT);

				/* take away carry bit from T[i] */
				*tmpc++ &= MP_MASK;
			}
		}

		/* add carry */
		*tmpc++ = u;

		/* clear digits above oldused */
		for (i = c->used; i < olduse; i++) {
			*tmpc++ = 0;
		}
	}

	mp_clamp(c);
	return MP_OKAY;
}

/* low level subtraction (assumes |a| > |b|), HAC pp.595 Algorithm 14.9 */
int s_mp_sub(mp_int * a, mp_int * b, mp_int * c)
{
	int     olduse, res, min, max;

	/* find sizes */
	min = b->used;
	max = a->used;

	/* init result */
	if (c->alloc < max) {
		if ((res = mp_grow(c, max)) != MP_OKAY) {
			return res;
		}
	}
	olduse = c->used;
	c->used = max;

	{
		register mp_digit u, *tmpa, *tmpb, *tmpc;
		register int i;

		/* alias for digit pointers */
		tmpa = a->dp;
		tmpb = b->dp;
		tmpc = c->dp;

		/* set carry to zero */
		u = 0;
		for (i = 0; i < min; i++) {
			/* T[i] = A[i] - B[i] - U */
			*tmpc = *tmpa++ - *tmpb++ - u;

			/* U = carry bit of T[i]
			* Note this saves performing an AND operation since
			* if a carry does occur it will propagate all the way to the
			* MSB.  As a result a single shift is enough to get the carry
			*/
			u = *tmpc >> ((mp_digit)(CHAR_BIT * sizeof(mp_digit) - 1));   // 

			/* Clear carry from T[i] */
			*tmpc++ &= MP_MASK;
		}

		/* now copy higher words if any, e.g. if A has more digits than B  */
		for (; i < max; i++) {
			/* T[i] = A[i] - U */
			*tmpc = *tmpa++ - u;

			/* U = carry bit of T[i] */
			u = *tmpc >> ((mp_digit)(CHAR_BIT * sizeof(mp_digit) - 1));

			/* Clear carry from T[i] */
			*tmpc++ &= MP_MASK;
		}

		/* clear digits above used (since we may not have grown result above) */
		for (i = c->used; i < olduse; i++) {
			*tmpc++ = 0;
		}
	}

	mp_clamp(c);
	return MP_OKAY;
}

int mp_init_multi(mp_int *mp, ...)
{
	mp_err res = MP_OKAY;      /* Assume ok until proven otherwise */
	int n = 0;                 /* Number of ok inits */
	mp_int* cur_arg = mp;
	va_list args;

	va_start(args, mp);        /* init args to next argument from caller */// ???????  va_start ??????ARM???????
	while (cur_arg != NULL) {
		if (mp_init(cur_arg) != MP_OKAY) {
			/* Oops - error! Back-track and mp_clear what we already
			succeeded in init-ing, then return error.
			*/
			va_list clean_args;

			/* end the current list */
			va_end(args);

			/* now start cleaning up */
			cur_arg = mp;
			va_start(clean_args, mp);
			while (n--) {
				mp_clear(cur_arg);
				cur_arg = va_arg(clean_args,mp_int*);
			}
			va_end(clean_args);
		
			res = MP_MEM;
			break;
		}
		n++;
		cur_arg = va_arg(args, mp_int*);
	}
	va_end(args);
	return res;                /* Assumed ok, if error flagged above. */
}


/* slower bit-bang division... also smaller */

int mp_abs(mp_int * a, mp_int * b)
{
	int     res;

	/* copy a to b */
	if (a != b) {
		if ((res = mp_copy(a, b)) != MP_OKAY) {
			return res;
		}
	}

	/* force the sign of b to positive */
	b->sign = MP_ZPOS;

	return MP_OKAY;
}

/* single digit addition */
int mp_add_d(mp_int * a, mp_digit b, mp_int * c)
{
	int     res, ix, oldused;
	mp_digit *tmpa, *tmpc, mu;

	/* grow c as required */
	if (c->alloc < a->used + 1) {
		if ((res = mp_grow(c, a->used + 1)) != MP_OKAY) {
			return res;
		}
	}

	/* if a is negative and |a| >= b, call c = |a| - b */
	if (a->sign == MP_NEG && (a->used > 1 || a->dp[0] >= b)) {
		/* temporarily fix sign of a */
		a->sign = MP_ZPOS;

		/* c = |a| - b */
		res = mp_sub_d(a, b, c);

		/* fix sign  */
		a->sign = c->sign = MP_NEG;

		return res;
	}

	/* old number of used digits in c */
	oldused = c->used;

	/* sign always positive */
	c->sign = MP_ZPOS;

	/* source alias */
	tmpa = a->dp;

	/* destination alias */
	tmpc = c->dp;

	/* if a is positive */
	if (a->sign == MP_ZPOS) {
		/* add digit, after this we're propagating
		* the carry.
		*/
		*tmpc = *tmpa++ + b;
		mu = *tmpc >> DIGIT_BIT;
		*tmpc++ &= MP_MASK;

		/* now handle rest of the digits */
		for (ix = 1; ix < a->used; ix++) {
			*tmpc = *tmpa++ + mu;
			mu = *tmpc >> DIGIT_BIT;
			*tmpc++ &= MP_MASK;
		}
		/* set final carry */
		ix++;
		*tmpc++ = mu;

		/* setup size */
		c->used = a->used + 1;
	}
	else {
		/* a was negative and |a| < b */
		c->used = 1;

		/* the result is a single digit */
		if (a->used == 1) {
			*tmpc++ = b - a->dp[0];
		}
		else {
			*tmpc++ = b;
		}

		/* setup count so the clearing of oldused
		* can fall through correctly
		*/
		ix = 1;
	}

	/* now zero to oldused */
	while (ix++ < oldused) {
		*tmpc++ = 0;
	}
	mp_clamp(c);

	return MP_OKAY;
}


/* set a 32-bit const */
int mp_set_int(mp_int * a, unsigned long b)
{
	int     x, res;

	mp_zero(a);

	/* set four bits at a time */
	for (x = 0; x < 8; x++) {
		/* shift the number up four bits */
		if ((res = mp_mul_2d(a, 4, a)) != MP_OKAY) {
			return res;
		}

		/* OR in the top four bits of the source */
		a->dp[0] |= (b >> 28) & 15;

		/* shift the source up to the next four bits */
		b <<= 4;

		/* ensure that digits are not clamped off */
		a->used += 1;
	}
	mp_clamp(a);
	return MP_OKAY;
}

	/* integer signed division.
	* c*b + d == a [e.g. a/b, c=quotient, d=remainder]
	* HAC pp.598 Algorithm 14.20
	*
	* Note that the description in HAC is horribly
	* incomplete.  For example, it doesn't consider
	* the case where digits are removed from 'x' in
	* the inner loop.  It also doesn't consider the
	* case that y has fewer than three digits, etc..
	*
	* The overall algorithm is as described as
	* 14.20 from HAC but fixed to treat these cases.
	*/
	//-----------------------------------------------------------
	//===========================================================
	//-----------------------------------------------------------
int mp_div(mp_int * a, mp_int * b, mp_int * c, mp_int * d)   // ????????????????
{
		mp_int  q, x, y, t1, t2;
		int     res, n, t, i, norm, neg;
	
		/* is divisor zero ? */
		if (mp_iszero(b) == 1) {
			return MP_VAL;
		}
	
		/* if a < b then q=0, r = a */
		if (mp_cmp_mag(a, b) == MP_LT) {
			if (d != NULL) {
				res = mp_copy(a, d);
			}
			else {
				res = MP_OKAY;
			}
			if (c != NULL) {
				mp_zero(c);
			}
			return res;
		}
	
		if ((res = mp_init_size(&q, a->used + 2)) != MP_OKAY) {
			return res;
		}
		q.used = a->used + 2;
	
		if ((res = mp_init(&t1)) != MP_OKAY) {
			goto LBL_Q;
		}
	
		if ((res = mp_init(&t2)) != MP_OKAY) {
			goto LBL_T1;
		}
	
		if ((res = mp_init_copy(&x, a)) != MP_OKAY) {
			goto LBL_T2;
		}
	
		if ((res = mp_init_copy(&y, b)) != MP_OKAY) {
			goto LBL_X;
		}
	
		/* fix the sign */
		neg = (a->sign == b->sign) ? MP_ZPOS : MP_NEG;
		x.sign = y.sign = MP_ZPOS;
	
		/* normalize both x and y, ensure that y >= b/2, [b == 2**DIGIT_BIT] */
		norm = mp_count_bits(&y) % DIGIT_BIT;
		if (norm < (int)(DIGIT_BIT - 1)) {
			norm = (DIGIT_BIT - 1) - norm;
			if ((res = mp_mul_2d(&x, norm, &x)) != MP_OKAY) {
				goto LBL_Y;
			}
			if ((res = mp_mul_2d(&y, norm, &y)) != MP_OKAY) {
				goto LBL_Y;
			}
		}
		else {
			norm = 0;
		}
	
		/* note hac does 0 based, so if used==5 then its 0,1,2,3,4, e.g. use 4 */
		n = x.used - 1;
		t = y.used - 1;
	
		/* while (x >= y*b**n-t) do { q[n-t] += 1; x -= y*b**{n-t} } */
		if ((res = mp_lshd(&y, n - t)) != MP_OKAY) { /* y = y*b**{n-t} */
			goto LBL_Y;
		}
	
		while (mp_cmp(&x, &y) != MP_LT) {
			++(q.dp[n - t]);
			if ((res = mp_sub(&x, &y, &x)) != MP_OKAY) {
				goto LBL_Y;
			}
		}
	
		/* reset y by shifting it back down */
		mp_rshd(&y, n - t);
	
		/* step 3. for i from n down to (t + 1) */
		for (i = n; i >= (t + 1); i--) {
			if (i > x.used) {
				continue;
			}
	
			/* step 3.1 if xi == yt then set q{i-t-1} to b-1,
			* otherwise set q{i-t-1} to (xi*b + x{i-1})/yt */
			if (x.dp[i] == y.dp[t]) {
				q.dp[i - t - 1] = ((((mp_digit)1) << DIGIT_BIT) - 1);
			}
			else {
				mp_word tmp;
				tmp = ((mp_word)x.dp[i]) << ((mp_word)DIGIT_BIT);
				tmp |= ((mp_word)x.dp[i - 1]);
				tmp /= ((mp_word)y.dp[t]);
				if (tmp > (mp_word)MP_MASK)
					tmp = MP_MASK;
				q.dp[i - t - 1] = (mp_digit)(tmp & (mp_word)(MP_MASK));
			}
	
			/* while (q{i-t-1} * (yt * b + y{t-1})) >
			xi * b**2 + xi-1 * b + xi-2
	
			do q{i-t-1} -= 1;
			*/
			q.dp[i - t - 1] = (q.dp[i - t - 1] + 1) & MP_MASK;
			do {
				q.dp[i - t - 1] = (q.dp[i - t - 1] - 1) & MP_MASK;
	
				/* find left hand */
				mp_zero(&t1);
				t1.dp[0] = (t - 1 < 0) ? 0 : y.dp[t - 1];
				t1.dp[1] = y.dp[t];
				t1.used = 2;
				if ((res = mp_mul_d(&t1, q.dp[i - t - 1], &t1)) != MP_OKAY) {
					goto LBL_Y;
				}
	
				/* find right hand */
				t2.dp[0] = (i - 2 < 0) ? 0 : x.dp[i - 2];
				t2.dp[1] = (i - 1 < 0) ? 0 : x.dp[i - 1];
				t2.dp[2] = x.dp[i];
				t2.used = 3;
			} while (mp_cmp_mag(&t1, &t2) == MP_GT);
	
			/* step 3.3 x = x - q{i-t-1} * y * b**{i-t-1} */
			if ((res = mp_mul_d(&y, q.dp[i - t - 1], &t1)) != MP_OKAY) {
				goto LBL_Y;
			}
	
			if ((res = mp_lshd(&t1, i - t - 1)) != MP_OKAY) {
				goto LBL_Y;
			}
	
			if ((res = mp_sub(&x, &t1, &x)) != MP_OKAY) {
				goto LBL_Y;
			}
	
			/* if x < 0 then { x = x + y*b**{i-t-1}; q{i-t-1} -= 1; } */
			if (x.sign == MP_NEG) {
				if ((res = mp_copy(&y, &t1)) != MP_OKAY) {
					goto LBL_Y;
				}
				if ((res = mp_lshd(&t1, i - t - 1)) != MP_OKAY) {
					goto LBL_Y;
				}
				if ((res = mp_add(&x, &t1, &x)) != MP_OKAY) {
					goto LBL_Y;
				}
	
				q.dp[i - t - 1] = (q.dp[i - t - 1] - 1UL) & MP_MASK;
			}
		}
	
		/* now q is the quotient and x is the remainder
		* [which we have to normalize]
		*/
	
		/* get sign before writing to c */
		x.sign = x.used == 0 ? MP_ZPOS : a->sign;
	
		if (c != NULL) {
			mp_clamp(&q);
			mp_exch(&q, c);
			c->sign = neg;
		}
	
		if (d != NULL) {
			mp_div_2d(&x, norm, &x, NULL);
			mp_exch(&x, d);
		}
	
		res = MP_OKAY;
	
	LBL_Y:mp_clear(&y);
	LBL_X:mp_clear(&x);
	LBL_T2:mp_clear(&t2);
	LBL_T1:mp_clear(&t1);
	LBL_Q:mp_clear(&q);
		return res;
	}



	/* d = a + b (mod c) */
	int mp_addmod(mp_int * a, mp_int * b, mp_int * c, mp_int * d)
	{
		int     res;
		mp_int  t;

		if ((res = mp_init(&t)) != MP_OKAY) {
			return res;
		}

		if ((res = mp_add(a, b, &t)) != MP_OKAY) {
			mp_clear(&t);
			return res;
		}
		res = mp_mod(&t, c, d);
		mp_clear(&t);
		return res;
	}

	int mp_cnt_lsb(mp_int *a)
	{
		int x;
		mp_digit q, qq;

		/* easy out */
		if (mp_iszero(a) == 1) {
			return 0;
		}

		/* scan lower digits until non-zero */
		for (x = 0; x < a->used && a->dp[x] == 0; x++);
		q = a->dp[x];
		x *= DIGIT_BIT;

		/* now scan this digit until a 1 is found */
		if ((q & 1) == 0) {
			do {
				qq = q & 15;
				x += lnz[qq];
				q >>= 4;
			} while (qq == 0);
		}
		return x;
	}
	

int mp_sqrmod(mp_int * a, mp_int * b, mp_int * c)
{
	int     res;
	mp_int  t;

	if ((res = mp_init(&t)) != MP_OKAY) {
		return res;
	}

	if ((res = mp_sqr(a, &t)) != MP_OKAY) {
		mp_clear(&t);
		return res;
	}
	res = mp_mod(&t, b, c);
	mp_clear(&t);
	return res;
}
	
int mp_prime_miller_rabin(mp_int * a, mp_int * b, int *result)
	{
		mp_int  n1, y, r;
		int s;
		int j;
		int err;

		/* default */
		*result = MP_NO;

		/* ensure b > 1 */
		if (mp_cmp_d(b, 1) != MP_GT) {
			return MP_VAL;
		}

		/* get n1 = a - 1 */
		if ((err = mp_init_copy(&n1, a)) != MP_OKAY) {
			return err;
		}
		if ((err = mp_sub_d(&n1, 1, &n1)) != MP_OKAY) {
			goto LBL_N1;
		}

		/* set 2**s * r = n1 */
		if ((err = mp_init_copy(&r, &n1)) != MP_OKAY) {
			goto LBL_N1;
		}

		/* count the number of least significant bits
		* which are zero
		*/
		s = mp_cnt_lsb(&r);

		/* now divide n - 1 by 2**s */
		if ((err = mp_div_2d(&r, s, &r, NULL)) != MP_OKAY) {
			goto LBL_R;
		}

		/* compute y = b**r mod a */
		if ((err = mp_init(&y)) != MP_OKAY) {
			goto LBL_R;
		}
		if ((err = mp_exptmod(b, &r, a, &y)) != MP_OKAY) {
			goto LBL_Y;
		}

		/* if y != 1 and y != n1 do */
		if (mp_cmp_d(&y, 1) != MP_EQ && mp_cmp(&y, &n1) != MP_EQ) {
			j = 1;
			/* while j <= s-1 and y != n1 */
			while ((j <= (s - 1)) && mp_cmp(&y, &n1) != MP_EQ) {
				if ((err = mp_sqrmod(&y, a, &y)) != MP_OKAY) {
					goto LBL_Y;
				}

				/* if y == 1 then composite */
				if (mp_cmp_d(&y, 1) == MP_EQ) {
					goto LBL_Y;
				}

				++j;
			}

			/* if y != n1 then composite */
			if (mp_cmp(&y, &n1) != MP_EQ) {
				goto LBL_Y;
			}
		}

		/* probably prime now */
		*result = MP_YES;
	LBL_Y:mp_clear(&y);
	LBL_R:mp_clear(&r);
	LBL_N1:mp_clear(&n1);
		return err;
	}
	/* single digit subtraction */
	int mp_sub_d(mp_int * a, mp_digit b, mp_int * c)
	{
		mp_digit *tmpa, *tmpc, mu;
		int       res, ix, oldused;

		/* grow c as required */
		if (c->alloc < a->used + 1) {
			if ((res = mp_grow(c, a->used + 1)) != MP_OKAY) {
				return res;
			}
		}

		/* if a is negative just do an unsigned
		* addition [with fudged signs]
		*/
		if (a->sign == MP_NEG) {
			a->sign = MP_ZPOS;
			res = mp_add_d(a, b, c);
			a->sign = c->sign = MP_NEG;
			return res;
		}

		/* setup regs */
		oldused = c->used;
		tmpa = a->dp;
		tmpc = c->dp;

		/* if a <= b simply fix the single digit */
		if ((a->used == 1 && a->dp[0] <= b) || a->used == 0) {
			if (a->used == 1) {
				*tmpc++ = b - *tmpa;
			}
			else {
				*tmpc++ = b;
			}
			ix = 1;

			/* negative/1digit */
			c->sign = MP_NEG;
			c->used = 1;
		}
		else {
			/* positive/size */
			c->sign = MP_ZPOS;
			c->used = a->used;

			/* subtract first digit */
			*tmpc = *tmpa++ - b;
			mu = *tmpc >> (sizeof(mp_digit) * CHAR_BIT - 1);
			*tmpc++ &= MP_MASK;

			/* handle rest of the digits */
			for (ix = 1; ix < a->used; ix++) {
				*tmpc = *tmpa++ - mu;
				mu = *tmpc >> (sizeof(mp_digit) * CHAR_BIT - 1);
				*tmpc++ &= MP_MASK;
			}
		}

		/* zero excess digits */
		while (ix++ < oldused) {
			*tmpc++ = 0;
		}
		mp_clamp(c);
		return MP_OKAY;
	}

static const char rem_128[128] = {
 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1
};

static const char rem_105[105] = {
 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1,
 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1,
 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1,
 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1
};

/* Store non-zero to ret if arg is square, and zero if not */
int mp_is_square(mp_int *arg,int *ret) 
{
  int           res;
  mp_digit      c;
  mp_int        t;
  unsigned long r;

  /* Default to Non-square :) */
  *ret = MP_NO; 

  if (arg->sign == MP_NEG) {
    return MP_VAL;
  }

  /* digits used?  (TSD) */
  if (arg->used == 0) {
     return MP_OKAY;
  }

  /* First check mod 128 (suppose that DIGIT_BIT is at least 7) */
  if (rem_128[127 & DIGIT(arg,0)] == 1) {
     return MP_OKAY;
  }

  /* Next check mod 105 (3*5*7) */
  if ((res = mp_mod_d(arg,105,&c)) != MP_OKAY) {
     return res;
  }
  if (rem_105[c] == 1) {
     return MP_OKAY;
  }


  if ((res = mp_init_set_int(&t,11L*13L*17L*19L*23L*29L*31L)) != MP_OKAY) {
     return res;
  }
  if ((res = mp_mod(arg,&t,&t)) != MP_OKAY) {
     goto ERR;
  }
  r = mp_get_int(&t);
  /* Check for other prime modules, note it's not an ERROR but we must
   * free "t" so the easiest way is to goto ERR.  We know that res
   * is already equal to MP_OKAY from the mp_mod call 
   */ 
  if ( (1L<<(r%11)) & 0x5C4L )             goto ERR;
  if ( (1L<<(r%13)) & 0x9E4L )             goto ERR;
  if ( (1L<<(r%17)) & 0x5CE8L )            goto ERR;
  if ( (1L<<(r%19)) & 0x4F50CL )           goto ERR;
  if ( (1L<<(r%23)) & 0x7ACCA0L )          goto ERR;
  if ( (1L<<(r%29)) & 0xC2EDD0CL )         goto ERR;
  if ( (1L<<(r%31)) & 0x6DE2B848L )        goto ERR;

  /* Final check - is sqr(sqrt(arg)) == arg ? */
  if ((res = mp_sqrt(arg,&t)) != MP_OKAY) {
     goto ERR;
  }
  if ((res = mp_sqr(&t,&t)) != MP_OKAY) {
     goto ERR;
  }

  *ret = (mp_cmp_mag(&t,arg) == MP_EQ) ? MP_YES : MP_NO;
ERR:mp_clear(&t);
  return res;
}
	
int mp_read_radix (mp_int * a, const char *str, int radix)
{
  int     y, res, neg;
  char    ch;

  /* make sure the radix is ok */
  if (radix < 2 || radix > 64) {
    return MP_VAL;
  }

  /* if the leading digit is a 
   * minus set the sign to negative. 
   */
  if (*str == '-') {
    ++str;
    neg = MP_NEG;
  } else 
	{
    neg = MP_ZPOS;
  }

  /* set the integer to the default of zero */
  mp_zero (a);
  
  /* process each digit of the string */
  while (*str) 
	{
    /* if the radix < 36 the conversion is case insensitive
     * this allows numbers like 1AB and 1ab to represent the same  value
     * [e.g. in hex]
     */
    ch = (char) ((radix < 36) ? toupper (*str) : *str);
    for (y = 0; y < 64; y++) 
		{
      if (ch == mp_s_rmap[y]) 
			{
         break;
      }
    }

    /* if the char was found in the map 
     * and is less than the given radix add it
     * to the number, otherwise exit the loop. 
     */
    if (y < radix) 
		{
      if ((res = mp_mul_d(a, (mp_digit) radix, a)) != MP_OKAY) 
			{
         return res;
      }
      if ((res = mp_add_d (a, (mp_digit) y, a)) != MP_OKAY) 
			{
         return res;
      }
    } 
			else 
		{
      break;
    }
    ++str;
  }
  
  /* set the sign only if a != 0 */
  if (mp_iszero(a) != 1) {
     a->sign = neg;
  }
  return MP_OKAY;
}
	
unsigned long mp_get_int(mp_int * a) 
{
  int i;
  unsigned long res;

  if (a->used == 0) {
     return 0;
  }

  /* get number of digits of the lsb we have to read */
  i = MIN(a->used,(int)((sizeof(unsigned long)*CHAR_BIT+DIGIT_BIT-1)/DIGIT_BIT))-1;

  /* get most significant digit of result */
  res = DIGIT(a,i);
   
  while (--i >= 0) {
    res = (res << DIGIT_BIT) | DIGIT(a,i);
  }

  /* force result to 32-bits always so it is consistent on non 32-bit platforms */
  return res & 0xFFFFFFFFUL;
}
	
















//-----------------------------------------------------
//-----------------------------------------------------



//得到lon比特长素数
int GetPrime(mp_int *m,int lon);
//得到B和G点X坐标G点Y坐标
void Get_B_X_Y(mp_int *x1,mp_int *y1,mp_int *b,  mp_int *a,  mp_int *p);
//点乘
//void Ecc_points_mul(mp_int *qx,mp_int *qy, mp_int *px, mp_int *py,mp_int *d,mp_int *a,mp_int *p);
//点加
int Two_points_add(mp_int *x1,mp_int *y1,mp_int *x2,mp_int *y2,mp_int *x3,mp_int *y3,mp_int *a,mp_int *p);

//ECC加密
void Ecc_encipher(mp_int *qx,mp_int *qy, mp_int *px, mp_int *py,mp_int *a,mp_int *p);
//ECC解密
void Ecc_decipher(mp_int *k, mp_int *a,mp_int *p);

//生成密钥对
void SetKey(mp_int*QX,mp_int*QY,mp_int*GX,mp_int*GY,mp_int*sk,mp_int*n,mp_int*A,mp_int*B,mp_int*P);
//检验密钥
bool CheckKey(mp_int*QX,mp_int*QY,mp_int*sk,mp_int*n,mp_int*A,mp_int*B,mp_int*P);

int myrng(unsigned char *dst, int len, void *dat) //产生指定长度的十六进制随机数
{
   int x;
   for (x = 0; x < len; x++) dst[x] = rand() & 0xFF;
   return len;
}

char temp[200]={0};

char *Check(mp_int *mp)
{	
		mp_toradix(mp,temp,16);
		return temp;	
}

int mp_init(mp_int * a)
{
	int i =0 ;
	/* allocate memory required and clear it */
	a->dp = OPT_CAST(mp_digit) XMALLOC(sizeof(mp_digit) * MP_PREC);


	if (!(a->dp))
	{  //if (a->dp == NULL)
		printf("mp_init error!!");
		return MP_MEM;
	}

	/* set the digits to zero */
	for (i = 0; i < MP_PREC; i++) {
		a->dp[i] = 0;
	}

	/* set the used to zero, allocated digits to the default precision
	* and sign to positive */
	a->used = 0;
	a->alloc = MP_PREC;
	a->sign = MP_ZPOS;

    extern int initflag;
	initflag++;
	return MP_OKAY;
}

void Ecc_points_mul(mp_int *qx,mp_int *qy, mp_int *px, mp_int *py,mp_int *d,mp_int *a,mp_int *p)
{
	mp_int X, Y;	
	unsigned int i;
	char Bt_array[800];

    mp_toradix(d,Bt_array,2); 
 	
    mp_zero(qx);
    mp_zero(qy);
	mp_init(&X);
	mp_init(&Y);

	for(i=0;i<strlen(Bt_array);i++)
	{	
		
	   Two_points_add(qx,qy,qx,qy,&X,&Y,a,p);  	   			
       mp_copy(&X, qx);
	   mp_copy(&Y, qy);
	  
	   if(Bt_array[i]=='1')
	   {	   	
		  Two_points_add(qx,qy,px,py,&X,&Y,a,p);
		  
		  mp_copy(&X, qx);
		  mp_copy(&Y, qy);
	   }	  		 
	}
    
   mp_clear(&X);
   mp_clear(&Y); 
}

bool Testy(mp_int *temp2, mp_int *y1,mp_int *p)
{
		int ret;
		mp_is_square(temp2,&ret);
		if(ret)
		{
			mp_sqrt(temp2,y1);
			return true;
		}
		mp_add(temp2,p,temp2);
		mp_is_square(temp2,&ret);
		if(ret)
		{
			mp_sqrt(temp2,y1);
			return true;
		}
		mp_add(temp2,p,temp2);
		mp_is_square(temp2,&ret);
		if(ret)
		{
			mp_sqrt(temp2,y1);
			return true;
		}
		//printf("ret=%d\n",ret);  
	return false;
}

//选择基点，并计算其阶  	BasePoint(&GX,&GY,&A, &B, &P, &n);
void BasePoint(mp_int *x1,mp_int *y1,mp_int *a,  mp_int *b,  mp_int *p, mp_int *n)
{

	mp_int temp1;
	mp_int temp2;
	mp_int temp3;
	mp_int two;
	mp_int four;
	mp_int quotient;
	mp_int pm1;
	mp_int intMax;
	mp_int x2,y2,x3,y3;

	mp_init(&temp1);
	mp_init(&temp2);
	mp_init(&temp3);
	mp_init(&quotient);
	mp_init(&pm1);
	mp_init(&intMax);
	mp_init(&x2);
	mp_init(&y2);
	mp_init(&x3);
	mp_init(&y3);


	mp_copy(p,&intMax);
	mp_add(&intMax,p,&intMax);	
	mp_init_set(&two,2);
	mp_init_set(&four,4);
	mp_sub_d(p,1,&pm1); 		  
	while (mp_cmp(x1,&pm1)==-1)
	{ 
	   //mp_toradix(x1,temp,10);
	   //printf("here x1=%s\n",temp);  

		mp_expt_d(x1, 3, &temp1);
		mp_mul(a, x1, &temp2);
		mp_add(&temp1, &temp2, &temp3);
		mp_add(&temp3, b, &temp1);
		mp_mod(&temp1,p,&temp2);

		if (mp_cmp_d(&temp2,0)==0)
		{
			mp_zero(y1);
		}		
		else
		{  		 		

			/* d = a**b (mod c) */
			//int mp_exptmod(mp_int *a, mp_int *b, mp_int *c, mp_int *d);
			/* a/b => cb + d == a */
			//int mp_div(mp_int *a, mp_int *b, mp_int *c, mp_int *d);			
			
			mp_div(&pm1,&two,&quotient,&temp3);
			mp_exptmod(&temp2,&quotient,p,&temp1);
			if (mp_cmp_d(&temp1,1)==0)
			{
				//A^{(p-1)/2}=1
			
				/* c = a mod b, 0 <= c < b  */
				//int mp_mod(mp_int *a, mp_int *b, mp_int *c);
				printf("0000000000000!\n");
				mp_mod(p,&four,&temp1);
				print_mp(&temp1);
				if (mp_cmp_d(&temp1,3)==0)
				{				
					// p = 4k + 3 for some positive integer k 
					// output z := g^{k + 1} mod p.
					mp_add_d(p,1,&temp1);
					mp_div(&temp1,&four,&quotient,&temp3);
					mp_exptmod(&temp2,&quotient,p,y1);
					print_mp(y1);
				}	
				else
				{
				//	if(!Testy(&temp2,y1,p))
					{
						mp_add_d(x1,1,x1);
						continue;
					}
				}			
			}
			else
			{
				//if(!Testy(&temp2,y1,p))
				{
					mp_add_d(x1,1,x1);
					continue;
				}
			}
		}
	
		mp_init_set(n,1);

		/* copy a to  b */
		//int mp_copy(mp_int *a, mp_int *b);

		mp_copy(x1,&x2);		  
		mp_copy(y1,&y2);
	
		int result;
		
	mp_int tempx;mp_init(&tempx);
	mp_int tempy;mp_init(&tempy);

	mp_expt_d(x1,3,&tempx);             //  tempx = x1**3  
	mp_mod(&tempx,p,&tempx);			//  tempx = tempx mod P
	mp_mulmod(x1,a,p,&tempy);			//  tempy = x1 * A (mod P)
	mp_addmod(&tempx,&tempy,p,&tempx);	//  tempx = tempx + tempy (mod P)
	mp_addmod(&tempx,b,p,&tempx);		//  tempx = tempx + B (mod P) 
	mp_expt_d(y1,2,&tempy);				//  tempy = y1**2
	mp_mod(&tempy,p,&tempy);			//  tempy = tempy mod P

	if (mp_cmp(&tempx,&tempy))
	{
		printf("here2  bad\n\n");
		return;
	}
	printf("eee!\n");
	//printf("here good\n\n");
	//return;
	  while(1)
	  {
				/*	
			 mp_toradix(x1,temp,10);
			printf("here2 x1=%s\n",temp);  
			 mp_toradix(y1,temp,10);
			printf("here2 y1=%s\n",temp);  
			*/

			Two_points_add(x1,y1,&x2,&y2,&x3,&y3,a,p);
			mp_copy(&x3,&x2);
			mp_copy(&y3,&y2);
			mp_add_d(n,1,n);
			if (mp_cmp(n,&intMax)==1) 
			{
				mp_add_d(x1,1,x1);					
				break;
			}
			if (mp_cmp_d(&x2,0)==0 && mp_cmp_d(&y2,0)==0)
			{			
				/* Sets result to 1 if probably prime, 0 otherwise */
				// int mp_prime_is_prime(mp_int *a, int t, int *result);				
				mp_prime_is_prime(n,10,&result);

				if (result==1)
				{	
					mp_clear(&temp1);
					mp_clear(&temp2);
					mp_clear(&temp3);
					mp_clear(&quotient);
					mp_clear(&pm1);
					mp_clear(&intMax);
					mp_clear(&x2);
					mp_clear(&y2);
					mp_clear(&x3);
					mp_clear(&y3);
					return;
				}
				else
				{
					//mp_add_d(x1,1,x1);
					//break;
					mp_zero(n);  

					mp_clear(&temp1);
					mp_clear(&temp2);
					mp_clear(&temp3);
					mp_clear(&quotient);
					mp_clear(&pm1);
					mp_clear(&intMax);
					mp_clear(&x2);
					mp_clear(&y2);
					mp_clear(&x3);
					mp_clear(&y3);
					return;
				}
			}
	  }
	}	
}

void BasePoint_test(mp_int *x1,mp_int *y1,mp_int *a,  mp_int *b,  mp_int *p, mp_int *n)
{

	mp_int temp1;
	mp_int temp2;
	mp_int temp3;
	mp_int two;
	mp_int four;
	mp_int quotient;
	mp_int pm1;
	mp_int intMax;
	mp_int x2,y2,x3,y3;

	mp_init(&temp1);
	mp_init(&temp2);
	mp_init(&temp3);
	mp_init(&quotient);
	mp_init(&pm1);
	mp_init(&intMax);
	mp_init(&x2);
	mp_init(&y2);
	mp_init(&x3);
	mp_init(&y3);


	mp_copy(p,&intMax);
	mp_add(&intMax,p,&intMax);	
	mp_init_set(&two,2);
	mp_init_set(&four,4);
	mp_sub_d(p,1,&pm1); 		  
	while (mp_cmp(x1,&pm1)==-1)
	{ 
	   //mp_toradix(x1,temp,10);
	   //printf("here x1=%s\n",temp);  

		mp_expt_d(x1, 3, &temp1);
		mp_mul(a, x1, &temp2);
		mp_add(&temp1, &temp2, &temp3);
		mp_add(&temp3, b, &temp1);
		mp_mod(&temp1,p,&temp2);

		if (mp_cmp_d(&temp2,0)==0)
		{
			mp_zero(y1);
		}		
		else
		{  		 		

			/* d = a**b (mod c) */
			//int mp_exptmod(mp_int *a, mp_int *b, mp_int *c, mp_int *d);
			/* a/b => cb + d == a */
			//int mp_div(mp_int *a, mp_int *b, mp_int *c, mp_int *d);			
			
			mp_div(&pm1,&two,&quotient,&temp3);
			mp_exptmod(&temp2,&quotient,p,&temp1);
			if (mp_cmp_d(&temp1,1)==0)
			{
				//A^{(p-1)/2}=1
			
				/* c = a mod b, 0 <= c < b  */
				//int mp_mod(mp_int *a, mp_int *b, mp_int *c);
				printf("0000000000000!\n");
				mp_mod(p,&four,&temp1);
				print_mp(&temp1);
				if (mp_cmp_d(&temp1,3)==0)
				{				
					// p = 4k + 3 for some positive integer k 
					// output z := g^{k + 1} mod p.
					mp_add_d(p,1,&temp1);
					mp_div(&temp1,&four,&quotient,&temp3);
					mp_exptmod(&temp2,&quotient,p,y1);
					print_mp(y1);
					point_is_on_curve(x1,y1,a,b,p);
				}	
				else
				{
				//	if(!Testy(&temp2,y1,p))
					{
						mp_add_d(x1,1,x1);
						continue;
					}
				}			
			}
			else
			{
				//if(!Testy(&temp2,y1,p))
				{
					mp_add_d(x1,1,x1);
					continue;
				}
			}
			mp_add_d(x1,1,x1);
	}

    }	
}

//生成密钥对
void SetKey(mp_int*QX,mp_int*QY,mp_int*GX,mp_int*GY,mp_int*sk,mp_int*n,mp_int*A, mp_int*B,mp_int*P)
{
	
	bool r=false;
	
	mp_int nm1;mp_init(&nm1);
	mp_sub_d(n,1,&nm1);  //nm1=n - 1



	while (mp_cmp(sk,&nm1)==-1)  // 确保 sk < n-1
	{
		Ecc_points_mul(QX,QY,GX,GY,sk,A,P);//Q=sk*G	123567676......................
		//printf("QX=%s\n",Check(QX));
		r=CheckKey(QX,QY,sk,n,A,B,P);
		if (r)
		{
			return;
		}
		else
		mp_add_d(sk,1,sk);
			
	}
	mp_set(sk,0);
}
//检验密钥
bool CheckKey(mp_int*QX,mp_int*QY,mp_int*sk,mp_int*n,mp_int*A,mp_int*B,mp_int*P)
{
	if (mp_cmp_d(QX,0)==0 && mp_cmp_d(QY,0)==0)	//Q=O
	{
		return false;
	}
	if (mp_cmp_d(QX,0)==-1 || mp_cmp(QX,P)>=0) //QX<0||QX>=q
	{
		return false;
	}	
	if (mp_cmp_d(QY,0)==-1 || mp_cmp(QY,P)>=0) //QY<0||QY>=q
	{
		return false;
	}	
	mp_int tempx;mp_init(&tempx);
	mp_int tempy;mp_init(&tempy);
	mp_expt_d(QX,3,&tempx);
	mp_mod(&tempx,P,&tempx);
	/* d = a * b (mod c) */
	//int mp_mulmod(mp_int *a, mp_int *b, mp_int *c, mp_int *d);
	mp_mulmod(QX,A,P,&tempy);
	mp_addmod(&tempx,&tempy,P,&tempx);
	mp_addmod(&tempx,B,P,&tempx);
	mp_expt_d(QY,2,&tempy);
	mp_mod(&tempy,P,&tempy);
	if (mp_cmp(&tempx,&tempy))
	{
		return false;
	}
	mp_int rx;mp_init(&rx);
	mp_int ry;mp_init(&ry);
	Ecc_points_mul(&rx,&ry,QX,QY,n,A,P);//R=nQ=ndG
	if (mp_cmp_d(&rx,0)==0 && mp_cmp_d(&ry,0)==0)	//R=O
	{
		printf("G OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK!\n");
		return true;
	}
	else
	{
		return false;
	}
}

int GetPrime(mp_int *m,int lon){
   mp_prime_random_ex(m, 10, lon, 
		(rand()&1)?LTM_PRIME_2MSB_OFF:LTM_PRIME_2MSB_ON, myrng, NULL);
   return MP_OKAY;
}

//两点加
//调用Two_points_add 时，x3不能等于x1或x2,y3不能等于y1或y2 !!!
// Two_points_add(&X1,&Y1,&X1,&Y1,&X2,&Y2,&A,&P);  is OK		 
// Two_points_add(&X1,&Y1,px,py,&X1,&Y1,&A,&P);	is dangerous	
int Two_points_add(mp_int *x1,mp_int *y1,mp_int *x2,mp_int *y2,mp_int *x3,mp_int *y3,mp_int *a,mp_int *p)
{
mp_int x2x1;
mp_int y2y1;
mp_int tempk;
mp_int temp1;

mp_init(&x2x1);
mp_init(&y2y1);
mp_init(&tempk);
mp_init(&temp1);



	 if(mp_cmp_d(x1,0)==0 && mp_cmp_d(y1,0)==0)  
	   {
		 mp_copy(x2,x3);
		 mp_copy(y2,y3);   
		 goto L;
	   }
	 
	 if(mp_cmp_d(x2,0)==0 && mp_cmp_d(y2,0)==0)  
	   {
		 mp_copy(x1,x3);
		 mp_copy(y1,y3);   
		 goto L;
	   }

	
   mp_submod(x2, x1,p, &x2x1);
  // mp_toradix(&x2x1,temp,10);
  //  printf("sub x2x1=%s\n",temp); 
   if(mp_cmp_d(&x2x1,0)==0)
   {	  
		// P-P=0
	   mp_addmod(y2, y1,p, &y2y1);	    
	   if(mp_cmp_d(&y2y1,0)==0)  
	   {
		   mp_zero(x3); 
		   mp_zero(y3);
		   goto L;
	   }
	   else
	   { 	   //P+P=2P	   
		   mp_invmod(&y2y1,p,&tempk);		   
		   mp_expt_d(x1,2,&temp1);
		   mp_mul_d(&temp1,3,&temp1);
		   mp_addmod(&temp1,a,p,&temp1);
		   mp_mulmod(&temp1,&tempk,p,&tempk);		   
	   }
	   /*else
	   {   // x1=x2
		   //y1 +y2 !=0;
		   //y1 !=y2
			printf("Here is a question\n");
			return 0;
		} */
   }
   else
   {
	   //  printf("Here 1\n");   	   
	   mp_invmod(&x2x1,p,&tempk);	    
	   mp_submod(y2, y1,p, &y2y1);
	   mp_mulmod(&y2y1,&tempk,p,&tempk);	   
   }
   
   //printf("Here\n");
   

   mp_sqr(&tempk, &temp1);
   mp_sub(&temp1, x1, &temp1);   
   mp_submod(&temp1, x2,p, x3);    
   mp_sub(x1, x3, &temp1);   
   mp_mul(&temp1, &tempk, &temp1); 
   mp_submod(&temp1, y1, p, y3);   

L:

   mp_clear(&x2x1);
   mp_clear(&y2y1);
   mp_clear(&tempk);   
   mp_clear(&temp1);
  
   return 1;
}

//加密  Encrypt(&GX,&GY,&QXa,&QYa,&A,&P, &n, tmp,&ux, &uy, &s);
void Encrypt(mp_int *GX,mp_int *GY, mp_int *QX,mp_int *QY, mp_int *a, mp_int *p, mp_int *n,unsigned char *mstr, mp_int *rx,  mp_int *ry, mp_int *c)
{
	mp_int tempk;
	mp_int tempr;
	mp_int m;
	mp_int ff;
	mp_int tempx, tempy;
	unsigned char *pm=mstr;
	int k;	
	
	mp_init(&tempk);
	mp_init(&tempr);
	mp_init_set(&m,0);
	mp_init_set(&ff,256);
	mp_init(&tempx);
	mp_init(&tempy);
	
	// strlen(mstr)<=20
	while(*pm)
	{
		mp_mulmod(&m,&ff,p,&m);
		mp_add_d(&m,*(pm++),&m);
	}
	printf("m= %s\n",Check(&m));

	while (1)
	{	
		k=rand();//srand() is in main
		k = 0x1234567;
		mp_init_set( &tempk, k);
		mp_mod(&tempk,n,&tempk);	
		print_mp(n);
		print_mp(&tempk);
		if (mp_cmp_d(&tempk,0)==0)mp_add_d(&tempk,1,&tempk);
		Ecc_points_mul(rx,ry, GX, GY,&tempk,a,p);
		Ecc_points_mul(&tempx,&tempy, QX, QY,&tempk,a,p);
		printf("Point\n");
		print_mp(&tempx);
		print_mp(&tempy);
		printf("point\n");


		//mp_mod(&tempx,n,&tempx);

		if (mp_cmp_d(&tempx,0)==0)continue;
				
		mp_mulmod(&m,&tempx,p,c);
		print_mp(&tempx);	
		return;
	}	

}

//解密
bool Decrypt(mp_int *GX,mp_int *GY, mp_int *a, mp_int *b, mp_int *p,  mp_int *sk, unsigned char *mstr, mp_int *rx,  mp_int *ry, mp_int *c)
{ //Decrypt(&GX,&GY,&A,&B,&P, &sk,tmp2,&ux, &uy, &s)
	mp_int tempx;
	mp_int tempy;
	
	mp_int two;
	mp_int three;
	mp_init(&two);
	mp_init(&three);
	mp_set(&two,2);
	mp_set(&three,3);	

	mp_int m;
	
	mp_int ff; 

	unsigned char str[20];
	int count=0;
	mp_int tempr;

	mp_init(&tempx);
	mp_init(&m);
	mp_init(&tempy);
	mp_init(&tempr);
	mp_init_set(&ff,256);
	char temp[20];
	unsigned char *pm=mstr;

	mp_expt_d(rx,3,&tempx);
	mp_mod(&tempx,p,&tempx);
	/* d = a * b (mod c) */
	//int mp_mulmod(mp_int *a, mp_int *b, mp_int *c, mp_int *d);
	mp_mulmod(rx,a,p,&tempy);
	mp_addmod(&tempx,&tempy,p,&tempx);
	mp_addmod(&tempx,b,p,&tempx);
	//mp_expt_d(ry,2,&tempy);
	//mp_mod(&tempy,p,&tempy);
	mp_exptmod(ry,&two,p,&tempy);
	if (mp_cmp(&tempx,&tempy)) //&& mp_cmp_d(&n,0))
	{
		printf("Bad  R\n\n\n"); 		
	    printf("right=%s\n",Check(&tempx)); 
	    printf("left=%s\n",Check(&tempy));	
		return false;
	}
	print_mp(&tempx);
	Ecc_points_mul(&tempx,&tempy, rx,ry,sk,a,p);

	if (mp_cmp_d(&tempx,0)==0)return false;
	mp_invmod(&tempx,p,&tempx);
	mp_mulmod(c,&tempx,p,&m);
	print_mp(&m);	
	while(count<20)
	{
		/* a/b => cb + d == a */
		//int mp_div(mp_int *a, mp_int *b, mp_int *c, mp_int *d);
		mp_div(&m, &ff,&m,&tempr);
		print_mp(&tempr);
		mp_toradix(&tempr,temp,10);
		str[count]=atoi(temp);
	        printf("%d!\n",str[count]);
		count++;
		if (mp_cmp_d(&m,0)==0)break;
	}
	

	while(count>0)
	{
		count--;
		*(pm++)=str[count];
	}
	*(pm)=0;

	mp_clear(&tempx);
	mp_clear(&m);
	mp_clear(&tempy);
	mp_clear(&tempr);
	mp_clear(&ff);

	return true;
}

//签名
void Sign(mp_int *GX,mp_int *GY, mp_int *a, mp_int *p, mp_int *n,mp_int *sk, unsigned char *m, mp_int *ux,  mp_int *uy, mp_int *s)
{
	mp_int tempk;
	mp_int tempr;
	mp_int Hash;
	unsigned char *pm=m;
	int k;
	int count;
	int tm;

	mp_init(&tempk);
	mp_init(&tempr);
	mp_init_set(&Hash,0);
	
	count=0;tm=0;
	while(*pm)
	{
			tm=tm*256+*(pm++);
			count++;
			if(count==4)
			{				
				mp_add_d(&Hash,tm,&Hash);
				mp_mod(&Hash,n,&Hash);
				count=0;tm=0;
			}
	}
	printf("hash: %s\n",Check(&Hash));

	while (1)
	{	
		k=rand();//srand() is in main
		mp_init_set( &tempk, k);
		mp_mod(&tempk,n,&tempk);
		if (mp_cmp_d(&tempk,0)==0)mp_add_d(&tempk,1,&tempk);
		Ecc_points_mul(ux,uy, GX, GY,&tempk,a,p);
		mp_mod(ux,n,&tempr);

		if (mp_cmp_d(&tempr,0)==0)continue;
		
		//hash(m)+sk*r(mod n)
		mp_mulmod(&tempr,sk,n,&tempr);
	
		mp_addmod(&Hash,&tempr,n,&tempr);
		if (mp_cmp_d(&tempr,0)==0)continue;

		mp_invmod(&tempk,n,&tempk);
		mp_mulmod(&tempk,&tempr,n,s);	
		return;
	}	
}
//认证

bool Authentication(mp_int *GX,mp_int *GY,mp_int *QX,mp_int *QY, mp_int *a, mp_int *p, mp_int *n, unsigned char *m, mp_int *ux,  mp_int *uy, mp_int *s)
{
	
	mp_int Hash;
	mp_int plusx;	
	mp_int plusy;
	mp_int tempgx;	
	mp_int tempgy;
	mp_int tempr;
	mp_int tempqx;	
	mp_int tempqy;

	mp_init(&Hash);
	mp_init(&tempr);
	mp_init(&tempgx);
	mp_init(&tempgy);
	mp_init(&tempqx);
	mp_init(&tempqy);
	mp_init(&plusx);
	mp_init(&plusy);

	unsigned char *pm=m;
	int count;
	int tm;
	count=0;tm=0;
	while(*pm)
	{
			tm=tm*256+*(pm++);
			count++;
			if(count==4)
			{
				mp_add_d(&Hash,tm,&Hash);
				mp_mod(&Hash,n,&Hash);
				count=0;tm=0;
			}
	}
	printf("h2 %s\n",Check(&Hash));


	if (mp_cmp_d(s,0)==0)return false;
	if (mp_cmp(s,n)>-1)return false;
	
	//sU=Hash[m]G+rQ
	mp_int tempux;	mp_int tempuy;
	mp_init(&tempux);mp_init(&tempuy);
	Ecc_points_mul(&tempux,&tempuy, ux, uy,s,a,p);

	Ecc_points_mul(&tempgx,&tempgy, GX, GY,&Hash,a,p);
	
	mp_mod(ux,n,&tempr);
	Ecc_points_mul(&tempqx,&tempqy, QX, QY,&tempr,a,p);

	Two_points_add(&tempgx,&tempgy,&tempqx,&tempqy,&plusx,&plusy,a,p);

	if (mp_cmp(&tempux,&plusx))return false;	

	if (mp_cmp(&tempuy,&plusy))return false;
	
	mp_clear(&Hash);
	mp_clear(&tempr);
	mp_clear(&tempgx);
	mp_clear(&tempgy);
	mp_clear(&tempqx);
	mp_clear(&tempqy);
	mp_clear(&plusx);
	mp_clear(&plusy);

	return true;
}


void FunctionCheck(mp_int GX,mp_int GY,mp_int P,mp_int A,mp_int B)
{
	mp_int tempx;mp_init(&tempx);
	mp_int tempy;mp_init(&tempy);
	mp_expt_d(&GX,3,&tempx);              //  tempx = GX**3  
	mp_mod(&tempx,&P,&tempx);			  //  tempx = tempx mod P
	mp_mulmod(&GX,&A,&P,&tempy);          //  tempy = GX * A (mod P)
	mp_addmod(&tempx,&tempy,&P,&tempx);   //  tempx = tempx + tempy (mod P) 
	mp_addmod(&tempx,&B,&P,&tempx);       //  tempx = tempx + B (mod P) 

	mp_expt_d(&GY,2,&tempy);              //  tempy = GY**2
	mp_mod(&tempy,&P,&tempy);             //  tempy = tempy mod P

	char temp[200]={0};
	if (mp_cmp(&tempx,&tempy))  //&& mp_cmp_d(&n,0))   验证 y**2 = x**3 + a*x + b 
	{  
		printf("after Bad\n\n\n"); 
		mp_toradix(&tempx,temp,16);
		printf("%s\n",temp); 
		mp_toradix(&tempy,temp,16);
		printf("%s\n",temp); 
	}

	printf("检查完成\n");

	mp_clear(&tempx);
	mp_clear(&tempy);
	
}





void print_mp(mp_int *mp)
{
         char *buf = NULL;
         static unsigned long i;
         
         buf = malloc(1024);
         if(!buf){
                 printf("Alloc mem failed!\n");
                 return;
         }
	 memset(buf,0,1024);
         i++;    
         
         mp_toradix(mp,buf,16);
         printf("Id : %d [%s]\n",i,buf);
         free(buf);                 
}

int initflag = 0 ;
int lon=0;

void main(){

	mp_int GX;
	mp_int GY;
	mp_int A;
	mp_int B;
	mp_int P;//Fp中的p(有限域P)
	mp_int n; 
	mp_int ska;//a 私有密钥
	mp_int skb;//b 私有密钥
	mp_init(&GX);// Here must be initialized !
	mp_init(&GY);
	mp_init(&A);
	mp_init(&B);
	mp_init(&P);
	mp_init(&n);
	mp_init(&ska);
	mp_init(&skb);
	
	int res;
	printf("sizeof mp_word is %d!\n",sizeof(mp_word));
	printf("sizeof mp_digit is %d!\n",sizeof(mp_digit));
	printf("DIGIT_BIT : %d!\n",DIGIT_BIT);
	printf("MP_MASK : %lx!\n",MP_MASK);
	/*char p="FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFF";
	char a=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFC;
	char b=28E9FA9E 9D9F5E34 4D5A9E4B CF6509A7 F39789F5 15AB8F92 DDBCBD41 4D940E93;
	char n=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF 7203DF6B 21C6052B 53BBF409 39D54123;
	char Gx=32C4AE2C 1F198119 5F990446 6A39C994 8FE30BBF F2660BE1 715A4589 334C74C7;
	char Gy=BC3736A2 F4F6779C 59BDCEE3 6B692153 D0A9877C C62A4740 02DF32E5 2139F0A0;
	char k=1649AB77 A00637BD 5E2EFE28 3FBF3535 34AA7F7C B89463F2 08DDBC29 20BB0DA0;*/

    time_t t;           
    srand( (unsigned) time( &t ) );

    printf("椭圆曲线的参数如下(以十六进制显示):\n");	
	//char temp[200]={0};
	//extern char temp[200];
	
   // GetPrime(&P,8); // 素数 

	int lon=0;
	printf("有限域 P 是:\n");			
	memcpy(temp,"FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF",64);	
	mp_read_radix(&P,temp,16);
	printf("%s\n",Check(&P));

	lon=lon+sizeof P;
	printf("累计长度%d\n",lon);

	//mp_read_radix(&P,temp,16);//保留

	memcpy(temp,"FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC",64);	
	mp_read_radix(&A,temp,16);
	printf("曲线参数 A 是:\n");	
    //mp_toradix(&A,temp,16);
    printf("%s\n",Check(&A));  
	lon=lon+sizeof A;
	printf("累计长度%d\n",lon);
	
	memcpy(temp,"28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93",64);
	mp_read_radix(&B,temp,16);
	printf("曲线参数 B 是:\n");	
    printf("%s\n",Check(&B));  
	lon=lon+sizeof B;
	printf("累计长度%d\n",lon);
	
	 
	printf("曲线G点X坐标是:\n");	
	memcpy(temp,"32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7",64);
	mp_read_radix(&GX,temp,16);
	//BasePoint(&GX,&GY,&A, &B, &P, &n);

	printf("曲线G点X坐标是:\n");	
	memcpy(temp,"32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7",64);
	mp_read_radix(&GX,temp,16);
	printf("%s\n",Check(&GX));  

	lon=lon+sizeof GX;
	printf("累计长度%d\n",lon);

	printf("曲线G点Y坐标是:\n");
	memcpy(temp,"BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0",64);
	mp_read_radix(&GY,temp,16);
	printf("%s\n",Check(&GY));  

	lon=lon+sizeof GY;
	printf("累计长度%d\n",lon);



	//printf("当前init函数被调用了%d次，占用大小为%d K\n",initflag,initflag*256/1024);
	//strcpy(temp,"9");
	//p_read_radix(&n,temp,10);
	//printf("%s\n",Check(&n));  
	//Ecc_points_mul(&QX,&QY,&GX,&GY,&n,&A,&P);
	 //printf("QX=%s\n",Check(&QX));  
	 //printf("QY=%s\n",Check(&QY));  


	mp_int tempx;mp_init(&tempx);
	mp_int tempy;mp_init(&tempy);
	res = mp_expt_d(&GX,3,&tempx);              //  tempx = GX**3
	if(res != MP_OKAY)	
		printf("MP_expt_d err:%d!\n",res);
	print_mp(&tempx); 
	res = mp_mod(&tempx,&P,&tempx);			  //  tempx = tempx mod P
	if(res != MP_OKAY)
		printf("MP_MOD err:%d!\n",res);
	print_mp(&tempx);
	res = mp_mulmod(&GX,&A,&P,&tempy);          //  tempy = GX * A (mod P)
	if(res != MP_OKAY)
		printf("MP_MULMOD err:%d!\n",res);
	print_mp(&tempy);
	res = mp_addmod(&tempx,&tempy,&P,&tempx);   //  tempx = tempx + tempy (mod P)
	if(res != MP_OKAY)
		printf("mp_addmod err:%d!\n",res);
	print_mp(&tempx); 
	res = mp_addmod(&tempx,&B,&P,&tempx);       //  tempx = tempx + B (mod P)
	if(res != MP_OKAY)
		printf("mp_addmod err:%d!\n",res); 
	print_mp(&tempx);

	res = mp_expt_d(&GY,2,&tempy);              //  tempy = GY**2
	if(res != MP_OKAY)
		printf("mp_expt_d err:%d!\n",res);
	print_mp(&tempy);
	res = mp_mod(&tempy,&P,&tempy);             //  tempy = tempy mod P
	if(res != MP_OKAY)
		printf("mp_mod err:%d!\n",res);
	print_mp(&tempy);

	//char temp[200]={0};
	if (mp_cmp(&tempx,&tempy))  //&& mp_cmp_d(&n,0))   验证 y**2 = x**3 + a*x + b 
	{  
		printf("after Bad\n\n\n"); 
		mp_toradix(&tempx,temp,16);
		printf("%s\n",temp); 
		mp_toradix(&tempy,temp,16);
		printf("%s\n",temp); 
	}

	printf("检查完成\n");
	
	printf("G点的阶n是:\n");
	memcpy(temp,"FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123",64);	
	mp_read_radix(&n,temp,16);

	//printf("GX %s\n",Check(&GX));
	//printf("GY %s\n",Check(&GY));

	printf("基点阶N%s\n",Check(&n));  

	int result;
	mp_prime_is_prime(&n,10,&result);	
	if (result==1)
	printf("n is prime\n"); 
	
//	GetPrime(&ska,256);
	memcpy(temp,"944E16955CE909023431EDEE27B5A936E48D3C858530F7BF2CE645ACF2C2A287",64);
	mp_read_radix(&ska,temp,16);

	printf("Alice私钥是：%s\n",Check(&ska));
	lon=lon+sizeof ska;
	printf("累计长度%d\n",lon);

	mp_int QAX;
	mp_int QAY;
	mp_init(&QAX);
	mp_init(&QAY);
	SetKey(&QAX,&QAY,&GX,&GY,&ska,&n,&A,&B,&P);

	//GetPrime(&skb,256);
	memcpy(temp,"944E16955CE909023431EDEE27B5A936E48D3C858530F7BF2CE645ACF2C2A287",64);
	mp_read_radix(&skb,temp,16);
	printf("Bob私钥是  ：%s\n",Check(&skb));
	lon=lon+sizeof skb;
	printf("累计长度%d\n",lon);

	mp_int QBX;
	mp_int QBY;
	mp_init(&QBX);
	mp_init(&QBY);
	SetKey(&QBX,&QBY,&GX,&GY,&skb,&n,&A,&B,&P);
	
	//交换后 Alice
	mp_int QXa;
	mp_int QYa;
	mp_init(&QXa);
	mp_init(&QYa);
	Ecc_points_mul(&QXa,&QYa,&QBX,&QBY,&ska,&A,&P);

	//交换后 bob
	mp_int QXb;
	mp_int QYb;
	mp_init(&QXb);
	mp_init(&QYb);
	Ecc_points_mul(&QXb,&QYb,&QAX,&QAY,&skb,&A,&P);

	printf("A公钥X坐标是:\n");
	printf("%s\n",Check(&QXa));  

	printf("B公钥X坐标是:\n");
	printf("%s\n",Check(&QXb));  

	printf("A公钥Y坐标是:\n");
	printf("%s\n",Check(&QYa));  


	printf("B公钥Y坐标是:\n");
	printf("%s\n",Check(&QYb));  


	mp_int ux;mp_init(&ux);
	mp_int uy;mp_init(&uy);
	mp_int s;mp_init(&s);

	unsigned char m[]="中天安泰";

	Sign(&GX,&GY,&A,&P, &n,&ska,m,&ux, &uy, &s);
	printf("ux=%s\n",Check(&ux));  
	
	printf("uy=%s\n",Check(&uy));  

	printf("s=%s\n",Check(&s));  
	if(Authentication(&GX,&GY,&QAX,&QAY,&A,&P, &n,m,&ux, &uy, &s))	  
		printf("Authentication OK\n"); 

	Ecc_points_mul(&QXa,&QYa,&GX,&GY,&ska,&A,&P);  //alice

	// Ecc_points_mul(&QXb,&QYb,&GX,&GY,&skb,&A,&P);  //bob

	//原为	Ecc_points_mul(&QX,&QY,&GX,&GY,&sk,&A,&P);



	unsigned  char teststr[1024]={"zhaoziyuantest"};
	unsigned char tmp[21], tmp2[21];
	unsigned char tmp3[200];
	unsigned char *pstr=teststr;
	unsigned char *pstr2=tmp3;
	int i;

	clock_t start,mid,end;
	int CLK_TCK;
	struct tms tms;
	
	CLK_TCK = sysconf(_SC_CLK_TCK);
	printf("CLK_TCK:%d!\n",CLK_TCK);
	
	start = times(&tms);

	printf("密文为：\n%s\n",teststr); 

	//printf("GetPrime后当前init函数被调用了%d次，占用大小为%d K\n",initflag,initflag*256/1024);
	tmp3[0]=0;
	while(*pstr)
	{
		for(i=0;i<20 && *pstr;i++)
			tmp[i]= *(pstr++);
		tmp[i]=0;
		Encrypt(&GX,&GY,&QXa,&QYa,&A,&P, &n, tmp,&ux, &uy, &s);   //使用 alice 密钥
			printf("rx=%s\n",Check(&ux));  
				printf("ry=%s\n",Check(&uy));  
					printf("c=%s\n",Check(&s));  //密文

	 

		mid = times(&tms);			
		
		if(Decrypt(&GX,&GY,&A,&B,&P, &ska,tmp2,&ux, &uy, &s))
			for(i=0;i<20 && tmp2[i];i++)*(pstr2++)= tmp2[i]; 
				printf("tmp2=%s\n",tmp2);  //解密后
			
		
	}
	*(pstr2)=0;
		printf("\n--- 解密后---\n");
	printf("Plain text=\n%s\n",tmp3); 


	end = times(&tms);
	printf("加密时间：%2.1f\n解密时间：%2.1f\n总用时:%2.1f\n",(double)(mid-start),(double)(end-mid),(double)(end-start));


	//getchar();
	//getchar();
	mp_clear(&GX);
	mp_clear(&GY);
	mp_clear(&ska);//私有密钥
	mp_clear(&skb);//私有密钥
	mp_clear(&A);
	mp_clear(&B);
	mp_clear(&QAX);
	mp_clear(&QAY);
	mp_clear(&QBX);
	mp_clear(&QBY);
	mp_clear(&P);//Fp中的p(有限域P)
	
	return;
	
}

bool point_is_on_curve(mp_int *rx,mp_int *ry,mp_int *a,mp_int *b,mp_int *p)
{
	bool res = true;
	mp_int tempx;
	mp_int tempy;

	mp_init(&tempx);
	mp_init(&tempy);

	mp_expt_d(rx,3,&tempx);
	mp_mod(&tempx,p,&tempx);
	mp_mulmod(rx,a,p,&tempy);
	mp_addmod(&tempx,&tempy,p,&tempx);
	mp_addmod(&tempx,b,p,&tempx);
	mp_expt_d(ry,2,&tempy);
	mp_mod(&tempy,p,&tempy);
	if (mp_cmp(&tempx,&tempy)) //&& mp_cmp_d(&n,0))
	{
		printf("Bad  R\n\n\n"); 		
		printf("right=%s\n",Check(&tempx)); 
	        printf("left=%s\n",Check(&tempy));
		res = false;
		goto ret;
	}

	res = true;
	printf("G point generarte ok!\n");
ret:
	mp_clear(&tempx);
	mp_clear(&tempy);

	return res;

}
/*-------------------------------------------------------------------------*/
