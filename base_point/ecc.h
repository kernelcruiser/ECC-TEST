#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

#undef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#undef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))


#define BN_MP_TOOM_SQR_C
#define BN_MP_KARATSUBA_SQR_C
#define BN_FAST_S_MP_SQR_C
#define BN_S_MP_SQR_C
#define BN_MP_TOOM_MUL_C
//#define BN_MP_KARATSUBA_MUL_C
#define BN_FAST_S_MP_MUL_DIGS_C
#define BN_S_MP_MUL_DIGS_C
#define BN_FAST_MP_INVMOD_C  //????
#define BN_MP_INVMOD_SLOW_C  //????
#define BN_MP_DIV_3_C
#define BN_MP_INVMOD_C  //????????
#define BN_MP_DR_IS_MODULUS_C
#define BN_MP_REDUCE_IS_2K_C
#define BN_MP_EXPTMOD_FAST_C
#define BN_S_MP_EXPTMOD_C
#define BN_MP_EXPTMOD_FAST_C

#define BN_MP_EXPTMOD_FAST_C
#define BN_S_MP_EXPTMOD_C
#define BN_MP_MONTGOMERY_SETUP_C
#define BN_MP_MONTGOMERY_REDUCE_C
#define BN_MP_DR_SETUP_C
#define BN_MP_DR_REDUCE_C
#define BN_MP_REDUCE_2K_SETUP_C
#define BN_MP_REDUCE_2K_C
#define BN_MP_MONTGOMERY_CALC_NORMLIZATION_C

//ÄÚ´æ¹ÜÀíº¯Êý¿ªÊ¼

typedef unsigned long  u32;
typedef unsigned short u16;
typedef unsigned char  u8;   
#ifndef NULL
#define NULL 0
#endif 

//#define MEM_BLOCK_SIZE   32          //??????32??
//#define MAX_MEM_SIZE   10*1024       //?????? 10K
//#define MEM_ALLOC_TABLE_SIZE MAX_MEM_SIZE/MEM_BLOCK_SIZE //?????
////???????
//struct _m_mallco_dev
//{
// void (*init)(void);     //???
// u8 (*perused)(void);         //?????
// u8  membase[MAX_MEM_SIZE];   //???
// u16 memmap[MEM_ALLOC_TABLE_SIZE];  //???????
// u8  memrdy;        //????????
//};
//extern struct _m_mallco_dev mallco_dev;  //?mallco.c????
//void mymemset(void *s,u8 c,u32 count);  //????
//void mymemcpy(void *des,void *src,u32 n);//???? 
//void mem_init(void);      //?????????(?/????)
//u32 mem_malloc(u32 size);     //????(????)
//u8 mem_free(u32 offset);     //????(????)
//u8 mem_perused(void);      //???????(?/????) 
//////////////////////////////////////////////////////////////////////////////////
////??????
//void myfree(void *ptr);       //????(????)
//void *mymalloc(u32 size);     //????(????)
//void *myrealloc(void *ptr,u32 size);  //??????(????)



//ÄÚ´æ¹ÜÀíº¯Êý½áÊø


#define  OPT_CAST(x)
/* detect 64-bit mode if possible */
#if defined(__x86_64__) 
#if !(defined(MP_64BIT) && defined(MP_16BIT) && defined(MP_8BIT))
#define MP_64BIT
#endif
#endif

#ifndef NULL
#ifdef __cplusplus
#define NULL    0
#else
#define NULL    ((void *)0)
#endif
#endif

#define bool char 
#define true 1
#define false 0

//#define CHAR_BIT      8         /* number of bits in a char */

#ifdef MP_8BIT
typedef unsigned char      mp_digit;
typedef unsigned short     mp_word;
#elif defined(MP_16BIT)
typedef unsigned short     mp_digit;
typedef unsigned long      mp_word;
#elif defined(MP_64BIT)
/* for GCC only on supported platforms */
#ifndef CRYPT
typedef unsigned long long ulong64;
typedef signed long long   long64;
#endif

typedef unsigned long      mp_digit;
//typedef unsigned long      mp_word __attribute__ ((mode(TI)));
typedef __uint128_t      mp_word;

#define DIGIT_BIT          60
#else
/* this is the default case, 28-bit digits */

/* this is to make porting into LibTomCrypt easier :-) */
#ifndef CRYPT
#if defined(_MSC_VER) || defined(__BORLANDC__) 
typedef unsigned __int64   ulong64;
typedef signed __int64     long64;
#else
typedef unsigned long long ulong64;
typedef signed long long   long64;
#endif
#endif

typedef unsigned long      mp_digit;
typedef ulong64            mp_word;

#ifdef MP_31BIT   
/* this is an extension that uses 31-bit digits */
#define DIGIT_BIT          31
#else
/* default case is 28-bit digits, defines MP_28BIT as a handy macro to test */
#define DIGIT_BIT          28
#define MP_28BIT
#endif   
#endif

/* number of primes */
#ifdef MP_8BIT
#define PRIME_SIZE      31
#else
#define PRIME_SIZE      256
#endif

/* define heap macros */
#ifndef CRYPT
/* default to libc stuff */
#ifndef XMALLOC 
#define XMALLOC  malloc
#define XFREE    free
#define XREALLOC realloc
#define XCALLOC  calloc
#else
/* prototypes for our heap functions */
extern void *XMALLOC(size_t n);
extern void *REALLOC(void *p, size_t n);
extern void *XCALLOC(size_t n, size_t s);
extern void XFREE(void *p);
#endif
#endif


/* otherwise the bits per digit is calculated automatically from the size of a mp_digit */
#ifndef DIGIT_BIT
#define DIGIT_BIT     ((int)((CHAR_BIT * sizeof(mp_digit) - 1)))  /* bits per digit */
#endif
//#define CHAR_BIT  8  //zzy 2016-11-4
#define MP_DIGIT_BIT     DIGIT_BIT
#define MP_MASK          ((((mp_digit)1)<<((mp_digit)DIGIT_BIT))-((mp_digit)1))
#define MP_DIGIT_MAX     MP_MASK

/* equalities */
#define MP_LT        -1   /* less than */
#define MP_EQ         0   /* equal to */
#define MP_GT         1   /* greater than */

#define MP_ZPOS       0   /* positive integer */
#define MP_NEG        1   /* negative */

#define MP_OKAY       0   /* ok result */
#define MP_MEM        -2  /* out of mem */
#define MP_VAL        -3  /* invalid input */
#define MP_RANGE      MP_VAL

#define MP_YES        1   /* yes response */
#define MP_NO         0   /* no response */

/* Primality generation flags */
#define LTM_PRIME_BBS      0x0001 /* BBS style prime */
#define LTM_PRIME_SAFE     0x0002 /* Safe prime (p-1)/2 == prime */
#define LTM_PRIME_2MSB_OFF 0x0004 /* force 2nd MSB to 0 */
#define LTM_PRIME_2MSB_ON  0x0008 /* force 2nd MSB to 1 */

typedef int           mp_err;

#define MP_MASK          ((((mp_digit)1)<<((mp_digit)DIGIT_BIT))-((mp_digit)1))

#define mp_iszero(a) (((a)->used == 0) ? MP_YES : MP_NO)
#define mp_iseven(a) (((a)->used > 0 && (((a)->dp[0] & 1) == 0)) ? MP_YES : MP_NO)
#define mp_isodd(a)  (((a)->used > 0 && (((a)->dp[0] & 1) == 1)) ? MP_YES : MP_NO)



extern int KARATSUBA_MUL_CUTOFF,
	KARATSUBA_SQR_CUTOFF,
	TOOM_MUL_CUTOFF,
	TOOM_SQR_CUTOFF;

/* define this to use lower memory usage routines (exptmods mostly) */
/* #define MP_LOW_MEM */

/* default precision */
#ifndef MP_PREC
#ifndef MP_LOW_MEM
#define MP_PREC                 64     /* default digits of precision */
#else
#define MP_PREC                 8      /* default digits of precision */
#endif   
#endif

/* size of comba arrays, should be at least 2 * 2**(BITS_PER_WORD - BITS_PER_DIGIT*2) */
#define MP_WARRAY               (1 << (sizeof(mp_word) * CHAR_BIT - 2 * DIGIT_BIT + 1))

/* the infamous mp_int structure */
typedef struct{
	int used, alloc, sign;
	mp_digit *dp;
}mp_int;
//mp_int mp_obj;


static const int lnz[16] = { 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0 };

/* callback for mp_prime_random, should fill dst with random bytes and return how many read [upto len] */
typedef int ltm_prime_callback(unsigned char *dst, int len, void *dat);


//_Check_return_ _CRT_JIT_INTRINSIC _CRTIMP int __cdecl toupper(_In_ int _C);





#define USED(m)    ((m)->used)
#define DIGIT(m,k) ((m)->dp[(k)])
#define SIGN(m)    ((m)->sign)



#define BIT_LEN 800 
#define KEY_LONG 128  //?????
#define P_LONG 200    //???P???
#define EN_LONG 40    //????????(x,20)(y,20)

#define  Max 0xfff

//??lon?????
int GetPrime(mp_int *m,int lon);
//??B?G?X??G?Y??
void Get_B_X_Y(mp_int *x1,mp_int *y1,mp_int *b,  mp_int *a,  mp_int *p);
//??
//void Ecc_points_mul(mp_int *qx,mp_int *qy, mp_int *px, mp_int *py,mp_int *d,mp_int *a,mp_int *p);
//??
int Two_points_add(mp_int *x1,mp_int *y1,mp_int *x2,mp_int *y2,mp_int *x3,mp_int *y3,mp_int *a,mp_int *p);

//ECC??
void Ecc_encipher(mp_int *qx,mp_int *qy, mp_int *px, mp_int *py,mp_int *a,mp_int *p);
//ECC??
void Ecc_decipher(mp_int *k, mp_int *a,mp_int *p);

//?????
void SetKey(mp_int*QX,mp_int*QY,mp_int*GX,mp_int*GY,mp_int*sk,mp_int*n,mp_int*A,mp_int*B,mp_int*P);
//????
bool CheckKey(mp_int*QX,mp_int*QY,mp_int*sk,mp_int*n,mp_int*A,mp_int*B,mp_int*P);
void mp_clamp(mp_int * a);

void mp_rshd(mp_int * a, int b);
void bn_reverse(unsigned char *s, int len);
void mp_zero(mp_int * a);
void mp_clear(mp_int * a);
int mp_init_copy(mp_int * a, mp_int * b);
void mp_exch(mp_int * a, mp_int * b);
int mp_read_unsigned_bin(mp_int * a, const unsigned char *b, int c);
int mp_prime_is_prime(mp_int * a, int t, int *result);
int mp_sub_d(mp_int * a, mp_digit b, mp_int * c);

unsigned long mp_get_int(mp_int * a) ;

void mp_set(mp_int * a, mp_digit b);
int mp_toradix(mp_int * a, char *str, int radix);
int mp_init(mp_int * a);
int mp_init_set_int(mp_int * a, unsigned long b);
int mp_init_set(mp_int * a, mp_digit b);
int mp_prime_random_ex(mp_int *a, int t, int size, int flags, ltm_prime_callback cb, void *dat);
int mp_expt_d(mp_int * a, mp_digit b, mp_int * c);
int mp_sqr(mp_int * a, mp_int * b);
int mp_sqrt(mp_int *arg, mp_int *ret);
int mp_mul_d(mp_int * a, mp_digit b, mp_int * c);
int mp_add(mp_int * a, mp_int * b, mp_int * c);
int mp_sub(mp_int * a, mp_int * b, mp_int * c);
int mp_mod(mp_int * a, mp_int * b, mp_int * c);
int mp_cmp(mp_int * a, mp_int * b);
void mp_clear(mp_int * a);
int mp_copy(mp_int * a, mp_int * b);
int mp_init_copy(mp_int * a, mp_int * b);
void mp_zero(mp_int * a);
int mp_invmod(mp_int * a, mp_int * b, mp_int * c);
int mp_mulmod(mp_int * a, mp_int * b, mp_int * c, mp_int * d);

int mp_neg(mp_int * a, mp_int * b);
int mp_submod(mp_int * a, mp_int * b, mp_int * c, mp_int * d);
int mp_lshd(mp_int * a, int b);

int mp_read_unsigned_bin(mp_int * a, const unsigned char *b, int c);

void bn_reverse(unsigned char *s, int len);
int mp_cmp_d(mp_int * a, mp_digit b);

int mp_mod_2d(mp_int * a, int b, mp_int * c);
int mp_mul_2d(mp_int * a, int b, mp_int * c);
static int s_is_power_of_two(mp_digit b, int *p);
void mp_exch(mp_int * a, mp_int * b);
int mp_grow(mp_int * a, int size);
int mp_init_size(mp_int * a, int size);
void mp_clamp(mp_int * a);
int mp_mod_d(mp_int * a, mp_digit b, mp_digit * c);                                  
int mp_div_2d(mp_int * a, int b, mp_int * c, mp_int * d);
void mp_rshd(mp_int * a, int b);
int mp_div_d(mp_int * a, mp_digit b, mp_int * c, mp_digit * d);
int mp_cmp_mag(mp_int * a, mp_int * b);
int mp_prime_is_divisible(mp_int * a, int *result);
int mp_prime_is_prime(mp_int * a, int t, int *result);
int mp_exptmod(mp_int * G, mp_int * X, mp_int * P, mp_int * Y);
int mp_count_bits(mp_int * a);

int mp_mul_2(mp_int * a, mp_int * b);
int mp_div_2(mp_int * a, mp_int * b);
int s_mp_add(mp_int * a, mp_int * b, mp_int * c);
int s_mp_sub(mp_int * a, mp_int * b, mp_int * c);
int mp_init_multi(mp_int *mp, ...);
int mp_abs(mp_int * a, mp_int * b);

int mp_add_d(mp_int * a, mp_digit b, mp_int * c);
int mp_set_int(mp_int * a, unsigned long b);
int mp_div(mp_int * a, mp_int * b, mp_int * c, mp_int * d);
int mp_addmod(mp_int * a, mp_int * b, mp_int * c, mp_int * d);
int mp_cnt_lsb(mp_int *a);
int mp_sqrmod(mp_int * a, mp_int * b, mp_int * c);
int mp_prime_miller_rabin(mp_int * a, mp_int * b, int *result);
int mp_sub_d(mp_int * a, mp_digit b, mp_int * c);
int mp_is_square(mp_int *arg,int *ret);
int mp_read_radix (mp_int * a, const char *str, int radix);
unsigned long mp_get_int(mp_int * a);
int mp_mul(mp_int * a, mp_int * b, mp_int * c);
void BasePoint_test(mp_int *x1,mp_int *y1,mp_int *a,  mp_int *b,  mp_int *p, mp_int *n);
bool point_is_on_curve(mp_int *rx,mp_int *ry,mp_int *a,mp_int *b,mp_int *p);


