/* bqf.h - binary quadratic forms = ideal class group Cl(D) of an
 * imaginary quadratic order, D < 0, D = 1 (mod 4).
 *
 * This is the algebraic platform of section 3.2/3.3 of the report:
 *   G = Cl(Delta),  F(X) = X^2,  F^r(X) = X^(2^r).
 */
#ifndef BQF_H
#define BQF_H

#include <gmp.h>
#include <stddef.h>

typedef struct { mpz_t a, b, c; } bqf_t;

typedef struct {
    mpz_t D;             /* discriminant (negative, = 1 mod 4)          */
    size_t abytes;       /* byte width of a and |b| in canonical encoding */
    size_t elt_bytes;    /* canonical serialised element size            */
    /* operation counters, matching the cost units of section 9.1 */
    unsigned long long n_sq;        /* CSQ  : class-group squarings      */
    unsigned long long n_mul;       /* CGM  : general multiplications    */
    unsigned long long n_red_steps; /* reduction (rho) steps             */
    /* scratch */
    mpz_t t[16];
    bqf_t tf;
} cg_ctx;

void   cg_ctx_init(cg_ctx *x, const mpz_t D);
void   cg_ctx_clear(cg_ctx *x);
void   cg_reset_counters(cg_ctx *x);

void   bqf_init(bqf_t *f);
void   bqf_clear(bqf_t *f);
void   bqf_set(bqf_t *r, const bqf_t *f);
int    bqf_equal(const bqf_t *f, const bqf_t *g);

void   bqf_identity(bqf_t *r, cg_ctx *x);
int    bqf_is_identity(const bqf_t *f);
void   bqf_inv(bqf_t *r, const bqf_t *f, cg_ctx *x);
int    bqf_check(const bqf_t *f, cg_ctx *x);   /* 1 if reduced & correct disc */

void   bqf_reduce(bqf_t *f, cg_ctx *x);
void   bqf_mul(bqf_t *r, const bqf_t *f, const bqf_t *g, cg_ctx *x);  /* CGM */
void   bqf_sqr(bqf_t *r, const bqf_t *f, cg_ctx *x);                  /* CSQ */
void   bqf_chain(bqf_t *r, const bqf_t *f, unsigned long steps, cg_ctx *x); /* F^steps */
void   bqf_pow(bqf_t *r, const bqf_t *f, const mpz_t e, cg_ctx *x);

void   bqf_random(bqf_t *r, cg_ctx *x, gmp_randstate_t st);
void   bqf_prime_form(bqf_t *r, cg_ctx *x, unsigned long q_bits, gmp_randstate_t st);

size_t bqf_serialize(unsigned char *out, size_t cap, const bqf_t *f, cg_ctx *x);

/* discriminant generation: D = -p, p prime, p = 3 (mod 4) */
void   cg_gen_discriminant(mpz_t D, unsigned long bits, gmp_randstate_t st, int mr_rounds);

#endif
