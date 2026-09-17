/* bqf.c - class group of an imaginary quadratic field via binary quadratic forms.
 *
 * Composition uses the classical Dirichlet/Gauss formula
 *     e   = gcd(a1, a2, s),  s = (b1+b2)/2,  n = (b1-b2)/2
 *     e   = U*a1 + V*a2 + W*s
 *     a3  = a1*a2 / e^2
 *     b3  = b2 + (2*a2/e) * (V*n - W*c2)      (mod 2*a3)
 *     c3  = (b3^2 - D) / (4*a3)
 * followed by full reduction.  Duplication is the specialised case
 * a1=a2=a, b1=b2=b, n=0.
 *
 * NOTE: this is a straightforward "compose-then-reduce" implementation.  It is
 * NOT NUCOMP/NUDUPL, so the measured cost of one squaring is an upper bound;
 * see the report for the implication.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "bqf.h"

void bqf_init(bqf_t *f){ mpz_inits(f->a, f->b, f->c, NULL); }
void bqf_clear(bqf_t *f){ mpz_clears(f->a, f->b, f->c, NULL); }
void bqf_set(bqf_t *r, const bqf_t *f){
    if (r == f) return;
    mpz_set(r->a, f->a); mpz_set(r->b, f->b); mpz_set(r->c, f->c);
}
int bqf_equal(const bqf_t *f, const bqf_t *g){
    return mpz_cmp(f->a,g->a)==0 && mpz_cmp(f->b,g->b)==0 && mpz_cmp(f->c,g->c)==0;
}

void cg_ctx_init(cg_ctx *x, const mpz_t D)
{
    int i;
    mpz_init_set(x->D, D);
    for (i = 0; i < 16; i++) mpz_init(x->t[i]);
    bqf_init(&x->tf);
    /* reduced forms satisfy a <= sqrt(|D|/3), |b| <= a */
    x->abytes = (mpz_sizeinbase(D, 2) / 2 + 8) / 8;
    x->elt_bytes = 1 + 2 * x->abytes;   /* sign byte + a + |b| */
    cg_reset_counters(x);
}
void cg_ctx_clear(cg_ctx *x)
{
    int i;
    mpz_clear(x->D);
    for (i = 0; i < 16; i++) mpz_clear(x->t[i]);
    bqf_clear(&x->tf);
}
void cg_reset_counters(cg_ctx *x){ x->n_sq = x->n_mul = x->n_red_steps = 0; }

/* --- reduction ---------------------------------------------------------- */

static void bqf_normalize(bqf_t *f, cg_ctx *x)
{
    mpz_neg(x->t[0], f->a);
    if (mpz_cmp(f->b, x->t[0]) > 0 && mpz_cmp(f->b, f->a) <= 0) return;
    mpz_sub(x->t[0], f->a, f->b);        /* a - b            */
    mpz_mul_2exp(x->t[1], f->a, 1);      /* 2a               */
    mpz_fdiv_q(x->t[2], x->t[0], x->t[1]);/* r               */
    mpz_mul(x->t[3], f->a, x->t[2]);     /* a r              */
    mpz_add(x->t[3], x->t[3], f->b);     /* a r + b          */
    mpz_mul(x->t[3], x->t[3], x->t[2]);  /* r(a r + b)       */
    mpz_add(f->c, f->c, x->t[3]);        /* c += r(ar+b)     */
    mpz_mul(x->t[3], x->t[2], x->t[1]);  /* 2 a r            */
    mpz_add(f->b, f->b, x->t[3]);
}

/* one rho step: (a,b,c) -> (c, -b + 2cr, a - br + cr^2), r = floor((c+b)/(2c)) */
static void bqf_rho(bqf_t *f, cg_ctx *x)
{
    mpz_add(x->t[0], f->c, f->b);
    mpz_mul_2exp(x->t[1], f->c, 1);
    mpz_fdiv_q(x->t[2], x->t[0], x->t[1]);   /* r      */
    mpz_mul(x->t[3], f->c, x->t[2]);         /* c r    */
    mpz_sub(x->t[4], x->t[3], f->b);         /* cr - b */
    mpz_mul(x->t[4], x->t[4], x->t[2]);
    mpz_add(x->t[4], x->t[4], f->a);         /* new c  */
    mpz_mul_2exp(x->t[5], x->t[3], 1);
    mpz_sub(x->t[5], x->t[5], f->b);         /* new b  */
    mpz_set(f->a, f->c);
    mpz_set(f->b, x->t[5]);
    mpz_set(f->c, x->t[4]);
}

void bqf_reduce(bqf_t *f, cg_ctx *x)
{
    bqf_normalize(f, x);
    while (mpz_cmp(f->a, f->c) > 0) { bqf_rho(f, x); x->n_red_steps++; }
    if (mpz_cmp(f->a, f->c) == 0 && mpz_sgn(f->b) < 0) mpz_neg(f->b, f->b);
}

/* --- basic elements ----------------------------------------------------- */

void bqf_identity(bqf_t *r, cg_ctx *x)
{
    mpz_set_ui(r->a, 1);
    if (mpz_odd_p(x->D)) {
        mpz_set_ui(r->b, 1);
        mpz_ui_sub(r->c, 1, x->D);        /* (1-D)/4          */
        mpz_fdiv_q_2exp(r->c, r->c, 2);
    } else {
        mpz_set_ui(r->b, 0);
        mpz_neg(r->c, x->D);
        mpz_fdiv_q_2exp(r->c, r->c, 2);   /* -D/4             */
    }
}
int bqf_is_identity(const bqf_t *f){ return mpz_cmp_ui(f->a, 1) == 0; }

void bqf_inv(bqf_t *r, const bqf_t *f, cg_ctx *x)
{
    bqf_set(r, f);
    mpz_neg(r->b, r->b);
    bqf_reduce(r, x);
}

int bqf_check(const bqf_t *f, cg_ctx *x)
{
    /* discriminant b^2 - 4ac == D and reduced */
    mpz_mul(x->t[0], f->b, f->b);
    mpz_mul(x->t[1], f->a, f->c);
    mpz_mul_2exp(x->t[1], x->t[1], 2);
    mpz_sub(x->t[0], x->t[0], x->t[1]);
    if (mpz_cmp(x->t[0], x->D) != 0) return 0;
    if (mpz_sgn(f->a) <= 0) return 0;
    if (mpz_cmpabs(f->b, f->a) > 0) return 0;
    if (mpz_cmp(f->a, f->c) > 0) return 0;
    if ((mpz_cmp(f->a, f->c) == 0 || mpz_cmpabs(f->b, f->a) == 0) && mpz_sgn(f->b) < 0)
        return 0;
    return 1;
}

/* --- composition -------------------------------------------------------- */

void bqf_mul(bqf_t *r, const bqf_t *f, const bqf_t *g, cg_ctx *x)
{
#define S   x->t[0]
#define N   x->t[1]
#define E1  x->t[2]
#define U1  x->t[3]
#define V1  x->t[4]
#define EE  x->t[5]
#define U2  x->t[6]
#define V2  x->t[7]
#define VV  x->t[8]
#define A3  x->t[9]
#define B3  x->t[10]
#define C3  x->t[11]
#define T0  x->t[12]
#define T1  x->t[13]
#define T2  x->t[14]

    x->n_mul++;

    mpz_add(S, f->b, g->b); mpz_divexact_ui(S, S, 2);   /* s = (b1+b2)/2 */
    mpz_sub(N, f->b, g->b); mpz_divexact_ui(N, N, 2);   /* n = (b1-b2)/2 */

    mpz_gcdext(E1, U1, V1, f->a, g->a);                 /* U1*a1 + V1*a2 = e1 */
    mpz_gcdext(EE, U2, V2, E1, S);                      /* U2*e1 + V2*s  = e  */
    mpz_mul(VV, U2, V1);                                /* V = U2*V1 (coeff of a2) */

    mpz_mul(A3, f->a, g->a);
    mpz_mul(T0, EE, EE);
    mpz_divexact(A3, A3, T0);                           /* a3 = a1 a2 / e^2 */

    mpz_divexact(T0, g->a, EE);
    mpz_mul_2exp(T0, T0, 1);                            /* 2 a2 / e */
    mpz_mul(T1, VV, N);
    mpz_mul(T2, V2, g->c);
    mpz_sub(T1, T1, T2);                                /* V*n - W*c2 */
    mpz_mul(T0, T0, T1);
    mpz_add(B3, g->b, T0);

    mpz_mul_2exp(T0, A3, 1);                            /* 2 a3 */
    mpz_mod(B3, B3, T0);                                /* [0, 2a3) */
    if (mpz_cmp(B3, A3) > 0) mpz_sub(B3, B3, T0);       /* (-a3, a3] */

    mpz_mul(T1, B3, B3);
    mpz_sub(T1, T1, x->D);
    mpz_mul_2exp(T2, A3, 2);
    mpz_divexact(C3, T1, T2);

    mpz_set(r->a, A3); mpz_set(r->b, B3); mpz_set(r->c, C3);
    bqf_reduce(r, x);

#undef S
#undef N
#undef E1
#undef U1
#undef V1
#undef EE
#undef U2
#undef V2
#undef VV
#undef A3
#undef B3
#undef C3
#undef T0
#undef T1
#undef T2
}

/* duplication: the case a1=a2, b1=b2, n=0 */
void bqf_sqr(bqf_t *r, const bqf_t *f, cg_ctx *x)
{
#define G   x->t[0]
#define U   x->t[1]
#define V   x->t[2]
#define AA  x->t[3]
#define A3  x->t[4]
#define B3  x->t[5]
#define C3  x->t[6]
#define T0  x->t[7]
#define T1  x->t[8]

    x->n_sq++;

    mpz_gcdext(G, U, V, f->a, f->b);      /* U*a + V*b = g = gcd(a,b) */
    mpz_divexact(AA, f->a, G);            /* A = a/g                  */
    if (mpz_cmp_ui(AA, 1) != 0) mpz_mod(V, V, AA);
    else mpz_set_ui(V, 0);

    mpz_mul(A3, AA, AA);                  /* a3 = A^2                 */

    mpz_mul(T0, AA, V);
    mpz_mul_2exp(T0, T0, 1);
    mpz_mul(T0, T0, f->c);
    mpz_sub(B3, f->b, T0);                /* b3 = b - 2*A*V*c         */

    mpz_mul_2exp(T0, A3, 1);
    mpz_mod(B3, B3, T0);
    if (mpz_cmp(B3, A3) > 0) mpz_sub(B3, B3, T0);

    mpz_mul(T1, B3, B3);
    mpz_sub(T1, T1, x->D);
    mpz_mul_2exp(T0, A3, 2);
    mpz_divexact(C3, T1, T0);

    mpz_set(r->a, A3); mpz_set(r->b, B3); mpz_set(r->c, C3);
    bqf_reduce(r, x);

#undef G
#undef U
#undef V
#undef AA
#undef A3
#undef B3
#undef C3
#undef T0
#undef T1
}

void bqf_chain(bqf_t *r, const bqf_t *f, unsigned long steps, cg_ctx *x)
{
    unsigned long i;
    bqf_set(r, f);
    for (i = 0; i < steps; i++) bqf_sqr(r, r, x);
}

void bqf_pow(bqf_t *r, const bqf_t *f, const mpz_t e, cg_ctx *x)
{
    bqf_t base, acc;
    long i, bits;
    if (mpz_sgn(e) == 0) { bqf_identity(r, x); return; }
    bqf_init(&base); bqf_init(&acc);
    bqf_set(&base, f);
    bqf_identity(&acc, x);
    bits = (long)mpz_sizeinbase(e, 2);
    for (i = bits - 1; i >= 0; i--) {
        bqf_sqr(&acc, &acc, x);
        if (mpz_tstbit(e, i)) bqf_mul(&acc, &acc, &base, x);
    }
    bqf_set(r, &acc);
    bqf_clear(&base); bqf_clear(&acc);
}

/* --- sampling ----------------------------------------------------------- */

/* prime form: a = q prime, q = 3 mod 4, (D|q) = 1, b = sqrt(D) mod q odd. */
void bqf_prime_form(bqf_t *r, cg_ctx *x, unsigned long q_bits, gmp_randstate_t st)
{
    mpz_t q, b, e, t;
    mpz_inits(q, b, e, t, NULL);
    for (;;) {
        mpz_urandomb(q, st, q_bits);
        mpz_setbit(q, q_bits - 1);
        /* force q = 3 mod 4 then walk to a prime keeping the residue */
        mpz_sub_ui(q, q, mpz_fdiv_ui(q, 4));
        mpz_add_ui(q, q, 3);
        while (!mpz_probab_prime_p(q, 12)) mpz_add_ui(q, q, 4);
        if (mpz_kronecker(x->D, q) != 1) continue;
        mpz_add_ui(e, q, 1);
        mpz_fdiv_q_2exp(e, e, 2);            /* (q+1)/4 */
        mpz_mod(t, x->D, q);
        mpz_powm(b, t, e, q);
        mpz_mul(t, b, b);
        mpz_sub(t, t, x->D);
        mpz_mod(t, t, q);
        if (mpz_sgn(t) != 0) continue;       /* should not happen */
        if (mpz_even_p(b)) mpz_sub(b, q, b); /* need b odd so 4q | b^2 - D */
        mpz_set(r->a, q);
        mpz_set(r->b, b);
        mpz_mul(t, b, b);
        mpz_sub(t, t, x->D);
        mpz_mul_2exp(e, q, 2);
        mpz_divexact(r->c, t, e);
        bqf_reduce(r, x);
        break;
    }
    mpz_clears(q, b, e, t, NULL);
}

/* Random group element: product of m independent random 64-bit prime forms,
 * with m chosen so that the unreduced product exceeds the reduction bound
 * sqrt(|D|/3).  Fewer factors would leave a form whose leading coefficient is
 * far smaller than a typical reduced form: such an element is cheap to square
 * for the first few steps and is not distributed like a random class.  This is
 * a heuristic sampler based on equidistribution of prime forms, not a proof of
 * uniformity; see the limitations section of the report. */
void bqf_random(bqf_t *r, cg_ctx *x, gmp_randstate_t st)
{
    bqf_t p;
    int i, m;
    m = (int)(mpz_sizeinbase(x->D, 2) / 2 / 64) + 4;
    bqf_init(&p);
    bqf_prime_form(r, x, 64, st);
    for (i = 1; i < m; i++) {
        bqf_prime_form(&p, x, 64, st);
        bqf_mul(r, r, &p, x);
    }
    bqf_clear(&p);
}

/* --- canonical encoding -------------------------------------------------- */

size_t bqf_serialize(unsigned char *out, size_t cap, const bqf_t *f, cg_ctx *x)
{
    size_t na, nb;
    if (cap < x->elt_bytes) return 0;
    memset(out, 0, x->elt_bytes);
    out[0] = (unsigned char)(mpz_sgn(f->b) < 0 ? 1 : 0);
    na = (mpz_sizeinbase(f->a, 2) + 7) / 8;
    nb = (mpz_sizeinbase(f->b, 2) + 7) / 8;
    if (na > x->abytes || nb > x->abytes) return 0;
    mpz_export(out + 1 + (x->abytes - na), NULL, 1, 1, 1, 0, f->a);
    if (mpz_sgn(f->b) != 0)
        mpz_export(out + 1 + x->abytes + (x->abytes - nb), NULL, 1, 1, 1, 0, f->b);
    return x->elt_bytes;
}

/* --- discriminant generation -------------------------------------------- */

void cg_gen_discriminant(mpz_t D, unsigned long bits, gmp_randstate_t st, int mr_rounds)
{
    mpz_t p;
    static const unsigned small[] = {
        3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
        101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
        191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
        281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
        389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
        491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,
        607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
        719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
        829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,
        953,967,971,977,983,991,997,0};
    mpz_init(p);
    mpz_urandomb(p, st, bits);
    mpz_setbit(p, bits - 1);
    mpz_sub_ui(p, p, mpz_fdiv_ui(p, 4));
    mpz_add_ui(p, p, 3);                  /* p = 3 mod 4 */
    for (;;) {
        int i, ok = 1;
        for (i = 0; small[i]; i++)
            if (mpz_divisible_ui_p(p, small[i])) { ok = 0; break; }
        if (ok && mpz_probab_prime_p(p, mr_rounds)) break;
        mpz_add_ui(p, p, 4);
    }
    mpz_neg(D, p);
    mpz_clear(p);
}
