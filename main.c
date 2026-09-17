/* main.c - CG-HWOTS reference implementation: self-tests and benchmarks. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <dlfcn.h>
#include "cghwots.h"

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int failures = 0;
static void check(int cond, const char *what)
{
    printf("  [%s] %s\n", cond ? " ok " : "FAIL", what);
    if (!cond) failures++;
}

/* ================= brute-force small class group ======================== */

typedef struct { long a, b, c; } sform;

static long enum_reduced(long D, sform *out, long cap)
{
    long a, b, c, n = 0;
    long amax = (long)floor(sqrt((double)(-D) / 3.0)) + 1;
    for (a = 1; a <= amax; a++) {
        for (b = -a; b <= a; b++) {
            long num;
            if (((b % 2) + 2) % 2 != ((D % 2) + 2) % 2) continue;
            num = b * b - D;
            if (num % (4 * a)) continue;
            c = num / (4 * a);
            if (c < a) continue;
            if ((a == c || b == -a) && b < 0) continue;
            if (n < cap) { out[n].a = a; out[n].b = b; out[n].c = c; }
            n++;
        }
    }
    return n;
}

static void small_group_test(long Dsmall)
{
    mpz_t D;
    cg_ctx x;
    sform reduced[4096];
    long h, i, j, bad = 0, sq_collision = 0;
    bqf_t *g, t1, t2, t3, id;
    char msg[128];

    mpz_init_set_si(D, Dsmall);
    cg_ctx_init(&x, D);
    h = enum_reduced(Dsmall, reduced, 4096);

    g = malloc(sizeof(bqf_t) * h);
    for (i = 0; i < h; i++) {
        bqf_init(&g[i]);
        mpz_set_si(g[i].a, reduced[i].a);
        mpz_set_si(g[i].b, reduced[i].b);
        mpz_set_si(g[i].c, reduced[i].c);
    }
    bqf_init(&t1); bqf_init(&t2); bqf_init(&t3); bqf_init(&id);
    bqf_identity(&id, &x);
    bqf_reduce(&id, &x);

    /* closure: every product is one of the enumerated reduced forms */
    for (i = 0; i < h; i++)
        for (j = 0; j < h; j++) {
            long m; int found = 0;
            bqf_mul(&t1, &g[i], &g[j], &x);
            if (!bqf_check(&t1, &x)) { bad++; continue; }
            for (m = 0; m < h; m++) if (bqf_equal(&t1, &g[m])) { found = 1; break; }
            if (!found) bad++;
        }
    sprintf(msg, "D=%ld  h=%ld  composition closed & reduced", Dsmall, h);
    check(bad == 0, msg);

    /* identity and inverse */
    bad = 0;
    for (i = 0; i < h; i++) {
        bqf_mul(&t1, &g[i], &id, &x);
        if (!bqf_equal(&t1, &g[i])) bad++;
        bqf_inv(&t2, &g[i], &x);
        bqf_mul(&t1, &g[i], &t2, &x);
        if (!bqf_equal(&t1, &id)) bad++;
    }
    check(bad == 0, "identity and inverses behave");

    /* associativity on all triples (h is small) */
    bad = 0;
    for (i = 0; i < h; i++)
        for (j = 0; j < h; j++) {
            long m;
            for (m = 0; m < h; m++) {
                bqf_mul(&t1, &g[i], &g[j], &x); bqf_mul(&t1, &t1, &g[m], &x);
                bqf_mul(&t2, &g[j], &g[m], &x); bqf_mul(&t2, &g[i], &t2, &x);
                if (!bqf_equal(&t1, &t2)) bad++;
            }
        }
    check(bad == 0, "associativity on all triples");

    /* Lagrange: f^h = identity */
    bad = 0;
    { mpz_t e; mpz_init_set_si(e, h);
      for (i = 0; i < h; i++) { bqf_pow(&t1, &g[i], e, &x); if (!bqf_equal(&t1, &id)) bad++; }
      mpz_clear(e); }
    sprintf(msg, "f^h = 1 for all %ld classes (h %s)", h, (h % 2) ? "odd" : "even");
    check(bad == 0, msg);

    /* squaring is a permutation iff h is odd (genus theory, section 3.3) */
    {
        bqf_t *sq = malloc(sizeof(bqf_t) * h);
        for (i = 0; i < h; i++) { bqf_init(&sq[i]); bqf_sqr(&sq[i], &g[i], &x); }
        for (i = 0; i < h; i++)
            for (j = i + 1; j < h; j++)
                if (bqf_equal(&sq[i], &sq[j])) sq_collision++;
        for (i = 0; i < h; i++) bqf_clear(&sq[i]);
        free(sq);
    }
    sprintf(msg, "squaring injective = %s (h %s)",
            sq_collision == 0 ? "yes" : "no", (h % 2) ? "odd" : "even");
    check((sq_collision == 0) == ((h % 2) == 1), msg);

    /* Dirichlet concordant-form cross-check on composite a */
    bad = 0;
    for (i = 0; i < h; i++) {
        long a = reduced[i].a, b = reduced[i].b, c = reduced[i].c, a1;
        for (a1 = 2; a1 * a1 <= a; a1++) {
            long a2;
            if (a % a1) continue;
            a2 = a / a1;
            {   long gg = a1 < a2 ? a1 : a2, gg2 = a1 > a2 ? a1 : a2, r;
                while (gg) { r = gg2 % gg; gg2 = gg; gg = r; }
                if (gg2 != 1) continue;             /* need gcd(a1,a2)=1 */
            }
            /* f1 = (a1, b, a2*c), f2 = (a2, b, a1*c) compose to (a1*a2, b, c) */
            mpz_set_si(t1.a, a1); mpz_set_si(t1.b, b); mpz_set_si(t1.c, a2 * c);
            mpz_set_si(t2.a, a2); mpz_set_si(t2.b, b); mpz_set_si(t2.c, a1 * c);
            bqf_reduce(&t1, &x); bqf_reduce(&t2, &x);
            bqf_mul(&t3, &t1, &t2, &x);
            if (!bqf_equal(&t3, &g[i])) bad++;
        }
    }
    check(bad == 0, "Dirichlet concordant composition cross-check");

    for (i = 0; i < h; i++) bqf_clear(&g[i]);
    free(g);
    bqf_clear(&t1); bqf_clear(&t2); bqf_clear(&t3); bqf_clear(&id);
    cg_ctx_clear(&x);
    mpz_clear(D);
}

/* ================= self tests =========================================== */

static void big_group_test(unsigned long bits, gmp_randstate_t st)
{
    mpz_t D, e1, e2, e3;
    cg_ctx x;
    bqf_t f, g, h, t1, t2, t3, id;
    int i, bad;
    char msg[160];

    mpz_inits(D, e1, e2, e3, NULL);
    cg_gen_discriminant(D, bits, st, 25);
    cg_ctx_init(&x, D);
    bqf_init(&f); bqf_init(&g); bqf_init(&h);
    bqf_init(&t1); bqf_init(&t2); bqf_init(&t3); bqf_init(&id);
    bqf_identity(&id, &x); bqf_reduce(&id, &x);

    bad = 0;
    for (i = 0; i < 20; i++) {
        bqf_random(&f, &x, st); bqf_random(&g, &x, st); bqf_random(&h, &x, st);
        if (!bqf_check(&f, &x)) bad++;
        bqf_mul(&t1, &f, &g, &x); bqf_mul(&t1, &t1, &h, &x);
        bqf_mul(&t2, &g, &h, &x); bqf_mul(&t2, &f, &t2, &x);
        if (!bqf_equal(&t1, &t2)) bad++;
        bqf_inv(&t3, &f, &x); bqf_mul(&t1, &f, &t3, &x);
        if (!bqf_equal(&t1, &id)) bad++;
    }
    sprintf(msg, "|D|=%lu bits: reduced forms, associativity, inverses", bits);
    check(bad == 0, msg);

    /* the homomorphism identity F^r(XY) = F^r(X) F^r(Y)  (section 3.2) */
    bad = 0;
    for (i = 0; i < 10; i++) {
        unsigned long r = 1 + (unsigned long)(rand() % 20);
        bqf_random(&f, &x, st); bqf_random(&g, &x, st);
        bqf_mul(&t1, &f, &g, &x); bqf_chain(&t1, &t1, r, &x);
        bqf_chain(&t2, &f, r, &x); bqf_chain(&t3, &g, r, &x);
        bqf_mul(&t2, &t2, &t3, &x);
        if (!bqf_equal(&t1, &t2)) bad++;
    }
    check(bad == 0, "F^r(XY) = F^r(X)F^r(Y)  (homomorphic chain)");

    /* F^r(X) = X^(2^r) */
    bad = 0;
    for (i = 0; i < 5; i++) {
        unsigned long r = 1 + (unsigned long)(rand() % 16);
        bqf_random(&f, &x, st);
        bqf_chain(&t1, &f, r, &x);
        mpz_set_ui(e1, 1); mpz_mul_2exp(e1, e1, r);
        bqf_pow(&t2, &f, e1, &x);
        if (!bqf_equal(&t1, &t2)) bad++;
    }
    check(bad == 0, "F^r(X) = X^(2^r)");

    bqf_clear(&f); bqf_clear(&g); bqf_clear(&h);
    bqf_clear(&t1); bqf_clear(&t2); bqf_clear(&t3); bqf_clear(&id);
    cg_ctx_clear(&x);
    mpz_clears(D, e1, e2, e3, NULL);
}

static void encoding_test(void)
{
    wparams p;
    int w, trial, bad = 0;
    unsigned d[600], d2[600];
    unsigned char m1[32], m2[32];
    for (w = 1; w <= 8; w *= 2) {
        wots_params(&p, w, 256);
        for (trial = 0; trial < 20000; trial++) {
            int j, ok = 0;
            for (j = 0; j < 32; j++) { m1[j] = rand() & 0xff; m2[j] = rand() & 0xff; }
            if (memcmp(m1, m2, 32) == 0) continue;
            wots_encode(d, m1, &p);
            wots_encode(d2, m2, &p);
            /* Lemma 1: some coordinate strictly decreases in either direction */
            for (j = 0; j < p.l; j++) if (d2[j] < d[j]) { ok = 1; break; }
            if (!ok) bad++;
        }
    }
    check(bad == 0, "Lemma 1: chain-order separation holds on 80k random pairs");
    { char msg[200];
      for (w = 1; w <= 8; w *= 2) {
        wots_params(&p, w, 256);
        sprintf(msg, "params  B=%3u : l1=%d l2=%d l=%d L=%u  chain budget lL=%u",
                p.B, p.l1, p.l2, p.l, p.L, (unsigned)(p.l * p.L));
        printf("        %s\n", msg);
      } }
}

static void scheme_test(unsigned long bits, gmp_randstate_t st)
{
    mpz_t D;
    cg_ctx x;
    wparams p;
    unsigned d[600];
    unsigned char msg[32], dg[32];
    int i, j;

    mpz_init(D);
    cg_gen_discriminant(D, bits, st, 25);
    cg_ctx_init(&x, D);
    wots_params(&p, 4, 256);           /* B = 16, l = 67, L = 15 */

    for (i = 0; i < 32; i++) msg[i] = (unsigned char)(i * 3 + 5);
    sha256(msg, 32, dg);
    wots_encode(d, dg, &p);

    /* --- base CG-WOTS --- */
    {
        cgwots_key k;
        bqf_t *sig = malloc(sizeof(bqf_t) * p.l);
        for (j = 0; j < p.l; j++) bqf_init(&sig[j]);
        cgwots_key_init(&k, &p);
        cgwots_keygen(&k, &x, st);
        cgwots_sign(sig, &k, d, &x);
        check(cgwots_verify(sig, k.pk, d, &p, &x), "CG-WOTS sign/verify");
        /* tamper */
        bqf_sqr(&sig[3], &sig[3], &x);
        check(!cgwots_verify(sig, k.pk, d, &p, &x), "CG-WOTS rejects a modified signature");
        for (j = 0; j < p.l; j++) bqf_clear(&sig[j]);
        free(sig);
        cgwots_key_clear(&k);
    }

    /* --- n-of-n --- */
    {
        int n = 3;
        hwots_nn s;
        bqf_t **parts = malloc(sizeof(bqf_t *) * n);
        bqf_t *sig = malloc(sizeof(bqf_t) * p.l);
        int okp = 1;
        for (j = 0; j < p.l; j++) bqf_init(&sig[j]);
        for (i = 0; i < n; i++) {
            parts[i] = malloc(sizeof(bqf_t) * p.l);
            for (j = 0; j < p.l; j++) bqf_init(&parts[i][j]);
        }
        hwots_nn_init(&s, n, &p);
        hwots_nn_keygen(&s, &x, st);
        for (i = 0; i < n; i++) {
            hwots_nn_partial(parts[i], &s, i, d, &x);
            if (!hwots_nn_verify_partial(parts[i], &s, i, d, &x)) okp = 0;
        }
        check(okp, "n/n partial signatures individually verifiable");
        hwots_nn_combine(sig, parts, n, &p, &x);
        check(cgwots_verify(sig, s.pk, d, &p, &x), "CG-HWOTS-3/3 aggregate verifies");

        /* a missing participant must break the aggregate */
        hwots_nn_combine(sig, parts, n - 1, &p, &x);
        check(!cgwots_verify(sig, s.pk, d, &p, &x), "CG-HWOTS-3/3 rejects a 2-party aggregate");

        for (i = 0; i < n; i++) { for (j = 0; j < p.l; j++) bqf_clear(&parts[i][j]); free(parts[i]); }
        free(parts);
        for (j = 0; j < p.l; j++) bqf_clear(&sig[j]);
        free(sig);
        hwots_nn_clear(&s);
    }

    /* --- rogue-key attack (section 6.1.1) --- */
    {
        int n = 3;
        hwots_nn s;
        bqf_t *b = malloc(sizeof(bqf_t) * p.l);
        bqf_t *P = malloc(sizeof(bqf_t) * p.l);
        bqf_t *sig = malloc(sizeof(bqf_t) * p.l);
        bqf_t tmp;
        unsigned char rogue_pk[32];
        bqf_init(&tmp);
        for (j = 0; j < p.l; j++) { bqf_init(&b[j]); bqf_init(&P[j]); bqf_init(&sig[j]); }
        hwots_nn_init(&s, n, &p);
        hwots_nn_keygen(&s, &x, st);
        /* participant 2 (last) sees p_0, p_1 and picks
           p*_j = F^L(b_j) * (p_0j p_1j)^{-1}                                  */
        for (j = 0; j < p.l; j++) {
            bqf_random(&b[j], &x, st);
            bqf_chain(&P[j], &b[j], p.L, &x);            /* F^L(b_j) = target endpoint */
            bqf_mul(&tmp, &s.e[0][j], &s.e[1][j], &x);
            bqf_inv(&tmp, &tmp, &x);
            bqf_mul(&s.e[2][j], &P[j], &tmp, &x);        /* rogue endpoint */
            bqf_set(&s.P[j], &P[j]);                     /* aggregate collapses to F^L(b_j) */
        }
        cgwots_pkhash(rogue_pk, s.P, p.l, &x);
        /* the rogue signs alone with b_j */
        for (j = 0; j < p.l; j++) bqf_chain(&sig[j], &b[j], d[j], &x);
        check(cgwots_verify(sig, rogue_pk, d, &p, &x),
              "rogue-key attack succeeds without registration (as predicted)");
        for (j = 0; j < p.l; j++) { bqf_clear(&b[j]); bqf_clear(&P[j]); bqf_clear(&sig[j]); }
        bqf_clear(&tmp);
        free(b); free(P); free(sig);
        hwots_nn_clear(&s);
    }

    /* --- k-of-n, all quorums --- */
    {
        int n = 5, k = 3;
        hwots_kn s;
        bqf_t **parts;
        bqf_t *sig = malloc(sizeof(bqf_t) * p.l);
        int T[8], ok_all = 1, ok_partial = 1, ok_small = 1;
        int i1, i2, i3;
        for (j = 0; j < p.l; j++) bqf_init(&sig[j]);
        parts = malloc(sizeof(bqf_t *) * k);
        for (i = 0; i < k; i++) {
            parts[i] = malloc(sizeof(bqf_t) * p.l);
            for (j = 0; j < p.l; j++) bqf_init(&parts[i][j]);
        }
        hwots_kn_init(&s, n, k, &p);
        hwots_kn_keygen(&s, &x, st);

        for (i1 = 0; i1 < n; i1++)
        for (i2 = i1 + 1; i2 < n; i2++)
        for (i3 = i2 + 1; i3 < n; i3++) {
            T[0] = i1; T[1] = i2; T[2] = i3;
            for (i = 0; i < k; i++) hwots_kn_partial(parts[i], &s, T[i], T, d, &x);
            hwots_kn_combine(sig, parts, k, &p, &x);
            if (!cgwots_verify(sig, s.pk, d, &p, &x)) ok_all = 0;
        }
        check(ok_all, "CG-HWOTS-3/5: all 10 quorums produce the same valid signature");

        /* partial verification against V^T_{i,j} */
        {
            bqf_t *V = malloc(sizeof(bqf_t) * p.l);
            bqf_t q;
            bqf_init(&q);
            for (j = 0; j < p.l; j++) bqf_init(&V[j]);
            T[0] = 0; T[1] = 2; T[2] = 4;
            for (i = 0; i < k; i++) {
                hwots_kn_partial(parts[i], &s, T[i], T, d, &x);
                hwots_kn_partial_endpoint(V, &s, T[i], T, &x);
                for (j = 0; j < p.l; j++) {
                    bqf_chain(&q, &parts[i][j], p.L - d[j], &x);
                    if (!bqf_equal(&q, &V[j])) ok_partial = 0;
                }
            }
            for (j = 0; j < p.l; j++) bqf_clear(&V[j]);
            bqf_clear(&q);
            free(V);
        }
        check(ok_partial, "CG-HWOTS-3/5 partial signatures individually verifiable");

        /* an unauthorised coalition of size k-1 cannot produce the signature */
        {
            T[0] = 0; T[1] = 1; T[2] = 2;
            for (i = 0; i < 2; i++) hwots_kn_partial(parts[i], &s, T[i], T, d, &x);
            hwots_kn_combine(sig, parts, 2, &p, &x);
            if (cgwots_verify(sig, s.pk, d, &p, &x)) ok_small = 0;
        }
        check(ok_small, "CG-HWOTS-3/5 rejects a 2-party coalition");

        for (i = 0; i < k; i++) { for (j = 0; j < p.l; j++) bqf_clear(&parts[i][j]); free(parts[i]); }
        free(parts);
        for (j = 0; j < p.l; j++) bqf_clear(&sig[j]);
        free(sig);
        hwots_kn_clear(&s);
    }

    /* --- conventional baseline --- */
    {
        cwots_key k;
        unsigned char (*sig)[32] = malloc(32 * (size_t)p.l);
        unsigned long long hc = 0;
        cwots_key_init(&k, &p);
        cwots_keygen(&k, &hc);
        cwots_sign(sig, &k, d, &hc);
        check(cwots_verify(sig, &k, d, &hc), "conventional SHA-256 WOTS sign/verify");
        free(sig);
        cwots_key_clear(&k);
    }

    cg_ctx_clear(&x);
    mpz_clear(D);
}

/* ================= benchmarks =========================================== */

static void load_disc(mpz_t D, const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    if (mpz_inp_str(D, f, 10) == 0) { fprintf(stderr, "bad discriminant file\n"); exit(1); }
    fclose(f);
}

/* optional: hardware-accelerated SHA-256 from libcrypto, loaded at run time,
   so that the conventional baseline is not penalised by a portable C hash. */
typedef unsigned char *(*ssl_sha256_fn)(const unsigned char *, size_t, unsigned char *);
static ssl_sha256_fn ssl_sha256 = NULL;
static void try_load_openssl(void)
{
    void *h = dlopen("libcrypto.so.3", RTLD_LAZY);
    if (!h) h = dlopen("libcrypto.so", RTLD_LAZY);
    if (!h) return;
    ssl_sha256 = (ssl_sha256_fn)dlsym(h, "SHA256");
    if (ssl_sha256) {   /* validate against the portable implementation */
        unsigned char a[32], b[32], m[55];
        memset(m, 0x42, sizeof m);
        sha256(m, 55, a);
        ssl_sha256(m, 55, b);
        if (memcmp(a, b, 32) != 0) ssl_sha256 = NULL;
    }
}

static void bench_primitives(const char *path, gmp_randstate_t st)
{
    mpz_t D;
    cg_ctx x;
    bqf_t f, g, r;
    double t0, t1;
    long iters, i;
    double csq, cgm;
    unsigned char h[32], buf[55];
    double hnat, hssl = 0, rho_per_sq = 0, tsample = 0;

    mpz_init(D);
    load_disc(D, path);
    cg_ctx_init(&x, D);
    bqf_init(&f); bqf_init(&g); bqf_init(&r);
    bqf_random(&f, &x, st);
    bqf_random(&g, &x, st);

    /* calibrate iteration count */
    iters = 64;
    for (;;) {
        t0 = now();
        for (i = 0; i < iters; i++) bqf_sqr(&r, &f, &x);
        t1 = now();
        if (t1 - t0 > 0.35 || iters > (1L<<22)) break;
        iters *= 4;
    }
    csq = (t1 - t0) / iters;
    { double redsteps;
      cg_reset_counters(&x);
      for (i = 0; i < 200; i++) bqf_sqr(&r, &f, &x);
      redsteps = (double)x.n_red_steps / 200.0;
      (void)redsteps; rho_per_sq = redsteps; }

    { long it2 = iters;
      t0 = now();
      for (i = 0; i < it2; i++) bqf_mul(&r, &f, &g, &x);
      t1 = now();
      cgm = (t1 - t0) / it2; }

    memset(buf, 0xa5, sizeof buf);
    { long it3 = 400000;
      t0 = now();
      for (i = 0; i < it3; i++) { buf[0] = (unsigned char)i; sha256(buf, 55, h); }
      t1 = now();
      hnat = (t1 - t0) / it3; }

    /* hardware SHA-256: bulk throughput on a 1 MiB buffer, converted to a
       per-compression cost.  This is the floor for a conventional chain step
       on a SHA-NI capable core, ignoring per-call overhead. */
    if (ssl_sha256) { long it3 = 400; size_t nb = 1u << 20;
      unsigned char *big = malloc(nb);
      memset(big, 0x37, nb);
      t0 = now();
      for (i = 0; i < it3; i++) ssl_sha256(big, nb, h);
      t1 = now();
      hssl = (t1 - t0) / it3 / (double)(nb / 64);
      free(big); }

    { long it4 = 200;
      t0 = now();
      for (i = 0; i < it4; i++) bqf_random(&r, &x, st);
      t1 = now();
      tsample = (t1 - t0) / it4; }

    printf("%7lu %10.3f %10.3f %8.1f %10.4f %8.1f %10.4f %8.1f %8.3f %6zu\n",
           (unsigned long)mpz_sizeinbase(D, 2),
           csq * 1e6, cgm * 1e6, rho_per_sq,
           hnat * 1e6, csq / hnat,
           hssl * 1e6, hssl > 0 ? csq / hssl : 0.0, tsample * 1e3, x.elt_bytes);
    fflush(stdout);

    bqf_clear(&f); bqf_clear(&g); bqf_clear(&r);
    cg_ctx_clear(&x);
    mpz_clear(D);
}

static void bench_scheme(const char *path, int w, int n, int k, int reps,
                         gmp_randstate_t st)
{
    mpz_t D;
    cg_ctx x;
    wparams p;
    unsigned d[600];
    unsigned char dg[32];
    double t0, tsign, tver, tcomb, tkeygen;
    int i, j, rep;
    unsigned long long sq_sign = 0, sq_ver = 0, mul_comb = 0, kg_sq = 0, kg_mul = 0;
    unsigned long Dsum = 0;

    mpz_init(D);
    load_disc(D, path);
    cg_ctx_init(&x, D);
    wots_params(&p, w, 256);

    if (k <= 0 || k == n) {
        hwots_nn s;
        bqf_t **parts, *sig;
        hwots_nn_init(&s, n, &p);
        cg_reset_counters(&x);
        t0 = now();
        hwots_nn_keygen(&s, &x, st);
        tkeygen = now() - t0;
        kg_sq = x.n_sq; kg_mul = x.n_mul;

        parts = malloc(sizeof(bqf_t *) * n);
        for (i = 0; i < n; i++) {
            parts[i] = malloc(sizeof(bqf_t) * p.l);
            for (j = 0; j < p.l; j++) bqf_init(&parts[i][j]);
        }
        sig = malloc(sizeof(bqf_t) * p.l);
        for (j = 0; j < p.l; j++) bqf_init(&sig[j]);

        tsign = tver = tcomb = 0;
        for (rep = 0; rep < reps; rep++) {
            unsigned char m[8];
            for (i = 0; i < 8; i++) m[i] = (unsigned char)(rep * 31 + i);
            sha256(m, 8, dg);
            wots_encode(d, dg, &p);
            Dsum += (unsigned long)wots_sum(d, &p);

            cg_reset_counters(&x);
            t0 = now();
            for (i = 0; i < n; i++) hwots_nn_partial(parts[i], &s, i, d, &x);
            tsign += now() - t0;
            sq_sign += x.n_sq;

            cg_reset_counters(&x);
            t0 = now();
            hwots_nn_combine(sig, parts, n, &p, &x);
            tcomb += now() - t0;
            mul_comb += x.n_mul;

            cg_reset_counters(&x);
            t0 = now();
            if (!cgwots_verify(sig, s.pk, d, &p, &x)) { printf("VERIFY FAILED\n"); exit(1); }
            tver += now() - t0;
            sq_ver += x.n_sq;
        }
        printf("%-8s B=%-4u n=%-3d k=%-3d keygen %8.2f s (%9llu CSQ %7llu CGM) | signers %8.3f s (%6.1f CSQ ea) | combine %8.4f s (%4llu CGM) | verify %7.3f s (%6.1f CSQ) | sig %7.1f KiB\n",
               "n-of-n", p.B, n, n, tkeygen, kg_sq, kg_mul, tsign / reps,
               (double)sq_sign / reps / n, tcomb / reps,
               (unsigned long long)(mul_comb / reps), tver / reps,
               (double)sq_ver / reps, p.l * x.elt_bytes / 1024.0);
        fflush(stdout);

        for (i = 0; i < n; i++) { for (j = 0; j < p.l; j++) bqf_clear(&parts[i][j]); free(parts[i]); }
        free(parts);
        for (j = 0; j < p.l; j++) bqf_clear(&sig[j]);
        free(sig);
        hwots_nn_clear(&s);
    } else {
        hwots_kn s;
        bqf_t **parts, *sig;
        int T[64];
        hwots_kn_init(&s, n, k, &p);
        cg_reset_counters(&x);
        t0 = now();
        hwots_kn_keygen(&s, &x, st);
        tkeygen = now() - t0;
        kg_sq = x.n_sq; kg_mul = x.n_mul;

        parts = malloc(sizeof(bqf_t *) * k);
        for (i = 0; i < k; i++) {
            parts[i] = malloc(sizeof(bqf_t) * p.l);
            for (j = 0; j < p.l; j++) bqf_init(&parts[i][j]);
        }
        sig = malloc(sizeof(bqf_t) * p.l);
        for (j = 0; j < p.l; j++) bqf_init(&sig[j]);
        for (i = 0; i < k; i++) T[i] = i;

        tsign = tver = tcomb = 0;
        for (rep = 0; rep < reps; rep++) {
            unsigned char m[8];
            for (i = 0; i < 8; i++) m[i] = (unsigned char)(rep * 31 + i);
            sha256(m, 8, dg);
            wots_encode(d, dg, &p);

            cg_reset_counters(&x);
            t0 = now();
            for (i = 0; i < k; i++) hwots_kn_partial(parts[i], &s, T[i], T, d, &x);
            tsign += now() - t0;
            sq_sign += x.n_sq; mul_comb += x.n_mul;

            t0 = now();
            hwots_kn_combine(sig, parts, k, &p, &x);
            tcomb += now() - t0;

            cg_reset_counters(&x);
            t0 = now();
            if (!cgwots_verify(sig, s.pk, d, &p, &x)) { printf("VERIFY FAILED\n"); exit(1); }
            tver += now() - t0;
            sq_ver += x.n_sq;
        }
        printf("%-8s B=%-4u n=%-3d k=%-3d keygen %8.2f s (%9llu CSQ %7llu CGM) | quorum  %8.3f s (%6.1f CSQ ea) | combine %8.4f s (%4llu CGM) | verify %7.3f s (%6.1f CSQ) | sig %7.1f KiB | N=%ld\n",
               "k-of-n", p.B, n, k, tkeygen, kg_sq, kg_mul, tsign / reps,
               (double)sq_sign / reps / k, tcomb / reps,
               (unsigned long long)(mul_comb / reps / k), tver / reps,
               (double)sq_ver / reps, p.l * x.elt_bytes / 1024.0, s.N);
        fflush(stdout);

        for (i = 0; i < k; i++) { for (j = 0; j < p.l; j++) bqf_clear(&parts[i][j]); free(parts[i]); }
        free(parts);
        for (j = 0; j < p.l; j++) bqf_clear(&sig[j]);
        free(sig);
        hwots_kn_clear(&s);
    }
    (void)Dsum;
    cg_ctx_clear(&x);
    mpz_clear(D);
}

static void bench_conventional(int w, int reps)
{
    wparams p;
    cwots_key k;
    unsigned char (*sig)[32];
    unsigned d[600];
    unsigned char dg[32], m[8];
    unsigned long long hc;
    double t0, tk, ts = 0, tv = 0;
    int i, rep;
    unsigned long long hs = 0, hv = 0;

    wots_params(&p, w, 256);
    cwots_key_init(&k, &p);
    sig = malloc(32 * (size_t)p.l);
    hc = 0;
    t0 = now(); cwots_keygen(&k, &hc); tk = now() - t0;

    for (rep = 0; rep < reps; rep++) {
        for (i = 0; i < 8; i++) m[i] = (unsigned char)(rep * 31 + i);
        sha256(m, 8, dg);
        wots_encode(d, dg, &p);
        hc = 0; t0 = now(); cwots_sign(sig, &k, d, &hc); ts += now() - t0; hs += hc;
        hc = 0; t0 = now();
        if (!cwots_verify(sig, &k, d, &hc)) { printf("CONV VERIFY FAILED\n"); exit(1); }
        tv += now() - t0; hv += hc;
    }
    printf("%-10s B=%-4u l=%-4d       keygen %9.6f s | sign %11.6f s (%6.1f hashes) "
           "| verify %11.6f s (%6.1f hashes) | sig %6.3f KiB\n",
           "SHA256-WOTS", p.B, p.l, tk, ts / reps, (double)hs / reps,
           tv / reps, (double)hv / reps, p.l * 32 / 1024.0);
    fflush(stdout);
    free(sig);
    cwots_key_clear(&k);
}

static void bench_digits(void)
{
    /* empirical mean of D(d) and V(d) over random digests, to confirm
       the expectations tabulated in section 9.3 */
    int w;
    printf("\n  Winternitz digit statistics over 200000 random digests:\n");
    printf("  %-5s %-5s %-5s %-6s %-12s %-12s\n", "B", "l", "L", "lL", "mean D", "mean V");
    for (w = 1; w <= 8; w *= 2) {
        wparams p;
        unsigned d[600];
        unsigned char dg[32];
        double sum = 0;
        int t, j;
        wots_params(&p, w, 256);
        for (t = 0; t < 200000; t++) {
            for (j = 0; j < 32; j++) dg[j] = rand() & 0xff;
            wots_encode(d, dg, &p);
            sum += (double)wots_sum(d, &p);
        }
        sum /= 200000;
        printf("  %-5u %-5d %-5u %-6u %-12.1f %-12.1f\n",
               p.B, p.l, p.L, (unsigned)(p.l * p.L), sum, p.l * p.L - sum);
    }
}

/* ================= main ================================================= */

int main(int argc, char **argv)
{
    gmp_randstate_t st;
    gmp_randinit_mt(st);
    gmp_randseed_ui(st, 20260905UL);
    srand(12345);

    if (argc >= 2 && strcmp(argv[1], "gendisc") == 0) {
        mpz_t D;
        unsigned long bits = strtoul(argv[2], NULL, 10);
        int mr = (argc >= 5) ? atoi(argv[4]) : 25;
        double t0 = now();
        mpz_init(D);
        gmp_randseed_ui(st, 20260905UL + bits);
        cg_gen_discriminant(D, bits, st, mr);
        { FILE *f = fopen(argv[3], "w");
          mpz_out_str(f, 10, D);
          fputc('\n', f);
          fclose(f); }
        fprintf(stderr, "discriminant %lu bits written to %s in %.1f s\n",
                bits, argv[3], now() - t0);
        mpz_clear(D);
        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "selftest") == 0) {
        printf("=== class group arithmetic, exhaustive small groups ===\n");
        small_group_test(-23);
        small_group_test(-47);
        small_group_test(-71);
        small_group_test(-199);
        small_group_test(-1019);
        small_group_test(-2003);
        small_group_test(-15);      /* even class number, squaring NOT injective */
        small_group_test(-231);     /* even class number                          */
        small_group_test(-2311);
        printf("\n=== class group arithmetic, cryptographic sizes ===\n");
        big_group_test(256, st);
        big_group_test(512, st);
        big_group_test(1024, st);
        printf("\n=== Winternitz encoding ===\n");
        encoding_test();
        printf("\n=== schemes (512-bit discriminant, B=16) ===\n");
        scheme_test(512, st);
        printf("\n%s: %d failure(s)\n", failures ? "RESULT" : "RESULT", failures);
        return failures ? 1 : 0;
    }

    if (argc >= 2 && strcmp(argv[1], "primitives") == 0) {
        int i;
        try_load_openssl();
        printf("%7s %10s %10s %8s %10s %8s %10s %8s %8s %6s\n",
               "|D|bits", "CSQ(us)", "CGM(us)", "rho/sq",
               "SHA-C(us)", "CSQ/H", "SHAni-blk", "CSQ/Hni", "smpl(ms)", "bytes");
        for (i = 2; i < argc; i++) bench_primitives(argv[i], st);
        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "scheme") == 0) {
        /* scheme <discfile> <w> <n> <k> <reps> */
        bench_scheme(argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), st);
        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "conv") == 0) {
        bench_conventional(atoi(argv[2]), atoi(argv[3]));
        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "digits") == 0) { bench_digits(); return 0; }

    fprintf(stderr,
      "usage:\n"
      "  %s gendisc <bits> <file> [mr_rounds]\n"
      "  %s selftest\n"
      "  %s primitives <discfile>...\n"
      "  %s scheme <discfile> <w> <n> <k> <reps>\n"
      "  %s conv <w> <reps>\n"
      "  %s digits\n", argv[0], argv[0], argv[0], argv[0], argv[0], argv[0]);
    return 1;
}
