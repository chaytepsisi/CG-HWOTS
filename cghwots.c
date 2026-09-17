/* cghwots.c */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "cghwots.h"

/* ---------------- Winternitz encoding ----------------------------------- */

void wots_params(wparams *p, int w, int digest_bits)
{
    unsigned long prod;
    int l2;
    double lg;
    p->w = w;
    p->B = 1u << w;
    p->L = p->B - 1;
    p->digest_bits = digest_bits;
    p->l1 = (digest_bits + w - 1) / w;
    prod = (unsigned long)p->l1 * p->L;
    /* l2 = floor(log_B(l1*L)) + 1 */
    l2 = 0;
    { unsigned long v = prod; while (v >= p->B) { v /= p->B; l2++; } }
    lg = 0; (void)lg;
    p->l2 = l2 + 1;
    p->l = p->l1 + p->l2;
}

void wots_encode(unsigned *d, const unsigned char *digest, const wparams *p)
{
    int i, j;
    unsigned long cs = 0;
    int w = p->w;
    /* message digits, MSB-first */
    for (i = 0; i < p->l1; i++) {
        int bitpos = i * w;
        int byte = bitpos >> 3;
        int shift = 8 - w - (bitpos & 7);
        d[i] = (digest[byte] >> shift) & p->L;
    }
    for (i = 0; i < p->l1; i++) cs += p->L - d[i];
    /* checksum digits, MSB-first, exactly l2 base-B digits */
    for (j = 0; j < p->l2; j++) {
        int sh = (p->l2 - 1 - j) * w;
        d[p->l1 + j] = (unsigned)((cs >> sh) & p->L);
    }
}

unsigned long wots_sum(const unsigned *d, const wparams *p)
{
    int i; unsigned long s = 0;
    for (i = 0; i < p->l; i++) s += d[i];
    return s;
}

/* ---------------- CG-WOTS ----------------------------------------------- */

static bqf_t *forms_alloc(int n)
{
    bqf_t *v = (bqf_t *)malloc(sizeof(bqf_t) * n);
    int i; for (i = 0; i < n; i++) bqf_init(&v[i]);
    return v;
}
static void forms_free(bqf_t *v, int n)
{
    int i; for (i = 0; i < n; i++) bqf_clear(&v[i]);
    free(v);
}

void cgwots_key_init(cgwots_key *k, const wparams *p)
{
    k->p = *p;
    k->A = forms_alloc(p->l);
    k->P = forms_alloc(p->l);
}
void cgwots_key_clear(cgwots_key *k)
{
    forms_free(k->A, k->p.l);
    forms_free(k->P, k->p.l);
}

void cgwots_pkhash(unsigned char out[32], bqf_t *P, int l, cg_ctx *x)
{
    sha256_ctx c;
    unsigned char *buf = (unsigned char *)malloc(x->elt_bytes);
    int j;
    sha256_init(&c);
    sha256_update(&c, "CG-HWOTS/pk/v1", 14);
    for (j = 0; j < l; j++) {
        size_t nb = bqf_serialize(buf, x->elt_bytes, &P[j], x);
        assert(nb == x->elt_bytes);
        sha256_update(&c, buf, nb);
    }
    sha256_final(&c, out);
    free(buf);
}

void cgwots_keygen(cgwots_key *k, cg_ctx *x, gmp_randstate_t st)
{
    int j;
    for (j = 0; j < k->p.l; j++) {
        bqf_random(&k->A[j], x, st);
        bqf_chain(&k->P[j], &k->A[j], k->p.L, x);
    }
    cgwots_pkhash(k->pk, k->P, k->p.l, x);
}

void cgwots_sign(bqf_t *sig, const cgwots_key *k, const unsigned *d, cg_ctx *x)
{
    int j;
    for (j = 0; j < k->p.l; j++)
        bqf_chain(&sig[j], &k->A[j], d[j], x);
}

int cgwots_verify(bqf_t *sig, const unsigned char pk[32], const unsigned *d,
                  const wparams *p, cg_ctx *x)
{
    bqf_t *Q = forms_alloc(p->l);
    unsigned char h[32];
    int j, ok;
    for (j = 0; j < p->l; j++)
        bqf_chain(&Q[j], &sig[j], p->L - d[j], x);
    cgwots_pkhash(h, Q, p->l, x);
    ok = memcmp(h, pk, 32) == 0;
    forms_free(Q, p->l);
    return ok;
}

/* ---------------- CG-HWOTS-n/n ------------------------------------------ */

void hwots_nn_init(hwots_nn *s, int n, const wparams *p)
{
    int i;
    s->n = n; s->p = *p;
    s->a = (bqf_t **)malloc(sizeof(bqf_t *) * n);
    s->e = (bqf_t **)malloc(sizeof(bqf_t *) * n);
    for (i = 0; i < n; i++) { s->a[i] = forms_alloc(p->l); s->e[i] = forms_alloc(p->l); }
    s->P = forms_alloc(p->l);
}
void hwots_nn_clear(hwots_nn *s)
{
    int i;
    for (i = 0; i < s->n; i++) { forms_free(s->a[i], s->p.l); forms_free(s->e[i], s->p.l); }
    free(s->a); free(s->e);
    forms_free(s->P, s->p.l);
}

void hwots_nn_keygen(hwots_nn *s, cg_ctx *x, gmp_randstate_t st)
{
    int i, j;
    for (i = 0; i < s->n; i++)
        for (j = 0; j < s->p.l; j++) {
            bqf_random(&s->a[i][j], x, st);
            bqf_chain(&s->e[i][j], &s->a[i][j], s->p.L, x);
        }
    for (j = 0; j < s->p.l; j++) {
        bqf_set(&s->P[j], &s->e[0][j]);
        for (i = 1; i < s->n; i++) bqf_mul(&s->P[j], &s->P[j], &s->e[i][j], x);
    }
    cgwots_pkhash(s->pk, s->P, s->p.l, x);
}

void hwots_nn_partial(bqf_t *out, const hwots_nn *s, int i, const unsigned *d, cg_ctx *x)
{
    int j;
    for (j = 0; j < s->p.l; j++) bqf_chain(&out[j], &s->a[i][j], d[j], x);
}

int hwots_nn_verify_partial(bqf_t *part, const hwots_nn *s, int i,
                            const unsigned *d, cg_ctx *x)
{
    bqf_t q; int j, ok = 1;
    bqf_init(&q);
    for (j = 0; j < s->p.l; j++) {
        bqf_chain(&q, &part[j], s->p.L - d[j], x);
        if (!bqf_equal(&q, &s->e[i][j])) { ok = 0; break; }
    }
    bqf_clear(&q);
    return ok;
}

void hwots_nn_combine(bqf_t *sig, bqf_t **parts, int n, const wparams *p, cg_ctx *x)
{
    int i, j;
    for (j = 0; j < p->l; j++) {
        bqf_set(&sig[j], &parts[0][j]);
        for (i = 1; i < n; i++) bqf_mul(&sig[j], &sig[j], &parts[i][j], x);
    }
}

/* ---------------- CG-HWOTS-k/n cumulative -------------------------------- */

long binom(int n, int r)
{
    long v = 1; int i;
    if (r < 0 || r > n) return 0;
    if (r > n - r) r = n - r;
    for (i = 0; i < r; i++) { v = v * (n - i) / (i + 1); }
    return v;
}

static void gen_subsets(int **U, int n, int r, long N)
{
    int *cur = (int *)malloc(sizeof(int) * (r > 0 ? r : 1));
    long idx = 0; int i;
    if (r == 0) { (void)cur; free(cur); return; }
    for (i = 0; i < r; i++) cur[i] = i;
    for (;;) {
        for (i = 0; i < r; i++) U[idx][i] = cur[i];
        idx++;
        if (idx >= N) break;
        i = r - 1;
        while (i >= 0 && cur[i] == n - r + i) i--;
        if (i < 0) break;
        cur[i]++;
        { int t; for (t = i + 1; t < r; t++) cur[t] = cur[t-1] + 1; }
    }
    free(cur);
}

void hwots_kn_init(hwots_kn *s, int n, int k, const wparams *p)
{
    long u;
    s->n = n; s->k = k; s->p = *p;
    s->N = binom(n, k - 1);
    s->U = (int **)malloc(sizeof(int *) * s->N);
    s->R = (bqf_t **)malloc(sizeof(bqf_t *) * s->N);
    s->C = (bqf_t **)malloc(sizeof(bqf_t *) * s->N);
    for (u = 0; u < s->N; u++) {
        s->U[u] = (int *)malloc(sizeof(int) * (k > 1 ? k - 1 : 1));
        s->R[u] = forms_alloc(p->l);
        s->C[u] = forms_alloc(p->l);
    }
    gen_subsets(s->U, n, k - 1, s->N);
    s->P = forms_alloc(p->l);
}

void hwots_kn_clear(hwots_kn *s)
{
    long u;
    for (u = 0; u < s->N; u++) {
        free(s->U[u]); forms_free(s->R[u], s->p.l); forms_free(s->C[u], s->p.l);
    }
    free(s->U); free(s->R); free(s->C);
    forms_free(s->P, s->p.l);
}

void hwots_kn_keygen(hwots_kn *s, cg_ctx *x, gmp_randstate_t st)
{
    long u; int j;
    for (u = 0; u < s->N; u++)
        for (j = 0; j < s->p.l; j++) {
            bqf_random(&s->R[u][j], x, st);
            bqf_chain(&s->C[u][j], &s->R[u][j], s->p.L, x);
        }
    for (j = 0; j < s->p.l; j++) {
        bqf_set(&s->P[j], &s->C[0][j]);
        for (u = 1; u < s->N; u++) bqf_mul(&s->P[j], &s->P[j], &s->C[u][j], x);
    }
    cgwots_pkhash(s->pk, s->P, s->p.l, x);
}

static int in_set(const int *set, int len, int v)
{
    int i; for (i = 0; i < len; i++) if (set[i] == v) return 1;
    return 0;
}

int hwots_kn_owner(const hwots_kn *s, long u, const int *T, int k)
{
    int t;
    for (t = 0; t < k; t++)
        if (!in_set(s->U[u], s->k - 1, T[t])) return T[t];
    return -1;   /* impossible for |T| = k */
}

void hwots_kn_partial(bqf_t *out, const hwots_kn *s, int i, const int *T,
                      const unsigned *d, cg_ctx *x)
{
    int j, first;
    long u;
    bqf_t B;
    bqf_init(&B);
    for (j = 0; j < s->p.l; j++) {
        first = 1;
        for (u = 0; u < s->N; u++) {
            if (hwots_kn_owner(s, u, T, s->k) != i) continue;
            if (first) { bqf_set(&B, &s->R[u][j]); first = 0; }
            else bqf_mul(&B, &B, &s->R[u][j], x);
        }
        if (first) bqf_identity(&B, x);
        bqf_chain(&out[j], &B, d[j], x);
    }
    bqf_clear(&B);
}

void hwots_kn_partial_endpoint(bqf_t *out, const hwots_kn *s, int i, const int *T,
                               cg_ctx *x)
{
    int j, first; long u;
    for (j = 0; j < s->p.l; j++) {
        first = 1;
        for (u = 0; u < s->N; u++) {
            if (hwots_kn_owner(s, u, T, s->k) != i) continue;
            if (first) { bqf_set(&out[j], &s->C[u][j]); first = 0; }
            else bqf_mul(&out[j], &out[j], &s->C[u][j], x);
        }
        if (first) bqf_identity(&out[j], x);
    }
}

void hwots_kn_combine(bqf_t *sig, bqf_t **parts, int k, const wparams *p, cg_ctx *x)
{
    int i, j;
    for (j = 0; j < p->l; j++) {
        bqf_set(&sig[j], &parts[0][j]);
        for (i = 1; i < k; i++) bqf_mul(&sig[j], &sig[j], &parts[i][j], x);
    }
}

/* ---------------- conventional SHA-256 WOTS (LMOTS-style) --------------- */

/* one chain step = SHA256(I(16) || q(4) || i(2) || j(1) || tmp(32)) = 55 bytes
 * = exactly one SHA-256 compression, as in RFC 8554 LMOTS. */
static void cw_step(unsigned char out[32], const unsigned char I[16],
                    unsigned q, unsigned i, unsigned j, const unsigned char in[32])
{
    unsigned char buf[55];
    memcpy(buf, I, 16);
    buf[16] = (unsigned char)(q >> 24); buf[17] = (unsigned char)(q >> 16);
    buf[18] = (unsigned char)(q >> 8);  buf[19] = (unsigned char)q;
    buf[20] = (unsigned char)(i >> 8);  buf[21] = (unsigned char)i;
    buf[22] = (unsigned char)j;
    memcpy(buf + 23, in, 32);
    sha256(buf, 55, out);
}

void cwots_key_init(cwots_key *k, const wparams *p)
{
    k->p = *p;
    k->sk = malloc(32 * (size_t)p->l);
    k->ep = malloc(32 * (size_t)p->l);
}
void cwots_key_clear(cwots_key *k){ free(k->sk); free(k->ep); }

void cwots_keygen(cwots_key *k, unsigned long long *hc)
{
    int i; unsigned j;
    unsigned char tmp[32];
    for (i = 0; i < 16; i++) k->I[i] = (unsigned char)(i * 7 + 1);
    for (i = 0; i < k->p.l; i++) {
        unsigned char seed[40];
        memcpy(seed, k->I, 16);
        seed[16] = (unsigned char)(i >> 8); seed[17] = (unsigned char)i;
        memset(seed + 18, 0x5a, 22);
        sha256(seed, 40, k->sk[i]);
        memcpy(tmp, k->sk[i], 32);
        for (j = 0; j < k->p.L; j++) { cw_step(tmp, k->I, 0, (unsigned)i, j, tmp); (*hc)++; }
        memcpy(k->ep[i], tmp, 32);
    }
    { sha256_ctx c; sha256_init(&c);
      sha256_update(&c, "CWOTS/pk/v1", 11);
      for (i = 0; i < k->p.l; i++) sha256_update(&c, k->ep[i], 32);
      sha256_final(&c, k->pk); }
}

void cwots_sign(unsigned char (*sig)[32], const cwots_key *k, const unsigned *d,
                unsigned long long *hc)
{
    int i; unsigned j;
    unsigned char tmp[32];
    for (i = 0; i < k->p.l; i++) {
        memcpy(tmp, k->sk[i], 32);
        for (j = 0; j < d[i]; j++) { cw_step(tmp, k->I, 0, (unsigned)i, j, tmp); (*hc)++; }
        memcpy(sig[i], tmp, 32);
    }
}

int cwots_verify(unsigned char (*sig)[32], const cwots_key *k, const unsigned *d,
                 unsigned long long *hc)
{
    int i; unsigned j;
    unsigned char tmp[32], pk[32];
    sha256_ctx c;
    sha256_init(&c);
    sha256_update(&c, "CWOTS/pk/v1", 11);
    for (i = 0; i < k->p.l; i++) {
        memcpy(tmp, sig[i], 32);
        for (j = d[i]; j < k->p.L; j++) { cw_step(tmp, k->I, 0, (unsigned)i, j, tmp); (*hc)++; }
        sha256_update(&c, tmp, 32);
    }
    sha256_final(&c, pk);
    return memcmp(pk, k->pk, 32) == 0;
}
