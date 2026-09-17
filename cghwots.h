/* cghwots.h - the schemes of sections 5, 6, 7 of the report, plus a
 * conventional LMOTS-style SHA-256 Winternitz OTS used as the baseline. */
#ifndef CGHWOTS_H
#define CGHWOTS_H

#include "bqf.h"
#include "sha256.h"

/* ---- Winternitz encoding (section 3.1) --------------------------------- */
typedef struct {
    int w;              /* log2(B)                     */
    unsigned B;         /* base                        */
    unsigned L;         /* B-1, max chain length       */
    int l1, l2, l;      /* message / checksum / total coordinates */
    int digest_bits;
} wparams;

void wots_params(wparams *p, int w, int digest_bits);
void wots_encode(unsigned *d, const unsigned char *digest, const wparams *p);
unsigned long wots_sum(const unsigned *d, const wparams *p);          /* D(d) */

/* ---- base one-time signature CG-WOTS (section 5) ----------------------- */
typedef struct {
    wparams p;
    bqf_t *A;                  /* secret chain starts, l elements */
    bqf_t *P;                  /* endpoints F^L(A_j)              */
    unsigned char pk[32];
} cgwots_key;

void cgwots_key_init(cgwots_key *k, const wparams *p);
void cgwots_key_clear(cgwots_key *k);
void cgwots_keygen(cgwots_key *k, cg_ctx *x, gmp_randstate_t st);
void cgwots_pkhash(unsigned char out[32], bqf_t *P, int l, cg_ctx *x);
void cgwots_sign(bqf_t *sig, const cgwots_key *k, const unsigned *d, cg_ctx *x);
int  cgwots_verify(bqf_t *sig, const unsigned char pk[32], const unsigned *d,
                   const wparams *p, cg_ctx *x);

/* ---- CG-HWOTS-n/n (section 6) ------------------------------------------ */
typedef struct {
    int n;
    wparams p;
    bqf_t **a;      /* a[i][j]  participant secrets   */
    bqf_t **e;      /* e[i][j]  participant endpoints */
    bqf_t *P;       /* aggregate endpoints            */
    unsigned char pk[32];
} hwots_nn;

void hwots_nn_init(hwots_nn *s, int n, const wparams *p);
void hwots_nn_clear(hwots_nn *s);
void hwots_nn_keygen(hwots_nn *s, cg_ctx *x, gmp_randstate_t st);
void hwots_nn_partial(bqf_t *out, const hwots_nn *s, int i, const unsigned *d, cg_ctx *x);
int  hwots_nn_verify_partial(bqf_t *part, const hwots_nn *s, int i,
                             const unsigned *d, cg_ctx *x);
void hwots_nn_combine(bqf_t *sig, bqf_t **parts, int n, const wparams *p, cg_ctx *x);

/* ---- CG-HWOTS-k/n cumulative (section 7) ------------------------------- */
typedef struct {
    int n, k;
    long N;          /* C(n, k-1) maximal unauthorised sets */
    wparams p;
    int **U;         /* U[u][0..k-2] the unauthorised sets  */
    bqf_t **R;       /* R[u][j] secret factors              */
    bqf_t **C;       /* C[u][j] = F^L(R[u][j])              */
    bqf_t *P;
    unsigned char pk[32];
} hwots_kn;

long binom(int n, int r);
void hwots_kn_init(hwots_kn *s, int n, int k, const wparams *p);
void hwots_kn_clear(hwots_kn *s);
void hwots_kn_keygen(hwots_kn *s, cg_ctx *x, gmp_randstate_t st);
/* owner of unauthorised set u for active quorum T (sorted, size k) */
int  hwots_kn_owner(const hwots_kn *s, long u, const int *T, int k);
void hwots_kn_partial(bqf_t *out, const hwots_kn *s, int i, const int *T,
                      const unsigned *d, cg_ctx *x);
void hwots_kn_partial_endpoint(bqf_t *out, const hwots_kn *s, int i, const int *T,
                               cg_ctx *x);
void hwots_kn_combine(bqf_t *sig, bqf_t **parts, int k, const wparams *p, cg_ctx *x);

/* ---- conventional SHA-256 Winternitz OTS (LMOTS-style baseline) -------- */
typedef struct {
    wparams p;
    unsigned char I[16];
    unsigned char (*sk)[32];
    unsigned char (*ep)[32];
    unsigned char pk[32];
} cwots_key;

void cwots_key_init(cwots_key *k, const wparams *p);
void cwots_key_clear(cwots_key *k);
void cwots_keygen(cwots_key *k, unsigned long long *hash_calls);
void cwots_sign(unsigned char (*sig)[32], const cwots_key *k, const unsigned *d,
                unsigned long long *hash_calls);
int  cwots_verify(unsigned char (*sig)[32], const cwots_key *k, const unsigned *d,
                  unsigned long long *hash_calls);

#endif
