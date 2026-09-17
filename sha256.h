/* sha256.h - compact portable SHA-256 */
#ifndef SHA256_H
#define SHA256_H
#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint32_t h[8];
    uint64_t len;
    unsigned char buf[64];
    size_t buflen;
} sha256_ctx;

void sha256_init(sha256_ctx *c);
void sha256_update(sha256_ctx *c, const void *data, size_t len);
void sha256_final(sha256_ctx *c, unsigned char out[32]);
void sha256(const void *data, size_t len, unsigned char out[32]);

/* single-block compression on an exactly 64-byte input; used for the
   conventional WOTS chain-step micro-benchmark (one compression call). */
void sha256_block64(const unsigned char in[64], unsigned char out[32]);

#endif
