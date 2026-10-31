#include "sampler.h"
#include "zeroize.h"
#include "fixpoint.h"
#include "gauss_q80_constants.h"
#include "symmetric.h"
#include <stdint.h>
#include <arm_neon.h>
/*************************************************
 * Name:        rej_uniform
 *
 * Description: Sample uniformly random coefficients in [0, Q-1] by
 *              performing rejection sampling on array of random bytes.
 *
 * Arguments:   - int32_t *a: pointer to output array (allocated)
 *              - unsigned int len: number of coefficients to be sampled
 *              - const uint8_t *buf: array of random bytes
 *              - unsigned int buflen: length of array of random bytes
 *
 * Returns number of sampled coefficients. Can be smaller than len 
 * if not enough random bytes were given.
 **************************************************/
unsigned int rej_uniform(int32_t *a, unsigned int len,
                         const uint8_t *buf, unsigned int buflen)
{
    unsigned int ctr = 0, pos = 0;

#if Q_BITS == 18
    while (ctr < len && pos + 9 <= buflen) {
        uint32_t t0 =  (uint32_t)buf[pos]
                      | ((uint32_t)buf[pos+1] << 8)
                      | (((uint32_t)buf[pos+2] & 0x03) << 16);

        uint32_t t1 =  ((uint32_t)buf[pos+2] >> 2)
                      | ((uint32_t)buf[pos+3] << 6)
                      | (((uint32_t)buf[pos+4] & 0x0F) << 14);

        uint32_t t2 =  ((uint32_t)buf[pos+4] >> 4)
                      | ((uint32_t)buf[pos+5] << 4)
                      | (((uint32_t)buf[pos+6] & 0x3F) << 12);

        uint32_t t3 =  ((uint32_t)buf[pos+6] >> 6)
                      | ((uint32_t)buf[pos+7] << 2)
                      | ((uint32_t)buf[pos+8] << 10);

        pos += 9;

        if (t0 < Q && ctr < len) a[ctr++] = (int32_t)t0;
        if (t1 < Q && ctr < len) a[ctr++] = (int32_t)t1;
        if (t2 < Q && ctr < len) a[ctr++] = (int32_t)t2;
        if (t3 < Q && ctr < len) a[ctr++] = (int32_t)t3;
    }

#elif Q_BITS == 17
    while (ctr < len && pos + 17 <= buflen) {
        const uint8_t *b = buf + pos;
        uint32_t t[8];

        t[0] =  (uint32_t)b[0]        | ((uint32_t)b[1]  << 8) | (((uint32_t)b[2]  & 0x01) << 16);
        t[1] = ((uint32_t)b[2]  >> 1) | ((uint32_t)b[3]  << 7) | (((uint32_t)b[4]  & 0x03) << 15);
        t[2] = ((uint32_t)b[4]  >> 2) | ((uint32_t)b[5]  << 6) | (((uint32_t)b[6]  & 0x07) << 14);
        t[3] = ((uint32_t)b[6]  >> 3) | ((uint32_t)b[7]  << 5) | (((uint32_t)b[8]  & 0x0F) << 13);
        t[4] = ((uint32_t)b[8]  >> 4) | ((uint32_t)b[9]  << 4) | (((uint32_t)b[10] & 0x1F) << 12);
        t[5] = ((uint32_t)b[10] >> 5) | ((uint32_t)b[11] << 3) | (((uint32_t)b[12] & 0x3F) << 11);
        t[6] = ((uint32_t)b[12] >> 6) | ((uint32_t)b[13] << 2) | (((uint32_t)b[14] & 0x7F) << 10);
        t[7] = ((uint32_t)b[14] >> 7) | ((uint32_t)b[15] << 1) | ((uint32_t)b[16]  << 9);

        pos += 17;

        for (unsigned int i = 0; i < 8 && ctr < len; i++) {
            if (t[i] < Q)
                a[ctr++] = (int32_t)t[i];
        }
    }

#else
    #error "Q_BITS must be 17 or 18"
#endif
    {
        uint32_t acc = 0;
        unsigned int acc_bits = 0;
        const uint32_t mask = (1u << Q_BITS) - 1u;

        while (ctr < len) {
            while (acc_bits < Q_BITS && pos < buflen) {
                acc |= ((uint32_t)buf[pos++]) << acc_bits;
                acc_bits += 8;
            }
            if (acc_bits < Q_BITS) break;

            uint32_t t = acc & mask;
            acc >>= Q_BITS;
            acc_bits -= Q_BITS;

            if (t < Q)
                a[ctr++] = (int32_t)t;
        }
    }

    return ctr;
}

/*************************************************
 * Name:        rej_p
 *
 * Description: Sample coefficients in {-1,0,1} with
 *              Pr[X = 1] = Pr[X = -1] = P, Pr[X = 0] = 1 - 2P
 *              by performing rejection sampling on array of random bytes.
 *
 * Arguments:   - int32_t *a: pointer to output array (allocated)
 *              - unsigned int len: number of coefficients to be sampled
 *              - const uint8_t *buf: array of random bytes
 *              - unsigned int buflen: length of array of random bytes
 *
 * Returns number of sampled coefficients. Can be smaller than len 
 * if not enough random bytes were given.
 **************************************************/
static inline int32_t floor10(uint16_t t) {
    return (int32_t)((t * 205u) >> 11);
}

static inline int32_t mod10(uint16_t t) {
    int32_t r;
    r = (t & 31) + ((t >> 5) << 1);
    r = (r & 31) + ((r >> 5) << 1);
    r = (r & 15) + ((r >> 4) * 6);
    r -= 10 * (r >= 10);
    return r - 10 * (r >= 10);
}

unsigned int rej_p(int32_t *a, unsigned int len,
                   const uint8_t *buf, unsigned int buflen)
{
    unsigned int ctr = 0, pos = 0, acc_bits = 0;
    uint32_t acc = 0;

    if (P == 0.15 || P == 0.2)
    {
        const int32_t thresh = (P == 0.15) ? 3 : 4;

        while (ctr < len) {
            while (acc_bits < 13 && pos < buflen) {
                acc |= ((uint32_t)buf[pos++]) << acc_bits;
                acc_bits += 8;
            }
            if (acc_bits < 13) break;

            uint16_t t = acc & 0x3FFu;
            acc >>= 10;
            acc_bits -= 10;

            if (t >= 1000u) continue;

            int32_t nz, s;

            nz = (int32_t)(mod10(t) < thresh);
            s  = (int32_t)(acc & 1u) * 2 - 1;
            a[ctr++] = nz * s;
            if (ctr >= len) break;
            acc >>= (nz & 1);
            acc_bits -= (nz & 1);
            t = floor10(t);

            nz = (int32_t)(mod10(t) < thresh);
            s  = (int32_t)(acc & 1u) * 2 - 1;
            a[ctr++] = nz * s;
            if (ctr >= len) break;
            acc >>= (nz & 1);
            acc_bits -= (nz & 1);
            t = floor10(t);

            nz = (int32_t)(t < thresh);
            s  = (int32_t)(acc & 1u) * 2 - 1;
            a[ctr++] = nz * s;
            if (ctr >= len) break;
            acc >>= (nz & 1);
            acc_bits -= (nz & 1);
        }
    }
    else if (P == 0.3125)
    {
        while (ctr < len) {
            while (acc_bits < 15 && pos < buflen) {
                acc |= ((uint32_t)buf[pos++]) << acc_bits;
                acc_bits += 8;
            }
            if (acc_bits < 15) break;

            uint32_t t = acc & 0xFFFu;
            acc >>= 12;
            acc_bits -= 12;

            int32_t nz, s;
            uint32_t d;

            d  = t & 0xFu;
            nz = (int32_t)(d < 10u);
            s  = (int32_t)(acc & 1u) * 2 - 1;
            a[ctr++] = nz * s;
            if (ctr >= len) break;
            acc >>= (nz & 1);
            acc_bits -= (nz & 1);

            d  = (t >> 4) & 0xFu;
            nz = (int32_t)(d < 10u);
            s  = (int32_t)(acc & 1u) * 2 - 1;
            a[ctr++] = nz * s;
            if (ctr >= len) break;
            acc >>= (nz & 1);
            acc_bits -= (nz & 1);

            d  = (t >> 8) & 0xFu;
            nz = (int32_t)(d < 10u);
            s  = (int32_t)(acc & 1u) * 2 - 1;
            a[ctr++] = nz * s;
            if (ctr >= len) break;
            acc >>= (nz & 1);
            acc_bits -= (nz & 1);
        }
    }

    return ctr;
}

/* Q63 one-sided degree-12 approximation of exp(-x). */
static inline __attribute__((always_inline)) uint64_t approx_exp(const uint64_t x) {
    int64_t result = EXP_C12;
    uint64_t threshold;

    result = smulh63_ceil(result, x) + EXP_C11;
    result = smulh63_ceil(result, x) + EXP_C10;
    result = smulh63_ceil(result, x) + EXP_C9;
    result = smulh63_ceil(result, x) + EXP_C8;
    result = smulh63_ceil(result, x) + EXP_C7;
    result = smulh63_ceil(result, x) + EXP_C6;
    result = smulh63_ceil(result, x) + EXP_C5;
    result = smulh63_ceil(result, x) + EXP_C4;
    result = smulh63_ceil(result, x) + EXP_C3;
    result = smulh63_ceil(result, x) + EXP_C2;
    result = smulh63_ceil(result, x) + EXP_C1;
    result = smulh63_ceil(result, x);

    /* Unsigned addition gives EXP_C0 + result even when result is negative. */
    threshold = EXP_C0 + (uint64_t)result;
    if (threshold > (UINT64_C(1) << EXP_FRAC_BITS))
        threshold = UINT64_C(1) << EXP_FRAC_BITS;
    return threshold;
}

/* Constant-time inverse-CDF sampling from the generated 112-bit CDT. */
static inline __attribute__((always_inline)) uint64_t sample_gauss112(const uint64_t rand_lo,
                                const uint64_t rand_hi) {
    unsigned int i;
    uint64_t r = 0;

    for (i = 0; i < CDTLEN; i++) {
        const uint64_t borrow = rand_lo > CDT112_LO[i];
        const uint64_t hi_diff = CDT112_HI[i] - rand_hi - borrow;
        r += hi_diff >> 63;
    }
    return r;
}

/* Fused CDF scan of four attempts in one NEON loop over the shared CDT
 * tables: each attempt owns a 2x64-bit-lane accumulator whose lanes cover
 * the even/odd CDT entries of the block, and the rand value is splatted
 * into both lanes. Per (entry, attempt) the lane arithmetic is identical
 * to the scalar loop above (same borrow, same subtract, same sign test);
 * the two lane sums are combined with one exact horizontal 64-bit add, so
 * each x[k] equals the scalar count bit-for-bit. CDTLEN is odd (195): the
 * last entry is handled by the scalar tail. */
static inline __attribute__((always_inline)) void sample_gauss112_4x(
    uint64_t x[4], const uint64_t rand_lo[4], const uint64_t rand_hi[4]) {
    uint64x2_t acc[4];
    uint64x2_t vrh[4], vrl[4];
    unsigned int i;
    unsigned int k;

    for (k = 0; k < 4; k++) {
        acc[k] = vdupq_n_u64(0);
        vrl[k] = vdupq_n_u64(rand_lo[k]);
        vrh[k] = vdupq_n_u64(rand_hi[k]);
    }

    for (i = 0; i + 2 <= CDTLEN; i += 2) {
        const uint64x2_t vlo = vld1q_u64(&CDT112_LO[i]);
        const uint64x2_t vhi = vld1q_u64(&CDT112_HI[i]);
        uint64x2_t d;

        for (k = 0; k < 4; k++) {
            d = vsubq_u64(vhi, vrh[k]);
            d = vsubq_u64(d, vshrq_n_u64(vcgtq_u64(vrl[k], vlo), 63));
            acc[k] = vaddq_u64(acc[k], vshrq_n_u64(d, 63));
        }
    }
    x[0] = vaddvq_u64(acc[0]);
    x[1] = vaddvq_u64(acc[1]);
    x[2] = vaddvq_u64(acc[2]);
    x[3] = vaddvq_u64(acc[3]);

    for (; i < CDTLEN; i++) {
        for (k = 0; k < 4; k++) {
            const uint64_t borrow = rand_lo[k] > CDT112_LO[i];
            const uint64_t hi_diff = CDT112_HI[i] - rand_hi[k] - borrow;
            x[k] += hi_diff >> 63;
        }
    }
}

static inline __attribute__((always_inline)) uint64_t load64_le(const uint8_t src[8]) {
    return (uint64_t)src[0] | ((uint64_t)src[1] << 8) |
           ((uint64_t)src[2] << 16) | ((uint64_t)src[3] << 24) |
           ((uint64_t)src[4] << 32) | ((uint64_t)src[5] << 40) |
           ((uint64_t)src[6] << 48) | ((uint64_t)src[7] << 56);
}

#define GAUSS_RAND (72 + CDT_BITS + 64)
#define GAUSS_RAND_BYTES ((GAUSS_RAND + 7) / 8) /* 31 */

/* File-private inlined copy of fixpoint_square (body verbatim from
 * fixpoint.c): the batched retry loop needs it inlined into gauss_finish,
 * but the header-inline version changed the compilation layout of every
 * including TU (128 keygen -2.5% via the fixpoint.c Newton loop). Only
 * this sampler TU sees the inline; all other callers keep the external
 * definition. */
static inline __attribute__((always_inline)) void fixpoint_square_sampler(
    fp96_76 *sqx, const fp96_76 *x) {
    uint64_t tmp[2];
    sq48(&sqx->limb48[0], x->limb48[0]);

    // shift right by 48, rounding
    //sqx->limb48[0] += 1ULL << 47;
    sqx->limb48[0] >>= 48;
    sqx->limb48[0] += sqx->limb48[1];

    // mul
    mul48(tmp, x->limb48[0], x->limb48[1]);
    sqx->limb48[0] += tmp[0] << 1;
    sqx->limb48[1] = tmp[1] << 1;

    // shift right by 28, rounding
    //sqx->limb48[0] += 1ULL << 27;
    sqx->limb48[0] >>= 28;
    sqx->limb48[0] += (sqx->limb48[1] << 20) & ((1ULL << 48) - 1);
    sqx->limb48[1] >>= 28;

    sq64(tmp, x->limb48[1]);
    sqx->limb48[0] += (tmp[0] << 20) & ((1ULL << 48) - 1);
    sqx->limb48[1] += (tmp[0] >> 28) + (tmp[1] << 36);

    renormalize(sqx);
}

/* Everything of sample_gauss_sigma76 after the CDT lookup: y construction,
 * rounding, y*y, exp input, rejection draw. Split out so the batched loop
 * can share one fused CDT scan across four attempts. Expressions are
 * identical to the scalar version. */
static inline __attribute__((always_inline)) int gauss_finish(
    uint64_t *r, fp96_76 *sqr, const uint8_t rand[GAUSS_RAND_BYTES],
    const uint64_t rejection_word, const uint64_t x) {
    const uint64_t rand_rej = rejection_word & INT64_MAX;
    const uint64_t zero_coin = rejection_word >> 63;
    uint64_t exp_in, exp_threshold, diff_hi;
    fp96_76 y;

    // y := append x to y
    // leave 16 bit for carries
    y.limb48[0] = load64_le(rand + 22) & ((UINT64_C(1) << 48) - 1);
    y.limb48[1] = rand[28] | ((uint64_t)rand[29] << 8) |
                  ((uint64_t)rand[30] << 16) | (x << 24);

    /* Round y without losing the high bit when x >= 128. */
    *r = ((y.limb48[0] >> 15) + 1) >> 1;
    *r += y.limb48[1] << 32;

    // sqr := y*y
    fixpoint_square_sampler(sqr, &y);

    /* If D=(sqr-x^2*2^68), Q63 exp_in is round(D/2^14). */
    diff_hi = sqr->limb48[1] - ((x * x) << 20);
    exp_in = (diff_hi << 34) | (sqr->limb48[0] >> 14);
    exp_in += (sqr->limb48[0] >> 13) & 1;

    exp_threshold = approx_exp(exp_in);
    return (int)(((rand_rej - exp_threshold) >> 63) &
                 (((*r | (UINT64_C(0) - *r)) >> 63) | zero_coin));
}

static inline __attribute__((always_inline)) int sample_gauss_sigma76(uint64_t *r, fp96_76 *sqr,
                                const uint8_t rand[GAUSS_RAND_BYTES]) {
    const uint64_t rand_gauss_lo = load64_le(rand);
    const uint64_t rand_gauss_hi = load64_le(rand + 8) &
                                   ((UINT64_C(1) << 48) - 1);
    const uint64_t rejection_word = load64_le(rand + 14);
    const uint64_t x = sample_gauss112(rand_gauss_lo, rand_gauss_hi);

    return gauss_finish(r, sqr, rand, rejection_word, x);
}

int sample_gauss(uint64_t *r, fp96_76 *sqsum, const uint8_t *buf, const size_t buflen, const size_t len, const int dont_write_last)
{
    const uint8_t *pos = buf;
    fp96_76 sqr;
    size_t bytecnt = buflen, coefcnt = 0;
    int accepted;
    uint64_t dummy;

    /* 4-way attempt batching. The guard coefcnt + 4 <= len keeps all
     * written slots coefcnt..coefcnt+3 strictly below len, and only the
     * attempt-3 slot can reach slot len-1 (prefix p_k <= k <= 3). */
    while (coefcnt + 4 <= len && bytecnt >= 4 * GAUSS_RAND_BYTES) {
        uint64_t r0, r1, r2, r3;
        fp96_76 s0, s1, s2, s3;
        uint64_t x[4];
        uint64_t rand_lo[4], rand_hi[4], rej[4];
        int a0, a1, a2, a3;
        uint64_t p1, p2, p3;

        for (unsigned int k = 0; k < 4; k++) {
            rand_lo[k] = load64_le(pos + k * GAUSS_RAND_BYTES);
            rand_hi[k] = load64_le(pos + k * GAUSS_RAND_BYTES + 8) &
                         ((UINT64_C(1) << 48) - 1);
            rej[k] = load64_le(pos + k * GAUSS_RAND_BYTES + 14);
        }
        sample_gauss112_4x(x, rand_lo, rand_hi);
        a0 = gauss_finish(&r0, &s0, pos, rej[0], x[0]);
        a1 = gauss_finish(&r1, &s1, pos + GAUSS_RAND_BYTES, rej[1], x[1]);
        a2 = gauss_finish(&r2, &s2, pos + 2 * GAUSS_RAND_BYTES, rej[2], x[2]);
        a3 = gauss_finish(&r3, &s3, pos + 3 * GAUSS_RAND_BYTES, rej[3], x[3]);
        p1 = (uint64_t)a0;
        p2 = p1 + (uint64_t)a1;
        p3 = p2 + (uint64_t)a2;
        /* Write slots = prefix-accept counts (a rejected attempt's dead
         * write is overwritten by the next accepted attempt, matching the
         * sequential loop; only slots < coefcnt+4 are ever read). */
        r[coefcnt] = r0;
        r[coefcnt + p1] = r1;
        r[coefcnt + p2] = r2;
        if (!(dont_write_last && coefcnt + p3 == len - 1)) {
            r[coefcnt + p3] = r3; /* replicate the dummy branch for slot len-1 */
        }
        coefcnt += p3 + (uint64_t)a3;
        pos += 4 * GAUSS_RAND_BYTES;
        bytecnt -= 4 * GAUSS_RAND_BYTES;

        sqsum->limb48[0] += (s0.limb48[0] & -(int64_t)a0) +
                            (s1.limb48[0] & -(int64_t)a1) +
                            (s2.limb48[0] & -(int64_t)a2) +
                            (s3.limb48[0] & -(int64_t)a3);
        sqsum->limb48[1] += (s0.limb48[1] & -(int64_t)a0) +
                            (s1.limb48[1] & -(int64_t)a1) +
                            (s2.limb48[1] & -(int64_t)a2) +
                            (s3.limb48[1] & -(int64_t)a3);
    }

    while (coefcnt < len) {
        if (bytecnt < GAUSS_RAND_BYTES) {
          renormalize(sqsum);
          return coefcnt;
        }

        if (dont_write_last && coefcnt == len-1)
        {
          accepted = sample_gauss_sigma76(&dummy, &sqr, pos);
        } else {
          accepted = sample_gauss_sigma76(&r[coefcnt], &sqr, pos);
        }
        coefcnt += accepted;
        pos += GAUSS_RAND_BYTES;
        bytecnt -= GAUSS_RAND_BYTES;

        sqsum->limb48[0] += sqr.limb48[0] & -(int64_t)accepted;
        sqsum->limb48[1] += sqr.limb48[1] & -(int64_t)accepted;
    }

    renormalize(sqsum);
    return len;
}

#define POLY_HYPERBALL_BUFLEN (GAUSS_RAND_BYTES * N)
#define POLY_HYPERBALL_NBLOCKS ((POLY_HYPERBALL_BUFLEN + SM3_XOF_BLOCKBYTES - 1) / SM3_XOF_BLOCKBYTES)

/*************************************************
 * Name:        sample_gauss_N
 *
 * Description: Sample len coefficients from discrete Gaussian with
 *              standard deviation sigma=76/2^13
 *              using output of SM3-XOF(seed|nonce).
 *              Also compute the squared norm of the sampled coefficients.
 *
 * Arguments:   - uint64_t *r: pointer to output array of sampled coefficients
 *              - uint8_t *signs: pointer to output array of sign bits
 *              - fp96_76 *sqsum: pointer to output squared norm
 *              - const uint8_t seed[CRHBYTES]: byte array with seed of length CRHBYTES
 *              - const uint16_t nonce: 2-byte nonce
 *              - const size_t len: number of coefficients to be sampled
 **************************************************/
void sample_gauss_N(uint64_t *r, uint8_t *signs, fp96_76 *sqsum,
                    const uint8_t seed[CRHBYTES], const uint16_t nonce,
                    const size_t len) {
    uint8_t buf[POLY_HYPERBALL_NBLOCKS * SM3_XOF_BLOCKBYTES];
    size_t bytecnt, coefcnt, firstflag = 1;
    sm3_xof_state state;
    sm3_xof_stream_init(&state, seed, CRHBYTES, nonce);

    sm3_xof_squeezeblocks(buf, POLY_HYPERBALL_NBLOCKS, &state);
    for (size_t i = 0; i < len / 8; i++) {
        signs[i] = buf[i];
    }
    bytecnt = POLY_HYPERBALL_NBLOCKS * SM3_XOF_BLOCKBYTES - len / 8;
    coefcnt = sample_gauss(r, sqsum, buf + len / 8, bytecnt, len, len%N);
    while (coefcnt < len) {
        size_t off = bytecnt % GAUSS_RAND_BYTES;
        for (size_t i = 0; i < off; i++) {
            buf[i] = buf[bytecnt + len/8*firstflag - off + i];
        }
        sm3_xof_squeezeblocks(buf + off, 1, &state);
        bytecnt = SM3_XOF_BLOCKBYTES + off;

        coefcnt += sample_gauss(r + coefcnt, sqsum, buf, bytecnt, len - coefcnt, len%N);
        firstflag = 0;
    }

    /* Erase the Gaussian sampling stream (secret y material). */
    darts_zeroize(buf, sizeof(buf));
}
