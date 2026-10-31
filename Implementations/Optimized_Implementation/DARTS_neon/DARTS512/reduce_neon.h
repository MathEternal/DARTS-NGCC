#ifndef SIGN_REDUCE_NEON_H
#define SIGN_REDUCE_NEON_H

#include "params.h"
#include <stdint.h>

#if !defined(__ARM_NEON) && !defined(__ARM_NEON__)
#error "ARM_NEON_Implementation/DARTS requires AArch64/ARM NEON. Build with ARCH=armv8 on an ARM target or use an AArch64 cross compiler."
#endif

#include <arm_neon.h>

#if DARTS_MODE == 128 || DARTS_MODE == 256

#define MONT 114369
#define MONTSQ 7148
#define QINV 4244570369
#define QREC 32831
#define DQREC 16415

#elif DARTS_MODE == 512

#define MONT 130976
#define MONTSQ 125151
#define QINV 2820670977
#define QREC 16480
#define DQREC 8240

#endif

#define montgomery_reduce DARTS_NAMESPACE(montgomery_reduce)
#define caddq DARTS_NAMESPACE(caddq)
#define freeze DARTS_NAMESPACE(freeze)
#define freeze_centered DARTS_NAMESPACE(freeze_centered)

static inline int32x4_t darts_neon_mulhi_s32(int32x4_t a, int32x4_t b) {
    int64x2_t lo = vmull_s32(vget_low_s32(a), vget_low_s32(b));
    int64x2_t hi = vmull_s32(vget_high_s32(a), vget_high_s32(b));
    return vcombine_s32(vshrn_n_s64(lo, 32), vshrn_n_s64(hi, 32));
}

/* DARTS v2 perf NTT 腿（DARTS 参数版推导，恒等式经 1.6M 对实证零失配）：
 * mr(a*b) = (SQRDMULH(a,b) − SQRDMULH(t,Q)) >> 1，t = lo32(a*b)*QINV。
 *
 * 推导：标量 fqmul(a,b) = (p − t*Q)>>32（p = a*b 精确 64 位积）。对任意
 * 64 位 x 取欧几里得分解 x = h*2^32 + l（h = x>>32，l ∈ [0,2^32)），ARM
 * SQRDMULH(u,v) = (2uv + 2^31)>>32 = 2*hi(uv) + g(l(uv))，其中
 * g(l) = [l ≥ 2^30] + [l ≥ 2^32−2^30]。因 QINV*Q ≡ 1 (mod 2^32)（q=
 * 130817: 4244570369；q=260609: 2820670977，数值核验），t*Q ≡ a*b
 * (mod 2^32) ⇒ l(t*Q) = l(a*b) ⇒ g 项在差中精确相消：
 *   SQRDMULH(a,b) − SQRDMULH(t,Q) = 2*(hi(p) − hi(t*Q)) = 2*fqmul(a,b)
 * 故右移 1 位（算术）即得，逐位等于 darts_scalar_fqmul。
 *
 * 前提（DARTS 参数重推）：调用点 |a*b| ≤ 2^41（ntt/invntt 中间值 ≤ 2^20 与
 * |zeta| < Q 相乘；basemul 输入 < 2^21 相乘），远低于 SQRDMULH 饱和界
 * |uv| < 2^62，无饱和；vmulq_s32 低 32 位与标量 (int32_t) 截断同义。
 * 常时：固定指令序列，无数据相关分支。
 * 实证：160 万随机对（|a|,|b| ∈ {2^20,2^21,2^24,2^27} 档，双 q），与
 * 标量 mr 零失配（perf_ntt_evidence 恒等式门）。
 */
static inline int32x4_t darts_neon_fqmul(int32x4_t a, int32x4_t b) {
    const int32x4_t v_qinv = vdupq_n_s32((int32_t)QINV);
    const int32x4_t v_q = vdupq_n_s32(Q);
    int32x4_t t = vmulq_s32(vmulq_s32(a, b), v_qinv);
    return vshrq_n_s32(vsubq_s32(vqrdmulhq_s32(a, b),
                                 vqrdmulhq_s32(t, v_q)), 1);
}

static inline int32x4_t darts_neon_montgomery_reduce(int32x4_t a) {
    return darts_neon_fqmul(a, vdupq_n_s32(1));
}

static inline int32_t darts_scalar_montgomery_reduce(int64_t a) {
    int32_t t = (int64_t)(int32_t)a * QINV;
    return (int32_t)((a - (int64_t)t * Q) >> 32);
}

static inline __attribute__((always_inline))
int32_t darts_scalar_fqmul(int32_t a, int32_t b) {
    return darts_scalar_montgomery_reduce((int64_t)a * b);
}

static inline int32_t darts_neon_fqmul_lane(int32_t a, int32_t b) {
    return darts_scalar_fqmul(a, b);
}

static inline __attribute__((always_inline))
int32_t darts_scalar_freeze(int32_t a) {
    int64_t t = (int64_t)a * QREC;
    t >>= 32;
    t = a - t * Q;
    t += (t >> 31) & DQ;
    t -= ~((t - Q) >> 31) & Q;
    return (int32_t)t;
}

static inline __attribute__((always_inline))
int32_t darts_scalar_fqinv(int32_t a) {
#if DARTS_MODE == 128 || DARTS_MODE == 256
    int32_t t2, t3;

    t2 = darts_scalar_fqmul(a, a);
    t2 = darts_scalar_fqmul(t2, a);
    t3 = darts_scalar_fqmul(t2, t2);
    t3 = darts_scalar_fqmul(t3, t3);
    t2 = darts_scalar_fqmul(t3, t2);

    t3 = darts_scalar_fqmul(t2, t2);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t2 = darts_scalar_fqmul(t3, t2);

    t3 = darts_scalar_fqmul(t2, t2);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t3);

    return darts_scalar_fqmul(t3, t2);
#elif DARTS_MODE == 512
    int32_t t2, t3, t4;

    t2 = darts_scalar_fqmul(a, a);
    t2 = darts_scalar_fqmul(t2, a);

    t3 = darts_scalar_fqmul(t2, t2);
    t3 = darts_scalar_fqmul(t3, a);

    t4 = darts_scalar_fqmul(t3, t3);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t3);

    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, a);

    t3 = darts_scalar_fqmul(t4, t4);
    t3 = darts_scalar_fqmul(t3, t3);
    t3 = darts_scalar_fqmul(t3, t2);

    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);
    t4 = darts_scalar_fqmul(t4, t4);

    return darts_scalar_fqmul(t4, t3);
#endif
}

static inline int32x4_t darts_neon_caddq(int32x4_t a) {
    const int32x4_t v_q = vdupq_n_s32(Q);
    return vaddq_s32(a, vandq_s32(vshrq_n_s32(a, 31), v_q));
}

static inline int32x4_t darts_neon_freeze(int32x4_t a) {
    const int32x4_t v_qrec = vdupq_n_s32(QREC);
    const int32x4_t v_q = vdupq_n_s32(Q);
    const int32x4_t v_dq = vdupq_n_s32(DQ);
    int32x4_t t = darts_neon_mulhi_s32(a, v_qrec);
    t = vsubq_s32(a, vmulq_s32(t, v_q));
    t = vaddq_s32(t, vandq_s32(vshrq_n_s32(t, 31), v_dq));
    int32x4_t geq_q = vmvnq_s32(vshrq_n_s32(vsubq_s32(t, v_q), 31));
    return vsubq_s32(t, vandq_s32(geq_q, v_q));
}

static inline int32x4_t darts_neon_freeze_centered(int32x4_t a) {
    const int32x4_t v_q = vdupq_n_s32(Q);
    const int32x4_t v_half_q = vdupq_n_s32(Q >> 1);
    int32x4_t t = darts_neon_freeze(a);
    uint32x4_t gt_half = vcgtq_s32(t, v_half_q);
    return vsubq_s32(t, vandq_s32(vreinterpretq_s32_u32(gt_half), v_q));
}

static inline int32_t montgomery_reduce(int64_t a) {
    return darts_scalar_montgomery_reduce(a);
}

static inline int32_t caddq(int32_t a) {
    return vgetq_lane_s32(darts_neon_caddq(vdupq_n_s32(a)), 0);
}

static inline int32_t freeze(int32_t a) {
    return darts_scalar_freeze(a);
}

static inline int32_t freeze_centered(int32_t a) {
    int32_t t = freeze(a);
    int32_t b = t - (Q >> 1);
    t -= Q & -(b > 0);
    return t;
}

#endif
