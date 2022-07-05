/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_FIELDS_H__
#define __BLS12_381_ASM_FIELDS_H__

#include "blst_consts.h"
#include "blst_vect.h"

#ifndef __CUDA_ARCHI__
/*
 * BLS12-381-specifc Fp shortcuts to assembly.
 */

static const limb_t blst_p0 = (limb_t)0x89f3fffcfffcfffd; /* -1/P*/

static inline void sqr_n_mul_fp(vec384 out, const vec384 a, size_t count,
                                const vec384 b)
{
  sqr_n_mul_mont_383(out, a, count, BLS12_381_P, blst_p0, b);
}

static inline void add_fp(vec384 ret, const vec384 a, const vec384 b)
{
  add_mod_384(ret, a, b, BLS12_381_P);
}

static inline void sub_fp(vec384 ret, const vec384 a, const vec384 b)
{
  sub_mod_384(ret, a, b, BLS12_381_P);
}

static inline void mul_by_3_fp(vec384 ret, const vec384 a)
{
  mul_by_3_mod_384(ret, a, BLS12_381_P);
}

static inline void mul_by_8_fp(vec384 ret, const vec384 a)
{
  mul_by_8_mod_384(ret, a, BLS12_381_P);
}

static inline void lshift_fp(vec384 ret, const vec384 a, size_t count)
{
  lshift_mod_384(ret, a, count, BLS12_381_P);
}

static inline void rshift_fp(vec384 ret, const vec384 a, size_t count)
{
  rshift_mod_384(ret, a, count, BLS12_381_P);
}

static inline void div_by_2_fp(vec384 ret, const vec384 a)
{
  div_by_2_mod_384(ret, a, BLS12_381_P);
}

static inline void mul_fp(vec384 ret, const vec384 a, const vec384 b)
{
  mul_mont_384(ret, a, b, BLS12_381_P, blst_p0);
}

static inline void sqr_fp(vec384 ret, const vec384 a)
{
  sqr_mont_384(ret, a, BLS12_381_P, blst_p0);
}

static inline void cneg_fp(vec384 ret, const vec384 a, bool_t flag)
{
  cneg_mod_384(ret, a, flag, BLS12_381_P);
}

static inline void from_fp(vec384 ret, const vec384 a)
{
  from_mont_384(ret, a, BLS12_381_P, blst_p0);
}

static inline void redc_fp(vec384 ret, const vec768 a)
{
  redc_mont_384(ret, a, BLS12_381_P, blst_p0);
}

/*
 * BLS12-381-specifc Fp2 shortcuts to assembly.
 */
static inline void add_fp2(vec384x ret, const vec384x a, const vec384x b)
{
  add_mod_384x(ret, a, b, BLS12_381_P);
}

static inline void sub_fp2(vec384x ret, const vec384x a, const vec384x b)
{
  sub_mod_384x(ret, a, b, BLS12_381_P);
}

static inline void mul_by_3_fp2(vec384x ret, const vec384x a)
{
  mul_by_3_mod_384x(ret, a, BLS12_381_P);
}

static inline void mul_by_8_fp2(vec384x ret, const vec384x a)
{
  mul_by_8_mod_384x(ret, a, BLS12_381_P);
}

static inline void lshift_fp2(vec384x ret, const vec384x a, size_t count)
{
  lshift_mod_384(ret[0], a[0], count, BLS12_381_P);
  lshift_mod_384(ret[1], a[1], count, BLS12_381_P);
}

static inline void mul_fp2(vec384x ret, const vec384x a, const vec384x b)
{
  mul_mont_384x(ret, a, b, BLS12_381_P, blst_p0);
}

static inline void sqr_fp2(vec384x ret, const vec384x a)
{
  sqr_mont_384x(ret, a, BLS12_381_P, blst_p0);
}

static inline void cneg_fp2(vec384x ret, const vec384x a, bool_t flag)
{
  cneg_mod_384(ret[0], a[0], flag, BLS12_381_P);
  cneg_mod_384(ret[1], a[1], flag, BLS12_381_P);
}

#define vec_load_global vec_copy

statik void reciprocal_fp(vec384 out, const vec384 inp);
statik void flt_reciprocal_fp(vec384 out, const vec384 inp);
statik bool_t recip_sqrt_fp(vec384 out, const vec384 inp);
statik bool_t sqrt_fp(vec384 out, const vec384 inp);

statik void reciprocal_fp2(vec384x out, const vec384x inp);
statik void flt_reciprocal_fp2(vec384x out, const vec384x inp);
statik bool_t recip_sqrt_fp2(vec384x out, const vec384x inp,
                             const vec384x recip_ZZZ, const vec384x magic_ZZZ);
statik bool_t sqrt_fp2(vec384x out, const vec384x inp);
statik bool_t sqrt_align_fp2(vec384x out, const vec384x ret, const vec384x sqrt,
                             const vec384x inp);

typedef vec384x vec384fp2;
typedef vec384fp2 vec384fp6[3];
typedef vec384fp6 vec384fp12[2];

statik void sqr_fp12(vec384fp12 ret, const vec384fp12 a);
statik void cyclotomic_sqr_fp12(vec384fp12 ret, const vec384fp12 a);
statik void mul_fp12(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
statik void mul_by_xy00z0_fp12(vec384fp12 ret, const vec384fp12 a,
                               const vec384fp6 xy00z0);
statik void conjugate_fp12(vec384fp12 a);
statik void inverse_fp12(vec384fp12 ret, const vec384fp12 a);
/* caveat lector! |n| has to be non-zero and not more than 3! */
statik void frobenius_map_fp12(vec384fp12 ret, const vec384fp12 a, size_t n);

#else

extern "C" {
__device__ void mul_fp(vec384 ret, const vec384 a, const vec384 b);
__device__ void sqr_fp(vec384 ret, const vec384 a);
__device__ void add_fp(vec384 ret, const vec384 a, const vec384 b);
__device__ void sub_fp(vec384 ret, const vec384 a, const vec384 b);
__device__ void cneg_fp(vec384 ret, const vec384 ap, unsigned int flag);
__device__ void rshift_fp(vec384 ret, const vec384 a, unsigned int cnt);
__device__ void lshift_fp(vec384 ret, const vec384 a, unsigned int cnt);
__device__ void mul_by_3_fp(vec384 ret, const vec384 a);
__device__ void from_fp(vec384 ret, const vec384 a);

#pragma diag_suppress 3151
__device__ void mul_384(vec768 ret, const vec384 a, const vec384 b);
__device__ void sqr_384(vec768 ret, const vec384 a);
#pragma diag_default 3151
__device__ void redc_fp(vec384 ret, const vec768 a);
__device__ void add_fpx2(vec768 ret, const vec768 a, const vec768 b);
__device__ void sub_fpx2(vec768 ret, const vec768 a, const vec768 b);

__device__ void vec_load_global(limb_t *ret, const limb_t *a,
                                unsigned int sz = 48);
}

statik inline void mul_by_8_fp(vec384 ret, const vec384 a)
{
  lshift_fp(ret, a, 3);
}

statik inline void add_fp2(vec384x ret, const vec384x a, const vec384x b)
{
  add_fp(ret[0], a[0], b[0]);
  add_fp(ret[1], a[1], b[1]);
}

statik inline void sub_fp2(vec384x ret, const vec384x a, const vec384x b)
{
  sub_fp(ret[0], a[0], b[0]);
  sub_fp(ret[1], a[1], b[1]);
}

statik inline void mul_by_3_fp2(vec384x ret, const vec384x a)
{
  mul_by_3_fp(ret[0], a[0]);
  mul_by_3_fp(ret[1], a[1]);
}

statik inline void mul_by_8_fp2(vec384x ret, const vec384x a)
{
  lshift_fp(ret[0], a[0], 3);
  lshift_fp(ret[1], a[1], 3);
}

statik inline void lshift_fp2(vec384x ret, const vec384x a, size_t count)
{
  lshift_fp(ret[0], a[0], count);
  lshift_fp(ret[1], a[1], count);
}

statik inline void cneg_fp2(vec384x ret, const vec384x a, limb_t flag)
{
  cneg_fp(ret[0], a[0], flag);
  cneg_fp(ret[1], a[1], flag);
}

statik inline void mul_fp2(vec384x ret, const vec384x a, const vec384x b)
{
  vec384 aa, bb, cc;

  add_fp(aa, a[0], a[1]);
  add_fp(bb, b[0], b[1]);
  mul_fp(bb, bb, aa);

  mul_fp(aa, a[0], b[0]);
  mul_fp(cc, a[1], b[1]);

  sub_fp(ret[0], aa, cc);
  sub_fp(ret[1], bb, aa);
  sub_fp(ret[1], ret[1], cc);
}

statik inline void sqr_fp2(vec384x ret, const vec384x a)
{
  vec384 t0, t1;

  add_fp(t0, a[0], a[1]);
  sub_fp(t1, a[0], a[1]);

  mul_fp(ret[1], a[0], a[1]);
  add_fp(ret[1], ret[1], ret[1]);

  mul_fp(ret[0], t0, t1);
}
#endif

#define neg_fp(r, a)  cneg_fp((r), (a), 1)
#define neg_fp2(r, a) cneg_fp2((r), (a), 1)

#endif /* __BLS12_381_ASM_FIELDS_H__ */
