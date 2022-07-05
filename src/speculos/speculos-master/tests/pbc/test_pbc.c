/*
 * Copyright Ledger
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <malloc.h>
#include <setjmp.h>
#include <stdbool.h>
#include <string.h>

#include <openssl/bn.h>
#include <openssl/ec.h>
#include <stdbool.h>
#include <stdio.h>

#include "bolos/cx_ec.h"
#include "bolos/cxlib.h"
#include "innovation/cy_pbc.h"

#define CX_TRUE  true
#define CX_FALSE false

#define _NB_TESTS 2
#define SIZEFIELD 384

#if (SIZEFIELD == 256)
#define CX_CURVE_ID CX_CURVE_256K1
#endif
#if (SIZEFIELD == 381)
#define CX_CURVE_ID CX_CURVE_BLS12_381_G1
#endif
#if (SIZEFIELD == 384)
#define CX_CURVE_ID CX_CURVE_SECP384R1
#endif

/* A loop of _NB_TESTS random generation of x and y extraction to ensure
 * correctness of cx_curve_ysquare*/
static cx_err_t test_curve_ysquare(cx_curve_t curve)
{
  /* I.Declarations*/
  cx_err_t error = CX_OK;
  cx_err_t is_square = false;

  cx_bn_t bn_a;
  cx_bn_t bn_b;
  cx_bn_t bn_x;
  cx_bn_t bn_p;
  cx_bn_t bn_y;
  cx_bn_t bn_y2;
  cx_ecpoint_t P;
  uint32_t i = 0;

  bool is_on_curve = true;

  const cx_curve_weierstrass_t *l_curve = cx_ecfp_get_domain(curve);
  if (l_curve == NULL)
    return CX_INVALID_PARAMETER_VALUE;

  size_t tu8_p = l_curve->length; /* prime field byte length*/
  /* II. Allocations*/;
  CX_CHECK(sys_cx_bn_alloc(&bn_y2, tu8_p));
  CX_CHECK(sys_cx_bn_alloc(&bn_y, tu8_p));
  CX_CHECK(sys_cx_bn_alloc(&bn_x, tu8_p));
  CX_CHECK(sys_cx_ecpoint_alloc(&P, curve));
  CX_CHECK(sys_cx_bn_alloc_init(
      &bn_p, tu8_p, l_curve->p,
      tu8_p)); /* allocating and initializing big number p*/

  CX_CHECK(sys_cx_bn_is_prime(
      bn_p, &is_on_curve)); /* primality test for prime field modulus*/

  if (is_on_curve == false) {
    printf("\n primality failure !");
    return CX_INVALID_PARAMETER_VALUE;
  }

  CX_CHECK(sys_cx_bn_alloc_init(
      &bn_a, tu8_p, l_curve->a,
      tu8_p)); /* allocating and initializing big number a*/
  CX_CHECK(sys_cx_bn_alloc_init(
      &bn_b, tu8_p, l_curve->b,
      tu8_p)); /* allocating and initializing big number b*/

  /* III. Computations*/
  for (i = 0; i < _NB_TESTS; i++) {
    CX_CHECK(sys_cx_bn_rand(bn_x));

    CX_CHECK(sys_cx_bn_reduce(bn_x, bn_x, bn_p));
    /* evaluate*/
    CX_CHECK(cy_curve_eval_ysquare(bn_y2, bn_a, bn_b, bn_x, bn_p, tu8_p));

    /* check curve compliance*/
    is_square =
        (sys_cx_bn_mod_sqrt(bn_y, bn_y2, bn_p, 1)); /* extract square root*/

    /* test only if y^2 is really a square*/
    if (is_square == CX_OK) {
      /*printf("\n square");*/
      CX_CHECK(sys_cx_ecpoint_init_bn(&P, bn_x, bn_y)); /* P=(x,y) */
      CX_CHECK(sys_cx_ecpoint_is_on_curve(&P, &is_on_curve));

      if (is_on_curve == false) {
        error = CX_EC_INVALID_POINT;
        printf("\n Check failure");
      }
      // else printf("\n Check OK");
    }
    /*else printf("\n not square, CX=%X",is_square);*/
  }
  /* IV. Free*/
  CX_CHECK(sys_cx_bn_destroy(&bn_y2));
  CX_CHECK(sys_cx_bn_destroy(&bn_y));
  CX_CHECK(sys_cx_bn_destroy(&bn_x));
  CX_CHECK(sys_cx_bn_destroy(&bn_p));
  CX_CHECK(sys_cx_bn_destroy(&bn_a));
  CX_CHECK(sys_cx_bn_destroy(&bn_b));
  CX_CHECK(sys_cx_ecpoint_destroy(&P));
/* V. Return errcode*/
end:
  return error;
}

int main()
{
  /* I. Declarations*/
  cx_err_t error = CX_INTERNAL_ERROR;
  cx_bn_t bn_x;
  cx_bn_t bn_y;
  const cx_curve_weierstrass_t *l_curve = cx_ecfp_get_domain(CX_CURVE_ID);
  if (l_curve == NULL) {
    printf("\n Curve not known !!");
    return CX_INVALID_PARAMETER_VALUE;
  }
  size_t tu8_p = l_curve->length; /* prime field byte length*/
  uint32_t i, j;                  /* loop variables*/

  /* II. Allocations*/
  CX_CHECK(sys_cx_bn_lock(48, 0)); /* Init BN*/
  CX_CHECK(sys_cx_bn_alloc(&bn_x, tu8_p));
  CX_CHECK(sys_cx_bn_alloc(&bn_y, tu8_p));
  printf("\n curve id=%d", CX_CURVE_ID);

  uint8_t u8_hPoint[1 + SIZEFIELD * 2]; // note : check size of allocated point
                                        // (0x0*+ 2 size ?)
  uint8_t digest[SIZEFIELD];
  printf("\n test_curve_ysquare: ");

  CX_CHECK(test_curve_ysquare(CX_CURVE_ID));
  printf(" OK");

  printf("\n test_cx_swu_point ");
  for (i = 0; i < _NB_TESTS; i++) {
    printf("\n ----------------------------\n Digest:\n h=");
    for (j = 0; j<SIZEFIELD>> 3; j++)
      printf(" %02x", digest[j]);
    sys_cx_get_random_bytes(digest, SIZEFIELD >> 3);
    CX_CHECK(cy_swu_hashpoint(CX_CURVE_ID, bn_x, bn_y, digest));
    printf("\n Point:\n x=");
    sys_cx_bn_export(bn_x, digest, SIZEFIELD >> 3);
    for (j = 0; j<SIZEFIELD>> 3; j++)
      printf(" %02x", digest[j]);
    printf("\n y=");
    sys_cx_bn_export(bn_y, digest, SIZEFIELD >> 3);
    for (j = 0; j<SIZEFIELD>> 3; j++)
      printf(" %02x", digest[j]);
  }
  if (SIZEFIELD == 384) {
    printf("\n--------------------\n A cool digest Sha384[Ledger: Secured By "
           "Design]:\n h=");
    uint8_t digestcool[SIZEFIELD >> 3] = {
      0xc4, 0x2f, 0x8a, 0x53, 0x4f, 0x8b, 0x0f, 0xf0, 0x2e, 0x7c, 0xba, 0x07,
      0x11, 0xe2, 0xdd, 0x13, 0x0e, 0xb8, 0x29, 0xef, 0x32, 0x9e, 0x8f, 0x68,
      0xf5, 0x8f, 0xfa, 0x75, 0xdb, 0x4c, 0xb8, 0xa5, 0x39, 0xac, 0x04, 0x98,
      0xee, 0x34, 0xac, 0x5e, 0x60, 0xeb, 0xef, 0xdd, 0x60, 0x28, 0x00, 0x1a
    };

    for (j = 0; j<SIZEFIELD >> 3; j++)
      printf(" %02x", digestcool[j]);
    CX_CHECK(cy_swu_hashpoint(CX_CURVE_ID, bn_x, bn_y, digestcool));
    printf("\n Point:\n x=");
    sys_cx_bn_export(bn_x, digestcool, SIZEFIELD >> 3);
    for (j = 0; j<SIZEFIELD >> 3; j++)
      printf(" %02x", digestcool[j]);
    printf("\n y=");
    sys_cx_bn_export(bn_y, digestcool, SIZEFIELD >> 3);
    for (j = 0; j<SIZEFIELD >> 3; j++)
      printf(" %02x", digestcool[j]);
  }

end:
  printf("\n Global return :");
  if (error == 0)
    printf(" OK\n ");
  else
    printf("\n KO !!!");

  return error;
}
