
#ifndef _CY_PBC_H
#define _CY_PBC_H

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

#define CY_UNUSED(x) x

/* Evaluation of X over elliptic curve of weierstrass form X^3+aX=b*/
cx_err_t cy_curve_eval_ysquare(cx_bn_t o_bn_y2, cx_bn_t i_bn_a, cx_bn_t i_bn_b,
                               cx_bn_t i_bn_x, cx_bn_t i_bn_p, size_t i_tu8_p);

/* hashing of a Point over a Weierstrass curve from a digest using
 * SWU/Brier-Coron and all technique (alternative to Pedersen Hash)
 * https://eprint.iacr.org/2009/340.pdf*/
cx_err_t cy_swu_hashpoint(cx_curve_t curve, cx_bn_t o_bn_x, cx_bn_t o_bn_y,
                          uint8_t *i_digest);

cx_err_t cy_pairing(cx_curve_t curve, cx_ecpoint_t *P1, cx_ecpoint_t *P2,
                    uint8_t *pairing, size_t size8_pairing);

cx_err_t cy_pairing_asn1(cx_curve_t curve, unsigned char *P1, size_t P1_len,
                         unsigned char *P2, size_t P2_len,
                         const unsigned char *k, unsigned int k_len);

#endif
