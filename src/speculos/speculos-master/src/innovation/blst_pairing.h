/*
 * blst_pairing.h
 *
 *  Created on: Jul 4, 2022
 *      Author: dubois
 */

#ifndef SRC_INNOVATION_BLST_PAIRING_H_
#define SRC_INNOVATION_BLST_PAIRING_H_

void miller_loop_n(vec384fp12 ret, const POINTonE2_affine Q[],
                   const POINTonE1_affine P[], size_t n);

void final_exp(vec384fp12 ret, const vec384fp12 f);

#endif /* SRC_INNOVATION_BLST_PAIRING_H_ */
