/*
 * blst_e2.h
 *
 *  Created on: Jul 4, 2022
 *      Author: dubois
 */

#ifndef SRC_INNOVATION_BLST_E2_H_
#define SRC_INNOVATION_BLST_E2_H_

extern const POINTonE2 BLS12_381_G2;

statik BLST_ERROR POINTonE2_Deserialize_BE(POINTonE2_affine *out,
                                           const unsigned char in[192]);

#endif /* SRC_INNOVATION_BLST_E2_H_ */
