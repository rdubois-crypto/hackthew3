/*
 * blst_e1.h
 *
 *  Created on: Jul 4, 2022
 *      Author: dubois
 */

#ifndef SRC_INNOVATION_BLST_E1_H_
#define SRC_INNOVATION_BLST_E1_H_

#include "blst_errors.h"
#include "blst_point.h"

extern const POINTonE1 BLS12_381_G1;

 BLST_ERROR POINTonE1_Deserialize_BE(POINTonE1_affine *out,
                                           const unsigned char in[96]);

#endif /* SRC_INNOVATION_BLST_E1_H_ */
