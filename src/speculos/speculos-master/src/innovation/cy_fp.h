/*
 * cy_fp.h
 *
 *  Created on: Jul 5, 2022
 *      Author: dubois
 */

#ifndef API_CY_FP_H_
#define API_CY_FP_H_

#include "cy_errors.h"

cy_error_t fp_init(fp_ctx_t *ctx, uint8_t *Mem, int argc, uint8_t **argv);

cy_error_t fp_destroy(fp_ctx_t *ctx);


#endif /* API_CY_FP_H_ */
