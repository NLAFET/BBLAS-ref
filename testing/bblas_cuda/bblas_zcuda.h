/**
 *
 * @file bblas_zcuda.h
 *
 * @brief BBLAS CuBLAS routines for double _Complex precision.
 *
 *  BBLAS is a software package provided by Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-03-31
 *
 * Contains wrappers to batch CuBLAS functions.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#ifndef BBLAS_ZCUBLAS_H
#define BBLAS_ZCUBLAS_H

/* Check whether we are compiling with CUBLAS before including headers */
#if defined(BBLAS_WITH_CUBLAS)
#include "cublas_v2.h"
#include "cuda_runtime.h"
#endif

#define COMPLEX

#include "bblas.h"
#include "bblas_testing.h"

void cublas_zgemm_batch( const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
        const int *M, const int *N, const int *K,
        const BBLAS_Complex64_t *alpha,
        const BBLAS_Complex64_t **A, const int *lda,
        const BBLAS_Complex64_t **B, const int *ldb,
        const BBLAS_Complex64_t *beta,
        BBLAS_Complex64_t **C, const int *ldc, 
        const int batch_count, enum BBLAS_OPTS batch_opts, int *info,
        bblas_ztest_t *test );

void cublas_ztrsm_batch( const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
	const int  *M, const int *N,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **A, const int *lda,
    BBLAS_Complex64_t **B, const int *ldb,
    const int batch_count, enum BBLAS_OPTS batch_opts, int *info, 
    bblas_ztest_t *test);


#undef COMPLEX
#endif
