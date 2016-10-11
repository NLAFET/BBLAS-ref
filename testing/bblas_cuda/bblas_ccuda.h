/**
 *
 * @file bblas_ccuda.h
 *
 * @brief BBLAS CuBLAS routines for float _Complex precision.
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
 * @generated from ./bblas_cuda/bblas_zcuda.h normal z -> c, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_CCUBLAS_H
#define BBLAS_CCUBLAS_H

/* Check whether we are compiling with CUBLAS before including headers */
#if defined(BBLAS_WITH_CUBLAS)
#include "cublas_v2.h"
#include "cuda_runtime.h"
#endif

#define COMPLEX

#include "bblas.h"
#include "bblas_testing.h"

void cublas_cgemm_batch( const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
        const int *M, const int *N, const int *K,
        const BBLAS_Complex32_t *alpha,
        const BBLAS_Complex32_t **A, const int *lda,
        const BBLAS_Complex32_t **B, const int *ldb,
        const BBLAS_Complex32_t *beta,
        BBLAS_Complex32_t **C, const int *ldc, 
        const int batch_count, enum BBLAS_OPTS batch_opts, int *info,
        bblas_ctest_t *test );

void cublas_ctrsm_batch( const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
	const int  *M, const int *N,
    const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **A, const int *lda,
    BBLAS_Complex32_t **B, const int *ldb,
    const int batch_count, enum BBLAS_OPTS batch_opts, int *info, 
    bblas_ctest_t *test);


#undef COMPLEX
#endif
