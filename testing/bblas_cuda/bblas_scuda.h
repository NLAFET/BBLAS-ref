/**
 *
 * @file bblas_scuda.h
 *
 * @brief BBLAS CuBLAS routines for float precision.
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
 * @generated from ./bblas_cuda/bblas_zcuda.h normal z -> s, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_SCUBLAS_H
#define BBLAS_SCUBLAS_H

/* Check whether we are compiling with CUBLAS before including headers */
#if defined(BBLAS_WITH_CUBLAS)
#include "cublas_v2.h"
#include "cuda_runtime.h"
#endif

#define REAL

#include "bblas.h"
#include "bblas_testing.h"

void cublas_sgemm_batch( const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
        const int *M, const int *N, const int *K,
        const float *alpha,
        const float **A, const int *lda,
        const float **B, const int *ldb,
        const float *beta,
        float **C, const int *ldc, 
        const int batch_count, enum BBLAS_OPTS batch_opts, int *info,
        bblas_stest_t *test );

void cublas_strsm_batch( const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
	const int  *M, const int *N,
    const float *alpha,
    const float **A, const int *lda,
    float **B, const int *ldb,
    const int batch_count, enum BBLAS_OPTS batch_opts, int *info, 
    bblas_stest_t *test);


#undef REAL
#endif
