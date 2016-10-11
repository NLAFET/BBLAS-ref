/**
 *
 * @file bblas_smagma.h
 *
 * @brief BBLAS Magma routines for float precision.
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
 * Contains wrappers to batch Magma functions.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from ./bblas_magma/bblas_zmagma.h normal z -> s, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_SMAGMA_H
#define BBLAS_SMAGMA_H

/* Check whether we are compiling with MAGMA before including headers */
#if defined(BBLAS_WITH_MAGMA)
#include "magma.h"
#endif

#define REAL
#include "bblas.h"

void magma_sgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M, const int *N, const int *K,
		       const float *alpha,
		       const float **A, const int *lda,
		       const float **B, const int *ldb,
		       const float *beta,
		       float **C, const int *ldc, 
		       const int batch_count, enum BBLAS_OPTS batch_opts, int *info);


#ifdef COMPLEX
void magma_ssymm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const float *alpha,
		       const float **arrayA, const int *lda,
		       const float **arrayB, const int *ldb,
		       const float *beta, float **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_ssymm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const float *alpha,
		       const float **arrayA, const int *lda,
		       const float **arrayB, const int *ldb,
		       const float *beta, float **arrayC,
		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void magma_ssyr2k_batch(
			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
			const int *N, const int *K, const float *alpha,
			const float **arrayA, const int *lda,
			const float **arrayB, const int *ldb,
			const float *beta, float **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_ssyr2k_batch(
			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
			const int *N, const int *K, const float *alpha,
			const float **arrayA, const int *lda,
			const float **arrayB, const int *ldb,
			const float  *beta, float **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_ssyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
		       const int *N, const int *K, const float *alpha,
		       const float **arrayA, const int *lda,
		       const float *beta, float **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_ssyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const float *alpha,
		       const float **arrayA, const int *lda,
		       const float  *beta, float **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_strmm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const float *alpha,
		       const float **arrayA, const int *lda,
		       float **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void magma_strsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const float *alpha,
		       const float **arrayA, const int *lda,
		       float **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef REAL
#endif //BBLAS_SMAGMA_H
