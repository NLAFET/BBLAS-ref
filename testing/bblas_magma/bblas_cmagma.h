/**
 *
 * @file bblas_cmagma.h
 *
 * @brief BBLAS Magma routines for float _Complex precision.
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
 * @generated from ./bblas_magma/bblas_zmagma.h normal z -> c, Mon Jun  6 09:44:14 2016
 **/
#endif

#ifndef BBLAS_CMAGMA_H
#define BBLAS_CMAGMA_H

/* Check whether we are compiling with MAGMA before including headers */
#if defined(BBLAS_WITH_MAGMA)
#include "magma.h"
#endif

#define COMPLEX
#include "bblas.h"

void magma_cgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M, const int *N, const int *K,
		       const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **A, const int *lda,
		       const BBLAS_Complex32_t **B, const int *ldb,
		       const BBLAS_Complex32_t *beta,
		       BBLAS_Complex32_t **C, const int *ldc, 
		       const int batch_count, enum BBLAS_OPTS batch_opts, int *info);


#ifdef COMPLEX
void magma_chemm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       const BBLAS_Complex32_t **arrayB, const int *ldb,
		       const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_csymm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       const BBLAS_Complex32_t **arrayB, const int *ldb,
		       const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void magma_csyr2k_batch(
			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
			const int *N, const int *K, const BBLAS_Complex32_t *alpha,
			const BBLAS_Complex32_t **arrayA, const int *lda,
			const BBLAS_Complex32_t **arrayB, const int *ldb,
			const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_cher2k_batch(
			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
			const int *N, const int *K, const BBLAS_Complex32_t *alpha,
			const BBLAS_Complex32_t **arrayA, const int *lda,
			const BBLAS_Complex32_t **arrayB, const int *ldb,
			const float  *beta, BBLAS_Complex32_t **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_csyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
		       const int *N, const int *K, const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_cherk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const float *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       const float  *beta, BBLAS_Complex32_t **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_ctrmm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       BBLAS_Complex32_t **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void magma_ctrsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       BBLAS_Complex32_t **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef COMPLEX
#endif //BBLAS_CMAGMA_H
