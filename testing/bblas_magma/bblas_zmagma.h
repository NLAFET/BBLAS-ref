/**
 *
 * @file bblas_zmagma.h
 *
 * @brief BBLAS Magma routines for double _Complex precision.
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
 * @precisions normal z -> c d s
 **/
#endif

#ifndef BBLAS_ZMAGMA_H
#define BBLAS_ZMAGMA_H

/* Check whether we are compiling with MAGMA before including headers */
#if defined(BBLAS_WITH_MAGMA)
#include "magma.h"
#endif

#define COMPLEX
#include "bblas.h"

void magma_zgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M, const int *N, const int *K,
		       const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **A, const int *lda,
		       const BBLAS_Complex64_t **B, const int *ldb,
		       const BBLAS_Complex64_t *beta,
		       BBLAS_Complex64_t **C, const int *ldc, 
		       const int batch_count, enum BBLAS_OPTS batch_opts, int *info);


#ifdef COMPLEX
void magma_zhemm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       const BBLAS_Complex64_t **arrayB, const int *ldb,
		       const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_zsymm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       const BBLAS_Complex64_t **arrayB, const int *ldb,
		       const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void magma_zsyr2k_batch(
			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
			const int *N, const int *K, const BBLAS_Complex64_t *alpha,
			const BBLAS_Complex64_t **arrayA, const int *lda,
			const BBLAS_Complex64_t **arrayB, const int *ldb,
			const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_zher2k_batch(
			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
			const int *N, const int *K, const BBLAS_Complex64_t *alpha,
			const BBLAS_Complex64_t **arrayA, const int *lda,
			const BBLAS_Complex64_t **arrayB, const int *ldb,
			const double  *beta, BBLAS_Complex64_t **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_zsyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
		       const int *N, const int *K, const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_zherk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const double *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       const double  *beta, BBLAS_Complex64_t **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_ztrmm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       BBLAS_Complex64_t **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void magma_ztrsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       BBLAS_Complex64_t **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef COMPLEX
#endif //BBLAS_ZMAGMA_H
