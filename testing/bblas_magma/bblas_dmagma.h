/**
 *
 * @file bblas_dmagma.h
 *
 * @brief BBLAS Magma routines for double precision.
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
 * @generated from ./bblas_magma/bblas_zmagma.h normal z -> d, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_DMAGMA_H
#define BBLAS_DMAGMA_H

/* Check whether we are compiling with MAGMA before including headers */
#if defined(BBLAS_WITH_MAGMA)
#include "magma.h"
#endif

#define REAL
#include "bblas.h"

void magma_dgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M, const int *N, const int *K,
		       const double *alpha,
		       const double **A, const int *lda,
		       const double **B, const int *ldb,
		       const double *beta,
		       double **C, const int *ldc, 
		       const int batch_count, enum BBLAS_OPTS batch_opts, int *info);


#ifdef COMPLEX
void magma_dsymm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const double *alpha,
		       const double **arrayA, const int *lda,
		       const double **arrayB, const int *ldb,
		       const double *beta, double **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_dsymm_batch(
		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
		       const int *M, const int *N, const double *alpha,
		       const double **arrayA, const int *lda,
		       const double **arrayB, const int *ldb,
		       const double *beta, double **arrayC,
		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void magma_dsyr2k_batch(
			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
			const int *N, const int *K, const double *alpha,
			const double **arrayA, const int *lda,
			const double **arrayB, const int *ldb,
			const double *beta, double **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_dsyr2k_batch(
			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
			const int *N, const int *K, const double *alpha,
			const double **arrayA, const int *lda,
			const double **arrayB, const int *ldb,
			const double  *beta, double **arrayC,
			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_dsyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
		       const int *N, const int *K, const double *alpha,
		       const double **arrayA, const int *lda,
		       const double *beta, double **arrayC,
		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void magma_dsyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const double *alpha,
		       const double **arrayA, const int *lda,
		       const double  *beta, double **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void magma_dtrmm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const double *alpha,
		       const double **arrayA, const int *lda,
		       double **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void magma_dtrsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const double *alpha,
		       const double **arrayA, const int *lda,
		       double **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef REAL
#endif //BBLAS_DMAGMA_H
