/**
 * @file bblas_zomp.h
 *
 * @brief BBLAS header file for double _Complex routines.
 *
 *  BBLAS is a software package provided by Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 * Reference implementations of BBLAS using OpenMP when looping over the subproblems.
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#ifndef OMP_BBLAS_Z_H
#define OMP_BBLAS_Z_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"
#define COMPLEX

void omp_zgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB,
    const int *ldb, const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_zhemm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_zsymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void omp_zsyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_zher2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const double  *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_zsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_zherk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const double  *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_ztrsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    BBLAS_Complex64_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void omp_ztrmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    BBLAS_Complex64_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef COMPLEX
#endif /* OMP_BBLAS_Z_H */
