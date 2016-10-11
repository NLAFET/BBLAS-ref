/**
 * @file bblas_comp.h
 *
 * @brief BBLAS header file for float _Complex routines.
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
 * @generated from bblas_omp/bblas_zomp.h normal z -> c, Mon Jun  6 09:44:14 2016
 **/
#endif

#ifndef OMP_BBLAS_C_H
#define OMP_BBLAS_C_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"
#define COMPLEX

void omp_cgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB,
    const int *ldb, const BBLAS_Complex32_t *beta,
    BBLAS_Complex32_t **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_chemm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_csymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void omp_csyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_cher2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const float  *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_csyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_cherk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const float *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const float  *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_ctrsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    BBLAS_Complex32_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void omp_ctrmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    BBLAS_Complex32_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef COMPLEX
#endif /* OMP_BBLAS_C_H */
