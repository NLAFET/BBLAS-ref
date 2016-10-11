/**
 * @file bblas_somp.h
 *
 * @brief BBLAS header file for float routines.
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
 * @generated from ./bblas_omp/bblas_zomp.h normal z -> s, Mon Jun  6 09:44:14 2016
 **/
#endif

#ifndef OMP_BBLAS_S_H
#define OMP_BBLAS_S_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"
#define REAL

void omp_sgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const float *alpha,
    const float **arrayA, const int *lda,
    const float **arrayB,
    const int *ldb, const float *beta,
    float **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_ssymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const float *alpha,
    const float **arrayA, const int *lda,
    const float **arrayB, const int *ldb,
    const float *beta, float **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_ssymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const float *alpha,
    const float **arrayA, const int *lda,
    const float **arrayB, const int *ldb,
    const float *beta, float **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void omp_ssyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const float *alpha,
    const float **arrayA, const int *lda,
    const float **arrayB, const int *ldb,
    const float *beta, float **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_ssyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const float *alpha,
    const float **arrayA, const int *lda,
    const float **arrayB, const int *ldb,
    const float  *beta, float **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_ssyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const float *alpha,
    const float **arrayA, const int *lda,
    const float *beta, float **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_ssyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const float *alpha,
    const float **arrayA, const int *lda,
    const float  *beta, float **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_strsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const float *alpha,
    const float **arrayA, const int *lda,
    float **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void omp_strmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const float *alpha,
    const float **arrayA, const int *lda,
    float **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef REAL
#endif /* OMP_BBLAS_S_H */
