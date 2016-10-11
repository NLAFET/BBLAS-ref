/**
 * @file bblas_domp.h
 *
 * @brief BBLAS header file for double routines.
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
 * @generated from ./bblas_omp/bblas_zomp.h normal z -> d, Mon Jun  6 09:44:14 2016
 **/
#endif

#ifndef OMP_BBLAS_D_H
#define OMP_BBLAS_D_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"
#define REAL

void omp_dgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB,
    const int *ldb, const double *beta,
    double **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_dsymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_dsymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double *beta, double **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void omp_dsyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_dsyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double  *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_dsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void omp_dsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double  *beta, double **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void omp_dtrsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    double **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void omp_dtrmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    double **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

#undef REAL
#endif /* OMP_BBLAS_D_H */
