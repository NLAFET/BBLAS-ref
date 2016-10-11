/**
 * @file bblas_dmagma.c
 *
 * @brief BBLAS MAGMA routines for double precision.
 *
 *  BBLAS is a software package provided by Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-03-17
 *
 * Contains wrappers to batch MAGMA functions.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from ./bblas_magma/bblas_zmagma.c normal z -> d, Mon Jun  6 09:44:13 2016
 **/
#endif
#include "bblas_magma.h"

/**
 * Wrapper for dgemm_batch in MAGMA.
 *
 * Calls cblas_dgemm_batch (linked to MAGMA) with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS.
 * The arguments are converted to an MAGMA function call within this function.
 **/

#define REAL
void magma_dgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M,  const int *N, const int *K,
		       const double *alpha,
		       const double **arrayA, const int *lda,
		       const double **arrayB,
		       const int *ldb, const double *beta,
		       double **arrayC, const int *ldc, const int batch_count,
		       enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_dgemm_batched(
			  transA[first_index], transB[first_index], M[first_index],
			  N[first_index], K[first_index], 
			  ((double *) alpha)[first_index],
			  (double const * const * )arrayA, lda[first_index],
			  (double const * const * )arrayB, ldb[first_index],
			  ((double *)beta)[first_index],
			  (double **)arrayC, ldc[first_index], batch_count, queue);
  magma_sync_wtime(queue );
#endif
}

void magma_dtrsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const double *alpha,
		       const double **arrayA, const int *lda,
		       double **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_dtrsm_batched(side[first_index], uplo[first_index], transA[first_index],
			  diag[first_index], M[first_index], N[first_index], 
			  ((double *) alpha)[first_index],
			  (double **)arrayA, lda[first_index],
			  (double **)arrayB, ldb[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}

#ifdef COMPLEX
void magma_dsyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const double *alpha,
		       const double **arrayA, const int *lda,
		       const double  *beta, double **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_dsyrk_batched(uplo[first_index], trans[first_index],
			  N[first_index], K[first_index], 
			  alpha[first_index],
			  (double const * const * )arrayA, lda[first_index],
			  beta[first_index],
			  (double **)arrayC, ldc[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}
#endif //#ifdef COMPLEX


/* #ifdef COMPLEX */
/* void magma_dsymm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const double *alpha, */
/* 		       const double **arrayA, const int *lda, */
/* 		       const double **arrayB, const int *ldb, */
/* 		       const double *beta, double **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_dsymm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((double *) alpha)[first_index], */
/* 			  (double const * const * )arrayA, lda[first_index], */
/* 			  (double const * const * )arrayB, ldb[first_index], */
/* 			  ((double *)beta)[first_index], */
/* 			  (double **)arrayC, ldc[first_index], batch_count, queue); */
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */
/* #endif // ifdef COMPLEX */


/* void magma_dsymm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const double *alpha, */
/* 		       const double **arrayA, const int *lda, */
/* 		       const double **arrayB, const int *ldb, */
/* 		       const double *beta, double **arrayC, */
/* 		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */

/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_dsymm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((double *) alpha)[first_index], */
/* 			  (double const * const * )arrayA, lda[first_index], */
/* 			  (double const * const * )arrayB, ldb[first_index], */
/* 			  ((double *)beta)[first_index], */
/* 			  (double **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */



/* void magma_dsyr2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const double *alpha, */
/* 			const double **arrayA, const int *lda, */
/* 			const double **arrayB, const int *ldb, */
/* 			const double *beta, double **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_dsyr2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((double *) alpha)[first_index], */
/* 			  (double const * const * )arrayA, lda[first_index], */
/* 			  (double const * const * )arrayB, ldb[first_index], */
/* 			  ((double *)beta)[first_index], */
/* 			  (double **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */


/* #ifdef COMPLEX */
/* void magma_dsyr2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const double *alpha, */
/* 			const double **arrayA, const int *lda, */
/* 			const double **arrayB, const int *ldb, */
/* 			const double  *beta, double **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_dsyr2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((double *) alpha)[first_index], */
/* 			  (double const * const * )arrayA, lda[first_index], */
/* 			  (double const * const * )arrayB, ldb[first_index], */
/* 			  ((double *)beta)[first_index], */
/* 			  (double **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */
/* #endif //#ifdef COMPLEX */


/* void magma_dsyrk_batch( */
/* 		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 		       const int *N, const int *K, const double *alpha, */
/* 		       const double **arrayA, const int *lda, */
/* 		       const double *beta, double **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_dsyrk_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((double *) alpha)[first_index], */
/* 			  (double const * const * )arrayA, lda[first_index], */
/* 			  ((double *)beta)[first_index], */
/* 			  (double **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */

/* void magma_dtrmm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo, */
/* 		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag, */
/* 		       const int *M, const int *N, const double *alpha, */
/* 		       const double **arrayA, const int *lda, */
/* 		       double **arrayB, const int *ldb, */
/* 		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_dtrmm_batched(side[first_index], uplo[first_index], transA[first_index], */
/* 			  diag[first_index], M[first_index], N[first_index],  */
/* 			  ((double *) alpha)[first_index], */
/* 			  (double const * const * )arrayA, lda[first_index], */
/* 			  (double **)arrayB, ldb[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */




#undef REAL
