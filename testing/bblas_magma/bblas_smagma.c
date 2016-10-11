/**
 * @file bblas_smagma.c
 *
 * @brief BBLAS MAGMA routines for float precision.
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
 * @generated from ./bblas_magma/bblas_zmagma.c normal z -> s, Mon Jun  6 09:44:13 2016
 **/
#endif
#include "bblas_magma.h"

/**
 * Wrapper for sgemm_batch in MAGMA.
 *
 * Calls cblas_sgemm_batch (linked to MAGMA) with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS.
 * The arguments are converted to an MAGMA function call within this function.
 **/

#define REAL
void magma_sgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M,  const int *N, const int *K,
		       const float *alpha,
		       const float **arrayA, const int *lda,
		       const float **arrayB,
		       const int *ldb, const float *beta,
		       float **arrayC, const int *ldc, const int batch_count,
		       enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_sgemm_batched(
			  transA[first_index], transB[first_index], M[first_index],
			  N[first_index], K[first_index], 
			  ((float *) alpha)[first_index],
			  (float const * const * )arrayA, lda[first_index],
			  (float const * const * )arrayB, ldb[first_index],
			  ((float *)beta)[first_index],
			  (float **)arrayC, ldc[first_index], batch_count, queue);
  magma_sync_wtime(queue );
#endif
}

void magma_strsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const float *alpha,
		       const float **arrayA, const int *lda,
		       float **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_strsm_batched(side[first_index], uplo[first_index], transA[first_index],
			  diag[first_index], M[first_index], N[first_index], 
			  ((float *) alpha)[first_index],
			  (float **)arrayA, lda[first_index],
			  (float **)arrayB, ldb[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}

#ifdef COMPLEX
void magma_ssyrk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const float *alpha,
		       const float **arrayA, const int *lda,
		       const float  *beta, float **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_ssyrk_batched(uplo[first_index], trans[first_index],
			  N[first_index], K[first_index], 
			  alpha[first_index],
			  (float const * const * )arrayA, lda[first_index],
			  beta[first_index],
			  (float **)arrayC, ldc[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}
#endif //#ifdef COMPLEX


/* #ifdef COMPLEX */
/* void magma_ssymm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const float *alpha, */
/* 		       const float **arrayA, const int *lda, */
/* 		       const float **arrayB, const int *ldb, */
/* 		       const float *beta, float **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ssymm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((float *) alpha)[first_index], */
/* 			  (float const * const * )arrayA, lda[first_index], */
/* 			  (float const * const * )arrayB, ldb[first_index], */
/* 			  ((float *)beta)[first_index], */
/* 			  (float **)arrayC, ldc[first_index], batch_count, queue); */
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */
/* #endif // ifdef COMPLEX */


/* void magma_ssymm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const float *alpha, */
/* 		       const float **arrayA, const int *lda, */
/* 		       const float **arrayB, const int *ldb, */
/* 		       const float *beta, float **arrayC, */
/* 		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */

/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ssymm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((float *) alpha)[first_index], */
/* 			  (float const * const * )arrayA, lda[first_index], */
/* 			  (float const * const * )arrayB, ldb[first_index], */
/* 			  ((float *)beta)[first_index], */
/* 			  (float **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */



/* void magma_ssyr2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const float *alpha, */
/* 			const float **arrayA, const int *lda, */
/* 			const float **arrayB, const int *ldb, */
/* 			const float *beta, float **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ssyr2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((float *) alpha)[first_index], */
/* 			  (float const * const * )arrayA, lda[first_index], */
/* 			  (float const * const * )arrayB, ldb[first_index], */
/* 			  ((float *)beta)[first_index], */
/* 			  (float **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */


/* #ifdef COMPLEX */
/* void magma_ssyr2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const float *alpha, */
/* 			const float **arrayA, const int *lda, */
/* 			const float **arrayB, const int *ldb, */
/* 			const float  *beta, float **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ssyr2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((float *) alpha)[first_index], */
/* 			  (float const * const * )arrayA, lda[first_index], */
/* 			  (float const * const * )arrayB, ldb[first_index], */
/* 			  ((float *)beta)[first_index], */
/* 			  (float **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */
/* #endif //#ifdef COMPLEX */


/* void magma_ssyrk_batch( */
/* 		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 		       const int *N, const int *K, const float *alpha, */
/* 		       const float **arrayA, const int *lda, */
/* 		       const float *beta, float **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ssyrk_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((float *) alpha)[first_index], */
/* 			  (float const * const * )arrayA, lda[first_index], */
/* 			  ((float *)beta)[first_index], */
/* 			  (float **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */

/* void magma_strmm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo, */
/* 		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag, */
/* 		       const int *M, const int *N, const float *alpha, */
/* 		       const float **arrayA, const int *lda, */
/* 		       float **arrayB, const int *ldb, */
/* 		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_strmm_batched(side[first_index], uplo[first_index], transA[first_index], */
/* 			  diag[first_index], M[first_index], N[first_index],  */
/* 			  ((float *) alpha)[first_index], */
/* 			  (float const * const * )arrayA, lda[first_index], */
/* 			  (float **)arrayB, ldb[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */




#undef REAL
