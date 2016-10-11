/**
 * @file bblas_zmagma.c
 *
 * @brief BBLAS MAGMA routines for double _Complex precision.
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
 * @precisions normal z -> c d s
 **/
#endif
#include "bblas_magma.h"

/**
 * Wrapper for zgemm_batch in MAGMA.
 *
 * Calls cblas_zgemm_batch (linked to MAGMA) with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS.
 * The arguments are converted to an MAGMA function call within this function.
 **/

#define COMPLEX
void magma_zgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M,  const int *N, const int *K,
		       const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       const BBLAS_Complex64_t **arrayB,
		       const int *ldb, const BBLAS_Complex64_t *beta,
		       BBLAS_Complex64_t **arrayC, const int *ldc, const int batch_count,
		       enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_zgemm_batched(
			  transA[first_index], transB[first_index], M[first_index],
			  N[first_index], K[first_index], 
			  ((magmaDoubleComplex *) alpha)[first_index],
			  (magmaDoubleComplex const * const * )arrayA, lda[first_index],
			  (magmaDoubleComplex const * const * )arrayB, ldb[first_index],
			  ((magmaDoubleComplex *)beta)[first_index],
			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue);
  magma_sync_wtime(queue );
#endif
}

void magma_ztrsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const BBLAS_Complex64_t *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       BBLAS_Complex64_t **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_ztrsm_batched(side[first_index], uplo[first_index], transA[first_index],
			  diag[first_index], M[first_index], N[first_index], 
			  ((magmaDoubleComplex *) alpha)[first_index],
			  (magmaDoubleComplex **)arrayA, lda[first_index],
			  (magmaDoubleComplex **)arrayB, ldb[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}

#ifdef COMPLEX
void magma_zherk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const double *alpha,
		       const BBLAS_Complex64_t **arrayA, const int *lda,
		       const double  *beta, BBLAS_Complex64_t **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_zherk_batched(uplo[first_index], trans[first_index],
			  N[first_index], K[first_index], 
			  alpha[first_index],
			  (magmaDoubleComplex const * const * )arrayA, lda[first_index],
			  beta[first_index],
			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}
#endif //#ifdef COMPLEX


/* #ifdef COMPLEX */
/* void magma_zhemm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const BBLAS_Complex64_t *alpha, */
/* 		       const BBLAS_Complex64_t **arrayA, const int *lda, */
/* 		       const BBLAS_Complex64_t **arrayB, const int *ldb, */
/* 		       const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_zhemm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((magmaDoubleComplex *) alpha)[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaDoubleComplex *)beta)[first_index], */
/* 			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue); */
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */
/* #endif // ifdef COMPLEX */


/* void magma_zsymm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const BBLAS_Complex64_t *alpha, */
/* 		       const BBLAS_Complex64_t **arrayA, const int *lda, */
/* 		       const BBLAS_Complex64_t **arrayB, const int *ldb, */
/* 		       const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC, */
/* 		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */

/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_zsymm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((magmaDoubleComplex *) alpha)[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaDoubleComplex *)beta)[first_index], */
/* 			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */



/* void magma_zsyr2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const BBLAS_Complex64_t *alpha, */
/* 			const BBLAS_Complex64_t **arrayA, const int *lda, */
/* 			const BBLAS_Complex64_t **arrayB, const int *ldb, */
/* 			const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_zsyr2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((magmaDoubleComplex *) alpha)[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaDoubleComplex *)beta)[first_index], */
/* 			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */


/* #ifdef COMPLEX */
/* void magma_zher2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const BBLAS_Complex64_t *alpha, */
/* 			const BBLAS_Complex64_t **arrayA, const int *lda, */
/* 			const BBLAS_Complex64_t **arrayB, const int *ldb, */
/* 			const double  *beta, BBLAS_Complex64_t **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_zher2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((magmaDoubleComplex *) alpha)[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaDoubleComplex *)beta)[first_index], */
/* 			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */
/* #endif //#ifdef COMPLEX */


/* void magma_zsyrk_batch( */
/* 		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 		       const int *N, const int *K, const BBLAS_Complex64_t *alpha, */
/* 		       const BBLAS_Complex64_t **arrayA, const int *lda, */
/* 		       const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_zsyrk_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((magmaDoubleComplex *) alpha)[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayA, lda[first_index], */
/* 			  ((magmaDoubleComplex *)beta)[first_index], */
/* 			  (magmaDoubleComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */

/* void magma_ztrmm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo, */
/* 		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag, */
/* 		       const int *M, const int *N, const BBLAS_Complex64_t *alpha, */
/* 		       const BBLAS_Complex64_t **arrayA, const int *lda, */
/* 		       BBLAS_Complex64_t **arrayB, const int *ldb, */
/* 		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ztrmm_batched(side[first_index], uplo[first_index], transA[first_index], */
/* 			  diag[first_index], M[first_index], N[first_index],  */
/* 			  ((magmaDoubleComplex *) alpha)[first_index], */
/* 			  (magmaDoubleComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaDoubleComplex **)arrayB, ldb[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */




#undef COMPLEX
