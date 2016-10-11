/**
 * @file bblas_cmagma.c
 *
 * @brief BBLAS MAGMA routines for float _Complex precision.
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
 * @generated from ./bblas_magma/bblas_zmagma.c normal z -> c, Mon Jun  6 09:44:13 2016
 **/
#endif
#include "bblas_magma.h"

/**
 * Wrapper for cgemm_batch in MAGMA.
 *
 * Calls cblas_cgemm_batch (linked to MAGMA) with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS.
 * The arguments are converted to an MAGMA function call within this function.
 **/

#define COMPLEX
void magma_cgemm_batch(
		       const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
		       const int *M,  const int *N, const int *K,
		       const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       const BBLAS_Complex32_t **arrayB,
		       const int *ldb, const BBLAS_Complex32_t *beta,
		       BBLAS_Complex32_t **arrayC, const int *ldc, const int batch_count,
		       enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_cgemm_batched(
			  transA[first_index], transB[first_index], M[first_index],
			  N[first_index], K[first_index], 
			  ((magmaFloatComplex *) alpha)[first_index],
			  (magmaFloatComplex const * const * )arrayA, lda[first_index],
			  (magmaFloatComplex const * const * )arrayB, ldb[first_index],
			  ((magmaFloatComplex *)beta)[first_index],
			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue);
  magma_sync_wtime(queue );
#endif
}

void magma_ctrsm_batch(
		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
		       const int *M, const int *N, const BBLAS_Complex32_t *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       BBLAS_Complex32_t **arrayB, const int *ldb,
		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_ctrsm_batched(side[first_index], uplo[first_index], transA[first_index],
			  diag[first_index], M[first_index], N[first_index], 
			  ((magmaFloatComplex *) alpha)[first_index],
			  (magmaFloatComplex **)arrayA, lda[first_index],
			  (magmaFloatComplex **)arrayB, ldb[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}

#ifdef COMPLEX
void magma_cherk_batch(
		       const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
		       const int *N, const int *K, const float *alpha,
		       const BBLAS_Complex32_t **arrayA, const int *lda,
		       const float  *beta, BBLAS_Complex32_t **arrayC,
		       const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MAGMA)
  int first_index = 0;
  magma_init();
  magma_queue_t queue;
  magma_queue_create(&queue);
  
  magmablas_cherk_batched(uplo[first_index], trans[first_index],
			  N[first_index], K[first_index], 
			  alpha[first_index],
			  (magmaFloatComplex const * const * )arrayA, lda[first_index],
			  beta[first_index],
			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue);
  
  magma_sync_wtime(queue );			  
#endif  
}
#endif //#ifdef COMPLEX


/* #ifdef COMPLEX */
/* void magma_chemm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const BBLAS_Complex32_t *alpha, */
/* 		       const BBLAS_Complex32_t **arrayA, const int *lda, */
/* 		       const BBLAS_Complex32_t **arrayB, const int *ldb, */
/* 		       const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_chemm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((magmaFloatComplex *) alpha)[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaFloatComplex *)beta)[first_index], */
/* 			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue); */
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */
/* #endif // ifdef COMPLEX */


/* void magma_csymm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo, */
/* 		       const int *M, const int *N, const BBLAS_Complex32_t *alpha, */
/* 		       const BBLAS_Complex32_t **arrayA, const int *lda, */
/* 		       const BBLAS_Complex32_t **arrayB, const int *ldb, */
/* 		       const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC, */
/* 		       const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */

/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_csymm_batched(side[first_index], uplo[first_index], */
/* 			  M[first_index], N[first_index],  */
/* 			  ((magmaFloatComplex *) alpha)[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaFloatComplex *)beta)[first_index], */
/* 			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* }   */



/* void magma_csyr2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const BBLAS_Complex32_t *alpha, */
/* 			const BBLAS_Complex32_t **arrayA, const int *lda, */
/* 			const BBLAS_Complex32_t **arrayB, const int *ldb, */
/* 			const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_csyr2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((magmaFloatComplex *) alpha)[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaFloatComplex *)beta)[first_index], */
/* 			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */


/* #ifdef COMPLEX */
/* void magma_cher2k_batch( */
/* 			const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans, */
/* 			const int *N, const int *K, const BBLAS_Complex32_t *alpha, */
/* 			const BBLAS_Complex32_t **arrayA, const int *lda, */
/* 			const BBLAS_Complex32_t **arrayB, const int *ldb, */
/* 			const float  *beta, BBLAS_Complex32_t **arrayC, */
/* 			const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_cher2k_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((magmaFloatComplex *) alpha)[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayB, ldb[first_index], */
/* 			  ((magmaFloatComplex *)beta)[first_index], */
/* 			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */
/* #endif //#ifdef COMPLEX */


/* void magma_csyrk_batch( */
/* 		       const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans, */
/* 		       const int *N, const int *K, const BBLAS_Complex32_t *alpha, */
/* 		       const BBLAS_Complex32_t **arrayA, const int *lda, */
/* 		       const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC, */
/* 		       const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_csyrk_batched(uplo[first_index], trans[first_index], */
/* 			  N[first_index], K[first_index],  */
/* 			  ((magmaFloatComplex *) alpha)[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayA, lda[first_index], */
/* 			  ((magmaFloatComplex *)beta)[first_index], */
/* 			  (magmaFloatComplex **)arrayC, ldc[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */

/* void magma_ctrmm_batch( */
/* 		       const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo, */
/* 		       const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag, */
/* 		       const int *M, const int *N, const BBLAS_Complex32_t *alpha, */
/* 		       const BBLAS_Complex32_t **arrayA, const int *lda, */
/* 		       BBLAS_Complex32_t **arrayB, const int *ldb, */
/* 		       const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info) */
/* { */
/* #if defined(BBLAS_WITH_MAGMA) */
/*   int first_index = 0; */
/*   magma_init(); */
/*   magma_queue_t queue; */
/*   magma_queue_create(&queue); */
  
/*   magmablas_ctrmm_batched(side[first_index], uplo[first_index], transA[first_index], */
/* 			  diag[first_index], M[first_index], N[first_index],  */
/* 			  ((magmaFloatComplex *) alpha)[first_index], */
/* 			  (magmaFloatComplex const * const * )arrayA, lda[first_index], */
/* 			  (magmaFloatComplex **)arrayB, ldb[first_index], batch_count, queue); */
  
/*   magma_sync_wtime(queue );			   */
/* #endif   */
/* } */




#undef COMPLEX
