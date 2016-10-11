/**
 *
 * @file bblas_cmkl.h
 *
 *  BBLAS MKL routines for float _Complex precision.
 *  BBLAS is a software package provided by Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-03-17
 *
 * Contains wrappers to batch MKL functions.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from ./bblas_mkl/bblas_zmkl.h normal z -> c, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_CMKL_H
#define BBLAS_CMKL_H


/* Check whether we are compiling with MKL before including headers */
#if defined(BBLAS_WITH_MKL)
#include <mkl.h>
#endif

#define COMPLEX

#include "bblas.h"
#include<stdlib.h>
#include<stdio.h>
/*
 * MKL TYPES
 */
#ifndef BBLAS_MKL_VOID_H
#define BBLAS_MKL_VOID_H
/**
 * MKL uses void for the return type of complex functions.
 * But not for the real arithmetic versions.
 * This typedef is used by our code generation software to avoid this issue.
 **/
typedef void void;
#endif

void mkl_cgemm_batch(const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
					 const int *M,  const int *N, const int *K,
					 const BBLAS_Complex32_t *alpha,
					 const BBLAS_Complex32_t **arrayA, const int *lda,
					 const BBLAS_Complex32_t **arrayB,
					 const int *ldb, const BBLAS_Complex32_t *beta,
					 BBLAS_Complex32_t **arrayC, const int *ldc, const int batch_count,
					 enum BBLAS_OPTS batch_opts, int *info);

#undef COMPLEX
#endif
