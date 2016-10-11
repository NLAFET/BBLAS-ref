/**
 *
 * @file bblas_dmkl.h
 *
 *  BBLAS MKL routines for double precision.
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
 * @generated from ./bblas_mkl/bblas_zmkl.h normal z -> d, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_DMKL_H
#define BBLAS_DMKL_H


/* Check whether we are compiling with MKL before including headers */
#if defined(BBLAS_WITH_MKL)
#include <mkl.h>
#endif

#define REAL

#include "bblas.h"
#include<stdlib.h>
#include<stdio.h>
/*
 * MKL TYPES
 */
#ifndef BBLAS_MKL_VOID_H
#define BBLAS_MKL_VOID_H
/**
 * MKL uses void for the return type of real functions.
 * But not for the real arithmetic versions.
 * This typedef is used by our code generation software to avoid this issue.
 **/
typedef void double;
#endif

void mkl_dgemm_batch(const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
					 const int *M,  const int *N, const int *K,
					 const double *alpha,
					 const double **arrayA, const int *lda,
					 const double **arrayB,
					 const int *ldb, const double *beta,
					 double **arrayC, const int *ldc, const int batch_count,
					 enum BBLAS_OPTS batch_opts, int *info);

#undef REAL
#endif
