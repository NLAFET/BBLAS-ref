/**
 * @file bblas_testing_common.h
 *
 * @brief BBLAS testing common header for routines that don't require precision.
 *
 *  BBLAS is a software package provided by Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-02-20
 *
 * Contains structures and functions related to testing.
 *
 **/

#ifndef BBLAS_TESTING_COMMON_H
/** Include guard **/
#define BBLAS_TESTING_COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bblas_types.h"
#include "bblas_macros.h"
#include<time.h>
#include "auxiliary.h"
#if defined(BBLAS_WITH_CUBLAS)
#include "cublas_v2.h"
#include "cuda_runtime.h"
#endif


/* Definition of macros */

/** Use oo-norm **/
#define BblasInfNorm      175

/** Set UPLO to lower **/
#define UPLO_LOWER        1
/** Set UPLO to upper **/
#define UPLO_UPPER        2
/** Set UPLO to a random value **/
#define UPLO_RAND         3

/** Set TRANS to no_trans **/
#define NO_TRANS          1
/** Set TRANS to trans **/
#define TRANS             2
/** Set TRANS to conjugate **/
#define CONJ              3
/** Set TRANS to a random value **/
#define TRANS_RAND        4

/** Set SIDE to left **/
#define SIDE_LEFT         1
/** Set SIDE to right **/
#define SIDE_RIGHT        2
/** Set SIDE to a random value **/
#define SIDE_RAND         3

/** Set DIAG to not unit **/
#define DIAG_NO_U         1
/** Set DIAG to unit **/
#define DIAG_U            2
/** Set DIAG to a random value **/
#define DIAG_RAND         3


/**
 * Enumerate to specify which BLAS routine we are testing
 **/
enum BBLAS_ROUTINE{BBLAS_GEMM=1, BBLAS_HEMM=2,
		   BBLAS_HER2K=3, BBLAS_HERK=4,
		   BBLAS_SYMM=5, BBLAS_SYR2K=6,
		   BBLAS_SYRK=7, BBLAS_TRMM=8,
		   BBLAS_TRSM=9};

/**
 * Enumerate to specify which other software packages we are
 * comparing against.
 * BBLAS_OTHER is, by default, an OpenMP implementation of the reference implementation.
 **/
enum BBLAS_TARGET{BBLAS_MKL    =1,
		  BBLAS_CUBLAS =2,
		  BBLAS_MAGMA  =3,
		  BBLAS_OTHER  =4,
		  BBLAS_CUMKL  =5};

/*******************************************************
 *Prototype of functions declared in bblas_ztestings.c
 *******************************************************/

int irandRange(int min_n, int max_n);
BBLAS_Double_t bblas_wtime( void );
int bblas_avgarrayI(int *myArray, int size);
int bblas_minarrayI(int *myArray, int size);
void bblas_print_platform_info();
void bblas_print_header(enum BBLAS_TARGET target);
void bblas_print_omp_info();
void bblas_print_mkl_info(int mkl_sequential);
void bblas_print_time();
char* bblas_getroutine(enum BBLAS_ROUTINE routine);
char* bblas_getoption(enum BBLAS_OPTS opts);
#if defined(BBLAS_WITH_CUDA)
cublasOperation_t op_bb2cu(unsigned int op);
#endif
void bblas_print_tolerance(enum BBLAS_ROUTINE routine, int new_accuracy);
void bblas_print_title(enum BBLAS_ROUTINE routine, enum BBLAS_OPTS batch_opts);
char* bblas_op2char(unsigned int op);

#endif /* BBLAS_TESTING_COMMON_H */
