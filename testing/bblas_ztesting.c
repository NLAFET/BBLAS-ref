/**
 * @file bblas_ztesting.c
 *
 * @brief BBLAS testing  for double _Complex routines.
 *
 * BBLAS is a software package provided by Univ. of Manchester,
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#include "bblas_common.h"
#if defined(BBLAS_WITH_MKL)
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <cblas.h>



#define COMPLEX

/**
 * Initialize test parameters to their default values.
 **/

void bblas_zinit_config (bblas_ztest_t *test)
{
    test->gen_uplo       =1;
    test->gen_transA     =1;
    test->gen_transB     =1;
    test->gen_trans      =1;
    test->gen_side       =1;
    test->gen_diag       =1;
    test->minM           = 0;
    test->minN           = 0;
    test->minK           = 0;
    test->maxM           = 0;
    test->maxN           = 0;
    test->maxK           = 0;
    test->minbatch_count = 1;
    test->maxbatch_count = 1;
    test->batch_opts     = 0;
    test->routine        = 1;
    test->nb_test        = 1;
    test->set_error      = 0;
    test->global_error   = 0;
    test->faulty_iter    = 0;
    test->mkl_sequential = 0;
    test->new_accuracy   = 1;
}

/**
 * Set the values of all the BBLAS parameters inside the test structure.
 **/

void bblas_zsettest(bblas_ztest_t *test)
{

    enum BBLAS_ROUTINE routine = test->routine;

    /*
     * Set the value of batch count
     */
    bblas_zset_batch_count(test);

    /*
     * Allocate memory and set values for uplo
     */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K)||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	bblas_zsetuplo(test);
    }

    /*
     * Allocate memory and set values for transA
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_TRMM) ||
	(routine == BBLAS_TRSM))
    {
	bblas_zsettransA(test);
    }

    /*
     * Allocate memory and set values for transB
     */
    if (routine == BBLAS_GEMM)
    {
	bblas_zsettransB(test);
    }

    /*
     * Allocate memory and set values for trans
     */
    if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
    {
	bblas_zsettrans(test);
    }

    /*
     * Allocate memory and set values for side
     */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	bblas_zsetside(test);
    }

    /*
     * Allocate memory and set values for diag
     */
    if  ((routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	bblas_zsetdiag(test);
    }

    /*
     * Allocate memory and set values for  M
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_TRMM) ||
	(routine == BBLAS_TRSM))
    {
	bblas_zsetM(test);
    }

    /*
     * Allocate memory and set values for N, all routines
     */
    bblas_zsetN(test);

    /*
     * Allocate memory and set values for K
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYRK)  ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K))
    {
	bblas_zsetK(test);
    }

    /*
     * Allocate memory and set values for lda, all routines
     */
    bblas_zsetlda(test);

    /*
     * Allocate memory and set values for ldb, all routines
     */
    if ((routine == BBLAS_GEMM)  || (routine == BBLAS_SYMM)  ||
	(routine == BBLAS_HEMM)  || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K) || (routine == BBLAS_TRMM)  ||
	(routine == BBLAS_TRSM))
    {
	bblas_zsetldb(test);
    }

    /*
     * Allocate memory and set values for ldc, all routines
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	(routine == BBLAS_HER2K))
    {
	bblas_zsetldc(test);
    }

    /*
     * Allocate memory and set values for alpha, all routines
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK)  ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K)||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	bblas_zsetalpha(test);
    }

    /*
     * Allocate memory and set values for alpha
     */
    if (routine == BBLAS_HERK)
    {
	bblas_zsetalpha_herk(test);
    }


    /*
     * Allocate memory and set values for beta
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
        (routine == BBLAS_SYR2K))
    {
	bblas_zsetbeta(test);
    }

    /*
     * Allocate memory and set values for beta_herk
     */
    if ((routine == BBLAS_HERK) || (routine == BBLAS_HER2K))
    {
	bblas_zsetbeta_herk(test);
    }

    /*
     * Allocate memory and set values for arrayA
     */
    bblas_zsetarrayA(test);

    /*
     * Allocate memory and set values for arrayB
     */
    if ((routine == BBLAS_GEMM)  || (routine == BBLAS_SYMM)  ||
	(routine == BBLAS_HEMM)  || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K) || (routine == BBLAS_TRMM)  ||
	(routine == BBLAS_TRSM))
    {
	bblas_zsetarrayB(test);
    }


    /*
     * Allocate memory and set values for arrayC
     */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	(routine == BBLAS_HER2K))
    {
	bblas_zsetarrayC(test);
    }

    /* Memory allocation for result and error variables */
    bblas_zmalloc_result_error(test);
}

/**
 * Allocate memory and set the values of uplo
 **/

void bblas_zsetuplo(bblas_ztest_t *test)
{
    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "uplo";
    int random_number;
    /*Initialize random number generation */
    srand ( time(NULL) ) ;

    /*Memory allocation */
    test->uplo =  (enum BBLAS_UPLO*) malloc(nb_data*sizeof(enum BBLAS_UPLO));

    /*Malloc checking */
    bblas_malloc_check(test->uplo, ptr_name);

    /*set UPLO values */
    switch (test->gen_uplo)
    {
	case UPLO_LOWER:

	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		test->uplo[batch_iter] = BblasLower;
	    }
	    break;

	case UPLO_UPPER:

	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		test->uplo[batch_iter] = BblasUpper;
	    }
	    break;

	default:
	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		/*Generate a random number */
		random_number = rand() % 100;

		if (random_number < 50 )
		{
		    test->uplo[batch_iter] = BblasLower;

		}else
		{
		    test->uplo[batch_iter] = BblasUpper;
		}
		break;
	    }
    }
}


/**
 * Allocate memory and set the values of transA
 **/

void bblas_zsettransA(bblas_ztest_t *test)
{
    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "transA";
    int random_number;

    /*Initialize random number generation */
    srand ( time(NULL) ) ;

    /*Memory allocation */
    test->transA =  (enum BBLAS_TRANS*) malloc(nb_data*sizeof(enum BBLAS_TRANS));

    /*Malloc checking */
    bblas_malloc_check(test->transA, ptr_name);

    /*set transA */
    switch (test->gen_transA)
    {
	case NO_TRANS:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->transA[batch_iter] = BblasNoTrans;
	    }
	    break;

	case TRANS:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->transA[batch_iter] = BblasTrans;
	    }
	    break;

	case CONJ:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->transA[batch_iter] = BblasConjTrans;
	    }
	    break;

	default:
	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		/*Generate a random number */
		random_number = rand() % 100;

		if (random_number < 50)
		{
		    test->transA[batch_iter] = BblasNoTrans;

		}else if (random_number < 80)
		{
		    test->transA[batch_iter] = BblasTrans;
		}else
		{
		    test->transA[batch_iter] = BblasConjTrans;
		}
	    }
	    break;
    }
}

/**
 * Allocate memory and set the values of transB
 **/

void bblas_zsettransB(bblas_ztest_t *test)
{
    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "transB";
    int random_number;

    /*Initialize random number generation */
    srand ( time(NULL) ) ;

    /*Memory allocation for transB */
    test->transB =  (enum BBLAS_TRANS*) malloc(nb_data*sizeof(enum BBLAS_TRANS));

    /*checking memory allocation */
    bblas_malloc_check(test->transB, ptr_name);

    /*set transB */
    switch (test->gen_transB)
    {
	case NO_TRANS:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->transB[batch_iter] = BblasNoTrans;
	    }
	    break;

	case TRANS:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->transB[batch_iter] = BblasTrans;
	    }
	    break;

	case CONJ:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->transB[batch_iter] = BblasConjTrans;
	    }
	    break;

	default:
	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		/*Generate a random number */
		random_number = rand() % 100;

		if (random_number < 50)
		{
		    test->transB[batch_iter] = BblasNoTrans;

		}else if (random_number < 80)
		{
		    test->transB[batch_iter] = BblasTrans;
		}else
		{
		    test->transB[batch_iter] = BblasConjTrans;
		}
	    }
	    break;
    }

}


/**
 * Allocate memory and set the values of trans
 **/

void bblas_zsettrans(bblas_ztest_t *test)
{
    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "trans";
    int random_number;

    /*Initialize random number generation */
    srand ( time(NULL) ) ;

    /*Memory allocation for trans */
    test->trans =  (enum BBLAS_TRANS*) malloc(nb_data*sizeof(enum BBLAS_TRANS));

    /*checking memory allocation */
    bblas_malloc_check(test->trans, ptr_name);

    /*Set the values of trans */
    switch (test->gen_trans)
    {
	case NO_TRANS:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->trans[batch_iter] = BblasNoTrans;
	    }
	    break;

	case TRANS:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->trans[batch_iter] = BblasTrans;
	    }
	    break;

	case CONJ:

	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		test->trans[batch_iter] = BblasConjTrans;
	    }
	    break;

	default:
	    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	    {
		/*Generate a random number */
		random_number = rand() % 100;

		if (random_number < 50)
		{
		    test->trans[batch_iter] = BblasNoTrans;

		}else if (random_number < 80)
		{
		    test->trans[batch_iter] = BblasTrans;
		}else
		{
		    test->trans[batch_iter] = BblasConjTrans;
		}
	    }
	    break;
    }
}



/**
 * Allocate memory and set the values of side
 **/

void bblas_zsetside(bblas_ztest_t *test)
{

    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "side";
    int random_number;

    /*Initialize random number generation */
    srand ( time(NULL) ) ;

    /*Memory allocation for side */
    test->side =  (enum BBLAS_SIDE*) malloc(nb_data*sizeof(enum BBLAS_SIDE));

    /*checking memory allocation */
    bblas_malloc_check(test->side, ptr_name);

    /*Set the values of side */
    switch (test->gen_side)
    {
	case SIDE_LEFT:

	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		test->side[batch_iter] = BblasLeft;
	    }
	    break;

	case SIDE_RIGHT:

	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		test->side[batch_iter] = BblasRight;
	    }
	    break;

	default:
	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		/*Generate a random number */
		random_number = rand() % 100;

		if (random_number < 50 )
		{
		    test->side[batch_iter] = BblasLeft;

		}else
		{
		    test->side[batch_iter] = BblasRight;
		}
		break;
	    }
    }
}



/**
 * Allocate memory and set the values of diag
 **/

void bblas_zsetdiag(bblas_ztest_t *test)
{

    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "diag";
    int random_number;

    /*Initialize random number generation */
    srand ( time(NULL) ) ;

    /*Memory allocation for diag */
    test->diag =  (enum BBLAS_DIAG*) malloc(nb_data*sizeof(enum BBLAS_DIAG));

    /*checking memory allocation */
    bblas_malloc_check(test->diag, ptr_name);

    /*Set the values of diag */
    switch (test->gen_diag)
    {
	case DIAG_NO_U:

	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		test->diag[batch_iter] = BblasNonUnit;
	    }
	    break;

	case DIAG_U:

	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		test->diag[batch_iter] = BblasUnit;
	    }
	    break;

	default:
	    for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	    {
		/*Generate a random number */
		random_number = rand() % 100;

		if (random_number < 50 )
		{
		    test->diag[batch_iter] = BblasNonUnit;

		}else
		{
		    test->diag[batch_iter] = BblasUnit;
		}
		break;
	    }
    }
}


/**
 * Allocate memory and set the values of M
 **/

void bblas_zsetM(bblas_ztest_t *test)
{

    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "M";

    /*Memory allocation for M */
    test->M = (int*) malloc(nb_data*sizeof(int));

    /*checking memory allocation */
    bblas_malloc_check(test->M, ptr_name);

    /*Set the values of M */
    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
    {
	test->M[batch_iter] = irandRange(test->minM, test->maxM);
    }
}


/**
 * Allocate memory and set the values of N
 **/

void bblas_zsetN(bblas_ztest_t *test)
{

    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "N";

    /*Memory allocation for N */
    test->N = (int*) malloc(nb_data*sizeof(int));

    /*checking memory allocation */
    bblas_malloc_check(test->N, ptr_name);

    /*Set the values of N */
    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
    {
	test->N[batch_iter] = irandRange(test->minN, test->maxN);
    }
}


/**
 * Allocate memory and set the values of K
 **/

void bblas_zsetK(bblas_ztest_t *test)
{
    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "K";

    /*Memory allocation for K */
    test->K = (int*) malloc(nb_data*sizeof(int));

    /*checking memory allocation */
    bblas_malloc_check(test->K, ptr_name);

    /*Set the values of K */
    for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
    {
	test->K[batch_iter] = irandRange(test->minK, test->maxK);
    }
}


/**
 * Allocate memory and set the values of lda
 **/

void bblas_zsetlda(bblas_ztest_t *test)
{
    int batch_iter;
    int nb_data                 = bblas_znbdata(test);
    int routine                 = test->routine;
    char ptr_name[NAME_LENGTH]  = "lda";

    /*Memory allocation for lda */
    test->lda = (int*) malloc(nb_data*sizeof(int));

    /*checking memory allocation */
    bblas_malloc_check(test->lda, ptr_name);

    /*LDA for GEMM */
    if (routine == BBLAS_GEMM)
    {
	for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	{
	    if (test->transA[batch_iter] == BblasNoTrans)
	    {
		test->lda[batch_iter] = test->M[batch_iter];
	    }else
	    {
		test->lda[batch_iter] =test->K[batch_iter] ;
	    }
	}
    }

    /*LDA for SYMM, HEMM, TRMM AND TRSM */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	{
	    if (test->side[batch_iter] == BblasLeft)
	    {
		test->lda[batch_iter] = test->M[batch_iter];
	    }else
	    {
		test->lda[batch_iter] =test->N[batch_iter] ;
	    }
	}
    }

    /*LDA for SYRK, HERK, SYR2K, HER2K */
    if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
    {
	for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	{
	    if (test->trans[batch_iter] == BblasNoTrans)
	    {
		test->lda[batch_iter] = test->N[batch_iter];
	    }else
	    {
		test->lda[batch_iter] =test->K[batch_iter] ;
	    }
	}
    }
}

/**
 * Allocate memory and set the values of ldb
 **/

void bblas_zsetldb(bblas_ztest_t *test)
{
    int batch_iter;
    int routine                 = test->routine;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "ldb";

    /*Memory allocation for ldb */
    test->ldb = (int*) malloc(nb_data*sizeof(int));

    /*checking memory allocation */
    bblas_malloc_check(test->ldb, ptr_name);

    /*LDB for GEMM */
    if ((routine == BBLAS_GEMM) )
    {
	for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	{
	    if (test->transB[batch_iter] == BblasNoTrans)
	    {
		test->ldb[batch_iter] = test->K[batch_iter];
	    }else
	    {
		test->ldb[batch_iter] =test->N[batch_iter] ;
	    }
	}
    }

    /*LDB SYMM, HEMM, TRMM AND TRSM */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	{
	    test->ldb[batch_iter] = test->M[batch_iter];
	}
    }

    /*LDB for SYR2K, HER2K */
    if ((routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
    {
	for( batch_iter =0; batch_iter < nb_data; batch_iter++)
	{
	    if (test->trans[batch_iter] == BblasNoTrans)
	    {
		test->ldb[batch_iter] = test->N[batch_iter];
	    }else
	    {
		test->ldb[batch_iter] =test->K[batch_iter] ;
	    }
	}
    }
}


/**
 * Allocate memory and set the values of ldc
 **/

void bblas_zsetldc(bblas_ztest_t *test)
{
    int batch_iter;
    int routine                 = test->routine;
    int nb_data                 = bblas_znbdata(test);
    char ptr_name[NAME_LENGTH]  = "ldc";

    /*Memory allocation for ldc */
    test->ldc = (int*) malloc(nb_data*sizeof(int));

    /*checking memory allocation */
    bblas_malloc_check(test->ldc, ptr_name);

    /*LDC for GEMM */
    if (routine == BBLAS_GEMM)
    {
	for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	{
	    test->ldc[batch_iter] = test->M[batch_iter];
	}
    }

    /*LDC for SYMM, HEMM */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM))
    {
	for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	{
	    test->ldc[batch_iter] = test->M[batch_iter];
	}
    }

    /*LDC for SYRK, HERK, SYR2K, HER2K */
    if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
    {
	for( batch_iter =0; batch_iter < nb_data ; batch_iter++)
	{
	    test->ldc[batch_iter] = test->N[batch_iter];
	}
    }
}


/**
 * Allocate memory and set the values of alpha
 **/

void bblas_zsetalpha(bblas_ztest_t *test)
{
    int batch_iter;
    int batch_count             = test->batch_count;
    char ptr_name[NAME_LENGTH]  = "alpha";

    /*Memory allocation for alpha */
    test->alpha = (BBLAS_Complex64_t*) malloc(batch_count*sizeof(BBLAS_Complex64_t));

    /*checking memory allocation */
    bblas_malloc_check(test->alpha, ptr_name);

    /* Set value of alpha */
    for( batch_iter =0; batch_iter < batch_count; batch_iter++)
    {
	test->alpha[batch_iter] = ((BBLAS_Complex64_t)rand()/(BBLAS_Complex64_t)RAND_MAX);
    }
}


/**
 * Allocate memory and set the values of alpha specifically for herk.
 **/

void bblas_zsetalpha_herk(bblas_ztest_t *test)
{

    int batch_iter;
    int batch_count             = test->batch_count;
    char ptr_name[NAME_LENGTH]  = "alpha_herk";

    /*Memory allocation for alpha_herk */
    test->alpha_herk = (double*) malloc(batch_count*sizeof(double));

    /*checking memory allocation */
    bblas_malloc_check(test->alpha_herk, ptr_name);

    /* Set value of alpha_herk */
    for( batch_iter =0; batch_iter < batch_count; batch_iter++)
    {
	test->alpha_herk[batch_iter] = ((double)rand()/(double)RAND_MAX);
    }
}


/**
 * Allocate memory and set the values of beta.
 **/

void bblas_zsetbeta(bblas_ztest_t *test)
{
    int batch_iter;
    int batch_count             = test->batch_count;
    char ptr_name[NAME_LENGTH]  = "beta";

    /*Memory allocation for beta */
    test->beta = (BBLAS_Complex64_t*) malloc(batch_count*sizeof(BBLAS_Complex64_t));

    /*checking memory allocation */
    bblas_malloc_check(test->beta, ptr_name);

    /* Set value of beta */
    for( batch_iter =0; batch_iter < batch_count; batch_iter++)
    {
	test->beta[batch_iter] = ((BBLAS_Complex64_t)rand()/(BBLAS_Complex64_t)RAND_MAX);
    }
}

/**
 * Allocate memory and set the values of beta specifically for herk.
 **/

void bblas_zsetbeta_herk(bblas_ztest_t *test)
{
     int batch_iter;
     int batch_count             = test->batch_count;
     char ptr_name[NAME_LENGTH]  = "beta_herk";

     /*Memory allocation for beta_herk */
     test->beta_herk = (double*) malloc(batch_count*sizeof(double));

     /*checking memory allocation */
     bblas_malloc_check(test->beta_herk, ptr_name);

     /* Set value of beta_herk */
     for( batch_iter =0; batch_iter < batch_count; batch_iter++)
     {
	 test->beta_herk[batch_iter] = ((double)rand()/(double)RAND_MAX);
     }
}


/**
 * Allocate memory and set the values of arrayA.
 **/

void bblas_zsetarrayA(bblas_ztest_t *test)
{
    int batch_iter, nb_row, nb_col;
     int first_index             = 0;
     int routine                 = test->routine;
     int batch_count             = test->batch_count;
     char ptr_name[NAME_LENGTH]  = "arrayA";

     int IONE     = 1;
     int ISEED[4] ={0,0,0,1};

     /*Memory allocation for **arrayA  */
     test->arrayA = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
     bblas_malloc_check(test->arrayA, ptr_name);

     if( test->batch_opts == BBLAS_VARIABLE )
     {
	 for( batch_iter =0; batch_iter < batch_count; batch_iter++)
	 {
	     nb_row = test->lda[batch_iter];

	     /* nb_col  for GEMM */
	     if (routine == BBLAS_GEMM)
	     {
		 if (test->transA[batch_iter] == BblasNoTrans)
		 {
		     nb_col =test->K[batch_iter];
		 }else
		 {
		     nb_col =test->M[batch_iter];
		 }
	     }

	     /*  nb_col  SYMM, HEMM, TRMM AND TRSM */
	     if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
		 (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	     {
		 if(test->side[batch_iter] == BblasLeft )
		 {
		     nb_col = test->M[batch_iter];
		 }else
		 {
		     nb_col = test->N[batch_iter];
		 }
	     }

	     /* nb_col for SYRK, HERK, SYR2K, HER2K */
	     if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
		 (routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
	     {
		 if (test->trans[batch_iter] == BblasNoTrans)
		 {
		     nb_col = test->K[batch_iter];
		 }else
		 {
		     nb_col = test->N[batch_iter];
		 }
	     }

	     /*Matrix filling */
	     test->arrayA[batch_iter] = (BBLAS_Complex64_t *) malloc(nb_row*nb_col* sizeof(BBLAS_Complex64_t ));
	     bblas_malloc_check(test->arrayA[batch_iter], ptr_name);

#if defined(BBLAS_WITH_MKL)
	     LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, (MKL_Complex16*) test->arrayA[batch_iter]);
#else
		 LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, test->arrayA[batch_iter]);
#endif

	     if( (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	       {
		 for(int i=0; i<max(nb_row,nb_col); i++)
		   {
		     test->arrayA[batch_iter][nb_col*i+i] = test->arrayA[batch_iter][nb_col*i+i] + 1.0;
		   }
	       }

	     if(routine == BBLAS_HEMM )
	       {
		 for(int i=0; i< nb_row; i++)
		   {
		     test->arrayA[batch_iter][nb_col*i+i] = creal(test->arrayA[batch_iter][nb_col*i+i]);
		   }
	       }
	 }

     }else if( test->batch_opts == BBLAS_FIXED )
     {
	 nb_row = test->lda[first_index];

	 /* nb_col  for GEMM */
	 if (routine == BBLAS_GEMM)
	 {
	     if (test->transA[first_index] == BblasNoTrans)
	     {
		 nb_col =test->K[first_index];
	     }else
	     {
		 nb_col =test->M[first_index];
	     }
	 }

	 /*  nb_col  SYMM, HEMM, TRMM AND TRSM */
	 if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	     (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	 {
	     if(test->side[first_index] == BblasLeft )
	     {
		 nb_col = test->M[first_index];
	     }else
	     {
		 nb_col = test->N[first_index];
	     }
	 }

	 /* nb_col for SYRK, HERK, SYR2K, HER2K */
	 if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	     (routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
	 {
	     if (test->trans[first_index] == BblasNoTrans)
	     {
		 nb_col = test->K[first_index];
	     }else
	     {
		 nb_col = test->N[first_index];
	     }
	 }

	 /*Matrix filling */
	 for( batch_iter =0; batch_iter < batch_count; batch_iter++)
	 {
	     test->arrayA[batch_iter] = (BBLAS_Complex64_t *) malloc(nb_row*nb_col* sizeof(BBLAS_Complex64_t ));
	     bblas_malloc_check(test->arrayA[batch_iter], ptr_name);

#if defined(BBLAS_WITH_MKL)
	     LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, (MKL_Complex16*) test->arrayA[batch_iter]);
#else
		 LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, test->arrayA[batch_iter]);
#endif

	     if( (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	       {
		 for(int i=0; i<max(nb_row,nb_col); i++)
		   {
		     test->arrayA[batch_iter][nb_col*i+i] = test->arrayA[batch_iter][nb_col*i+i] + 1.0;
		   }
	       }

	     if(routine == BBLAS_HEMM )
	       {
		 for(int i=0; i<max(nb_row,nb_col); i++)
		   {
		     test->arrayA[batch_iter][nb_col*i+i] = creal(test->arrayA[batch_iter][nb_col*i+i]);
		   }
	       }
	 }
     }else
     {
	 bblas_error("bblas_ztesting.c", "wrong batch_opts value");
     }
}

/**
 * Allocate memory and set the values of arrayB.
 **/

void bblas_zsetarrayB(bblas_ztest_t *test)
{
    int batch_iter, nb_row, nb_col;
    int first_index             = 0;
    int routine                 = test->routine;
    int batch_count             = test->batch_count;
    char ptr_name[NAME_LENGTH]  = "arrayB";
    int max_work_size = max(test->maxK, max(test->maxM, test->maxN));
    int IONE     = 1;
    int ISEED[4] ={0,0,0,1};

    /*Memory allocation for arrayB */
    test->arrayB = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
    bblas_malloc_check(test->arrayB, ptr_name);

    test->Binitnorm = (double *)malloc(batch_count*sizeof(double));
    bblas_malloc_check(test->Binitnorm, "Binitnorm");

    double *work = (double *)malloc(max_work_size*sizeof(double));
    bblas_malloc_check(work, "work");

    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    nb_row =test->ldb[batch_iter];

	    /* nb_col  for GEMM */
	    if (routine == BBLAS_GEMM)
	    {
		if (test->transB[batch_iter] == BblasNoTrans)
		{
		    nb_col = test->N[batch_iter];
		}else
		{
		    nb_col = test->K[batch_iter];
		}
	    }

	    /* nb_col for SYMM, HEMM, TRMM AND TRSM */
	    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
		(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	    {
		nb_col =test->N[batch_iter];
	    }

	    /* nb_col for SYR2K, HER2K */
	    if ((routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
	    {
		if (test->trans[batch_iter] == BblasNoTrans)
		{
		    nb_col = test->K[batch_iter];
		}else
		{
		    nb_col = test->N[batch_iter];
		}
	    }

	    test->arrayB[batch_iter] = (BBLAS_Complex64_t *) malloc(nb_row*nb_col* sizeof(BBLAS_Complex64_t ));
	    bblas_malloc_check(test->arrayB[batch_iter], ptr_name);

#if defined(BBLAS_WITH_MKL)
	    LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, (MKL_Complex16*) test->arrayB[batch_iter]);
#else
		LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, test->arrayB[batch_iter]);
#endif

	    /*Compute the infinity norm of B  */
	    if( (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	    {
#if defined(BBLAS_WITH_MKL)
		test->Binitnorm[batch_iter] = (BBLAS_Complex64_t) LAPACKE_zlange_work(LAPACK_COL_MAJOR,
								  'I', nb_row, nb_col,
								  (MKL_Complex16*) test->arrayB[batch_iter],
								  test->ldb[batch_iter], work);
#else
		test->Binitnorm[batch_iter] = LAPACKE_zlange_work(LAPACK_COL_MAJOR,
								  'I', nb_row, nb_col,
								  test->arrayB[batch_iter],
								  test->ldb[batch_iter], work);
#endif
	    }
	}

    }else if( test->batch_opts == BBLAS_FIXED )
    {
	nb_row =test->ldb[first_index];

	/* nb_col  for GEMM */
	if (routine == BBLAS_GEMM)
	{
	    if (test->transB[first_index] == BblasNoTrans)
	    {
		nb_col = test->N[first_index];
	    }else
	    {
		nb_col = test->K[first_index];
	    }
	}

	/* nb_col for SYMM, HEMM, TRMM AND TRSM */
	if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	    (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	{
	    nb_col =test->N[first_index];
	}

	/* nb_col for SYR2K, HER2K */
	if ((routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
	{
	    if (test->trans[first_index] == BblasNoTrans)
	    {
		nb_col = test->K[first_index];
	    }else
	    {
		nb_col = test->N[first_index];
	    }
	}
	/*Matrix filling */
	for( batch_iter =0; batch_iter < batch_count; batch_iter++)
	{
	     test->arrayB[batch_iter] = (BBLAS_Complex64_t *) malloc(nb_row*nb_col* sizeof(BBLAS_Complex64_t ));
#if defined(BBLAS_WITH_MKL)
	     LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, (MKL_Complex16*) test->arrayB[batch_iter]);
#else
		 LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, test->arrayB[batch_iter]);
#endif

	     /*Compute the infinity norm of B  */
	     if( (routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
	     {
#if defined(BBLAS_WITH_MKL)
		 test->Binitnorm[batch_iter] = (BBLAS_Complex64_t) LAPACKE_zlange_work(LAPACK_COL_MAJOR,
								   'I', nb_row, nb_col,
								   (MKL_Complex16*) test->arrayB[batch_iter],
								   test->ldb[first_index], work);
#else
		 test->Binitnorm[batch_iter] = LAPACKE_zlange_work(LAPACK_COL_MAJOR,
								   'I', nb_row, nb_col,
								   test->arrayB[batch_iter],
								   test->ldb[first_index], work);
#endif
	     }
	 }
    }else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value");
    }
    /*Free work */
    free(work);
}


/**
 * Allocate memory and set the values of arrayC.
 **/

void bblas_zsetarrayC(bblas_ztest_t *test)
{
    int batch_iter, nb_row, nb_col;
    int first_index             = 0;
    int batch_count             = test->batch_count;
    int max_work_size = max(test->maxK, max(test->maxM, test->maxN));
    char ptr_name[NAME_LENGTH]  = "arrayC";

    int IONE     = 1;
    int ISEED[4] ={0,0,0,1};

    test->arrayC = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
    bblas_malloc_check(test->arrayC, ptr_name);

    test->Cinitnorm = (double *)malloc(batch_count*sizeof(double));
    bblas_malloc_check(test->Cinitnorm, "Cinitnorm");

    double *work = (double *)malloc(max_work_size*sizeof(double));
    bblas_malloc_check(work, "work");


    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    nb_row =test->ldc[batch_iter];
	    nb_col =test->N[batch_iter];

	    test->arrayC[batch_iter] = (BBLAS_Complex64_t *) malloc(nb_row*nb_col* sizeof(BBLAS_Complex64_t ));
	    bblas_malloc_check(test->arrayC[batch_iter], ptr_name);
#if defined(BBLAS_WITH_MKL)
	    LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, (MKL_Complex16*) test->arrayC[batch_iter]);
#else
		LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, test->arrayC[batch_iter]);
#endif

	    if( (test->routine == BBLAS_HERK) || (test->routine == BBLAS_HER2K))
	      {
		for(int i=0; i<max(nb_row,nb_col); i++)
		  {
		    test->arrayC[batch_iter][nb_col*i+i] = creal(test->arrayC[batch_iter][nb_col*i+i]);
		  }
	      }

	    /*Compuptation of the norm of C */
#if defined(BBLAS_WITH_MKL)
	    test->Cinitnorm[batch_iter] = (BBLAS_Complex64_t) LAPACKE_zlange_work(LAPACK_COL_MAJOR,
								  'I', nb_row, nb_col, (MKL_Complex16*) test->arrayC[batch_iter],
								  test->ldc[batch_iter], work);
#else
		test->Cinitnorm[batch_iter] = LAPACKE_zlange_work(LAPACK_COL_MAJOR,
							      'I', nb_row, nb_col, test->arrayC[batch_iter],
								  test->ldc[batch_iter], work);
#endif
	}

    }else if( test->batch_opts == BBLAS_FIXED )
    {
	nb_row =test->ldc[first_index];
	nb_col =test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    test->arrayC[batch_iter] = (BBLAS_Complex64_t *) malloc(nb_row*nb_col* sizeof(BBLAS_Complex64_t ));
	    bblas_malloc_check(test->arrayC[batch_iter], ptr_name);

#if defined(BBLAS_WITH_MKL)
	    LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, (MKL_Complex16*) test->arrayC[batch_iter]);
#else
		LAPACKE_zlarnv_work(IONE, ISEED, nb_row*nb_col, test->arrayC[batch_iter]);
#endif

	    if( (test->routine == BBLAS_HERK) || (test->routine == BBLAS_HER2K))
	      {
		for(int i=0; i<max(nb_row,nb_col); i++)
		  {
		    test->arrayC[batch_iter][nb_col*i+i] = creal(test->arrayC[batch_iter][nb_col*i+i]);
		  }
	      }

	    /*Compuptation of the norm of C */
#if defined(BBLAS_WITH_MKL)
	    test->Cinitnorm[batch_iter] = (BBLAS_Complex64_t) LAPACKE_zlange_work(LAPACK_COL_MAJOR,
									  'I', nb_row, nb_col,
									  (MKL_Complex16*) test->arrayC[batch_iter],
									  test->ldc[first_index], work);
		#else
		test->Cinitnorm[batch_iter] = LAPACKE_zlange_work(LAPACK_COL_MAJOR,
							      'I',nb_row, nb_col, test->arrayC[batch_iter],
								  test->ldc[first_index], work);
		#endif
	}
    } else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value");
    }

    /*Free work */
    free(work);
}

/**
 * Allocate memory for error checking.
 **/

void bblas_zmalloc_result_error(bblas_ztest_t *test)
{

    int batch_count             = test->batch_count;

    /*Memory for error computation  */
    switch(test->target)
    {
	case BBLAS_MKL:
	    test->mkl_result  = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
	    bblas_malloc_check(test->mkl_result, "mkl_result");

	    test->mkl_error = (double*) malloc(batch_count*sizeof(double));
	    bblas_malloc_check(test->mkl_error, "mkl_error");
	    break;

	case BBLAS_CUBLAS:
	case BBLAS_MAGMA:
	    test->device_result = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
	    bblas_malloc_check(test->device_result, "device_result");

	    test->device_error = (double*) malloc(batch_count*sizeof(double));
	    bblas_malloc_check(test->device_error, "device_error");
	    break;

	case BBLAS_OTHER:
	    test->other_result  = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
	    bblas_malloc_check(test->other_result, "other_result");

	    test->other_error = (double*) malloc(batch_count*sizeof(double));
	    bblas_malloc_check(test->other_error, "other_error");
	    break;

	case BBLAS_CUMKL:
	    test->mkl_result    = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
	    bblas_malloc_check(test->mkl_result, "mkl_result");

	    test->device_result = (BBLAS_Complex64_t**) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
	    bblas_malloc_check(test->device_result, "device_result");

	    test->mkl_error = (double*) malloc(batch_count*sizeof(double));
	    bblas_malloc_check(test->mkl_error, "mkl_error");

	    test->device_error = (double*) malloc(batch_count*sizeof(double));
	    bblas_malloc_check(test->device_error, "device_error");
	    break;

	default:
	    printf("Memory alloation for error: Target no defined\n");
	    exit(EXIT_FAILURE);
    }

    /*Memory for info  */
    test->info = (int*) malloc(batch_count*sizeof(int));
}


/**
 * Allocate memory and copy the initial arrayC values.
 * This is required to take norms after the computation is complete.
 **/

void bblas_zcopy_Cinit(bblas_ztest_t *test, BBLAS_Complex64_t **C_copy)
{
    /*Local variables */
    enum BBLAS_ROUTINE routine  = test->routine;
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int ldc, N;

    if (!((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	  (routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	  (routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	  (routine == BBLAS_HER2K)))
    {
	printf("BBLAS FATAL ERROR: bblas_ztesting.c():\n");
	printf("\t bblas_zcopy_Cinit not defined for %s \n", bblas_getroutine(test->routine));
	exit(EXIT_FAILURE);
    }

    if( test->batch_opts == BBLAS_VARIABLE ) // Varible size
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    ldc = test->ldc[batch_iter];
	    N   = test->N[batch_iter];
	    C_copy[batch_iter] = (BBLAS_Complex64_t *) malloc(ldc*N* sizeof(BBLAS_Complex64_t ));

	    /*Copy the matrix {C}_i */
	    cblas_zcopy (ldc*N, test->arrayC[batch_iter], 1, C_copy[batch_iter], 1);
	}

    }else if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	ldc = test->ldc[first_index];
	N = test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    C_copy[batch_iter] = (BBLAS_Complex64_t *) malloc(ldc*N* sizeof(BBLAS_Complex64_t ));

	    /*Copy the matrix {C}_i */
	    cblas_zcopy (ldc*N, test->arrayC[batch_iter], 1,  C_copy[batch_iter], 1);
	}
    }else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value\n");
    }
}

/**
 * Allocate memory and copy the initial arrayB values.
 * This is required to take norms after the computation is complete.
 **/

void bblas_zcopy_Binit(bblas_ztest_t *test, BBLAS_Complex64_t **B_copy)
{
    /*Local variables */
    enum BBLAS_ROUTINE routine  = test->routine;
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int ldb, N;

    if (!((routine == BBLAS_TRMM) || (routine == BBLAS_TRSM)))
    {
	printf("BBLAS FATAL ERROR: bblas_ztesting.c():\n");
	printf("\t bblas_zcopy_Binit not defined for %s \n", bblas_getroutine(test->routine));
	exit(EXIT_FAILURE);
    }

    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    ldb = test->ldb[batch_iter];
	    N   = test->N[batch_iter];

	    B_copy[batch_iter] = (BBLAS_Complex64_t *) malloc(ldb*N* sizeof(BBLAS_Complex64_t ));

	    /*Copy the matrix {B}_i */
	    cblas_zcopy (ldb*N, test->arrayB[batch_iter], 1, B_copy[batch_iter], 1);
	}

    }else if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	ldb = test->ldb[first_index];
	N = test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    B_copy[batch_iter] = (BBLAS_Complex64_t *) malloc(ldb*N* sizeof(BBLAS_Complex64_t ));

	    /*Copy the matrix {B}_i */
	    cblas_zcopy (ldb*N, test->arrayB[batch_iter], 1, B_copy[batch_iter], 1);
	}
    }else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value\n");
    }
}

/**
 * Compute the relative error of each batch operation.
 **/

void bblas_zcheck_Cfinal(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{

    /*Local variables */
    enum BBLAS_ROUTINE routine  = test->routine;
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int ldc, N;

    /* Error calculation variables */
    double Cnorm, Error_norm;
    BBLAS_Complex64_t alpha =-1;

    /*Check if the call has been made by the correct routine */
    if (!((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	  (routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	  (routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	  (routine == BBLAS_HER2K)))
    {
	printf("BBLAS FATAL ERROR: bblas_ztesting.c():\n");
	printf("\t bblas_zcheck_Cfinal not defined for %s \n", bblas_getroutine(test->routine));
	exit(EXIT_FAILURE);
    }

    /*Temporary buffer to save (test-arrayC -C_final) */
    BBLAS_Complex64_t **C_diff;
    C_diff = (BBLAS_Complex64_t **) malloc(batch_count*sizeof(BBLAS_Complex64_t *));

    /*Make a copy of test->arrayC in C_diff */
    bblas_zcopy_Cinit(test, C_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    ldc = test->ldc[batch_iter];
	    N   = test->N[batch_iter];

	    /*Computation of the Frobenus norm of  {C}_batch_iter */
	    Cnorm = cblas_dznrm2(ldc*N, test->arrayC[batch_iter], 1);

	    /*Compute the error */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the norm assoicated  with the error */
	    Error_norm = cblas_dznrm2(ldc*N, C_diff[batch_iter], 1);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = Error_norm/Cnorm;
		    break;

		case BBLAS_CUBLAS:
                case BBLAS_MAGMA:
  		    test->device_error[batch_iter] = Error_norm/Cnorm;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = Error_norm/Cnorm;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	ldc = test->ldc[first_index];
	N   = test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    /*Computation of the Frobenus norm of  {C}_batch_iter */
	    Cnorm = cblas_dznrm2(ldc*N, test->arrayC[batch_iter], 1);

	    /*Compute the error */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the norm assoicated  with the error */
	    Error_norm = cblas_dznrm2(ldc*N, C_diff[batch_iter], 1);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = Error_norm/Cnorm;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = Error_norm/Cnorm;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = Error_norm/Cnorm;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal: Fixed, Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value");
    }

    /*Free  C_diff */
    for(batch_iter=0; batch_iter < batch_count ; batch_iter++)
    {
	free(C_diff[batch_iter]);
    }
    free(C_diff);
}

/**
 * Compute the relative error of each batch operation.
 **/

void bblas_zcheck_Bfinal(bblas_ztest_t *test, BBLAS_Complex64_t **B_final)
{

    /*Local variables */
    enum BBLAS_ROUTINE routine  = test->routine;
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int ldb, N;

    /* Error calculation variables */
    double Bnorm, Error_norm;
    BBLAS_Complex64_t alpha =-1;

    /*Check if the call has been made by the correct routine */
    if (!((routine == BBLAS_TRSM) || (routine == BBLAS_TRMM)))
    {
	printf("BBLAS FATAL ERROR: bblas_ztesting.c():\n");
	printf("\t bblas_zcheck_Bfinal not defined for %s \n", bblas_getroutine(test->routine));
	exit(EXIT_FAILURE);
    }

    /*Temporary buffer to save (test-arrayB -B_final) */
    BBLAS_Complex64_t **B_diff;
    B_diff = (BBLAS_Complex64_t **) malloc(batch_count*sizeof(BBLAS_Complex64_t *));

    /*Make a copy of test->arrayB in B_diff */
    bblas_zcopy_Binit(test, B_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    ldb = test->ldb[batch_iter];
	    N   = test->N[batch_iter];

	    /*Computation of the Frobenus norm of  {B}_batch_iter */
	    Bnorm = cblas_dznrm2(ldb*N, test->arrayB[batch_iter], 1);

	    /*Compute the error */
	    cblas_zaxpy (ldb*N, CBLAS_SADDR(alpha), B_final[batch_iter], 1, B_diff[batch_iter], 1);

	    /*Compute  the norm assoicated  with the error */
	    Error_norm = cblas_dznrm2(ldb*N, B_diff[batch_iter], 1);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = Error_norm/Bnorm;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = Error_norm/(Bnorm);
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = Error_norm/(Bnorm);
		    break;
		default:
		    printf("In bblas_zcheck_Bfinal(): Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED )
    {
	ldb = test->ldb[first_index];
	N = test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    /*Computation of the Frobenus norm of  {B}_batch_iter */
	    Bnorm = cblas_dznrm2(ldb*N, test->arrayB[batch_iter], 1);

	    /*Compute the error */
	    cblas_zaxpy (ldb*N, CBLAS_SADDR(alpha), B_final[batch_iter], 1, B_diff[batch_iter], 1);

	    /*Compute  the norm assoicated  with the error */
	    Error_norm = cblas_dznrm2(ldb*N, B_diff[batch_iter], 1);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = Error_norm/(Bnorm);
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = Error_norm/(Bnorm);
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = Error_norm/(Bnorm);
		    break;
		default:
		    printf("In bblas_zcheck_Bfinal(): Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value");
    }

    /*Free  B_diff */
    for(batch_iter=0; batch_iter < batch_count ; batch_iter++)
    {
	free(B_diff[batch_iter]);
    }
    free(B_diff);
}

/**
 * Free the memory associated with the test structure.
 **/

void bblas_zfreetest( bblas_ztest_t *test )
{
    enum BBLAS_ROUTINE routine  = test->routine;
    /*Free uplo */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K)||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	free(test->uplo);
    }

    /*Free TRANSA */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_TRMM) ||
	(routine == BBLAS_TRSM))
    {
	free( test->transA );
    }

    /*Free TRANSB */
    if (routine == BBLAS_GEMM)
    {
	free( test->transB );
    }

    /*Free TRANS */
    if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
    {
	free( test->trans );
    }

    /*Free SIDE*/
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	free( test->side );
    }

    /*Free DIAG*/
    if  ((routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	free(test->diag);
    }


    /* Free M  memory*/
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_TRMM) ||
	(routine == BBLAS_TRSM))
    {
	free( test->M );
    }

    /* Free N memory*/
    free( test->N );

    /* Free K memory*/
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYRK)  ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K))
    {
	free( test->K );
    }

    /*Free LDA memory */
    free( test->lda );

    /*Free LDB memory */
    if ((routine == BBLAS_GEMM)  || (routine == BBLAS_SYMM)  ||
	(routine == BBLAS_HEMM)  || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K) || (routine == BBLAS_TRMM)  ||
	(routine == BBLAS_TRSM))
    {
	free( test->ldb );
    }

    /*Free LDC memory */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	(routine == BBLAS_HER2K))
    {
	free( test->ldc );
    }

    /*Free ALPHA memory */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) ||(routine == BBLAS_SYRK)  ||
	(routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K)||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
    {
	free(test->alpha);

    }else if (routine == BBLAS_HERK)
    {
	free(test->alpha_herk);
    }

    /*Free ALPHA BETA memory */
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	(routine == BBLAS_SYR2K))
    {
	free( test->beta  );

    }else if ((routine == BBLAS_HERK) || (routine == BBLAS_HER2K))
    {
	free(test->beta_herk);
    }

    /*Free matrices {A}_i  */
    for ( int batch_iter = 0; batch_iter < test->batch_count; batch_iter++ )
    {
	free( test->arrayA[batch_iter] );
    }
    free( test->arrayA );

    /*Free matrices  {B}_i*/
    if ((routine == BBLAS_GEMM)  || (routine == BBLAS_SYMM)  ||
	(routine == BBLAS_HEMM)  || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K) || (routine == BBLAS_TRMM)  ||
	(routine == BBLAS_TRSM))
    {
	for ( int batch_iter = 0; batch_iter < test->batch_count; batch_iter++ )
	{
	    free( test->arrayB[batch_iter] );
	}
	free( test->arrayB );
    }


    /*Free matrices  {C}_i*/
    if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	(routine == BBLAS_HER2K))
    {
	for ( int batch_iter = 0; batch_iter < test->batch_count; batch_iter++ )
	{
	    free( test->arrayC[batch_iter] );
	}
	free( test->arrayC );
    }

    /*Free memory allocated for error computation */
    switch(test->target)
    {
	case BBLAS_MKL:
	    free(test->mkl_error);
	    //free(test->group_size);
	    for ( int batch_iter = 0; batch_iter < test->batch_count; batch_iter++ )
	    {
		free( test->mkl_result[batch_iter] );
	    }
	    free( test->mkl_result );

	    break;

	case BBLAS_CUBLAS:
	case BBLAS_MAGMA:
	    free(test->device_error);

	    for (int batch_iter =0; batch_iter < test->batch_count; batch_iter++)
	    {
		free(test->device_result[batch_iter]);
	    }

	    free(test->device_result);

	    break;

	case BBLAS_OTHER:

	    free(test->other_error);

	    for (int batch_iter =0; batch_iter < test->batch_count; batch_iter++)
	    {
		free(test->other_result[batch_iter]);
	    }
	    free(test->other_result);

	    break;

	case BBLAS_CUMKL:
	    free(test->mkl_error);
	    free(test->device_error);
	    free(test->group_size);

	    for (int batch_iter =0; batch_iter < test->batch_count; batch_iter++)
	    {
		free(test->device_result[batch_iter]);
		free(test->mkl_result[batch_iter]);
	    }
	    free(test->device_result);
	    free(test->mkl_result);

	    break;

	default:
	    printf("In bblas_zfreetest(): Target no defined\n");
	    exit(EXIT_FAILURE);
    }
    /*Free INFO */
    free( test->info );

    /*Free cuda memory */
    if( (test->target == BBLAS_CUBLAS) || (test->target == BBLAS_MAGMA))
    {
	bblas_zcudaFree(test);
    }

	if ((routine == BBLAS_GEMM)  || (routine == BBLAS_SYMM)  ||
	(routine == BBLAS_HEMM)  || (routine == BBLAS_SYR2K) ||
	(routine == BBLAS_HER2K) || (routine == BBLAS_TRMM)  ||
	(routine == BBLAS_TRSM))
	{
		free(test->Binitnorm);
    }

    /*Free memory allocated for norm computing */
	if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	(routine == BBLAS_HEMM) || (routine == BBLAS_SYRK) ||
	(routine == BBLAS_HERK) || (routine == BBLAS_SYR2K)||
	(routine == BBLAS_HER2K))
    {
		free(test->Cinitnorm);
    }
}

#undef COMPLEX
