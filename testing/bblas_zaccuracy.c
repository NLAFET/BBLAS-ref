/**
 * @file bblas_zaccuracy.c
 *
 * @brief BBLAS accuracy tests for double_Complex precision.
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
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#include "bblas_common.h"
#include <lapacke.h>
#include <cblas.h>



#define COMPLEX

/**
 * bblas_zgemm_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|B\|_{\infty} +
 *           |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_zgemm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{

    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, Bm, Bn, M, K, N, ldc;
    int max_work_size = max(test->maxK,max(test->maxM, test->maxN));
    double Anorm, Bnorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    BBLAS_Complex64_t alpha =-1;
    double *work = (double *)malloc(max_work_size*sizeof(double));

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
	    M = test->M[batch_iter];
	    K   = test->K[batch_iter];
	    N   = test->N[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if (test->transA[batch_iter] == BblasNoTrans) {
		Am = M; An = K;
	    } else {
		Am = K; An = M;
	    }
	    if (test->transB[batch_iter] == BblasNoTrans) {
		Bm = K; Bn = N;
	    } else {
		Bm = N; Bn = K;
	    }

	    /*Compute the error C - C_finial */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    result = Rnorm / ((K*cabs(test->alpha[batch_iter])*Anorm*Bnorm + cabs(test->beta[batch_iter])*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	M   = test->M[first_index];
	K   = test->K[first_index];
	ldc = test->ldc[first_index];
	N   = test->N[first_index];

	if (test->transA[first_index] == BblasNoTrans) {
	    Am = M; An = K;
	} else {
	    Am = K; An = M;
	}
	if (test->transB[first_index] == BblasNoTrans) {
	    Bm = K; Bn = N;
	} else {
	    Bm = N; Bn = K;
	}

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
		/*Compute the error C - C_finial */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);
	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', M, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    /*Compute the relative error */
	    result = Rnorm / ((K*cabs(test->alpha[first_index])*Anorm*Bnorm + cabs(test->beta[first_index])*test->Cinitnorm[batch_iter])*eps);

	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
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

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_zstarmm_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|B\|_{\infty} +
 * |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_zstarmm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, Bm, Bn, N, ldc;
    int max_work_size = max(test->maxM, test->maxN);
    double Anorm, Bnorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    BBLAS_Complex64_t alpha =-1;
    double *work = (double *)malloc(max_work_size*sizeof(double));

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
	    N = test->N[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if(test->side[batch_iter] == BblasLeft )
	    {
		Am = test->M[batch_iter];
	    } else {
		Am = test->N[batch_iter];
	    }
	    Bm = test->ldb[batch_iter];
	    Bn = test->N[batch_iter];
	    /*Compute the error C - C_finial */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, Am, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    result = Rnorm / ((N*cabs(test->alpha[batch_iter])*Anorm*Bnorm + cabs(test->beta[batch_iter])*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	N = test->N[first_index];
	ldc = test->ldc[first_index];

	if(test->side[first_index] == BblasLeft )
	{
	    Am = test->M[first_index];
	} else {
	    Am = test->N[first_index];
	}
	Bm = test->ldb[first_index];
	Bn = test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    /*Compute the error C - C_finial */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, Am, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    /*Compute the relative error */
	    result = Rnorm / ((N*cabs(test->alpha[first_index])*Anorm*Bnorm + cabs(test->beta[first_index])*test->Cinitnorm[batch_iter])*eps);

	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
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

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_zstark_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|A\|_1 +
 * |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_zstark_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, N, K, ldc;
    char u;
    int max_work_size = max(test->maxK, test->maxN);
    double Anorm, ATnorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    double alpha_norm, beta_norm;
    BBLAS_Complex64_t alpha =-1;
    double *work = (double *)malloc(max_work_size*sizeof(double));

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
	    N = test->N[batch_iter];
	    K = test->K[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if (test->trans[batch_iter] == BblasNoTrans) {
		Am = N; An = K;
	    } else {
		Am = K; An = N;
	    }

	    if(test->uplo[batch_iter] == BblasLower )
	      {
		u = 'L';
	      }else
	      {
		u = 'U';
	      }
	    /*Compute the error C - C_finial */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR, 'I', u, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the one norm of A */
	    ATnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

	    if (test->routine == BBLAS_HERK)
	    {
		beta_norm   = cabs(test->beta_herk[batch_iter]);
		alpha_norm  = cabs(test->alpha_herk[batch_iter]);
	    }else{
		beta_norm  = cabs(test->beta[batch_iter]);
		alpha_norm  = cabs(test->alpha[batch_iter]);
	    }

	    result = Rnorm / ((K*alpha_norm*Anorm*ATnorm + beta_norm*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	    N = test->N[first_index];
	    K = test->K[first_index];
	    ldc = test->ldc[first_index];

	    if (test->trans[first_index] == BblasNoTrans) {
		Am = N; An = K;
	    } else {
		Am = K; An = N;
	    }
	    if(test->uplo[first_index] == BblasLower )
	      {
		u = 'L';
	      }else
	      {
		u = 'U';
	      }

	    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	      {
		/*Compute the error C - C_finial */
		cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		Rnorm = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR, 'I', u, N, C_diff[batch_iter], ldc, work);

		/*Compute the infinity norm of A */
		Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

		/*Compute the one norm of A */
		ATnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

		if (test->routine == BBLAS_HERK)
		{
		    beta_norm   = cabs(test->beta_herk[first_index]);
		    alpha_norm  = cabs(test->alpha_herk[first_index]);
		}else{
		    beta_norm   = cabs(test->beta[first_index]);
		    alpha_norm  = cabs(test->alpha[first_index]);
		}

		/*Compute the relative error */
		result = Rnorm / ((K*alpha_norm*Anorm*ATnorm +
				   beta_norm*test->Cinitnorm[batch_iter])*eps);

		switch(test->target)
		{
		    case BBLAS_MKL:
			test->mkl_error[batch_iter] = result;
			break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		  test->device_error[batch_iter] = result;
		  break;

		    case BBLAS_OTHER:
			test->other_error[batch_iter] = result;
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

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_zstar2k_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|B\|_1 +
 * K|\alpha|\|B\|_{\infty}\|A\|_1  +
 * |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_zstar2k_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, Bn, Bm, N, K, ldc;
    int max_work_size = max(test->maxK, test->maxN);
    double Anorm, Bnorm, ATnorm, BTnorm,  Rnorm, result;
    double beta_norm;
    double eps = LAPACKE_dlamch_work('e');
    BBLAS_Complex64_t alpha =-1;
    double *work = (double *)malloc(max_work_size*sizeof(double));

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
	    N = test->N[batch_iter];
	    K = test->K[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if (test->trans[batch_iter] == BblasNoTrans) {
		Am = N; An = K;
		Bn = K;
	    } else {
		Am = K; An = N;
		Bn = N;
	    }
	    Bm = test->ldb[batch_iter];

	    /*Compute the error C - C_finial */
	    cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the one norm of A */
	    ATnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    /*Compute the one norm of B */
	    BTnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'O', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    if(test->routine == BBLAS_HER2K)
	    {
		beta_norm = cabs(test->beta_herk[batch_iter]);
	    }else
	    {
		beta_norm = cabs(test->beta[batch_iter]);
	    }

	    result = Rnorm / ((K*cabs(test->alpha[batch_iter])*Anorm*BTnorm +
			       K*cabs(test->alpha[batch_iter])*Bnorm*ATnorm +
			       beta_norm*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	N = test->N[first_index];
	K = test->K[first_index];
	ldc = test->ldc[first_index];

	if (test->trans[first_index] == BblasNoTrans) {
		Am = N; An = K;
		Bn = K;
	} else {
	    Am = K; An = N;
	    Bn = N;
	}
	Bm = test->ldb[first_index];


	    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	    {
		/*Compute the error C - C_finial */
		cblas_zaxpy (ldc*N, CBLAS_SADDR(alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
		Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

		/*Compute the infinity norm of A */
		Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

		/*Compute the one norm of A */
		ATnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

		/*Compute the infinity norm of B */
		Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

		/*Compute the one norm of B */
		BTnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'O', Bm, Bn, test->arrayB[batch_iter], Bm, work);

		if(test->routine == BBLAS_HER2K)
		{
		    beta_norm = cabs(test->beta_herk[first_index]);
		}else
		{
		    beta_norm = cabs(test->beta[first_index]);
		}

		/*Compute the relative error */
		result = Rnorm / ((K*cabs(test->alpha[first_index])*Anorm*BTnorm +
				   K*cabs(test->alpha[first_index])*Bnorm*ATnorm +
				   beta_norm*test->Cinitnorm[batch_iter])*eps);


		switch(test->target)
		{
		case BBLAS_MKL:
		  test->mkl_error[batch_iter] = result;
		  break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		  test->device_error[batch_iter] = result;
		  break;

		case BBLAS_OTHER:
		  test->other_error[batch_iter] = result;
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

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_ztrmm_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((|\alpha|\|A\|_{\infty}\|B\|_{\infty}N\epsilon).@f]
 *
 **/

void bblas_ztrmm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **B_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, ldb, N, M;
    int max_work_size = max(test->maxK, test->maxN);
    char u, diag;
    double Anorm, Rnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    BBLAS_Complex64_t alpha =-1;
    double *work = (double *)malloc(max_work_size*sizeof(double));



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
	    N = test->N[batch_iter];
	    M = test->M[batch_iter];
	    ldb = test->ldb[batch_iter];

	    Am = test->lda[batch_iter];

	    if(test->side[batch_iter] == BblasLeft )
	    {
		An = M;
	    }else
	    {
		An = N;
	    }

	    if(test->uplo[batch_iter] == BblasLower )
	    {
	       u = 'L';
	    }else
	    {
	       u = 'U';
	    }

	    if(test->diag[batch_iter] == BblasUnit )
	    {
		diag = 'U';
	    }else
	    {
		diag = 'N';
	    }

	    /*Compute the error B - B_finial */
	    cblas_zaxpy (ldb*N, CBLAS_SADDR(alpha), B_final[batch_iter], 1, B_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldb, N, B_diff[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);
	    
		result = Rnorm/(N*cabs(test->alpha[batch_iter])*Anorm*test->Binitnorm[batch_iter]*eps);
	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }

	}
    }else  if( test->batch_opts == BBLAS_FIXED )
    {

	N = test->N[first_index];
	M = test->M[first_index];
	ldb = test->ldb[first_index];

	Am = test->lda[first_index];

	if(test->side[first_index] == BblasLeft )
	{
	    An = M;
	}else
	{
	    An = N;
	}

	if(test->uplo[first_index] == BblasLower )
	{
	    u = 'L';
	}else
	{
	    u = 'U';
	}

	if(test->diag[first_index] == BblasUnit )
	{
	    diag = 'U';
	}else
	{
	    diag = 'N';
	}

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    /*Compute the error B - B_finial */
	    cblas_zaxpy (ldb*N, CBLAS_SADDR(alpha), B_final[batch_iter], 1, B_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldb, N, B_diff[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);

	    result = Rnorm/(N*cabs(test->alpha[first_index])*Anorm*test->Binitnorm[batch_iter]*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }
    else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value");
    }

    /*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(B_diff[index]);
    }
    free(B_diff);
    free(work);
}

/**
 * bblas_ztrsm_tolerance computes the backward error
 * @f[\|\alpha B - AX \|_{\infty}/((\|A\|_{\infty}\|X\| + |\alpha|\|B\|_{\infty})N\epsilon).@f]
 *
 **/

void bblas_ztrsm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **B_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, ldb, N, M;
    int max_work_size = max(test->maxK, test->maxN);
    char u, diag;
    double Anorm, Rnorm, result, Xnorm;
    double eps = LAPACKE_dlamch_work('e');
    BBLAS_Complex64_t alpha;
    double *work = (double *)malloc(max_work_size*sizeof(double));



    /*Variable to save the residual alphaB -AX */
    BBLAS_Complex64_t **residual;
    residual = (BBLAS_Complex64_t **) malloc(batch_count*sizeof(BBLAS_Complex64_t *));

    /*Make a copy of test->arrayB (X) in residual */
    bblas_zcopy_Binit(test, residual);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    N = test->N[batch_iter];
	    M = test->M[batch_iter];
	    ldb = test->ldb[batch_iter];

	    Am = test->lda[batch_iter];

	    if(test->side[batch_iter] == BblasLeft )
	    {
		An = M;
	    }else
	    {
		An = N;
	    }

	    if(test->uplo[batch_iter] == BblasLower )
	    {
	       u = 'L';
	    }else
	    {
	       u = 'U';
	    }

	    if(test->diag[batch_iter] == BblasUnit )
	    {
		diag = 'U';
	    }else
	    {
		diag = 'N';
	    }

	    /*Compute residual <- AX  (ztrmm) */
	    alpha = 1;
	    cblas_ztrmm(BblasColMajor, test->side[batch_iter], test->uplo[batch_iter], test->transA[batch_iter],
			test->diag[batch_iter], M, N, CBLAS_SADDR(alpha), test->arrayA[batch_iter],
			test->lda[batch_iter], residual[batch_iter], ldb);

	    /*Compute Binit = alpha Binit */
	    for (int index=0;  index < ldb*N; index++)
	    {
		test->Binit[batch_iter][index] = test->alpha[batch_iter]*test->Binit[batch_iter][index];
	    }

	    /*Compute the  residual <- alpha B - Ax  */
	    alpha= -1;
	    cblas_zaxpy (ldb*N, CBLAS_SADDR(alpha), test->Binit[batch_iter], 1, residual[batch_iter], 1);
		/*Compute  the infinity norm associated  with the residual */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldb, N, residual[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);

	    /*Compute  the infinity norm associated  with X (B) */
	    Xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldb, N, test->arrayB[batch_iter], ldb, work);

	    result = Rnorm/(Anorm*Xnorm*N*eps);
	    
	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }

	}
    }else  if( test->batch_opts == BBLAS_FIXED )
    {

	N = test->N[first_index];
	M = test->M[first_index];
	ldb = test->ldb[first_index];

	Am = test->lda[first_index];

	if(test->side[first_index] == BblasLeft )
	{
	    An = M;
	}else
	{
	    An = N;
	}

	if(test->uplo[first_index] == BblasLower )
	{
	    u = 'L';
	}else
	{
	    u = 'U';
	}

	if(test->diag[first_index] == BblasUnit )
	{
	    diag = 'U';
	}else
	{
	    diag = 'N';
	}

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{

	    /*Compute residual <- AX  (ztrmm) */
	    alpha = 1;
	    cblas_ztrmm(BblasColMajor, test->side[first_index], test->uplo[first_index], test->transA[first_index],
			test->diag[first_index], M, N, CBLAS_SADDR(alpha), (const BBLAS_Complex64_t*) test->arrayA[batch_iter],
			test->lda[first_index], residual[batch_iter], ldb);

	    /*Compute Binit = alpha Binit */
	    for (int index=0;  index < ldb*N; index++)
	    {
		test->Binit[batch_iter][index] = test->alpha[first_index]*test->Binit[batch_iter][index];
	    }

	    /*Compute the  residual <- alpha B - Ax  */
	    alpha = -1;
	    cblas_zaxpy (ldb*N, CBLAS_SADDR(alpha), test->Binit[batch_iter], 1, residual[batch_iter], 1);

		/*Compute  the infinity norm associated  with the residual */
	    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldb, N, residual[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);

	    /*Compute  the infinity norm associated  with X (B) */
	    Xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'I', ldb, N, test->arrayB[batch_iter], ldb, work);

	    result = Rnorm/(Anorm*Xnorm*N*eps);
	    
	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_zcheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }
    else
    {
	bblas_error("bblas_ztesting.c", "wrong batch_opts value");
    }

    /*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(residual[index]);
    }
    free(residual);
    free(work);
}


/**
 * Calls bblas_zstar2k_tolerance.
 **/
#ifdef COMPLEX
void bblas_zher2k_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    bblas_zstar2k_tolerance(test, C_final);

}
#endif
/**
 * Calls bblas_zstar2k_tolerance.
 **/

void bblas_zsyr2k_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    bblas_zstar2k_tolerance(test, C_final);
}

/**
 * Calls bblas_zstark_tolerance.
 **/

#ifdef COMPLEX
void bblas_zherk_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    bblas_zstark_tolerance(test, C_final);
}
#endif
/**
 * Calls bblas_zstark_tolerance.
 **/

void bblas_zsyrk_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    bblas_zstark_tolerance(test, C_final);
}

/**
 * Calls bblas_zstarmm_tolerance.
 **/

#ifdef COMPLEX
void bblas_zhemm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    bblas_zstarmm_tolerance(test, C_final);
}
#endif
/**
 * Calls bblas_zstamm_tolerance.
 **/

void bblas_zsymm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final)
{
    bblas_zstarmm_tolerance(test, C_final);
}
