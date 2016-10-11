/**
 * @file testing_dtrsm_batch.c
 *
 * @brief Testing dtrsm_batch double routine.
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
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from testing_ztrsm_batch.c normal z -> d, Mon Jun  6 09:44:14 2016
 **/
#endif

#define REAL

#include "bblas_common.h"
#include "bblas_mkl.h"
#include "bblas_magma.h"
#include "bblas_cuda.h"
#include "bblas.h"
#include "bblas_omp.h"

/**
 * Makes calls to the reference implementation, CuBLAS, MKL, etc.
 * gets the statistics of the associated errors and prints the results.
 **/

void testing_dtrsm_batch(bblas_dtest_t *test){

#if defined(BBLAS_WITH_CUBLAS)
    /*Creation of cublas handle for CuBLAS calls*/
	if(test->target == BBLAS_CUBLAS || test->target == BBLAS_CUMKL)
		cublasCreate(&test->cublas_handle);
#endif

    for(test->current_iter = 0; test->current_iter < test->nb_test; ++test->current_iter )
    {
	/* set test parameters */
	bblas_dsettest(test);

	/* Copy the arrayB for backward error computation */
	test->Binit  = (double**) malloc(test->batch_count*sizeof(double*));
	bblas_dcopy_Binit(test, test->Binit);


	/* Copy the arrayB in temporay buffer for accuracy testing */
	switch (test->target)
	{
	    case BBLAS_MKL:
		bblas_dcopy_Binit(test, test->mkl_result);
		break;

	    case BBLAS_CUBLAS:
	    case BBLAS_MAGMA:
		bblas_dsettest_cuda(test);
		break;

	    case BBLAS_OTHER:
		bblas_dcopy_Binit(test, test->other_result);
		break;

	    case BBLAS_CUMKL:
		bblas_dsettest_cuda(test);
		break;

	    default:
		printf("In testing_dtrsm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Compute flops */
	bblas_dflops(test);

	/* Error injection */
	if (test->current_iter == test->faulty_iter)
	{
	    bblas_dset_error(test);
	}


	/* =====================================================================
	   Performs operation using SEC
	   =================================================================== */
	test->ref_time = bblas_wtime();

	bblas_dtrsm_batch(test->side, test->uplo, test->transA, test->diag,
			  test->M, test->N, test->alpha, (const double**)test->arrayA,
			  test->lda, test->arrayB, test->ldb,	test->batch_count,
			  test->batch_opts, test->info);

	test->ref_time = bblas_wtime() - test->ref_time;

	/*Compute seq_perf */
	test->ref_perf = (test->flops/1e9)/test->ref_time;

	switch(test->target)
	{
	    case BBLAS_MKL:
		/* =====================================================================
		   Performs operation using  MKL BBLAS
		   =================================================================== */

		test->mkl_time = bblas_wtime();
		bblas_dtrsm_batch(test->side, test->uplo, test->transA, test->diag,
				  test->M, test->N, test->alpha, (const double**)test->arrayA,
				  test->lda, test->mkl_result, test->ldb, test->batch_count,
				  test->batch_opts, test->info);

		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_dtrsm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_dcheck_Bfinal(test, test->mkl_result);
		}

		/*Perform a statistic */
		bblas_dstatistic(test);

		/*Print results*/
		bblas_dprint_result(test);

		break;

	    case BBLAS_CUBLAS:
		/* =====================================================================
		   Performs operation using CUBLAS BBLAS
		   =================================================================== */

		test->device_time = bblas_wtime();

		cublas_dtrsm_batch(test->side, test->uplo,
					  test->transA, test->diag,
					  test->M, test->N,
					  test->alpha,
					  (const double**) test->arrayA_device, test->lda,
					  test->arrayB_device, test->ldb,
					  test->batch_count, test->batch_opts,
                      test->info,
                      test);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Copy B matrix from device to host*/
		bblas_dget_Bfinal(test);

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_dtrsm_tolerance(test, test->device_result);
		} else
		{
		    bblas_dcheck_Bfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_dstatistic(test);

		/*Print results*/
		bblas_dprint_result(test);

		break;


	    case BBLAS_MAGMA:
		/* =====================================================================
		   Performs operation using MAGMA BBLAS
		   =================================================================== */

		test->device_time = bblas_wtime();

		magma_dtrsm_batch(test->side, test->uplo,
					  test->transA, test->diag,
					  test->M, test->N,
					  test->alpha,
					  (const double**) test->arrayA_device, test->lda,
					  test->arrayB_device, test->ldb,
					  test->batch_count, test->batch_opts,
				  test->info);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Copy B matrix from device to host*/
		bblas_dget_Bfinal(test);

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_dtrsm_tolerance(test, test->device_result);
		} else
		{
		    bblas_dcheck_Bfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_dstatistic(test);

		/*Print results*/
		bblas_dprint_result(test);

		break;

	    case BBLAS_CUMKL:

		/*******************************
		 *        [1] CUBLAS CASE      *
		 ********************************/

		/*Change target to particular case of BBLAS_CUBLAS */
		test->target = BBLAS_CUBLAS;

		test->device_time = bblas_wtime();

        cublas_dtrsm_batch(test->side, test->uplo,
					  test->transA, test->diag,
					  test->M, test->N,
					  test->alpha,
					  (const double**) test->arrayA_device, test->lda,
					  test->arrayB_device, test->ldb,
					  test->batch_count, test->batch_opts,
                      test->info,
                      test);

	    test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Copy B matrix from device to host*/
		bblas_dget_Bfinal(test);

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_dtrsm_tolerance(test, test->device_result);
		} else
		{
		    bblas_dcheck_Bfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_dstatistic(test);

		/*******************************
		 *        [2] MKL CASE      *
		 ********************************/

		/*Change target to particular case of BBLAS_MKL */
		test->target = BBLAS_MKL;

		test->mkl_time = bblas_wtime();
		bblas_dtrsm_batch(test->side, test->uplo, test->transA, test->diag,
				  test->M, test->N, test->alpha, (const double**)test->arrayA,
				  test->lda, test->mkl_result, test->ldb, test->batch_count,
				  test->batch_opts, test->info);

		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_dtrsm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_dcheck_Bfinal(test, test->mkl_result);
		}


		/*Perform a statistic */
		bblas_dstatistic(test);

		/*******************************
		 *[3] BOTH CUBLAS & MKL CASE   *
		 ********************************/

		/*Change target to BBLAS_CUMKL */
		test->target = BBLAS_CUMKL;

		/*Print results*/
		bblas_dprint_result(test);
		break;

	    case BBLAS_OTHER:
		/* =====================================================================
		   Performs operation using OMP CPU BLAS
		   =================================================================== */

		test->other_time = bblas_wtime();
		omp_dtrsm_batch(test->side, test->uplo, test->transA, test->diag,
				  test->M, test->N, test->alpha, (const double**)test->arrayA,
				  test->lda, test->other_result, test->ldb,	test->batch_count,
				  test->batch_opts, test->info);

		test->other_time = bblas_wtime() - test->other_time;


		/*Compute perf */
		test->other_perf = (test->flops/1e9)/test->other_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_dtrsm_tolerance(test, test->other_result);
		} else
		{
		    bblas_dcheck_Bfinal(test, test->other_result);
		}

		/*Perform a statistic */
		bblas_dstatistic(test);

		/*Print results*/
		bblas_dprint_result(test);
		break;

	    default:
		printf("In testing_dtrsm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Free batch configuration */
	for (int i = 0; i < test->batch_count; i++)
	{
		free(test->Binit[i]);
	}
	free(test->Binit);
	bblas_dfreetest(test);
	printf("\n");
    }

#if defined(BBLAS_WITH_CUBLAS)
    /*Destroy cublas_handle*/
	if(test->target == BBLAS_CUBLAS || test->target == BBLAS_CUMKL)
        cublasDestroy(test->cublas_handle);
#endif

}


#undef REAL
