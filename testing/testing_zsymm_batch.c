/**
 * @file testing_zsymm_batch.c
 *
 * @brief Testing zsymm_batch double _Complex routine.
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
 * @precisions normal z -> c d s
 **/
#endif

#define COMPLEX

#include "bblas_common.h"
#include "bblas_mkl.h"
#include "bblas.h"
#include "bblas_omp.h"

/**
 * Makes calls to the reference implementation, CuBLAS, MKL, etc.
 * gets the statistics of the associated errors and prints the results.
 **/

void testing_zsymm_batch(bblas_ztest_t *test){


    for(test->current_iter = 0; test->current_iter < test->nb_test; ++test->current_iter )
    {
	/* set test parameters */
	bblas_zsettest(test);

	/* Copy the arrayC in temporay buffer for accuracy testing */
	switch (test->target)
	{
	    case BBLAS_MKL:
		bblas_zcopy_Cinit(test, test->mkl_result);
		break;

	    case BBLAS_CUBLAS:
		bblas_zcopy_Cinit(test, test->device_result);
		break;

	    case BBLAS_OTHER:
		bblas_zcopy_Cinit(test, test->other_result);
		break;

	    case BBLAS_CUMKL:
		bblas_zcopy_Cinit(test, test->mkl_result);
		bblas_zcopy_Cinit(test, test->device_result);
		break;

	    default:
		printf("In testing_zsymm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Compute flops */
	bblas_zflops(test);

	/* Error injection */
	if (test->current_iter == test->faulty_iter)
	{
	    bblas_zset_error(test);
	}

	/* =====================================================================
	   Performs operation using SEC
	   =================================================================== */
	test->ref_time = bblas_wtime();

	bblas_zsymm_batch(test->side, test->uplo, test->M, test->N,
		    test->alpha, (const BBLAS_Complex64_t**)test->arrayA,
		    test->lda, (const BBLAS_Complex64_t**)test->arrayB,
		    test->ldb, test->beta, test->arrayC, test->ldc,
		    test->batch_count, test->batch_opts, test->info);

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

		bblas_zsymm_batch(test->side, test->uplo, test->M, test->N,
			    test->alpha, (const BBLAS_Complex64_t**)test->arrayA,
			    test->lda, (const BBLAS_Complex64_t**)test->arrayB,
			    test->ldb, test->beta, test->mkl_result, test->ldc,
			    test->batch_count, test->batch_opts, test->info);

		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_zsymm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_zcheck_Cfinal(test, test->mkl_result);
		}

		/*Perform a statistic */
		bblas_zstatistic(test);

		/*Print results*/
		bblas_zprint_result(test);

		break;

	    case BBLAS_CUBLAS:
		/* =====================================================================
		   Performs operation using CUBLAS BBLAS
		   =================================================================== */

		test->device_time = bblas_wtime();

		bblas_zsymm_batch(test->side, test->uplo, test->M, test->N,
			    test->alpha, (const BBLAS_Complex64_t**)test->arrayA,
			    test->lda, (const BBLAS_Complex64_t**)test->arrayB,
			    test->ldb, test->beta, test->device_result, test->ldc,
			    test->batch_count, test->batch_opts, test->info);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_zsymm_tolerance(test, test->device_result);
		} else
		{
		    bblas_zcheck_Cfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_zstatistic(test);

		/*Print results*/
		bblas_zprint_result(test);

		break;

	    case BBLAS_CUMKL:

		/*******************************
		 *        [1] CUBLAS CASE      *
		********************************/

		/*Change target to particular case of BBLAS_CUBLAS */
		test->target = BBLAS_CUBLAS;

		test->device_time = bblas_wtime();

		bblas_zsymm_batch(test->side, test->uplo, test->M, test->N,
			    test->alpha, (const BBLAS_Complex64_t**)test->arrayA,
			    test->lda, (const BBLAS_Complex64_t**)test->arrayB,
			    test->ldb, test->beta, test->device_result, test->ldc,
			    test->batch_count, test->batch_opts, test->info);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_zsymm_tolerance(test, test->device_result);
		} else
		{
		    bblas_zcheck_Cfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_zstatistic(test);

		/*******************************
		 *        [2] MKL CASE      *
		********************************/

		/*Change target to particular case of BBLAS_MKL */
		test->target = BBLAS_MKL;

		test->mkl_time = bblas_wtime();

		bblas_zsymm_batch(test->side, test->uplo, test->M, test->N,
			    test->alpha, (const BBLAS_Complex64_t**)test->arrayA,
			    test->lda, (const BBLAS_Complex64_t**)test->arrayB,
			    test->ldb, test->beta, test->mkl_result, test->ldc,
			    test->batch_count, test->batch_opts, test->info);

		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_zsymm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_zcheck_Cfinal(test, test->mkl_result);
		}

		/*Perform a statistic */
		bblas_zstatistic(test);

		/*******************************
		 *[3] BOTH CUBLAS & MKL CASE   *
		********************************/

		/*Change target to BBLAS_CUMKL */
		test->target = BBLAS_CUMKL;

		/*Print results*/
		bblas_zprint_result(test);
		break;

	    case BBLAS_OTHER:
		/* =====================================================================
		   Performs operation using OMP CPU BLAS
		   =================================================================== */

		test->other_time = bblas_wtime();

		omp_zsymm_batch(test->side, test->uplo, test->M, test->N,
				test->alpha, (const BBLAS_Complex64_t**)test->arrayA,
				test->lda, (const BBLAS_Complex64_t**)test->arrayB,
				test->ldb, test->beta, test->other_result, test->ldc,
				test->batch_count, test->batch_opts, test->info);

		test->other_time = bblas_wtime() - test->other_time;


		/*Compute perf */
		test->other_perf = (test->flops/1e9)/test->other_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_zsymm_tolerance(test, test->other_result);
		} else
		{
		    bblas_zcheck_Cfinal(test, test->other_result);
		}

		/*Perform a statistic */
		bblas_zstatistic(test);

		/*Print results*/
		bblas_zprint_result(test);
		break;

	    default:
		printf("In testing_zsymm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Free batch configuration */

	bblas_zfreetest(test);
	printf("\n");
    }
}


#undef COMPLEX
