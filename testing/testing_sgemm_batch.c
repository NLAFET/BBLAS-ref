/**
 * @file testing_sgemm_batch.c
 *
 * @brief Testing sgemm_batch float routine.
 *
 * BBLAS is a software package provided by Univ. of Manschester,
 * Univ. of Tennessee.
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
 * @generated from testing_zgemm_batch.c normal z -> s, Mon Jun  6 09:44:14 2016
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

void testing_sgemm_batch(bblas_stest_t *test){

#if defined(BBLAS_WITH_CUBLAS)
    /*Creation of cublas handle for CuBLAS calls*/
	if(test->target == BBLAS_CUBLAS || test->target == BBLAS_CUMKL)
    cublasCreate(&test->cublas_handle);
#endif

    for(test->current_iter = 0; test->current_iter < test->nb_test; ++test->current_iter )
    {
	/* set test parameters */
	bblas_ssettest(test);


	/* Copy the arrayC in temporay buffer for accuracy testing */
	switch (test->target)
	{
	case BBLAS_MKL:
	  bblas_scopy_Cinit(test, test->mkl_result);
	  break;

	case BBLAS_CUBLAS:
	case BBLAS_MAGMA:
	  bblas_ssettest_cuda(test);
	  break;

	case BBLAS_OTHER:
	  bblas_scopy_Cinit(test, test->other_result);
	  break;

	case BBLAS_CUMKL:
	  bblas_ssettest_cuda(test);
	  bblas_scopy_Cinit(test, test->mkl_result);
	  break;

	default:
	  printf("In testing_sgemm_batch(): Target no defined\n");
	  exit(EXIT_FAILURE);
	}

	/*Compute flops */
	bblas_sflops(test);

	/* Error injection */
	if (test->current_iter == test->faulty_iter)
	{
	    bblas_sset_error(test);
	}

	/* =====================================================================
	   Performs operation using SEQ
	   =================================================================== */
	test->ref_time = bblas_wtime();

	bblas_sgemm_batch(test->transA, test->transB, test->M, test->N, test->K,
			  test->alpha, (const float**)test->arrayA,
			  test->lda, (const float**)test->arrayB,
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
			mkl_sgemm_batch(test->transA, test->transB,
					 test->M,  test->N, test->K,
					 test->alpha,
					 (const float**) test->arrayA, test->lda,
					 (const float**) test->arrayB, test->ldb,
					 test->beta,
					 test->mkl_result, test->ldc,
					 test->batch_count,	test->batch_opts,
					 test->info);
			test->mkl_time = bblas_wtime() - test->mkl_time;
			/*Compute perf */
			test->mkl_perf = (test->flops/1e9)/test->mkl_time;

			/*Compute relative errors */
			if(test->new_accuracy)
			{
			    bblas_sgemm_tolerance(test, test->mkl_result);
			} else
			{
			    bblas_scheck_Cfinal(test, test->mkl_result);
			}

			/*Perform a statistic */
			bblas_sstatistic(test);

			/*Print results*/
			bblas_sprint_result(test);

		break;

	    case BBLAS_CUBLAS:
		/* =====================================================================
		   Performs operation using CUBLAS BBLAS
		   =================================================================== */

	      test->device_time = bblas_wtime();

	      cublas_sgemm_batch(test->transA, test->transB,
				 test->M, test->N, test->K,
				 test->alpha,
				 (const float**) test->arrayA_device, test->lda,
				 (const float**) test->arrayB_device, test->ldb,
				 test->beta,
				 test->arrayC_device, test->ldc,
				 test->batch_count, test->batch_opts,
				 test->info,
				 test);

	      test->device_time = bblas_wtime() - test->device_time;

	      /*Compute perf */
	      test->device_perf = (test->flops/1e9)/test->device_time;

	      /*Copy C matrix from device to host*/
	      bblas_sget_Cfinal(test);

		/*Compute relative errors */
	      if(test->new_accuracy)
		{
		    bblas_sgemm_tolerance(test, test->device_result);
		}
	      else
		{
		  bblas_scheck_Cfinal(test, test->device_result);
		}

	      /*Perform a statistic */
	      bblas_sstatistic(test);

	      /*Print results*/
	      bblas_sprint_result(test);

	      break;

	    case BBLAS_MAGMA:
		/* =====================================================================
		   Performs operation using MAGMA BBLAS
		   =================================================================== */

	      test->device_time = bblas_wtime();

	      magma_sgemm_batch(test->transA, test->transB,
				test->M, test->N, test->K,
				test->alpha,
				(const float**) test->arrayA_device, test->lda,
				(const float**) test->arrayB_device, test->ldb,
				test->beta,
				test->arrayC_device, test->ldc,
				test->batch_count, test->batch_opts,
				test->info);

	      test->device_time = bblas_wtime() - test->device_time;

	      /*Compute perf */
	      test->device_perf = (test->flops/1e9)/test->device_time;

	      /*Copy C matrix from device to host*/
	      bblas_sget_Cfinal(test);

		/*Compute relative errors */
	      if(test->new_accuracy)
		{
		    bblas_sgemm_tolerance(test, test->device_result);
		}
	      else
		{
		  bblas_scheck_Cfinal(test, test->device_result);
		}

	      /*Perform a statistic */
	      bblas_sstatistic(test);

	      /*Print results*/
	      bblas_sprint_result(test);

	      break;

	case BBLAS_CUMKL:

		/*******************************
		 *        [1] CUBLAS CASE      *
		 ********************************/

		/*Change target to particular case of BBLAS_CUBLAS */
		test->target = BBLAS_CUBLAS;

		test->device_time = bblas_wtime();

        cublas_sgemm_batch(test->transA, test->transB,
					  test->M, test->N, test->K,
					  test->alpha,
					  (const float**) test->arrayA_device, test->lda,
					  (const float**) test->arrayB_device, test->ldb,
					  test->beta,
					  test->arrayC_device, test->ldc,
					  test->batch_count, test->batch_opts,
                      test->info,
                      test);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Copy C matrix from device to host*/
		bblas_sget_Cfinal(test);

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_sgemm_tolerance(test, test->device_result);
		} else
		{
		    bblas_scheck_Cfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*Print results*/
		bblas_sprint_result(test);

		/*******************************
		 *        [2] MKL CASE      *
		 ********************************/

		/*Change target to particular case of BBLAS_MKL */
		test->target = BBLAS_MKL;

		test->mkl_time = bblas_wtime();
		mkl_sgemm_batch(test->transA, test->transB,
					 test->M,  test->N, test->K,
					 test->alpha,
					 (const float**) test->arrayA, test->lda,
					 (const float**) test->arrayB, test->ldb,
					 test->beta,
					 test->mkl_result, test->ldc,
					 test->batch_count,	test->batch_opts,
					 test->info);
		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_sgemm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_scheck_Cfinal(test, test->mkl_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*******************************
		 *[3] BOTH CUBLAS & MKL CASE   *
		 ********************************/

		/*Change target to BBLAS_CUMKL */
		test->target = BBLAS_CUMKL;

		/*Print results*/
		bblas_sprint_result(test);
		break;

	    case BBLAS_OTHER:
		/* =====================================================================
		   Performs operation using OMP CPU BLAS
		   =================================================================== */

		test->other_time = bblas_wtime();

		omp_sgemm_batch(test->transA, test->transB, test->M, test->N, test->K,
				test->alpha, (const float**)test->arrayA,
				test->lda, (const float**)test->arrayB,
				test->ldb, test->beta, test->other_result, test->ldc,
				test->batch_count, test->batch_opts, test->info);

		test->other_time = bblas_wtime() - test->other_time;

		/*Compute perf */
		test->other_perf = (test->flops/1e9)/test->other_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_sgemm_tolerance(test, test->other_result);
		} else
		{
		    bblas_scheck_Cfinal(test, test->other_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*Print results*/
		bblas_sprint_result(test);
		break;

	    default:
		printf("In testing_sgemm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Free batch configuration */
	bblas_sfreetest(test);
	printf("\n");
    }

#if defined(BBLAS_WITH_CUBLAS)
    /*Destroy cublas_handle*/
	if(test->target == BBLAS_CUBLAS || test->target == BBLAS_CUMKL)
        cublasDestroy(test->cublas_handle);
#endif


}


#undef REAL
