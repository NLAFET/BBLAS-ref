/**
 * @file bblas_dprint.c
 *
 * @brief BBLAS printing routines for double routines.
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
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from bblas_zprint.c normal z -> d, Mon Jun  6 09:44:14 2016
 **/
#endif


/** Define bool **/
typedef int bool;
/** Boolean TRUE **/
#define TRUE  1
/** Boolean FALSE **/
#define FALSE 0


#include "bblas_common.h"
#include <sys/utsname.h>
#define REAL

/**
 * Print out system information, information about the current test,
 * the results table header from bblas_printheader, and each line of the results table.
 **/

void bblas_dprint_result(bblas_dtest_t *test)
{
    char test_passed_failed[6];
    int error_in_info = 0;
    static bool  first_time = TRUE;

    if ( (test->batch_count < 0) ||
	 (bblas_minarrayI(test->info, test->batch_count) < 0))
    {
	error_in_info = 1;
    }

    /*Print parameters & header once  */
    if (first_time)
    {
	/*Print system information */
	bblas_print_time();
	bblas_dprint_sys_info(test);
	bblas_print_tolerance(test->routine, test->new_accuracy);
	#ifdef COMPLEX
	printf("Note: Gflop/s is reported in real flops: one real addition is equal to two real additions etc.\n");
	#endif
	/*Print routine name and option */
	bblas_print_title(test->routine, test->batch_opts);
	printf("\n");

	/*Print parameters */
	bblas_dprint_parameters(test);
	/*Print the header for results */
	printf("\n");
	bblas_print_header(test->target);
	first_time = FALSE;
    }

    switch(test->target)
    {

	case BBLAS_MKL:
	    if (error_in_info)
	    {
		test->mkl_perf      = INFINITY;
		test->ref_perf      = INFINITY;
		test->mkl_min_error = INFINITY;
		test->mkl_max_error = INFINITY;
		test->mkl_avg_error = INFINITY;
		test->mkl_std_error = INFINITY;

	    }
	    bblas_dpassed_failed(test, test->mkl_max_error, test_passed_failed, error_in_info);

	    printf("%3d %10d  %14.2f %9s%6.2f)  %10.2f%5s%6.2f)  %15.2e %12.2e %12.2e %12.2e  %15s |\n",
		   test->current_iter, test->batch_count,
		   test->ref_perf,  "(",1000.*test->ref_time,
		   test->mkl_perf, "(",1000.*test->mkl_time,
		   test->mkl_min_error, test->mkl_avg_error,
		   test->mkl_max_error, test->mkl_std_error, test_passed_failed);
	    break;

    case BBLAS_CUBLAS:
    case BBLAS_MAGMA:
	    if (error_in_info)
	    {
		test->ref_perf         = INFINITY;
		test->device_perf      = INFINITY;
		test->device_min_error = INFINITY;
		test->device_max_error = INFINITY;
		test->device_avg_error = INFINITY;
		test->device_std_error = INFINITY;
	    }
	    bblas_dpassed_failed(test, test->device_max_error, test_passed_failed, error_in_info);

	    printf("%3d %10d  %14.2f %9s%6.2f)  %10.2f%8s%6.2f)  %12.2e %12.2e %12.2e %12.2e  %15s |\n",
		   test->current_iter, test->batch_count,
		   test->ref_perf,  "(",1000.*test->ref_time,
		   test->device_perf, "(",1000.*test->device_time,
		   test->device_min_error, test->device_avg_error,
		   test->device_max_error, test->device_std_error, test_passed_failed);
	    break;

	case BBLAS_OTHER:
	    if (error_in_info)
	    {
		test->ref_perf         = INFINITY;
		test->other_perf      = INFINITY;
		test->other_min_error = INFINITY;
		test->other_max_error = INFINITY;
		test->other_avg_error = INFINITY;
		test->other_std_error = INFINITY;

	    }
	    bblas_dpassed_failed(test, test->other_max_error, test_passed_failed, error_in_info);

	    printf("%3d %10d  %14.2f %7s%6.2f)  %10.2f%7s%6.2f)  %15.2e %12.2e %12.2e %12.2e  %15s |\n",
		   test->current_iter, test->batch_count,
		   test->ref_perf,  "(",1000.*test->ref_time,
		   test->other_perf, "(",1000.*test->other_time,
		   test->other_min_error, test->other_avg_error,
		   test->other_max_error, test->other_std_error, test_passed_failed);
	    break;

	case BBLAS_CUMKL:

	    if (error_in_info)
	    {
		test->ref_perf         = INFINITY;
		test->mkl_perf      = INFINITY;
		test->mkl_min_error = INFINITY;
		test->mkl_max_error = INFINITY;
		test->mkl_avg_error = INFINITY;
		test->mkl_std_error = INFINITY;

		test->device_perf      = INFINITY;
		test->device_min_error = INFINITY;
		test->device_max_error = INFINITY;
		test->device_avg_error = INFINITY;
		test->device_std_error = INFINITY;
	    }
	    /*To be implemented */
	    break;
	default:
	    printf("In bblas_dprint_result(): Target not defined\n");
	    exit(EXIT_FAILURE);
    }
}


/**
 * Print out system information depending on what we are comparing against, i.e. MKL or CuBLAS etc.
 **/

void bblas_dprint_sys_info(bblas_dtest_t *test)
{
    bblas_print_platform_info();

    switch(test->target)
	{
        case BBLAS_MKL:
            bblas_print_mkl_info(test->mkl_sequential);
            break;
        case BBLAS_CUBLAS:
	case BBLAS_MAGMA:
            bblas_dprint_cuda_info(test);
            break;
        case BBLAS_CUMKL:
            bblas_print_mkl_info(test->mkl_sequential);
            bblas_dprint_cuda_info(test);
            break;
        case BBLAS_OTHER:
            bblas_print_omp_info();
            break;
        default:
	    printf("In testing_dgemm_batch(): Target no defined\n");
	    exit(EXIT_FAILURE);
	}

}



/**
 * Print out CuBLAS information.
 **/

void bblas_dprint_cuda_info(bblas_dtest_t *test)
{

#if defined(BBLAS_WITH_CUBLAS)
    printf("\nCUDA Hardware Information:\n");
    int ndevices = 0;
    int cublas_v;
    struct cudaDeviceProp prop;
    cudaGetDeviceCount( &ndevices );

    /* Add cuda error check */
    for( int dev = 0; dev < ndevices; dev++ )
    {
         cudaGetDeviceProperties( &prop, dev );
        /* Add cuda error check */
        printf( "   device %d: %s, %.1f MHz clock, %.1f MB memory, capability %d.%d\n",
	        dev,
	        prop.name,
	        prop.clockRate / 1000.,
	        prop.totalGlobalMem / (1024.*1024.),
	        prop.major,
	        prop.minor );
    }

    printf("\nCuBLAS Information:\n");
    /*Add cuda error check */

    cublasGetVersion(test->cublas_handle, &cublas_v);
    printf("    Version: %d.\n", cublas_v);
#endif
}


/**
 * Print out the parameters used for the current test.
 **/

void bblas_dprint_parameters(bblas_dtest_t *test)
{
    /*Local variable */
    enum BBLAS_ROUTINE routine = test->routine;
    char funcname[] = "bblas_dprint_parameters";
    char separator[] = "--";
    char  col[]="|";
    int first_index =0;
	char* transA;
	char* transB;
	char* trans;
	char* uplo;
	char* side;
	char* diag;

    if (test->batch_opts == BBLAS_FIXED )
    {

	switch(routine)
	{
	    case BBLAS_GEMM:
		transA = bblas_op2char(test->transA[first_index]);
		transB = bblas_op2char(test->transB[first_index]);
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%15s|%15s|%10s|%10s|%5s\n","TRANSA","TRANSB","M","N","K");
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%15s|%15s|%10d|%10d|%5d\n",
			   transA, transB,
		       test->M[first_index],test->N[first_index],test->K[first_index]);
		printf("\t\t\t\t============================================================\n");
		free(transA);
		free(transB);
		break;

	    case BBLAS_HEMM:
	    case BBLAS_SYMM:
		uplo = bblas_op2char(test->uplo[first_index]);
		side = bblas_op2char(test->side[first_index]);
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%15s|%15s|%10s|%10s\n","UPLO","SIDE","M","N");
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%15s|%15s|%10d|%10d\n",
			   uplo, side,
		       test->M[first_index],test->N[first_index]);
		printf("\t\t\t\t============================================================\n");
		free(uplo);
		free(side);
		break;

	    case BBLAS_HER2K:
	    case BBLAS_HERK:
	    case BBLAS_SYR2K:
	    case BBLAS_SYRK:
		uplo = bblas_op2char(test->uplo[first_index]);
		trans = bblas_op2char(test->trans[first_index]);
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%15s|%15s|%10s|%10s\n","UPLO","TRANS","N","K");
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%15s|%15s|%10d|%10d\n",
			   uplo, trans,
		       test->N[first_index],test->K[first_index]);
		printf("\t\t\t\t============================================================\n");
		free(uplo);
		free(trans);
		break;

	    case BBLAS_TRMM:
	    case BBLAS_TRSM:
		uplo = bblas_op2char(test->uplo[first_index]);
		transA = bblas_op2char(test->transA[first_index]);
		side = bblas_op2char(test->side[first_index]);
		diag = bblas_op2char(test->diag[first_index]);
		printf("\t\t\t==========================================================================\n");
		printf("\t\t\t%10s|%15s|%15s|%15s|%5s|%5s\n","UPLO","TRANSA","SIDE","DIAG","M", "N");
		printf("\t\t\t==========================================================================\n");
		printf("\t\t\t%10s|%15s|%15s|%15s|%5d|%5d\n",
			   uplo, transA, side, diag,
		       test->M[first_index],test->N[first_index]);
		printf("\t\t\t==========================================================================\n");
		free(uplo);
		free(transA);
		free(side);
		free(diag);
		break;

	    default:
		printf("ERROR in %s, undefined bblas routine name\n",funcname);
		exit(EXIT_FAILURE);
	}
    }else if (test->batch_opts == BBLAS_VARIABLE )
    {
	switch(routine)
	{
	    case BBLAS_GEMM:
		transA = bblas_op2char(test->transA[first_index]);
		transB = bblas_op2char(test->transB[first_index]);
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%13s|%13s|%5s      |%5s      |%5s    \n","TRANSA","TRANSB","M","N","K");
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%13s|%13s|%4d%s%3d%3s%4d%s%3d%3s%4d%s%3d\n",
		       transA, transB,
		       test->minM, separator, test->maxM, col, test->minN,separator, test->maxN, col,test->minK,separator, test->maxK);
		printf("\t\t\t\t============================================================\n");
		free(transA);
		free(transB);
		break;

	    case BBLAS_HEMM:
	    case BBLAS_SYMM:
		uplo = bblas_op2char(test->uplo[first_index]);
		side = bblas_op2char(test->side[first_index]);
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%13s|%13s|%5s      |%5s    \n","UPLO","SIDE","M","N");
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%13s|%13s|%4d%s%3d%3s%4d%s%3d\n",
		       uplo, side,
		       test->minM, separator, test->maxM, col, test->minN,separator, test->maxN);
		printf("\t\t\t\t============================================================\n");
		free(uplo);
		free(side);
		break;

	    case BBLAS_HER2K:
	    case BBLAS_HERK:
	    case BBLAS_SYR2K:
	    case BBLAS_SYRK:
		uplo = bblas_op2char(test->uplo[first_index]);
		trans = bblas_op2char(test->trans[first_index]);
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%13s|%13s|%5s      |%5s    \n","UPLO", "TRANS","N","K");
		printf("\t\t\t\t============================================================\n");
		printf("\t\t\t\t%13s|%13s|%4d%s%3d%3s%4d%s%3d\n",
			   uplo, trans,
		       test->minN, separator, test->maxN, col, test->minK,separator, test->maxK);
		printf("\t\t\t\t============================================================\n");
		free(uplo);
		free(trans);
		break;


	    case BBLAS_TRMM:
	    case BBLAS_TRSM:
		uplo = bblas_op2char(test->uplo[first_index]);
		transA = bblas_op2char(test->transA[first_index]);
		side = bblas_op2char(test->side[first_index]);
		diag = bblas_op2char(test->diag[first_index]);
		printf("\t\t\t\t=============================================================================\n");
		printf("\t\t\t\t%13s|%13s|%13s|%13s|%5s      |%5s    \n","UPLO", "TRANSA", "SIDE","DIAG","M","N");
		printf("\t\t\t\t=============================================================================\n");
		printf("\t\t\t\t%13s|%13s|%13s|%13s|%4d%s%3d%3s%4d%s%3d\n",
			   uplo, transA, side, diag,
		       test->minM, separator, test->maxM, col, test->minN,separator, test->maxN);
		printf("\t\t\t\t=============================================================================\n");
		free(uplo);
		free(transA);
		free(side);
		free(diag);
		break;
	    default:
		printf("ERROR in %s, undefined bblas routine name\n",funcname);
		exit(EXIT_FAILURE);
	}

    }

}
