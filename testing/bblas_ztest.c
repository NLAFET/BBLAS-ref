/**
 * @file bblas_ztest.c
 *
 * @brief Runs tests for double _Complex routines.
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

#if defined(BBLAS_WITH_MKL)
#include <mkl.h>
#endif

/**
 * Runs the tests, calling all other necessary functions.
 **/

int main(int argc, char **argv){

    /* Local variables */
    bblas_ztest_t *test       = ( bblas_ztest_t *) malloc (sizeof(bblas_ztest_t));

    /* Check call  */
    if ( argc != 2 )
    {
	printf( "usage: %s parameters_file\n", argv[0] );
	return 0;
    }

    /* Open the file */
	char* filename = argv[1];
    FILE * file = fopen( filename, "r" );

    /* fopen returns 0, the NULL pointer, on failure */
    if ( file == 0 )
    {
	printf( "Could not open file\n" );
	return 0;
    }

    /* initialize parameters */
    bblas_zinit_config (test);

    /*Read parameters  */
    parse_zconfigfile(file, test);
    /* Close the file */
    fclose( file );

	/* Check that compiling with MKL if using target == BBLAS_MKL */
#if !defined(BBLAS_WITH_MKL)
	if (test->target == BBLAS_MKL)
	{
		printf("ERROR: Cannot compare against MKL when MKL is not linked: change make.inc or the input file to continue.\n");
		return 0;
	}
#endif

    /* Check that compiling with CuBLAS if using target == BBLAS_CUBLAS */
#if !defined(BBLAS_WITH_CUBLAS)
	if (test->target == BBLAS_CUBLAS)
	{
		printf("ERROR: Cannot compare against CuBLAS when CuBLAS is not linked: change make.inc or the input file to continue.\n");
		return 0;
	}
#endif

    /* Check that compiling with MAGMA if using target == BBLAS_MAGMA */
#if !defined(BBLAS_WITH_MAGMA)
	if (test->target == BBLAS_MAGMA)
	{
		printf("ERROR: Cannot compare against MAGMA when MAGMA is not linked: change make.inc or the input file to continue.\n");
		return 0;
	}
#endif


    /* Write parameters, help to debug */
    /*bblas_zprintconfig(test); */

    /* set number of mkl threads  */
    if (test->mkl_sequential)
      {
#if defined(BBLAS_WITH_MKL)
	mkl_set_num_threads(1);
#endif
      }

    /*Set faulty iter */
    while(test->faulty_iter < 2)
    {
	test->faulty_iter = irandRange(0, test->nb_test);
    }

    switch (test->routine)
    {
	case BBLAS_GEMM:

	    testing_zgemm_batch(test);
	    break;

	case BBLAS_HEMM:
	    testing_zhemm_batch(test);
	    break;

	case BBLAS_HER2K:
	    testing_zher2k_batch(test);
	    break;

	case BBLAS_HERK:
	    testing_zherk_batch(test);
	    break;

	case BBLAS_SYMM:
	    testing_zsymm_batch(test);
	    break;

	case BBLAS_SYR2K:
	    testing_zsyr2k_batch(test);
	    break;

	case BBLAS_SYRK:
	    testing_zsyrk_batch(test);
	    break;

	case BBLAS_TRMM:
	    testing_ztrmm_batch(test);
	    break;

	case BBLAS_TRSM:
	    testing_ztrsm_batch(test);
	    break;

	default:
	    printf("In main: wrong routine value\n");
	    return 0;
    }


    /* Free test memory */
    free(test);


}

#undef COMPLEX
