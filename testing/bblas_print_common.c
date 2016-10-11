/**
 * @file bblas_print_common.c
 *
 * @brief BBLAS printing routines for double _Complex routines.
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

/** Define bool **/
typedef int bool;
/** Boolean TRUE **/
#define TRUE  1
/** Boolean FALSE **/
#define FALSE 0


#include "bblas_common.h"

#if defined(BBLAS_WITH_MKL)
#include <mkl_service.h>
#endif
#include <omp.h>

#include <sys/utsname.h>

/**
 * Print the header of the results table.
 * This depends upon the value of target within the input file (whether we compare against MKL etc).
 **/

void bblas_print_header(enum BBLAS_TARGET target)
{

    switch(target)
    {
	case BBLAS_MKL:

	    printf("\n");
	    printf("%% Id    BatchCount      BBLAS-SEQ Gflop/s ( ms)       MKL Gflop/s (ms)          min(error)   avg(error)   max(error)   std(error)   PASSED/FAILED\n");
	    printf("%%================================================================================================================================================\n");
	    break;

	case BBLAS_CUBLAS:
	    printf("\n");
	    printf("%% Id    BatchCount   BBLAS-SEQ Gflop/s ( ms)          CUBLAS Gflop/s (ms)       min(error)   avg(error)   max(error)   std(error)   PASSED/FAILED\n");
	    printf("%%================================================================================================================================================\n");
	    break;

	case BBLAS_MAGMA:
	    printf("\n");
	    printf("%% Id    BatchCount   BBLAS-SEQ Gflop/s ( ms)          MAGMA  Gflop/s (ms)       min(error)   avg(error)   max(error)   std(error)   PASSED/FAILED\n");
	    printf("%%================================================================================================================================================\n");
	    break;


	case BBLAS_OTHER:
	    printf("\n");
	    printf("%% Id    BatchCount    BBLAS-SEQ Gflop/s ( ms)         OMP Gflop/s ( ms)         min(error)   avg(error)   max(error)   std(error)  PASSED/FAILED\n");
	    printf("%%===============================================================================================================================================\n");
	    break;

	case BBLAS_CUMKL:
	    /*To be implemented */
	    break;

	default:
	    printf("In bblas_print_header(): Target no defined\n");
	    exit(EXIT_FAILURE);
    }
}


/**
 * Print out OpenMP information.
 **/

void bblas_print_omp_info()
{
    printf("\nOpenMP Information:\n");
    int omp_threads = omp_get_max_threads();
    printf( "OpenMP threads %d. ", omp_threads );
}

/**
 * Print out MKL information.
 **/

void bblas_print_mkl_info(int mkl_sequential)
{
#if defined(BBLAS_WITH_MKL)
    printf("\nMKL Information:\n");
    int max_threads;

    if(mkl_sequential)
      {
	max_threads = 1;
      }else
      {
	max_threads = mkl_get_max_threads();
      }

    MKLVersion mkl_version;
    mkl_get_version( &mkl_version );
    printf( "   MKL version: %d.%d.%d, \n   MKL threads: %d.\n",
	    mkl_version.MajorVersion,
	    mkl_version.MinorVersion,
	    mkl_version.UpdateVersion,
	    max_threads);
#endif
}


/**
 * Print out general information for all systems.
 **/

void bblas_print_platform_info()
{
    struct utsname info;
    FILE* fp;
    char buffer[1024];
    int bytes_read;
    char* match;
    char model_name[100];
    int cpu_cores;
    char cache_size[20];
    char mem_size_char[20];
    int mem_size_int;

    uname(&info);

    printf("Operating System: %s\n",info.sysname);
    printf("    Release: %s\n",info.release);
    printf("    Version: %s\n\n",info.version);
    printf("Hardware:\n");
    //printf("    Architecture: %s\n",uname.machine);

    /* Read the entire contents of /proc/cpuinfo into the buffer.  */
    fp = fopen ("/proc/cpuinfo", "r");

    bytes_read = fread (buffer, 1, sizeof (buffer), fp);
    fclose (fp);

    /* NULL-terminate the text.*/
    buffer[bytes_read] = '\0';

    /* Locate the line that starts with "model name". */
    match = strstr (buffer, "model name");

    if (match != NULL)
    {
        sscanf (match, "model name  :  %[^\n]s", model_name);
        printf ("   CPU Model: %s\n", model_name);
    }
    match = strstr(buffer, "cpu cores");
    if (match != NULL)
    {
        sscanf (match, "cpu cores  :  %d", &cpu_cores);
        printf ("   CPU cores: %d\n", cpu_cores);
    }
    match = strstr(buffer, "cache size");
    if (match != NULL)
    {
        sscanf (match, "cache size  :  %[^\n]s", cache_size);
        printf ("   Cache size: %s\n", cache_size);
    }
    if (match == NULL)
    {
        printf("Could not find cpu info.\n");
    }

    /* Read the entire contents of /proc/meminfo into the buffer.  */
    fp = fopen ("/proc/meminfo", "r");

    bytes_read = fread (buffer, 1, sizeof (buffer), fp);
    fclose (fp);

    /* NULL-terminate the text.*/
    buffer[bytes_read] = '\0';

    /* Locate the line that starts with "MemTotal". */
    match = strstr (buffer, "MemTotal");

    if (match != NULL)
    {
        sscanf (match, "MemTotal  :  %s", mem_size_char);
        mem_size_int=atoi(mem_size_char);
        printf ("   Memory (RAM) size: %d MB\n", mem_size_int/1024);
    }
    if (match == NULL)
    {
        printf("Could not find mem info.\n");
    }

}

/**
 * Print out the current time.
 **/

void bblas_print_time()
{
    int size_buffer = 26;
    time_t timer;
    char buffer[size_buffer];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(buffer, size_buffer, "\n%Y:%m:%d %H:%M:%S\n", tm_info);
    puts(buffer);

}


/**
 * Print out the error check used to define whether a computation passes or fails our accuracy test.
 **/

void bblas_print_tolerance(enum BBLAS_ROUTINE routine, int new_accuracy)
{

    if (new_accuracy)
    {
	switch(routine)
	{
	    case BBLAS_GEMM:
	    case BBLAS_HEMM:
	    case BBLAS_SYMM:
		printf("\n\t                      ||C - Cseq||_oo\n");
		printf("SUCCESS IF ___________________________________________________  <  1\n\n");
		printf("\t (N.|alpha|.||A||_oo.||B||_oo + |beta|.||Cinit||_oo).eps\n\n\n\n");
		break;

	    case BBLAS_SYRK:
	    case BBLAS_HERK:
		printf("\n\t                      ||C - Cseq||_oo\n");
		printf("SUCCESS IF ___________________________________________________  <  1\n\n");
		printf("\t(N.|alpha|.||A||_oo.||B||_one + |beta|.||Cinit||_oo).eps\n\n\n\n");
		break;

	    case BBLAS_SYR2K:
	    case BBLAS_HER2K:
		printf("\n\t                                     ||C - Cseq||_oo\n");
		printf("SUCCESS IF _______________________________________________________________________________  <  1\n\n");
		printf("\t(N.|alpha|.||A||_oo.||B||_one +  N.|alpha|.||B||_oo.||A||_one   |beta|.||Cinit||_oo).eps\n\n\n\n");
		break;


	    case BBLAS_TRSM:
		printf("\n\t                ||alpha*B - A*X||_oo\n");
		printf("SUCCESS IF _________________________________________  <  1\n\n");
		printf("\t\t(|alpha|.||B|| + ||A||.||X||).N.eps\n\n\n\n");
		break;

	    case BBLAS_TRMM:
		printf("\n\t             ||B - Bseq||_oo\n");
		printf("SUCCESS IF ________________________________  <  1\n\n");
		printf("\t (N.|alpha|.||A||_oo.||Binit||_oo).eps\n\n\n\n");
		break;
	}
    }else
     {
	 printf("\n\t     ||C - Cseq||_oo\n");
	 printf("SUCCESS IF ____________________  < TOLERANCE * EPS \n\n");
	 printf("\t     ||Cseq||_oo * N\n\n\n\n");

     }
}

/**
 * Print out a title for the table of parameters.
 **/

void bblas_print_title(enum BBLAS_ROUTINE routine, enum BBLAS_OPTS batch_opts)
{

    printf("\n\t\t\t\t========================================================\n");
    printf("\t\t\t\t     TIMING AND ACCURRACY TESTING of  BATCHED %s\n", bblas_getroutine(routine));
    printf("\t\t\t\t%38s  \n",bblas_getoption(batch_opts));
    printf("\t\t\t\t========================================================\n");
}
