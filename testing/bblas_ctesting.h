/**
 * @file bblas_ctesting.h
 *
 * @brief BBLAS testing header for float _Complex routines.
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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from bblas_ztesting.h normal z -> c, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_CTESTING_H
/** Include guard **/
#define BBLAS_CTESTING_H

#include "bblas_testing_common.h"

/**
 * Main test structure used throughout all of our tests.
 *
 * The structure holds all the input parameters,
 * the matrices which we want to operate on,
 * the results of the computation,
 * and other items such as the computation time and flop count etc.
 **/
typedef struct bblas_ctest
{

    /** Store the input parameters used to generate test matrices **/
	/** Store UPLO input parameter **/
    enum BBLAS_UPLO    *uplo;
	/** Store TRANSA input parameter **/
    enum BBLAS_TRANS   *transA;
	/** Store TRANSB input parameter **/
    enum BBLAS_TRANS   *transB;
	/** Store TRANS input parameter **/
    enum BBLAS_TRANS   *trans;
	/** Store SIDE input parameter **/
    enum BBLAS_SIDE    *side;
	/** Store DIAG input parameter **/
    enum BBLAS_DIAG    *diag;


    /** Min and max of M, N, K, BATCH_COUNT used to generate test matrices **/
	/** Minumum M **/
    int                minM;
	/** Maximum M **/
    int                maxM;
	/** Minumum N **/
    int                minN;
	/** Maximum N **/
    int                maxN;
	/** Minumum K **/
    int                minK;
	/** Maximum K **/
    int                maxK;
	/** Minumum BATCH_COUNT **/
    int                minbatch_count;
	/** Maximum BATCH_COUNT **/
    int                maxbatch_count;

    /** Matrix sizes and dimensions used for each test matrix **/
	/** M used for test matrices **/
    int               *M;
	/** N used for test matrices **/
    int               *N;
	/** K used for test matrices **/
    int               *K;
	/** LDA used for test matrices **/
    int               *lda;
	/** LDB used for test matrices **/
    int               *ldb;
	/** LDC used for test matrices **/
    int               *ldc;

    /** Pointer to scalars alpha & beta  **/
	/** Pointer to alpha values **/
    BBLAS_Complex32_t *alpha;
	/** Pointer to beta values **/
    BBLAS_Complex32_t *beta;

	/** Pointer to alpha values for herk **/
    float *alpha_herk;
	/** Pointer to beta values for herk **/
    float *beta_herk;

    /** Pointer to test matrices **/
	/** Pointer to A matrices **/
    BBLAS_Complex32_t **arrayA;
	/** Pointer to B matrices **/
    BBLAS_Complex32_t **arrayB;
	/** Pointer to C matrices **/
    BBLAS_Complex32_t **arrayC;

    /** Pointer to device (GPU etc.) matrices **/
	/** Pointer to A_device matrices **/
    BBLAS_Complex32_t **arrayA_device;
	/** Pointer to B_device matrices **/
    BBLAS_Complex32_t **arrayB_device;
	/** Pointer to C_device matrices **/
    BBLAS_Complex32_t **arrayC_device;

    /** Pointer on host to device matrices **/
	/** Pointer to A_device matrices on host **/
    BBLAS_Complex32_t **arrayA_h_d;
	/** Pointer to B_device matrices on host **/
    BBLAS_Complex32_t **arrayB_h_d;
	/** Pointer to C_device matrices on host **/
    BBLAS_Complex32_t **arrayC_h_d;

    /** Number of problems in the batch **/
    int               batch_count;

    /** Batch opts (fixed, variable, etc.) **/
    enum BBLAS_OPTS   batch_opts;

    /** Pointer to info for error handling **/
    int               *info;

    /** Batch BLAS routine to test **/
    enum BBLAS_ROUTINE routine;

    /** Number of tests to run (size increases each time) **/
    int               nb_test;

    /** Number of implementations to be compared to **/
    int               nb_comparison;

    /** Number of OpenMP threads **/
    int               nb_threads;

    /** Accuracy test tolerance **/
    BBLAS_Double_t    tolerance;

    /** Norms of input matrices for accuracy computation **/
	/** Initial norm of C matrices **/
    float  *Cinitnorm;
	/** Initial norm of B matrices **/
    float  *Binitnorm;

    /*
     * Timings of differents BBLAS, Reference BBLAS
     * batch Device, batch MKL, etc.
     */

    /** Timing of Batch BLAS implementations **/
	/** Reference BBLAS time **/
    BBLAS_Double_t           ref_time;
	/** MKL BBLAS time **/
    BBLAS_Double_t           mkl_time;
	/** Device BBLAS time **/
    BBLAS_Double_t           device_time;
	/** Other implementation (default OpenMP) BBLAS time **/
    BBLAS_Double_t           other_time;

    /*
     * Flop rates of differents BBLAS, Reference BBLAS
     * batch Device, batch MKL, etc.
     */
	/** Theoretical flops required **/
    BBLAS_Double_t          flops;
	/** Flop rate for reference BBLAS **/
    BBLAS_Double_t          ref_perf;
	/** Flop rate for MKL BBLAS **/
    BBLAS_Double_t          mkl_perf;
	/** Flop rate for Device BBLAS **/
    BBLAS_Double_t          device_perf;
	/** Flop rate for other (defulat OpenMP) BBLAS **/
    BBLAS_Double_t          other_perf;

    /*
     * Accurary against differents BBLAS,
     * batch Device, batch MKL, etc.
     */
    /** Target bblas to be compared to **/
    enum BBLAS_TARGET       target;

    /** MKL error variables **/
	/** MKL error list **/
    float          *mkl_error;
	/** MKL min error **/
	float          mkl_min_error;
	/** MKL max error **/
	float          mkl_max_error;
	/** MKL error statistics **/
	/** MKL average error **/
    float          mkl_avg_error;
	/** MKL standard deviation of error **/
    float          mkl_std_error;

	/** Device error variables **/
	/** Device error list **/
    float          *device_error;
	/** Device min error **/
	float          device_min_error;
	/** Device max error **/
    float          device_max_error;
	/** Device error statistics **/
	/** Device average error **/
    float          device_avg_error;
	/** Device standard deviation of error **/
	float          device_std_error;

    /** Other error variables **/
	/** Other error list **/
    float          *other_error;
	/** Other min error **/
	float          other_min_error;
	/** Other max error **/
	float          other_max_error;
	/** Other error statistics **/
	/** Other average error **/
    float          other_avg_error;
	/** Other standard deviation of error **/
    float          other_std_error;


    /** Some variables required for statistics **/
	/** Average value of M **/
    int                avgM;
	/** Average value of N **/
    int                avgN;
	/** Average value of K **/
    int                avgK;

    /** Array of pointers to copied results for comparison **/
	/** MKL results **/
    BBLAS_Complex32_t **mkl_result;
	/** Device results **/
    BBLAS_Complex32_t **device_result;
	/** Other results **/
    BBLAS_Complex32_t **other_result;
	/** Initial B values **/
    BBLAS_Complex32_t **Binit;

#if defined(BBLAS_WITH_CUBLAS)
    /** cublasHandle for the CuBLAS calls **/
    cublasHandle_t cublas_handle;

    /** cublasStream to synchronize the CuBLAS calls **/
    cudaStream_t stream;
#endif

    /** Variables to specify values of transA, transB ...  **/
	/** Value of uplo **/
    int gen_uplo;
	/** Value of transA **/
    int gen_transA;
	/** Value of transB **/
    int gen_transB;
	/** Value of trans **/
    int gen_trans;
	/** Value of side **/
    int gen_side;
	/** Value of diag **/
    int gen_diag;

    /** MKL group_size array for BATCH MKL calls **/
    int     *group_size;

    /** Variables for unit test **/
	/** Do we want to set an error in an iteration? **/
    int     set_error;
	/** Do we set a global error? **/
    int     global_error;
	/** Which iteration should be faulty? **/
    int     faulty_iter;
	/** Is it the current iteration? **/
    int     current_iter;

	/** Do we use sequential MKL or multithreaded? **/
    int     mkl_sequential;
	/** Do we use new accuracy tests (based on error analysis) or old (based on relative error) **/
    int     new_accuracy;
}bblas_ctest_t;





/************************************************************
 *Prototype of functions declared main BBLAS testing routines
 *************************************************************/
void testing_cgemm_batch(bblas_ctest_t *test);
void testing_chemm_batch(bblas_ctest_t *test);
void testing_cher2k_batch(bblas_ctest_t *test);
void testing_cherk_batch(bblas_ctest_t *test);
void testing_csymm_batch(bblas_ctest_t *test);
void testing_csyr2k_batch(bblas_ctest_t *test);
void testing_csyrk_batch(bblas_ctest_t *test);
void testing_ctrmm_batch(bblas_ctest_t *test);
void testing_ctrsm_batch(bblas_ctest_t *test);

/*******************************************************
 *Prototype of functions declared in bblas_ctestings.c
 *******************************************************/

void bblas_cinit_config (bblas_ctest_t *test);
void bblas_csettest(bblas_ctest_t *test);
void bblas_cfreetest(bblas_ctest_t *test);
void bblas_ccopy_Binit(bblas_ctest_t *test, BBLAS_Complex32_t **B_copy);
void bblas_ccopy_Cinit(bblas_ctest_t *test, BBLAS_Complex32_t **C_copy);
void bblas_ccheck_Cfinal(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_ccheck_Bfinal(bblas_ctest_t *test, BBLAS_Complex32_t **B_final);
void parse_cconfigfile (FILE *file,  bblas_ctest_t *test);
void bblas_cprintmatrix(BBLAS_Complex32_t *matrix, int row, int col);
void bblas_cprintconfig(bblas_ctest_t *test);
void bblas_cpassed_failed(bblas_ctest_t  *test, float error, char *result, int info);
BBLAS_Complex32_t bblas_crandRange( int max_n);
void bblas_cflops(bblas_ctest_t *test);
void bblas_cmake_hermitian(int lda, int N, BBLAS_Complex32_t* A );
void bblas_cmake_symmetric(int lda, int N, BBLAS_Complex32_t* A);
void bblas_cstatistic(bblas_ctest_t *test);
void bblas_cprint_result(bblas_ctest_t *test);
int  bblas_cnbdata(bblas_ctest_t *test);
void bblas_cprint_parameters(bblas_ctest_t *test);
void bblas_cset_error(bblas_ctest_t *test);

/*Accuracy computing functions */
void bblas_cgemm_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_cstark_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_cstar2k_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_csyr2k_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_csyrk_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_cher2k_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_cherk_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_chemm_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_csymm_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **C_final);
void bblas_ctrmm_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **B_final);
void bblas_ctrsm_tolerance(bblas_ctest_t *test, BBLAS_Complex32_t **B_final);

/* Splited functions for allocation */
void bblas_csetuplo(bblas_ctest_t *test);
void bblas_csettransA(bblas_ctest_t *test);
void bblas_csettransB(bblas_ctest_t *test);
void bblas_csettrans(bblas_ctest_t *test);
void bblas_csetside(bblas_ctest_t *test);
void bblas_csetdiag(bblas_ctest_t *test);
void bblas_csetM(bblas_ctest_t *test);
void bblas_csetK(bblas_ctest_t *test);
void bblas_csetN(bblas_ctest_t *test);
void bblas_csetlda(bblas_ctest_t *test);
void bblas_csetldb(bblas_ctest_t *test);
void bblas_csetldc(bblas_ctest_t *test);
void bblas_csetalpha(bblas_ctest_t *test);
void bblas_csetalpha_herk(bblas_ctest_t *test);
void bblas_csetbeta(bblas_ctest_t *test);
void bblas_csetbeta_herk(bblas_ctest_t *test);
void bblas_csetarrayA(bblas_ctest_t *test);
void bblas_csetarrayB(bblas_ctest_t *test);
void bblas_csetarrayC(bblas_ctest_t *test);
void bblas_cmalloc_result_error(bblas_ctest_t *test);

/*Auxiliary functions */
void bblas_cset_batch_count(bblas_ctest_t *test);
void bblas_cprint_sys_info(bblas_ctest_t *test);

float bblas_cmaxarrayD(float *myArray, int size);
float bblas_cminarrayD(float *myArray, int size);
float bblas_cavgarrayD(float *myArray, int size);
float bblas_cstdarrayD(float *myArray, int size);

/*Prototype  of cuda functions */
void bblas_csettest_cuda(bblas_ctest_t *test);
void bblas_csetarrayA_device(bblas_ctest_t *test);
void bblas_csetarrayB_device(bblas_ctest_t *test);
void bblas_csetarrayC_device(bblas_ctest_t *test);
void bblas_ccudaFree( bblas_ctest_t *test );

void bblas_cget_Cfinal(bblas_ctest_t *test);
void bblas_cget_Bfinal(bblas_ctest_t *test);
void bblas_cprint_cuda_info(bblas_ctest_t *test);


#endif /* BBLAS_CTESTINGS_H */
