/**
 * @file bblas_stesting.h
 *
 * @brief BBLAS testing header for float routines.
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
 * @generated from bblas_ztesting.h normal z -> s, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_STESTING_H
/** Include guard **/
#define BBLAS_STESTING_H

#include "bblas_testing_common.h"

/**
 * Main test structure used throughout all of our tests.
 *
 * The structure holds all the input parameters,
 * the matrices which we want to operate on,
 * the results of the computation,
 * and other items such as the computation time and flop count etc.
 **/
typedef struct bblas_stest
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
    float *alpha;
	/** Pointer to beta values **/
    float *beta;

	/** Pointer to alpha values for herk **/
    float *alpha_herk;
	/** Pointer to beta values for herk **/
    float *beta_herk;

    /** Pointer to test matrices **/
	/** Pointer to A matrices **/
    float **arrayA;
	/** Pointer to B matrices **/
    float **arrayB;
	/** Pointer to C matrices **/
    float **arrayC;

    /** Pointer to device (GPU etc.) matrices **/
	/** Pointer to A_device matrices **/
    float **arrayA_device;
	/** Pointer to B_device matrices **/
    float **arrayB_device;
	/** Pointer to C_device matrices **/
    float **arrayC_device;

    /** Pointer on host to device matrices **/
	/** Pointer to A_device matrices on host **/
    float **arrayA_h_d;
	/** Pointer to B_device matrices on host **/
    float **arrayB_h_d;
	/** Pointer to C_device matrices on host **/
    float **arrayC_h_d;

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
    float **mkl_result;
	/** Device results **/
    float **device_result;
	/** Other results **/
    float **other_result;
	/** Initial B values **/
    float **Binit;

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
}bblas_stest_t;





/************************************************************
 *Prototype of functions declared main BBLAS testing routines
 *************************************************************/
void testing_sgemm_batch(bblas_stest_t *test);
void testing_ssymm_batch(bblas_stest_t *test);
void testing_ssyr2k_batch(bblas_stest_t *test);
void testing_ssyrk_batch(bblas_stest_t *test);
void testing_ssymm_batch(bblas_stest_t *test);
void testing_ssyr2k_batch(bblas_stest_t *test);
void testing_ssyrk_batch(bblas_stest_t *test);
void testing_strmm_batch(bblas_stest_t *test);
void testing_strsm_batch(bblas_stest_t *test);

/*******************************************************
 *Prototype of functions declared in bblas_stestings.c
 *******************************************************/

void bblas_sinit_config (bblas_stest_t *test);
void bblas_ssettest(bblas_stest_t *test);
void bblas_sfreetest(bblas_stest_t *test);
void bblas_scopy_Binit(bblas_stest_t *test, float **B_copy);
void bblas_scopy_Cinit(bblas_stest_t *test, float **C_copy);
void bblas_scheck_Cfinal(bblas_stest_t *test, float **C_final);
void bblas_scheck_Bfinal(bblas_stest_t *test, float **B_final);
void parse_sconfigfile (FILE *file,  bblas_stest_t *test);
void bblas_sprintmatrix(float *matrix, int row, int col);
void bblas_sprintconfig(bblas_stest_t *test);
void bblas_spassed_failed(bblas_stest_t  *test, float error, char *result, int info);
float bblas_srandRange( int max_n);
void bblas_sflops(bblas_stest_t *test);
void bblas_smake_symmetric(int lda, int N, float* A );
void bblas_smake_symmetric(int lda, int N, float* A);
void bblas_sstatistic(bblas_stest_t *test);
void bblas_sprint_result(bblas_stest_t *test);
int  bblas_snbdata(bblas_stest_t *test);
void bblas_sprint_parameters(bblas_stest_t *test);
void bblas_sset_error(bblas_stest_t *test);

/*Accuracy computing functions */
void bblas_sgemm_tolerance(bblas_stest_t *test, float **C_final);
void bblas_sstark_tolerance(bblas_stest_t *test, float **C_final);
void bblas_sstar2k_tolerance(bblas_stest_t *test, float **C_final);
void bblas_ssyr2k_tolerance(bblas_stest_t *test, float **C_final);
void bblas_ssyrk_tolerance(bblas_stest_t *test, float **C_final);
void bblas_ssyr2k_tolerance(bblas_stest_t *test, float **C_final);
void bblas_ssyrk_tolerance(bblas_stest_t *test, float **C_final);
void bblas_ssymm_tolerance(bblas_stest_t *test, float **C_final);
void bblas_ssymm_tolerance(bblas_stest_t *test, float **C_final);
void bblas_strmm_tolerance(bblas_stest_t *test, float **B_final);
void bblas_strsm_tolerance(bblas_stest_t *test, float **B_final);

/* Splited functions for allocation */
void bblas_ssetuplo(bblas_stest_t *test);
void bblas_ssettransA(bblas_stest_t *test);
void bblas_ssettransB(bblas_stest_t *test);
void bblas_ssettrans(bblas_stest_t *test);
void bblas_ssetside(bblas_stest_t *test);
void bblas_ssetdiag(bblas_stest_t *test);
void bblas_ssetM(bblas_stest_t *test);
void bblas_ssetK(bblas_stest_t *test);
void bblas_ssetN(bblas_stest_t *test);
void bblas_ssetlda(bblas_stest_t *test);
void bblas_ssetldb(bblas_stest_t *test);
void bblas_ssetldc(bblas_stest_t *test);
void bblas_ssetalpha(bblas_stest_t *test);
void bblas_ssetalpha_herk(bblas_stest_t *test);
void bblas_ssetbeta(bblas_stest_t *test);
void bblas_ssetbeta_herk(bblas_stest_t *test);
void bblas_ssetarrayA(bblas_stest_t *test);
void bblas_ssetarrayB(bblas_stest_t *test);
void bblas_ssetarrayC(bblas_stest_t *test);
void bblas_smalloc_result_error(bblas_stest_t *test);

/*Auxiliary functions */
void bblas_sset_batch_count(bblas_stest_t *test);
void bblas_sprint_sys_info(bblas_stest_t *test);

float bblas_smaxarrayD(float *myArray, int size);
float bblas_sminarrayD(float *myArray, int size);
float bblas_savgarrayD(float *myArray, int size);
float bblas_sstdarrayD(float *myArray, int size);

/*Prototype  of cuda functions */
void bblas_ssettest_cuda(bblas_stest_t *test);
void bblas_ssetarrayA_device(bblas_stest_t *test);
void bblas_ssetarrayB_device(bblas_stest_t *test);
void bblas_ssetarrayC_device(bblas_stest_t *test);
void bblas_scudaFree( bblas_stest_t *test );

void bblas_sget_Cfinal(bblas_stest_t *test);
void bblas_sget_Bfinal(bblas_stest_t *test);
void bblas_sprint_cuda_info(bblas_stest_t *test);


#endif /* BBLAS_STESTINGS_H */
