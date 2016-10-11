/**
 * @file bblas_dtesting.h
 *
 * @brief BBLAS testing header for double routines.
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
 * @generated from bblas_ztesting.h normal z -> d, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_DTESTING_H
/** Include guard **/
#define BBLAS_DTESTING_H

#include "bblas_testing_common.h"

/**
 * Main test structure used throughout all of our tests.
 *
 * The structure holds all the input parameters,
 * the matrices which we want to operate on,
 * the results of the computation,
 * and other items such as the computation time and flop count etc.
 **/
typedef struct bblas_dtest
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
    double *alpha;
	/** Pointer to beta values **/
    double *beta;

	/** Pointer to alpha values for herk **/
    double *alpha_herk;
	/** Pointer to beta values for herk **/
    double *beta_herk;

    /** Pointer to test matrices **/
	/** Pointer to A matrices **/
    double **arrayA;
	/** Pointer to B matrices **/
    double **arrayB;
	/** Pointer to C matrices **/
    double **arrayC;

    /** Pointer to device (GPU etc.) matrices **/
	/** Pointer to A_device matrices **/
    double **arrayA_device;
	/** Pointer to B_device matrices **/
    double **arrayB_device;
	/** Pointer to C_device matrices **/
    double **arrayC_device;

    /** Pointer on host to device matrices **/
	/** Pointer to A_device matrices on host **/
    double **arrayA_h_d;
	/** Pointer to B_device matrices on host **/
    double **arrayB_h_d;
	/** Pointer to C_device matrices on host **/
    double **arrayC_h_d;

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
    double  *Cinitnorm;
	/** Initial norm of B matrices **/
    double  *Binitnorm;

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
    double          *mkl_error;
	/** MKL min error **/
	double          mkl_min_error;
	/** MKL max error **/
	double          mkl_max_error;
	/** MKL error statistics **/
	/** MKL average error **/
    double          mkl_avg_error;
	/** MKL standard deviation of error **/
    double          mkl_std_error;

	/** Device error variables **/
	/** Device error list **/
    double          *device_error;
	/** Device min error **/
	double          device_min_error;
	/** Device max error **/
    double          device_max_error;
	/** Device error statistics **/
	/** Device average error **/
    double          device_avg_error;
	/** Device standard deviation of error **/
	double          device_std_error;

    /** Other error variables **/
	/** Other error list **/
    double          *other_error;
	/** Other min error **/
	double          other_min_error;
	/** Other max error **/
	double          other_max_error;
	/** Other error statistics **/
	/** Other average error **/
    double          other_avg_error;
	/** Other standard deviation of error **/
    double          other_std_error;


    /** Some variables required for statistics **/
	/** Average value of M **/
    int                avgM;
	/** Average value of N **/
    int                avgN;
	/** Average value of K **/
    int                avgK;

    /** Array of pointers to copied results for comparison **/
	/** MKL results **/
    double **mkl_result;
	/** Device results **/
    double **device_result;
	/** Other results **/
    double **other_result;
	/** Initial B values **/
    double **Binit;

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
}bblas_dtest_t;





/************************************************************
 *Prototype of functions declared main BBLAS testing routines
 *************************************************************/
void testing_dgemm_batch(bblas_dtest_t *test);
void testing_dsymm_batch(bblas_dtest_t *test);
void testing_dsyr2k_batch(bblas_dtest_t *test);
void testing_dsyrk_batch(bblas_dtest_t *test);
void testing_dsymm_batch(bblas_dtest_t *test);
void testing_dsyr2k_batch(bblas_dtest_t *test);
void testing_dsyrk_batch(bblas_dtest_t *test);
void testing_dtrmm_batch(bblas_dtest_t *test);
void testing_dtrsm_batch(bblas_dtest_t *test);

/*******************************************************
 *Prototype of functions declared in bblas_dtestings.c
 *******************************************************/

void bblas_dinit_config (bblas_dtest_t *test);
void bblas_dsettest(bblas_dtest_t *test);
void bblas_dfreetest(bblas_dtest_t *test);
void bblas_dcopy_Binit(bblas_dtest_t *test, double **B_copy);
void bblas_dcopy_Cinit(bblas_dtest_t *test, double **C_copy);
void bblas_dcheck_Cfinal(bblas_dtest_t *test, double **C_final);
void bblas_dcheck_Bfinal(bblas_dtest_t *test, double **B_final);
void parse_dconfigfile (FILE *file,  bblas_dtest_t *test);
void bblas_dprintmatrix(double *matrix, int row, int col);
void bblas_dprintconfig(bblas_dtest_t *test);
void bblas_dpassed_failed(bblas_dtest_t  *test, double error, char *result, int info);
double bblas_drandRange( int max_n);
void bblas_dflops(bblas_dtest_t *test);
void bblas_dmake_symmetric(int lda, int N, double* A );
void bblas_dmake_symmetric(int lda, int N, double* A);
void bblas_dstatistic(bblas_dtest_t *test);
void bblas_dprint_result(bblas_dtest_t *test);
int  bblas_dnbdata(bblas_dtest_t *test);
void bblas_dprint_parameters(bblas_dtest_t *test);
void bblas_dset_error(bblas_dtest_t *test);

/*Accuracy computing functions */
void bblas_dgemm_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dstark_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dstar2k_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dsyr2k_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dsyrk_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dsyr2k_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dsyrk_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dsymm_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dsymm_tolerance(bblas_dtest_t *test, double **C_final);
void bblas_dtrmm_tolerance(bblas_dtest_t *test, double **B_final);
void bblas_dtrsm_tolerance(bblas_dtest_t *test, double **B_final);

/* Splited functions for allocation */
void bblas_dsetuplo(bblas_dtest_t *test);
void bblas_dsettransA(bblas_dtest_t *test);
void bblas_dsettransB(bblas_dtest_t *test);
void bblas_dsettrans(bblas_dtest_t *test);
void bblas_dsetside(bblas_dtest_t *test);
void bblas_dsetdiag(bblas_dtest_t *test);
void bblas_dsetM(bblas_dtest_t *test);
void bblas_dsetK(bblas_dtest_t *test);
void bblas_dsetN(bblas_dtest_t *test);
void bblas_dsetlda(bblas_dtest_t *test);
void bblas_dsetldb(bblas_dtest_t *test);
void bblas_dsetldc(bblas_dtest_t *test);
void bblas_dsetalpha(bblas_dtest_t *test);
void bblas_dsetalpha_herk(bblas_dtest_t *test);
void bblas_dsetbeta(bblas_dtest_t *test);
void bblas_dsetbeta_herk(bblas_dtest_t *test);
void bblas_dsetarrayA(bblas_dtest_t *test);
void bblas_dsetarrayB(bblas_dtest_t *test);
void bblas_dsetarrayC(bblas_dtest_t *test);
void bblas_dmalloc_result_error(bblas_dtest_t *test);

/*Auxiliary functions */
void bblas_dset_batch_count(bblas_dtest_t *test);
void bblas_dprint_sys_info(bblas_dtest_t *test);

double bblas_dmaxarrayD(double *myArray, int size);
double bblas_dminarrayD(double *myArray, int size);
double bblas_davgarrayD(double *myArray, int size);
double bblas_dstdarrayD(double *myArray, int size);

/*Prototype  of cuda functions */
void bblas_dsettest_cuda(bblas_dtest_t *test);
void bblas_dsetarrayA_device(bblas_dtest_t *test);
void bblas_dsetarrayB_device(bblas_dtest_t *test);
void bblas_dsetarrayC_device(bblas_dtest_t *test);
void bblas_dcudaFree( bblas_dtest_t *test );

void bblas_dget_Cfinal(bblas_dtest_t *test);
void bblas_dget_Bfinal(bblas_dtest_t *test);
void bblas_dprint_cuda_info(bblas_dtest_t *test);


#endif /* BBLAS_DTESTINGS_H */
