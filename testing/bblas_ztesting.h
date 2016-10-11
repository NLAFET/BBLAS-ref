/**
 * @file bblas_ztesting.h
 *
 * @brief BBLAS testing header for double _Complex routines.
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
 * @precisions normal z -> c d s
 **/
#endif

#ifndef BBLAS_ZTESTING_H
/** Include guard **/
#define BBLAS_ZTESTING_H

#include "bblas_testing_common.h"

/**
 * Main test structure used throughout all of our tests.
 *
 * The structure holds all the input parameters,
 * the matrices which we want to operate on,
 * the results of the computation,
 * and other items such as the computation time and flop count etc.
 **/
typedef struct bblas_ztest
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
    BBLAS_Complex64_t *alpha;
	/** Pointer to beta values **/
    BBLAS_Complex64_t *beta;

	/** Pointer to alpha values for herk **/
    double *alpha_herk;
	/** Pointer to beta values for herk **/
    double *beta_herk;

    /** Pointer to test matrices **/
	/** Pointer to A matrices **/
    BBLAS_Complex64_t **arrayA;
	/** Pointer to B matrices **/
    BBLAS_Complex64_t **arrayB;
	/** Pointer to C matrices **/
    BBLAS_Complex64_t **arrayC;

    /** Pointer to device (GPU etc.) matrices **/
	/** Pointer to A_device matrices **/
    BBLAS_Complex64_t **arrayA_device;
	/** Pointer to B_device matrices **/
    BBLAS_Complex64_t **arrayB_device;
	/** Pointer to C_device matrices **/
    BBLAS_Complex64_t **arrayC_device;

    /** Pointer on host to device matrices **/
	/** Pointer to A_device matrices on host **/
    BBLAS_Complex64_t **arrayA_h_d;
	/** Pointer to B_device matrices on host **/
    BBLAS_Complex64_t **arrayB_h_d;
	/** Pointer to C_device matrices on host **/
    BBLAS_Complex64_t **arrayC_h_d;

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
    BBLAS_Complex64_t **mkl_result;
	/** Device results **/
    BBLAS_Complex64_t **device_result;
	/** Other results **/
    BBLAS_Complex64_t **other_result;
	/** Initial B values **/
    BBLAS_Complex64_t **Binit;

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
}bblas_ztest_t;





/************************************************************
 *Prototype of functions declared main BBLAS testing routines
 *************************************************************/
void testing_zgemm_batch(bblas_ztest_t *test);
void testing_zhemm_batch(bblas_ztest_t *test);
void testing_zher2k_batch(bblas_ztest_t *test);
void testing_zherk_batch(bblas_ztest_t *test);
void testing_zsymm_batch(bblas_ztest_t *test);
void testing_zsyr2k_batch(bblas_ztest_t *test);
void testing_zsyrk_batch(bblas_ztest_t *test);
void testing_ztrmm_batch(bblas_ztest_t *test);
void testing_ztrsm_batch(bblas_ztest_t *test);

/*******************************************************
 *Prototype of functions declared in bblas_ztestings.c
 *******************************************************/

void bblas_zinit_config (bblas_ztest_t *test);
void bblas_zsettest(bblas_ztest_t *test);
void bblas_zfreetest(bblas_ztest_t *test);
void bblas_zcopy_Binit(bblas_ztest_t *test, BBLAS_Complex64_t **B_copy);
void bblas_zcopy_Cinit(bblas_ztest_t *test, BBLAS_Complex64_t **C_copy);
void bblas_zcheck_Cfinal(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zcheck_Bfinal(bblas_ztest_t *test, BBLAS_Complex64_t **B_final);
void parse_zconfigfile (FILE *file,  bblas_ztest_t *test);
void bblas_zprintmatrix(BBLAS_Complex64_t *matrix, int row, int col);
void bblas_zprintconfig(bblas_ztest_t *test);
void bblas_zpassed_failed(bblas_ztest_t  *test, double error, char *result, int info);
BBLAS_Complex64_t bblas_zrandRange( int max_n);
void bblas_zflops(bblas_ztest_t *test);
void bblas_zmake_hermitian(int lda, int N, BBLAS_Complex64_t* A );
void bblas_zmake_symmetric(int lda, int N, BBLAS_Complex64_t* A);
void bblas_zstatistic(bblas_ztest_t *test);
void bblas_zprint_result(bblas_ztest_t *test);
int  bblas_znbdata(bblas_ztest_t *test);
void bblas_zprint_parameters(bblas_ztest_t *test);
void bblas_zset_error(bblas_ztest_t *test);

/*Accuracy computing functions */
void bblas_zgemm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zstark_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zstar2k_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zsyr2k_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zsyrk_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zher2k_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zherk_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zhemm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_zsymm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **C_final);
void bblas_ztrmm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **B_final);
void bblas_ztrsm_tolerance(bblas_ztest_t *test, BBLAS_Complex64_t **B_final);

/* Splited functions for allocation */
void bblas_zsetuplo(bblas_ztest_t *test);
void bblas_zsettransA(bblas_ztest_t *test);
void bblas_zsettransB(bblas_ztest_t *test);
void bblas_zsettrans(bblas_ztest_t *test);
void bblas_zsetside(bblas_ztest_t *test);
void bblas_zsetdiag(bblas_ztest_t *test);
void bblas_zsetM(bblas_ztest_t *test);
void bblas_zsetK(bblas_ztest_t *test);
void bblas_zsetN(bblas_ztest_t *test);
void bblas_zsetlda(bblas_ztest_t *test);
void bblas_zsetldb(bblas_ztest_t *test);
void bblas_zsetldc(bblas_ztest_t *test);
void bblas_zsetalpha(bblas_ztest_t *test);
void bblas_zsetalpha_herk(bblas_ztest_t *test);
void bblas_zsetbeta(bblas_ztest_t *test);
void bblas_zsetbeta_herk(bblas_ztest_t *test);
void bblas_zsetarrayA(bblas_ztest_t *test);
void bblas_zsetarrayB(bblas_ztest_t *test);
void bblas_zsetarrayC(bblas_ztest_t *test);
void bblas_zmalloc_result_error(bblas_ztest_t *test);

/*Auxiliary functions */
void bblas_zset_batch_count(bblas_ztest_t *test);
void bblas_zprint_sys_info(bblas_ztest_t *test);

double bblas_zmaxarrayD(double *myArray, int size);
double bblas_zminarrayD(double *myArray, int size);
double bblas_zavgarrayD(double *myArray, int size);
double bblas_zstdarrayD(double *myArray, int size);

/*Prototype  of cuda functions */
void bblas_zsettest_cuda(bblas_ztest_t *test);
void bblas_zsetarrayA_device(bblas_ztest_t *test);
void bblas_zsetarrayB_device(bblas_ztest_t *test);
void bblas_zsetarrayC_device(bblas_ztest_t *test);
void bblas_zcudaFree( bblas_ztest_t *test );

void bblas_zget_Cfinal(bblas_ztest_t *test);
void bblas_zget_Bfinal(bblas_ztest_t *test);
void bblas_zprint_cuda_info(bblas_ztest_t *test);


#endif /* BBLAS_ZTESTINGS_H */
