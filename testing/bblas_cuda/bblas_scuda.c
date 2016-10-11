/**
 * @file bblas_scuda.c
 *
 * @brief BBLAS CuBLAS routines for float precision.
 *
 *  BBLAS is a software package provided by Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-03-31
 *
 * Contains wrappers to batch CuBLAS functions and all other CUDA calls required by the
 * testing routines.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from ./bblas_cuda/bblas_zcuda.c normal z -> s, Mon Jun  6 09:44:14 2016
 **/
#endif

#include "bblas_common.h"


#ifndef BBLAS_GPU_CHECK_H
/** Include guard for bblas_gpu_check **/
#define BBLAS_GPU_CHECK_H
/** Call GPU assert to check command **/
#define bblas_gpu_check(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#if defined(BBLAS_WITH_CUBLAS)
/**
 * Wrapper for checking success of CUDA calls.
 **/

void gpuAssert(cudaError_t code, const char *file, int line)
{

    if (code != cudaSuccess)
    {
	fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
	exit(EXIT_FAILURE);
    }

}
#endif
#endif
/** Include real and float real precision routines **/
#define REAL

#if defined(BBLAS_WITH_CUBLAS)

/**
 * Transform BBLAS enumerates to CuBLAS variants.
 **/

cublasOperation_t op_bb2cu(unsigned int op)
{
    switch(op)
    {
        case BblasNoTrans:
            return CUBLAS_OP_N;
            break;

        case BblasTrans:
            return CUBLAS_OP_T;
            break;

        case BblasConjTrans:
            return CUBLAS_OP_C;
            break;

        case BblasLower:
            return CUBLAS_FILL_MODE_LOWER;
            break;

        case BblasUpper:
            return CUBLAS_FILL_MODE_UPPER;
            break;

        case BblasNonUnit:
            return CUBLAS_DIAG_NON_UNIT;
            break;

        case BblasUnit:
            return CUBLAS_DIAG_UNIT;
            break;

        case BblasLeft:
            return CUBLAS_SIDE_LEFT;
            break;

        case BblasRight:
            return CUBLAS_SIDE_RIGHT;
            break;

        default:
            return 0;
            exit(EXIT_FAILURE);
    }

}
#endif

/**
 * Move arrays to the GPU for computation.
 *
 * Calls bblas_ssetarrayX_device where X is chosen from A, B, or C depending upon
 * the current routine that is being tested.
 *
 **/

void bblas_ssettest_cuda(bblas_stest_t *test)
{
  enum BBLAS_ROUTINE routine  = test->routine;

  /**
   * Set arrayA on GPU
   */
  bblas_ssetarrayA_device(test);

  if ((routine == BBLAS_GEMM) || (routine == BBLAS_TRSM))
    {
      /**
       * Set arrayB on GPU
       */
      bblas_ssetarrayB_device(test);
    }

  if ((routine == BBLAS_GEMM) || (routine == BBLAS_HERK))
    {
      /**
       * Set arrayC on GPU
       */
      bblas_ssetarrayC_device(test);
    }
}

/**
 * Copy the array A into GPU memory ready for computation.
 **/
void bblas_ssetarrayA_device(bblas_stest_t *test)
{
#if defined(BBLAS_WITH_CUBLAS)
    /*Set the scalar values batch_count */
    int batch_iter, nb_row, nb_col;
    int first_index                 = 0;
    enum BBLAS_ROUTINE routine  = test->routine;
    int batch_count             = test->batch_count;

    /* Error return for cuda calls */
    cudaError_t cuda_err;

    /*Memory allocation for host to device matrices*/
    test->arrayA_h_d = (float**) malloc(batch_count*sizeof(float*));
    bblas_malloc_check(test->arrayA_h_d, "arrayA_h_d");

    /*Memory allocation for device matrices*/
    cuda_err = cudaMalloc((void**)&test->arrayA_device,
			  batch_count*sizeof(float*));
    bblas_gpu_check(cuda_err);

    nb_row = test->lda[first_index];

    /* nb_col for GEMM */
    if (routine == BBLAS_GEMM)
    {
	if (test->transA[first_index] == BblasNoTrans)
	{
	    nb_col =test->K[first_index];
	}else
	{
	    nb_col =test->M[first_index];
	}
    }

    /* nb_col for HERK */
    if (routine == BBLAS_HERK)
      {
	if (test->trans[first_index] == BblasNoTrans)
	  {
	    nb_col = test->K[first_index];
	  }else
	  {
	    nb_col = test->N[first_index];
	  }
      }

    /*  nb_col for TRSM */
    if 	(routine == BBLAS_TRSM)
      {
	if(test->side[first_index] == BblasLeft )
	  {
	    nb_col = test->M[first_index];
	  }else
	  {
	    nb_col = test->N[first_index];
	  }
      }

    for( batch_iter =0; batch_iter < batch_count; batch_iter++)
    {
	cuda_err = cudaMalloc((void**)&test->arrayA_h_d[batch_iter],
			      nb_row*nb_col*sizeof(float));
	bblas_gpu_check(cuda_err);
    }

    /*Copy the host array of device pointers to the device*/
    cuda_err = cudaMemcpy(
	test->arrayA_device, test->arrayA_h_d,
	batch_count*sizeof(float*),
	cudaMemcpyHostToDevice);
    bblas_gpu_check(cuda_err);

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	cuda_err = cublasSetMatrix(
	    nb_row, nb_col, sizeof(float),
	    test->arrayA[batch_iter], test->lda[first_index],
	    test->arrayA_h_d[batch_iter], test->lda[first_index]);
	bblas_gpu_check(cuda_err);
    }

#endif
}


/**
 * Copy the array B into GPU memory ready for computation.
 **/
void bblas_ssetarrayB_device(bblas_stest_t *test)
{
#if defined(BBLAS_WITH_CUBLAS)
    /*Set the scalar values batch_count */
    int batch_iter, nb_row, nb_col;
    int first_index                 = 0;
    enum BBLAS_ROUTINE routine  = test->routine;
    int batch_count             = test->batch_count;

    /* Error return for cuda calls */
    cudaError_t cuda_err;

    /*Memory allocation for host to device matrices*/
    test->arrayB_h_d = (float**) malloc(batch_count*sizeof(float*));

    /*Memory allocation for device matrices*/
    cuda_err = cudaMalloc((void**)&test->arrayB_device, batch_count*sizeof(float*));
    bblas_gpu_check(cuda_err);

    nb_row =test->ldb[first_index];

    /*nb_col for BBLAS_GEMM */
    if (routine == BBLAS_GEMM)
    {
	if (test->transB[first_index] == BblasNoTrans)
	{
	    nb_col = test->N[first_index];
	}else
	{
	    nb_col = test->K[first_index];
	}
    }

    /* nb_col for TRSM */
    if ((routine == BBLAS_SYMM) || (routine == BBLAS_HEMM) ||
	(routine == BBLAS_TRMM) || (routine == BBLAS_TRSM))
      {
	nb_col =test->N[first_index];
      }

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	cuda_err = cudaMalloc((void**)&test->arrayB_h_d[batch_iter],
			      nb_row*nb_col*sizeof(float));
	bblas_gpu_check(cuda_err);
    }

    /*Copy the host array of device pointers to the device*/
    cuda_err = cudaMemcpy(test->arrayB_device, test->arrayB_h_d,
			  batch_count*sizeof(float*),
			  cudaMemcpyHostToDevice);

    bblas_gpu_check(cuda_err);

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	cuda_err = cublasSetMatrix(
	    nb_row, nb_col, sizeof(float),
	    test->arrayB[batch_iter], test->ldb[first_index],
	    test->arrayB_h_d[batch_iter], test->ldb[first_index]);
	bblas_gpu_check(cuda_err);

    }
#endif

}

/**
 * Copy the array C into GPU memory ready for computation.
 **/
void bblas_ssetarrayC_device(bblas_stest_t *test)
{
#if defined(BBLAS_WITH_CUBLAS)
    /*Set the scalar values batch_count */
    int batch_iter, nb_row, nb_col;
    int first_index             = 0;
    int batch_count             = test->batch_count;

    /* Error return for cuda calls */
    cudaError_t cuda_err;

    /*Memory allocation for host to device matrices*/
    test->arrayC_h_d = (float**) malloc(batch_count*sizeof(float*));

    /*Memory allocation for device matrices*/
    cuda_err = cudaMalloc((void**)&test->arrayC_device, batch_count*sizeof(float*));

    nb_row =test->ldc[first_index];
    nb_col =test->N[first_index];

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	cuda_err = cudaMalloc((void**)&test->arrayC_h_d[batch_iter], nb_row*nb_col*sizeof(float));
    }

    /* Copy the host array of device pointers to the device */
    cuda_err = cudaMemcpy(test->arrayC_device, test->arrayC_h_d,
			  batch_count*sizeof(float*),
			  cudaMemcpyHostToDevice);

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	cuda_err = cublasSetMatrix(nb_row, nb_col,
				   sizeof(float),
				   test->arrayC[batch_iter],
				   test->ldc[first_index],
				   test->arrayC_h_d[batch_iter],
				   test->ldc[first_index]);
	bblas_gpu_check(cuda_err);
    }
#endif
}

/**
 * Copy the array C from GPU memory after computation to check the relative error etc.
 **/
void bblas_sget_Cfinal(bblas_stest_t *test)
{
    /*Local variables */
    int ldc, N, batch_count;
    int batch_iter, first_index;

    first_index = 0;
    ldc = test->ldc[first_index];
    N   = test->N[first_index];
    batch_count = test->batch_count;

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	test->device_result[batch_iter] = (float *) malloc(ldc*N* sizeof(float ));
#if defined(BBLAS_WITH_CUBLAS)
        cublasGetMatrix(ldc, N, sizeof(float),
			test->arrayC_h_d[batch_iter],
			ldc, test->device_result[batch_iter], ldc);
#endif
    }
}

/**
 * Copy the array B from GPU memory after computation to check the relative error etc.
 **/
void bblas_sget_Bfinal(bblas_stest_t *test)
{
#if defined(BBLAS_WITH_CUBLAS)
    /*Local variables */
    int ldb, N, batch_count;
    int batch_iter, first_index;

    first_index = 0;
    ldb = test->ldb[first_index];
    N   = test->N[first_index];
    batch_count = test->batch_count;

    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
    {
	test->device_result[batch_iter] = (float *) malloc(ldb*N* sizeof(float ));
	cublasGetMatrix(ldb, N, sizeof(float),
			test->arrayB_h_d[batch_iter], ldb,
			test->device_result[batch_iter], ldb);
    }
#endif
}

/**
 * Free all GPU memory ready for the next computation.
 **/
void bblas_scudaFree( bblas_stest_t *test )
{
#if defined(BBLAS_WITH_CUBLAS)

  if ((test->routine == BBLAS_GEMM) || (test->routine == BBLAS_HERK) ||
	(test->routine == BBLAS_TRSM))
    {
	for (int batch_iter =0; batch_iter < test->batch_count; batch_iter++)
	{
	    cudaFree(test->arrayA_h_d[batch_iter]);
	}
	free(test->arrayA_h_d);
	cudaFree(test->arrayA_device);
    }

    if ((test->routine == BBLAS_GEMM) || (test->routine == BBLAS_TRSM))
    {
	for (int batch_iter =0; batch_iter < test->batch_count; batch_iter++)
	{
	    cudaFree(test->arrayB_h_d[batch_iter]);
	}
	free(test->arrayB_h_d);
	cudaFree(test->arrayB_device);
    }

    if ((test->routine == BBLAS_GEMM) ||(test->routine == BBLAS_HERK ))
    {
	for (int batch_iter =0; batch_iter < test->batch_count; batch_iter++)
	{
	    cudaFree(test->arrayC_h_d[batch_iter]);
	}
	free(test->arrayC_h_d);
	cudaFree(test->arrayC_device);
    }
#endif
}


/**
 * Wrapper for cublasSgemmBatched.
 *
 * Calls cublasSgemmBatched with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS
 * with the additional first parameter <tt>bblas_stest_t *test</tt> which is required
 * to obtain the CUDA handle to the GPU.
 * The arguments are converted to a cublas function call within this function.
 **/

void cublas_sgemm_batch( const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
        const int *M, const int *N, const int *K,
        const float *alpha,
        const float **A, const int *lda,
        const float **B, const int *ldb,
        const float *beta,
        float **C, const int *ldc,
        const int batch_count, enum BBLAS_OPTS batch_opts, int *info,
        bblas_stest_t *test )
{
#if defined(BBLAS_WITH_CUBLAS)

    cublasStatus_t stat;

    stat = cublasSgemmBatched(test->cublas_handle,
                      op_bb2cu(transA[0]), op_bb2cu(transB[0]),
					  M[0], N[0], K[0],
					  (const float*)alpha,
					  (const float**)A, lda[0],
					  (const float**)B, ldb[0],
					  (const float*)beta,
					  (float**)C, ldc[0],
					  batch_count);

    cublasGetStream(test->cublas_handle, &test->stream);
	cudaStreamSynchronize(test->stream);

    if(stat != CUBLAS_STATUS_SUCCESS)
    {
	    printf("FAILURE IN CUBLAS SGEMM BATCHED CALL:");

        if(stat == CUBLAS_STATUS_NOT_INITIALIZED)
	    {
			printf("NOT INITIALIZED\n");
	    }
	    else if(stat == CUBLAS_STATUS_INVALID_VALUE)
	    {
    		printf("INVALID VALUE\n");
	    }
	    else if(stat == CUBLAS_STATUS_MAPPING_ERROR)
	    {
	    	printf("MAPPING_ERROR\n");
		}
    }

#endif
}


/**
 * Wrapper for cublasStrsmBatched.
 *
 * Calls cublasStrsmBatched with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS
 * with the additional first parameter <tt>bblas_stest_t *test</tt> which is required
 * to obtain the CUDA handle to the GPU.
 * The arguments are converted to a cublas function call within this function.
 **/

void cublas_strsm_batch( const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
	const int  *M, const int *N,
    const float *alpha,
    const float **A, const int *lda,
    float **B, const int *ldb,
    const int batch_count, enum BBLAS_OPTS batch_opts, int *info,
    bblas_stest_t *test)
{
#if defined(BBLAS_WITH_CUBLAS)

  cublasStatus_t stat;

  stat = cublasStrsmBatched(test->cublas_handle,
			    op_bb2cu(side[0]), op_bb2cu(uplo[0]),
			    op_bb2cu(transA[0]), op_bb2cu(diag[0]),
			    M[0], N[0],
			    (const float*)alpha,
			    (const float**)A, lda[0],
			    (float**)B, ldb[0],
			    batch_count);

  cublasGetStream(test->cublas_handle, &test->stream);
  cudaStreamSynchronize(test->stream);

  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("FAILURE IN CUBLAS STRSM BATCHED CALL:");

      if(stat == CUBLAS_STATUS_NOT_INITIALIZED)
	{
	  printf("NOT INITIALIZED\n");
	}
      else if(stat == CUBLAS_STATUS_INVALID_VALUE)
	{
	  printf("INVALID VALUE\n");
	}
      else if(stat == CUBLAS_STATUS_MAPPING_ERROR)
	{
	  printf("MAPPING_ERROR\n");
	}
    }
#endif
}

#undef REAL
