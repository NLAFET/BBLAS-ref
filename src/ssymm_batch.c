/**
 * @file ssymm_batch.c
 *
 * @brief BBLAS ssymm_batch float routine.
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
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from ../src/zsymm_batch.c normal z -> s, Mon Jun  6 09:44:13 2016
 **/
#endif

#include<cblas.h>
#include "bblas_s.h"

#define REAL

/**
    Purpose
    -------
    <b>ssymm_batch</b> is a batch version of ssymm.
    It performs one of the matrix-matrix operations

	arrayC[i] = alpha[i]*arrayA[i]*arrayB[i] + beta[i]*arrayC[i], or
    arrayC[i] = alpha[i]*arrayB[i]*arrayA[i] + beta[i]*arrayC[i],

    where alpha[i] and beta[i] are scalars, arrayA[i] is a symmetric matrix
    and arrayB[i] and arrayC[i] are M[i] by N[i] matrices.


	Fixed and Variable Batch Operations
	-----------------------------------
	Two types of batch operation are supported depending upon the value of batch_opts.

	When <tt>batch_opts = BBLAS_VARIABLE</tt>
	- all parameters that are arrays must have length at least batch_count.
	- all parameters that are arrays must have all values set.

	When <tt>batch_opts = BBLAS_FIXED</tt>
	- all parameters that are arrays (except for arrayA, arrayB, arrayC, and info)
	must have length at least one.
	- all parameters that are arrays (except for arrayA, arrayB, arrayC, and info)
	need only to have their first value set.

	This means that for a <tt>BBLAS_FIXED</tt> batch,
	the values of side[0], uplo[0], M[0], N[0],
	alpha[0], beta[0], lda[0], ldb[0], and ldc[0] are used for all computations.


    Parameters
    ----------
    @param[in]
    side    Array of <tt>enum BBLAS_SIDE</tt>.
            Each element side[i] specifies whether the symmetric
            matrix arrayA[i] appears on the left or right side of the
            operation as follows:
      -     = 'BblasLeft'  arrayC[i] = alpha[i]*arrayA[i]*arrayB[i] + beta[i]*arrayC[i].
      -     = 'BblasRight' arrayC[i] = alpha[i]*arrayB[i]*arrayA[i] + beta[i]*arrayC[i].

    @param[in]
    uplo    Array of <tt>enum BBLAS_UPLO</tt>.
            On entry, uplo[i] specifies whether the upper or
            lower triangular part of the symmetric matrix
            arrayA[i] is to be referenced as follows:
      -     = 'BblasUpper' Only the upper triangular part of
              arrayA[i] is to be referenced.
      -     = 'BblasLower' Only the lower triangular part of
	           arrayA[i] is to be referenced.

    @param[in]
    M       Array of <tt>int</tt>.
			Each element M[i] specifies the number of rows of the matrix arrayC[i].
			M[i] must be greater than zero.

    @param[in]
    N       Array of <tt>int</tt>.
            Each element N[i] specifies the number of columns of the matrix arrayC[i].
			N[i] must be greater than zero.

    @param[in]
    alpha   Array of <tt>real_16</tt>.

    @param[in]
    arrayA   Array of pointers.
	         Each element arrayA[i] is a pointer to a REAL matrix of
			 dimension lda[i] by Ka[i],
			 where Ka[i] = M[i] when side[i] = BblasLeft and is N[i] otherwise.
			 When using side[i] = BblasLeft the M[i] by M[i] part of arrayA[i]
			 must contain the symmetric matrix:
			 when uplo[i] = BblasUpper, the upper triangular part of arrayA[i]
			 must contain the upper triangular part of the symmetric matrix whilst
			 the strictly lower triangular part is not used;
			 similarly when uplo[i] = BblasLower, the lower triangular part of arrayA[i]
			 must contain the lower triangular part of the symmetric matrix
			 whilst the strictly upper triangular part is not used.
			 When using side[i] = BblasRight the N[i] by N[i] part of arrayA[i] must
			 contain the symmetric matrix:
			 when uplo[i] = BblasUpper, the upper triangular part of arrayA[i]
			 must contain the upper triangular part of the symmetric matrix whilst
			 the strictly lower triangular part is not used;
			 similarly when uplo[i] = BblasLower, the lower triangular part of arrayA[i]
			 must contain the lower triangular part of the symmetric matrix
			 whilst the strictly upper triangular part is not used.

    @param[in]
    lda     Array of <tt>int</tt>.
            On entry, lda[i] specifies the first dimension of arrayA[i] as declared
            in the calling (sub) program. When  side[i] = BblasLeft
            then lda[i] must be at least  max( 1, M[i] ),
            otherwise lda[i] must be at least max( 1, N[i] ).

    @param[in]
    arrayB   Array of pointers.
             Each element arrayB[i] is a pointer to a REAL matrix of
			 dimension ldb[i] by N[i].
			 The leading M[i] by N[i] part of arrayB[i] must contain the matrix elements.

    @param[in]
	ldb     Array of <tt>int</tt>.
			Each element ldb[i] specifies the first dimension of arrayB[i] as declared
            in the calling (sub) program. Each element ldb[i] must be at least max( 1, M[i] ).

	@param[in]
	beta    Array of <tt>real_16</tt>.
	        When beta[i] is set to zero arrayC[i] need not be set on input.

    @param[in,out]
	arrayC  Array of pointers.
			Each element arrayC[i] is a pointer to a REAL matrix of
  		    dimension ldc[i] by N[i].
  		    Before entry, the leading M[i] by N[i] part of the arrayC[i] must
            contain a matrix C, except when beta is zero, in which
            case C need not be set on entry.
            On exit, the matrix arrayC[i] is overwritten by the M[i] by N[i] matrix output.

	@param[in]
    ldc     Array of <tt>int</tt>.
  		    Each element ldc[i] specifies the first dimension of arrayC[i] as declared
  		    in the calling (sub) program. The value ldc[i] must be at least
  		    max( 1, M[i] ).

    @param[in]
    batch_count  <tt>int</tt>
                 The number of matrices to operate on.

    @param[in]
    batch_opts   <tt>enum BBLAS_OPTS</tt>
                 One of BBLAS_FIXED or BBLAS_VARIABLE depending upon the type of
				 batch operation required.

    @param[out]
    info   Array of <tt>int</tt>.
	        Each element info[i] is the error return code of the ith zymm in the batch,
			these need not be set on entry.
  		    The error codes can be found in bblas_macros.h.
**/

void bblas_ssymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M,  const int *N, const float *alpha,
    const float **arrayA, const int *lda,
    const float **arrayB, const int *ldb,
    const float *beta, float **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info)
{
/* Local variables  */
    int first_index = 0;
    int batch_iter;
    int LDA;
    char func_name[15] = "ssymm_batch";
    /* Check input arguments */
    if (batch_count < 0)
    {
	xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
    }
    if (batch_opts == BBLAS_FIXED)
    {
	if ((side[first_index] != BblasLeft) &&
	    (side[first_index] != BblasRight))
	{
	    xerbla_batch(func_name, BBLAS_ERR_SIDE, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter]  = BBLAS_ERR_SIDE;
            }
	    return;
	}
	if ((uplo[first_index] != BblasUpper) &&
	    (uplo[first_index] != BblasLower))
	{
	    xerbla_batch(func_name, BBLAS_ERR_UPLO, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_UPLO;
            }
	    return;
	}
	if (M[first_index] < 0)
	{
	    xerbla_batch(func_name, BBLAS_ERR_M, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_M;
            }
	    return;
	}
	if (N[first_index] < 0)
        {
	    xerbla_batch(func_name, BBLAS_ERR_N, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_N;
            }
	    return;
	}
	if (side[first_index] == BblasLeft)
        {
	    LDA = M[first_index];
	} else
        {
	    LDA = N[first_index];
	}
	if (lda[first_index] < LDA)
        {
	    xerbla_batch(func_name, BBLAS_ERR_LDA, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_LDA;
            }
	    return;
	}
	if (ldb[first_index] < max(1, M[first_index]))
        {
	    xerbla_batch(func_name, BBLAS_ERR_LDB, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_LDB;
            }
	    return;
	}
	if (ldc[first_index] < max(1, M[first_index]))
        {
	    xerbla_batch(func_name, BBLAS_ERR_LDC, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_LDC;
            }
	    return;
	}
	/* Skip subproblems where nothing needs to be done */
	if (M[first_index] == 0 || N[first_index] == 0 ||
	    (alpha[first_index] == (float)0.0 &&
	     beta[first_index] == (float)1.0))
	{
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] =  BBLAS_SUCCESS;
            }
	    return;
	}
	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
        {
            /* Call to cblas_ssymm */
        cblas_ssymm(
		BblasColMajor,
		side[first_index],
		uplo[first_index],
		M[first_index],
		N[first_index],
		(alpha[first_index]),
		arrayA[batch_iter],
		lda[first_index],
		arrayB[batch_iter],
		ldb[first_index],
		(beta[first_index]),
		arrayC[batch_iter],
		ldc[first_index]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
	} /* END FIXED SIZE FOR LOOP */
    }else if (batch_opts ==  BBLAS_VARIABLE)
    {
        for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
        {
            /* Check input arguments */
            if ((side[batch_iter] != BblasLeft) &&
                (side[batch_iter] != BblasRight))
            {
                xerbla_batch(func_name, BBLAS_ERR_SIDE, batch_iter);
                info[batch_iter]  = BBLAS_ERR_SIDE;
                continue;
            }
            if ((uplo[batch_iter] != BblasUpper) &&
                (uplo[batch_iter] != BblasLower))
            {
                xerbla_batch(func_name, BBLAS_ERR_UPLO, batch_iter);
                info[batch_iter] = BBLAS_ERR_UPLO;
                continue;
            }
            if (M[batch_iter] < 0)
            {
                xerbla_batch(func_name, BBLAS_ERR_M, batch_iter);
                info[batch_iter] = BBLAS_ERR_M;
                continue;
            }
            if (N[batch_iter] < 0)
            {
                xerbla_batch(func_name, BBLAS_ERR_N, batch_iter);
                info[batch_iter] = BBLAS_ERR_N;
                continue;
            }
            if (side[batch_iter] == BblasLeft)
            {
                LDA = M[batch_iter];
            } else
            {
                LDA = N[batch_iter];
            }
            if (lda[batch_iter] < LDA)
            {
                xerbla_batch(func_name, BBLAS_ERR_LDA, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDA;
                continue;
            }
            if (ldb[batch_iter] < max(1, M[batch_iter]))
            {
                xerbla_batch(func_name, BBLAS_ERR_LDB, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDB;
                continue;
            }
            if (ldc[batch_iter] < max(1, M[batch_iter]))
            {
                xerbla_batch(func_name, BBLAS_ERR_LDC, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDC;
                continue;
            }
            /* Skip subproblems where nothing needs to be done */
            if (M[batch_iter] == 0 || N[batch_iter] == 0 ||
                (alpha[batch_iter] == (float)0.0 &&
		 beta[batch_iter] == (float)1.0))
            {
                info[batch_iter] =  BBLAS_SUCCESS;
                continue;
            }
            cblas_ssymm(
				BblasColMajor,
				side[batch_iter],
				uplo[batch_iter],
				M[batch_iter],
				N[batch_iter],
				(alpha[batch_iter]),
				arrayA[batch_iter],
				lda[batch_iter],
				arrayB[batch_iter],
				ldb[batch_iter],
				(beta[batch_iter]),
				arrayC[batch_iter],
				ldc[batch_iter]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
        }
    }else
    {
		/* Error in batch_opts */
		xerbla_batch(func_name, BBLAS_ERR_BATCH_OPTS, -1);
    }
}
#undef REAL
