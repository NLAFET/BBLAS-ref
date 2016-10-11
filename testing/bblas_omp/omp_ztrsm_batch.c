/**
 * @file omp_ztrsm_batch.c
 *
 * @brief BBLAS omp_ztrsm_batch double _Complex routine.
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
 * @precisions normal z -> c d s
 **/
#endif

#include<cblas.h>
#include "bblas_omp.h"
#include "bblas.h"
#include <omp.h>
#define COMPLEX

/**
    Purpose
    -------
    <b>ztrsm_batch</b> is an OpenMP version of ztrsm_batch.
    It solves for X in one of the matrix equations

    op( arrayA[i] )*X =  alpha*arrayB[i], or
	X*op( arrayA[i] ) = alpha[i]*arrayB[i],

	where op( X ) is one of
    -    op( X ) = X
    or
    -    op( X ) = X**T
    or
    -    op( X ) = X**H,

    alpha[i] is a scalar, X and B are M[i] by N[i] matrices,
	and arrayA[i] is a unit or non-unit, upper or lower triangular matrix.
    The solution matrix X overwrites arrayB[i] on exit.


	Fixed and Variable Batch Operations
	-----------------------------------
	Two types of batch operation are supported depending upon the value of batch_opts.

	When <tt>batch_opts = BBLAS_VARIABLE</tt>
	- all parameters that are arrays must have length at least batch_count.
	- all parameters that are arrays must have all values set.

	When <tt>batch_opts = BBLAS_FIXED</tt>
	- all parameters that are arrays (except for arrayA, arrayB, and info)
	must have length at least one.
	- all parameters that are arrays (except for arrayA, arrayB, and info)
	need only to have their first value set.

	This means that for a <tt>BBLAS_FIXED</tt> batch,
	the values of side[0], uplo[0], transA[0], diag[0], M[0], N[0],
	alpha[0], lda[0], and ldb[0] are used for all computations.


    Parameters
    ----------
    @param[in]
    side    Array of <tt>enum BBLAS_SIDE</tt>.
            Each element side[i] specifies whether op( arrayA[i] )
            appears on the left or right side of the operation as follows:
      -     = 'BblasLeft'  op( arrayA[i] )*X = alpha[i]*arrayB[i].
      -     = 'BblasRight' X*op( arrayA[i] ) = alpha[i]*arrayB[i].

    @param[in]
    uplo    Array of <tt>enum BBLAS_UPLO</tt>.
	        On entry, uplo[i] specifies whether the matrix arrayA[i]
            is upper or lower triangular as follows:
      -     = 'BblasUpper' arrayA[i] is an upper triangular matrix.
      -     = 'BblasLower' arrayA[i] is a lower triangular matrix.

    @param[in]
    transA  Array of <tt>enum BBLAS_TRANS</tt>.
            On entry, trans[i] specifies the form of op( arrayA[i] ) to be
            used in the operation as follows:
      -     = 'BblasNoTrans'     op( arrayA[i] ) = arrayA[i].
      -     = 'BblasTrans'       op( arrayA[i] ) = arrayA[i]**T.
      -     = 'BblasConjTrans'   op( arrayA[i] ) = arrayA'[i]**H.

    @param[in]
	diag   - Array of <tt>enum BBLAS_DIAG</tt>.
             On entry, diag[i] specifies whether or not arrayA[i] is unit
             triangular as follows:
      -      = 'BblasUnit'    arrayA[i] is assumed to be unit triangular.
      -      = 'BblasNonUnit' arrayA[i] is not assumed to be unit triangular.

    @param[in]
	M       Array of <tt>int</tt>.
			Each element M[i] specifies the number of rows of the matrix arrayB[i].
			M[i] must be greater than zero.

    @param[in]
	N       Array of <tt>int</tt>.
            Each element N[i] specifies the number of columns of the matrix arrayB[i].
			N[i] must be greater than zero.

    @param[in]
    alpha   Array of COMPLEX_16
	        When alpha[i] is set to zero arrayA[i] is not referenced and arrayB[i] need
            not be set before entry.

    @param[in]
    arrayA   Array of pointers.
	         Each element arrayA[i] is a pointer to a COMPLEX_16 matrix of
			 dimension lda[i] by Ka[i],
			 where Ka[i] = M[i] when side[i] = BblasLeft and is N[i] otherwise.
			 When using side[i] = BblasLeft the M[i] by M[i] part of arrayA[i]
			 must contain the triangular matrix:
			 when uplo[i] = BblasUpper, the upper triangular part of arrayA[i]
			 must contain the matrix whilst the strictly lower triangular part is not used;
			 similarly when uplo[i] = BblasLower, the lower triangular part of arrayA[i]
			 must contain the matrix whilst the strictly upper triangular part is not used.
			 When using side[i] = BblasRight the N[i] by N[i] part of arrayA[i] must
			 contain the symmetric matrix:
			 when uplo[i] = BblasUpper, the upper triangular part of arrayA[i]
			 must contain the matrix whilst the strictly lower triangular part is not used;
			 similarly when uplo[i] = BblasLower, the lower triangular part of arrayA[i]
			 must contain the matrix whilst the strictly upper triangular part is not used.
			 Note that when diag = BblasUnit the diagonal elements of arrayA[i] are
			 not used either, they are assumed to be equal to one.

    @param[in]
    lda     Array of <tt>int</tt>.
            On entry, lda[i] specifies the first dimension of arrayA[i] as declared
            in the calling (sub) program. When  side[i] = BblasLeft
            then lda[i] must be at least  max( 1, M[i] ),
            otherwise lda[i] must be at least max( 1, N[i] ).

    @param[in,out]
    arrayB   Array of pointers.
             Each element arrayB[i] is a pointer to a COMPLEX_16 matrix of
			 dimension ldb[i] by N[i].
			 The leading M[i] by N[i] part of arrayB[i] must contain the matrix elements.
			 On exit is arrayB[i] overwritten by the solution matrix X.

    @param[in]
	ldb     Array of <tt>int</tt>.
			Each element ldb[i] specifies the first dimension of arrayB[i] as declared
            in the calling (sub) program. Each element ldb[i] must be at least max( 1, M[i] ).

    @param[in]
    batch_count  <tt>int</tt>
                 The number of matrices to operate on.

	@param[in]
    batch_opts   <tt>enum BBLAS_OPTS</tt>
                 One of BBLAS_FIXED or BBLAS_VARIABLE depending upon the type of
				 batch operation required.

    @param[out]
    info   Array of <tt>int</tt>.
	       Each element info[i] is the error return code of the ith ztrsm in the batch,
		   these need not be set on entry.
		   The error codes can be found in bblas_macros.h.
**/

void omp_ztrsm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA,  const enum BBLAS_DIAG *diag,
    const int *M,  const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    BBLAS_Complex64_t **arrayB, const int *ldb,
    const int batch_count, enum BBLAS_OPTS batch_opts, int *info)
{
/*Local variables  */
    int first_index = 0;
    int batch_iter;
    int LDA;
    char func_name[15] = "ztrsm_batch";

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
        if ((transA[first_index] != BblasNoTrans) &&
            (transA[first_index] != BblasTrans) &&
            (transA[first_index] != BblasConjTrans))
        {
            xerbla_batch(func_name, BBLAS_ERR_TRANSA, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
                info[batch_iter]  = BBLAS_ERR_TRANSA;
            }
            return;
        }
        if ((diag[first_index] != BblasNonUnit) &&
            (diag[first_index] != BblasUnit))
        {
            xerbla_batch(func_name, BBLAS_ERR_DIAG, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
                info[batch_iter]  = BBLAS_ERR_DIAG;
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
	if (lda[first_index] < max(1, LDA))
        {
	    xerbla_batch(func_name, BBLAS_ERR_LDA, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] =  BBLAS_ERR_LDA;
            }
	    return;
	}
	if (ldb[first_index] < max(1, M[first_index])) {
	    xerbla_batch(func_name, BBLAS_ERR_LDB, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_LDB;
            }
	    return;
	}
	/* particular case */
	if (min(M[first_index], N[first_index]) == 0)
	{
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] =  BBLAS_SUCCESS;
            }
	    return;
	}

#pragma omp parallel for private(batch_iter)
	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	{
            /*Call to cblas_ztrsm */
            cblas_ztrsm(
		BblasColMajor,
		side[first_index],
		uplo[first_index],
		transA[first_index],
		diag[first_index],
		M[first_index],
		N[first_index],
		CBLAS_SADDR(alpha[first_index]),
		arrayA[batch_iter],
		lda[first_index],
		arrayB[batch_iter],
		ldb[first_index]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
	} /*END FIXED SIZE FOR LOOP */
    }else if (batch_opts ==  BBLAS_VARIABLE)
    {
#pragma omp parallel for private(batch_iter,LDA)
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
            if ((transA[batch_iter] != BblasNoTrans) &&
                (transA[batch_iter] != BblasTrans) &&
                (transA[batch_iter] != BblasConjTrans))
            {
                xerbla_batch(func_name, BBLAS_ERR_TRANSA, batch_iter);
                info[batch_iter]  = BBLAS_ERR_TRANSA;
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
            if (lda[batch_iter] < max(1, LDA))
            {
                xerbla_batch(func_name, BBLAS_ERR_LDA, batch_iter);
                info[batch_iter] =  BBLAS_ERR_LDA;
                continue;
            }
            if (ldb[batch_iter] < max(1, M[batch_iter]))
            {
                xerbla_batch(func_name, BBLAS_ERR_LDC, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDC;
                continue;
            }
            /* particular case */
	    if (min(M[batch_iter], N[batch_iter]) == 0)
            {
                info[batch_iter] =  BBLAS_SUCCESS;
                continue;
            }
            cblas_ztrsm(
		BblasColMajor,
		side[batch_iter],
		uplo[batch_iter],
		transA[batch_iter],
		diag[batch_iter],
		M[batch_iter],
		N[batch_iter],
		CBLAS_SADDR(alpha[batch_iter]),
		arrayA[batch_iter],
		lda[batch_iter],
		arrayB[batch_iter],
		ldb[batch_iter]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
        }
    } else
    {
	xerbla_batch(func_name, BBLAS_ERR_BATCH_OPTS, -1);
    }
}
#undef COMPLEX
