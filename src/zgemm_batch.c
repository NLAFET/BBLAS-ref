/**
 * @file zgemm_batch.c
 *
 *  @brief BBLAS zgemm_batch double _Complex routine.
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

#include <cblas.h>
#include "bblas_z.h"

#define COMPLEX

/**
    Purpose
    -------
    <b>zgemm_batch</b> is a batch version of zgemm.
	It performs the matrix-matrix operations

	arrayC[i] = alpha[i]*op( arrayA[i] )*op( arrayB[i] ) + beta[i]*arrayC[i],

	where op( X ) is one of

	op( X ) = X      or
	op( X ) = X**T   or
	op( X ) = X**H,

	alpha[i] and beta[i] are scalars, and arrayA[i], arrayB[i] and C are matrices, with
	op( arrayA[i] ) an m by k matrix, op( arrayB[i] ) a k by n matrix and arrayC[i]
	an m by n matrix.


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
	the values of transA[0], transB[0], M[0], N[0],
	K[0], alpha[0], beta[0], lda[0], ldb[0], and ldc[0] are used for all computations.


	Parameters
	----------
	@param[in]
	transA   Array of <tt>enum BBLAS_TRANS</tt>.
	         On entry, transA[i] specifies the form of op( arrayA[i] ) to be used in
			 the matrix multiplication as follows:
      -     = BblasNoTrans:    op( arrayA[i] ) = arrayA[i].
      -     = BblasTrans:      op( arrayA[i] ) = arrayA[i]**T.
      -     = BblasConjTrans:  op( arrayA[i] ) = arrayA[i]**H.

	@param[in]
	transB   Array of <tt>enum BBLAS_TRANS</tt>.
	         On entry, transB[i] specifies the form of op( arrayB[i] ) to be used in
             the matrix multiplication as follows:
      -     = BblasNoTrans:    op( arrayB[i] ) = arrayB[i].
      -     = BblasTrans:      op( arrayB[i] ) = arrayB[i]**T.
      -     = BblasConjTrans:  op( arrayB[i] ) = arrayB[i]**H.

	@param[in]
	M       Array of <tt>int</tt>.
            Each element M[i] specifies the number of rows of the matrix
            op( arrayA[i] ) and of the matrix arrayC[i]. M[i] must be greater than zero.

	@param[in]
	N       Array of <tt>int</tt>.
            Each element N[i] specifies the number of columns of the matrix
            op( arrayB[i] ) and the number of columns of the matrix arrayC[i].
 		    N[i] must be greater than zero.

    @param[in]
    K      Array of <tt>int</tt>.
	       Each element K[i] specifies the number of columns of the matrix
           op( arrayA[i] ) and the number of rows of the matrix op( arrayB[i] ).
		   K[i] must be greater than zero.

    @param[in]
    alpha   Array of <tt>complex_16</tt>.

    @param[in]
    arrayA  Array of pointers.
            Each element arrayA[i] is a pointer to a COMPLEX_16 matrix of
 		    dimension lda[i] by Ka[i], where Ka[i] is K[i] when transA[i] = BblasNoTrans,
 		    and is M[i] otherwise.
            When using transA[i] = BblasNoTrans the leading M[i] by K[i]
            part of arrayA[i] must contain the matrix elements, otherwise
            the leading  K[i] by M[i] part of arrayA[i] must contain the
            matrix elements.

    @param[in]
    lda     Array of <tt>int</tt>.
            Each element lda[i] specifies the first dimension of arrayA[i] as declared
            in the calling (sub) program. When transA[i] = BblasNoTrans then
            lda[i] must be at least max( 1, M[i] ), otherwise lda[i] must be at
            least max( 1, K[i] ).

    @param[in]
    arrayB  Array of pointers.
  		    Each element arrayB[i] is a pointer to a COMPLEX_16 matrix of
  		    dimension ldb[i] by Kb[i], where Kb[i] is N[i] when transB[i] = BblasNoTrans,
  		    and is K[i] otherwise.
            When using transB[i] = BblasNoTrans the leading K[i] by N[i]
            part of arrayB[i] must contain the matrix elements, otherwise
            the leading N[i] by K[i] part of arrayB[i] must contain the
            matrix elements.


    @param[in]
    ldb     Array of <tt>int</tt>.
  		    Each element ldb[i] specifies the first dimension of arrayB[i] as declared
            in the calling (sub) program. When transB[i] = BblasNoTrans then
            ldb[i] must be at least max( 1, K[i] ), otherwise ldb[i] must be at
            least max( 1, N[i] ).

    @param[in]
    beta    Array of <tt>complex_16</tt>.
            When beta[i] is set to zero arrayC[i] need not be set on input.

    @param[in,out]
    arrayC  Array of pointers.
	    Each element arrayC[i] is a pointer to a COMPLEX_16 matrix of
      	    dimension ldc[i] by N[i].
  	    Before entry, the leading M[i] by N[i] part of the arrayC[i] must
            contain a matrix C, except when beta is zero, in which
            case C need not be set on entry.
            On exit, the matrix arrayC[i] is overwritten by the M[i] by N[i] matrix
            ( alpha[i]*op( arrayA[i] )*op( arrayB[i] ) + beta[i]*arrayC[i] ).

    @param[in]
    ldc     Array of <tt>int</tt>.
  	    Each element ldc[i] specifies the first dimension of arrayC[i] as declared
  	    in the calling (sub) program. The value ldc[i] must be at least
  	    max( 1, M[i] )

    @param[in]
    batch_count  <tt>int</tt>
                 The number of matrices to operate on.

    @param[in]
    batch_opts   <tt>enum BBLAS_OPTS</tt>
                 One of BBLAS_FIXED or BBLAS_VARIABLE depending upon the type of
				 batch operation required.

    @param[in,out]
    info    Array of <tt>int</tt>.
            Each element info[i] is the error return code of the ith zgemm in the batch,
            these need not be set on entry.
  	    The error codes can be found in bblas_macros.h.
 **/
void bblas_zgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB,
    const int *ldb, const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc, const int batch_count,
    enum BBLAS_OPTS batch_opts, int *info)
{
    /* Local variables */
    int first_index = 0;
    int batch_iter;
    int LDA,  LDB;
    char func_name[15] = "zgemm_batch";

	/* Check input arguments */
    if (batch_count < 0)
    {
		xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
    }
    if (batch_opts == BBLAS_FIXED)
    {
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
        if ((transB[first_index] != BblasNoTrans) &&
            (transB[first_index] != BblasTrans) &&
            (transB[first_index] != BblasConjTrans))
        {
            xerbla_batch(func_name, BBLAS_ERR_TRANSB, first_index);
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] = BBLAS_ERR_TRANSB;
            }
            return;
        }
        if ( transA[first_index] == BblasNoTrans )
        {
            LDA = M[first_index];
        } else
        {
            LDA = K[first_index];
        }
        if ( transB[first_index] == BblasNoTrans )
        {
            LDB = K[first_index];
        } else
        {
            LDB = N[first_index];
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
        if (N[first_index] < 0) {
            xerbla_batch(func_name, BBLAS_ERR_N, first_index);
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] = BBLAS_ERR_N;
            }
            return;
        }
        if (K[first_index] < 0) {
            xerbla_batch(func_name, BBLAS_ERR_K, first_index);
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] = BBLAS_ERR_K;
            }
            return;
        }
        if (lda[first_index] < max(1, LDA)) {
            xerbla_batch(func_name, BBLAS_ERR_LDA, first_index);
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] =  BBLAS_ERR_LDA;
            }
            return;
        }
        if (ldb[first_index] < max(1, LDB)) {
            xerbla_batch(func_name, BBLAS_ERR_LDB, first_index);
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] = BBLAS_ERR_LDB;
            }
            return;
        }
        if (ldc[first_index] < max(1, M[first_index])) {
            xerbla_batch(func_name, BBLAS_ERR_LDC, first_index);
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] = BBLAS_ERR_LDC;
            }
            return;
        }
        /* Skip subproblems where nothing needs to be done */
        if (M[first_index] == 0 || N[first_index] == 0     ||
	        ((alpha[first_index] == (BBLAS_Complex64_t)0.0 || K[first_index] == 0) &&
	        beta[first_index] == (BBLAS_Complex64_t)1.0 ))
        {
            for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
            {
                info[batch_iter] =  BBLAS_SUCCESS;
            }
            return;
        }
        for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
        {
            /* Call to cblas_zgemm */
            cblas_zgemm(
            BblasColMajor,
            transA[first_index],
            transB[first_index],
            M[first_index],
            N[first_index],
            K[first_index],
            CBLAS_SADDR(alpha[first_index]),
            arrayA[batch_iter],
            lda[first_index],
            arrayB[batch_iter],
            ldb[first_index],
            CBLAS_SADDR(beta[first_index]),
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
            if ((transA[batch_iter] != BblasNoTrans) &&
                (transA[batch_iter] != BblasTrans) &&
                (transA[batch_iter] != BblasConjTrans))
            {
                xerbla_batch(func_name, BBLAS_ERR_TRANSA, batch_iter);
                info[batch_iter]  = BBLAS_ERR_TRANSA;
                continue;
            }
            if ((transB[batch_iter] != BblasNoTrans) &&
                (transB[batch_iter] != BblasTrans) &&
                (transB[batch_iter] != BblasConjTrans))
            {
                xerbla_batch(func_name, BBLAS_ERR_TRANSB, batch_iter);
                info[batch_iter] = BBLAS_ERR_TRANSB;
                continue;
            }
            if (transA[batch_iter] == BblasNoTrans)
            {
                LDA = M[batch_iter];
            } else
            {
                LDA = K[batch_iter];
            }
            if (transB[batch_iter] == BblasNoTrans)
            {
                LDB = K[batch_iter];
            } else
            {
                LDB = N[batch_iter];
            }
            if (M[batch_iter] < 0)
            {
                xerbla_batch(func_name, BBLAS_ERR_M, batch_iter);
                info[batch_iter] = BBLAS_ERR_M;
                continue;
            }
            if (N[batch_iter] < 0) {
                xerbla_batch(func_name, BBLAS_ERR_N, batch_iter);
                info[batch_iter] = BBLAS_ERR_N;
                continue;
            }
            if (K[batch_iter] < 0) {
                xerbla_batch(func_name, BBLAS_ERR_K, batch_iter);
                info[batch_iter] = BBLAS_ERR_K;
                continue;
            }
            if (lda[batch_iter] < max(1, LDA)) {
                xerbla_batch(func_name, BBLAS_ERR_LDA, batch_iter);
                info[batch_iter] =  BBLAS_ERR_LDA;
                continue;
            }
            if (ldb[batch_iter] < max(1, LDB)) {
                xerbla_batch(func_name, BBLAS_ERR_LDB, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDB;
                continue;
            }
            if (ldc[batch_iter] < max(1, M[batch_iter])) {
                xerbla_batch(func_name, BBLAS_ERR_LDC, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDC;
                continue;
            }
			/* Skip subproblems where nothing needs to be done */
            if (M[batch_iter] == 0 || N[batch_iter] == 0 ||
                ((alpha[batch_iter] == (BBLAS_Complex64_t)0.0 || K[batch_iter] == 0) &&
                beta[batch_iter] == (BBLAS_Complex64_t)1.0))
            {
                info[batch_iter] =  BBLAS_SUCCESS;
                continue;
            }
			/* Call to cblas_zgemm */
            cblas_zgemm(
            BblasColMajor,
            transA[batch_iter],
            transB[batch_iter],
            M[batch_iter],
            N[batch_iter],
            K[batch_iter],
            CBLAS_SADDR(alpha[batch_iter]),
            arrayA[batch_iter],
            lda[batch_iter],
            arrayB[batch_iter],
            ldb[batch_iter],
            CBLAS_SADDR(beta[batch_iter]),
            arrayC[batch_iter],
            ldc[batch_iter]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
        }
    } else
    {
		/* Error in batch_opts */
        xerbla_batch(func_name, BBLAS_ERR_BATCH_OPTS, -1);
    }
}
#undef COMPLEX
