/**
 * @file flops.h
 *
 * @brief Flop counts of all batch operations.
 *
 * BBLAS is a software package provided by Univ. of Manchester,
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-04-05
 *
 * Contains numerous macros to return the flop count of all level-3 BLAS operations.
 **/

#ifndef BBLAS_FLOPS_H
#define BBLAS_FLOPS_H

/*
 * Level 3 BLAS flop counts - Intermediate quantities
 */

/************* GEMM ***************/
/** Number of multiplies in GEMM **/
#define FMULS_GEMM(m_, n_, k_) ((m_) * (n_) * (k_))
/** Number of additions in GEMM **/
#define FADDS_GEMM(m_, n_, k_) ((m_) * (n_) * (k_))

/********** SYMM/HEMM *************/
/** Number of multiplies in SYMM **/
#define FMULS_SYMM(side_, m_, n_) ( ( (side_) == BblasLeft ) ? FMULS_GEMM((m_), (m_), (n_)) : FMULS_GEMM((m_), (n_), (n_)) )
/** Number of additions in SYMM **/
#define FADDS_SYMM(side_, m_, n_) ( ( (side_) == BblasLeft ) ? FADDS_GEMM((m_), (m_), (n_)) : FADDS_GEMM((m_), (n_), (n_)) )
/** Number of multiplies in HEMM **/
#define FMULS_HEMM FMULS_SYMM
/** Number of additions in HEMM **/
#define FADDS_HEMM FADDS_SYMM

/********** SYRK/HERM *************/
/** Number of multiplies in SYRK **/
#define FMULS_SYRK(k_, n_) (0.5 * (k_) * (n_) * ((n_)+1))
/** Number of additions in SYRK **/
#define FADDS_SYRK(k_, n_) (0.5 * (k_) * (n_) * ((n_)+1))
/** Number of multiplies in HERK **/
#define FMULS_HERK FMULS_SYRK
/** Number of additions in HERK **/
#define FADDS_HERK FADDS_SYRK

/********** SYR2K/HER2K ***********/
/** Number of multiplies in SYR2K **/
#define FMULS_SYR2K(k_, n_) ((k_) * (n_) * (n_)        )
/** Number of additions in SYR2K **/
#define FADDS_SYR2K(k_, n_) ((k_) * (n_) * (n_) + (n_))
/** Number of multiplies in HER2K **/
#define FMULS_HER2K FMULS_SYR2K
/** Number of additions in HER2K **/
#define FADDS_HER2K FADDS_SYR2K

/************* TRMM ***************/
/** Subfunction: Number of multiplies in TRMM **/
#define FMULS_TRMM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)+1))
/** Subfunction: Number of additions in TRMM **/
#define FADDS_TRMM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)-1))
/** Number of multiplies in TRMM **/
#define FMULS_TRMM(side_, m_, n_) ( ( (side_) == BblasLeft ) ? FMULS_TRMM_2((m_), (n_)) : FMULS_TRMM_2((n_), (m_)) )
/** Number of additions in TRMM **/
#define FADDS_TRMM(side_, m_, n_) ( ( (side_) == BblasLeft ) ? FADDS_TRMM_2((m_), (n_)) : FADDS_TRMM_2((n_), (m_)) )

/************* TRSM ***************/
/** Number of multiplies in TRSM **/
#define FMULS_TRSM FMULS_TRMM
/** Number of additions in TRSM **/
#define FADDS_TRSM FADDS_TRMM


/*
 * Level 3 BLAS flop counts - Full quantities
 */

/** Flops in ZGEMM **/
#define FLOPS_ZGEMM(m_, n_, k_) (6. * FMULS_GEMM((double)(m_), (double)(n_), (double)(k_)) + 2.0 * FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) )
/** Flops in CGEMM **/
#define FLOPS_CGEMM(m_, n_, k_) (6. * FMULS_GEMM((double)(m_), (double)(n_), (double)(k_)) + 2.0 * FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) )
/** Flops in DGEMM **/
#define FLOPS_DGEMM(m_, n_, k_) (     FMULS_GEMM((double)(m_), (double)(n_), (double)(k_)) +       FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) )
/** Flops in SGEMM **/
#define FLOPS_SGEMM(m_, n_, k_) (     FMULS_GEMM((double)(m_), (double)(n_), (double)(k_)) +       FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) )

/** Flops in ZHEMM **/
#define FLOPS_ZHEMM(side_, m_, n_) (6. * FMULS_HEMM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_HEMM(side_, (double)(m_), (double)(n_)) )
/** Flops in CHEMM **/
#define FLOPS_CHEMM(side_, m_, n_) (6. * FMULS_HEMM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_HEMM(side_, (double)(m_), (double)(n_)) )

/** Flops in ZSYMM **/
#define FLOPS_ZSYMM(side_, m_, n_) (6. * FMULS_SYMM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_SYMM(side_, (double)(m_), (double)(n_)) )
/** Flops in CSYMM **/
#define FLOPS_CSYMM(side_, m_, n_) (6. * FMULS_SYMM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_SYMM(side_, (double)(m_), (double)(n_)) )
/** Flops in DSYMM **/
#define FLOPS_DSYMM(side_, m_, n_) (     FMULS_SYMM(side_, (double)(m_), (double)(n_)) +       FADDS_SYMM(side_, (double)(m_), (double)(n_)) )
/** Flops in SSYMM **/
#define FLOPS_SSYMM(side_, m_, n_) (     FMULS_SYMM(side_, (double)(m_), (double)(n_)) +       FADDS_SYMM(side_, (double)(m_), (double)(n_)) )

/** Flops in ZHERK **/
#define FLOPS_ZHERK(k_, n_) (6. * FMULS_HERK((double)(k_), (double)(n_)) + 2.0 * FADDS_HERK((double)(k_), (double)(n_)) )
/** Flops in CHERK **/
#define FLOPS_CHERK(k_, n_) (6. * FMULS_HERK((double)(k_), (double)(n_)) + 2.0 * FADDS_HERK((double)(k_), (double)(n_)) )

/** Flops in ZSYRK **/
#define FLOPS_ZSYRK(k_, n_) (6. * FMULS_SYRK((double)(k_), (double)(n_)) + 2.0 * FADDS_SYRK((double)(k_), (double)(n_)) )
/** Flops in CSYRK **/
#define FLOPS_CSYRK(k_, n_) (6. * FMULS_SYRK((double)(k_), (double)(n_)) + 2.0 * FADDS_SYRK((double)(k_), (double)(n_)) )
/** Flops in DSYRK **/
#define FLOPS_DSYRK(k_, n_) (     FMULS_SYRK((double)(k_), (double)(n_)) +       FADDS_SYRK((double)(k_), (double)(n_)) )
/** Flops in SSYRK **/
#define FLOPS_SSYRK(k_, n_) (     FMULS_SYRK((double)(k_), (double)(n_)) +       FADDS_SYRK((double)(k_), (double)(n_)) )

/** Flops in ZHER2K **/
#define FLOPS_ZHER2K(k_, n_) (6. * FMULS_HER2K((double)(k_), (double)(n_)) + 2.0 * FADDS_HER2K((double)(k_), (double)(n_)) )
/** Flops in CHER2K **/
#define FLOPS_CHER2K(k_, n_) (6. * FMULS_HER2K((double)(k_), (double)(n_)) + 2.0 * FADDS_HER2K((double)(k_), (double)(n_)) )

/** Flops in ZSYR2K **/
#define FLOPS_ZSYR2K(k_, n_) (6. * FMULS_SYR2K((double)(k_), (double)(n_)) + 2.0 * FADDS_SYR2K((double)(k_), (double)(n_)) )
/** Flops in CSYR2K **/
#define FLOPS_CSYR2K(k_, n_) (6. * FMULS_SYR2K((double)(k_), (double)(n_)) + 2.0 * FADDS_SYR2K((double)(k_), (double)(n_)) )
/** Flops in DSYR2K **/
#define FLOPS_DSYR2K(k_, n_) (     FMULS_SYR2K((double)(k_), (double)(n_)) +       FADDS_SYR2K((double)(k_), (double)(n_)) )
/** Flops in SSYR2K **/
#define FLOPS_SSYR2K(k_, n_) (     FMULS_SYR2K((double)(k_), (double)(n_)) +       FADDS_SYR2K((double)(k_), (double)(n_)) )

/** Flops in ZTRMM **/
#define FLOPS_ZTRMM(side_, m_, n_) (6. * FMULS_TRMM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_TRMM(side_, (double)(m_), (double)(n_)) )
/** Flops in CTRMM **/
#define FLOPS_CTRMM(side_, m_, n_) (6. * FMULS_TRMM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_TRMM(side_, (double)(m_), (double)(n_)) )
/** Flops in DTRMM **/
#define FLOPS_DTRMM(side_, m_, n_) (     FMULS_TRMM(side_, (double)(m_), (double)(n_)) +       FADDS_TRMM(side_, (double)(m_), (double)(n_)) )
/** Flops in STRMM **/
#define FLOPS_STRMM(side_, m_, n_) (     FMULS_TRMM(side_, (double)(m_), (double)(n_)) +       FADDS_TRMM(side_, (double)(m_), (double)(n_)) )

/** Flops in ZTRSM **/
#define FLOPS_ZTRSM(side_, m_, n_) (6. * FMULS_TRSM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_TRSM(side_, (double)(m_), (double)(n_)) )
/** Flops in CTRSM **/
#define FLOPS_CTRSM(side_, m_, n_) (6. * FMULS_TRSM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_TRSM(side_, (double)(m_), (double)(n_)) )
/** Flops in DTRSM **/
#define FLOPS_DTRSM(side_, m_, n_) (     FMULS_TRSM(side_, (double)(m_), (double)(n_)) +       FADDS_TRSM(side_, (double)(m_), (double)(n_)) )
/** Flops in STRSM **/
#define FLOPS_STRSM(side_, m_, n_) (     FMULS_TRSM(side_, (double)(m_), (double)(n_)) +       FADDS_TRSM(side_, (double)(m_), (double)(n_)) )

#endif /* BBLAS_FLOPS_H */
