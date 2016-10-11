include ../make.inc
include ../Makefile.gen

BBLAS_SRC_LIST   = zgemm_batch.c zhemm_batch.c  \
		 zher2k_batch.c zherk_batch.c \
		 zsymm_batch.c zsyr2k_batch.c \
		 zsyrk_batch.c ztrmm_batch.c  \
		 ztrsm_batch.c 

BBLAS_COMMON_LIST = auxiliary.c xerbla_batch.c \

BBLAS_SRC=$(addprefix $(BBLAS_SRC_DIR)/, $(BBLAS_SRC_LIST))
BBLAS_COMMON=$(addprefix $(BBLAS_SRC_DIR)/, $(BBLAS_COMMON_LIST))

OMP_SRC_LIST	= omp_zgemm_batch.c omp_zhemm_batch.c\
		  omp_zher2k_batch.c omp_zherk_batch.c \
		  omp_zsymm_batch.c  omp_zsyr2k_batch.c\
		  omp_zsyrk_batch.c omp_ztrsm_batch.c\
		  omp_ztrmm_batch.c

OMP_SRC        =$(addprefix $(BBLAS_OMP_DIR)/, $(OMP_SRC_LIST)) 

MKL_SRC_LIST     = bblas_zmkl.c
MKL_SRC        =$(addprefix $(BBLAS_MKL_DIR)/, $(MKL_SRC_LIST)) 


CUDA_SRC_LIST    = bblas_zcuda.c

CUDA_SRC         =$(addprefix $(BBLAS_CUDA_DIR)/, $(CUDA_SRC_LIST)) 

MAGMA_SRC_LIST    = bblas_zmagma.c

MAGMA_SRC         =$(addprefix $(BBLAS_MAGMA_DIR)/, $(MAGMA_SRC_LIST))

TEST_SRC      =  bblas_zaccuracy.c bblas_zflops.c \
		 bblas_zutil.c bblas_zprint.c \
		 bblas_ztesting.c  bblas_ztest.c bblas_zparsefile.c \
		 testing_zgemm_batch.c  testing_zhemm_batch.c \
		 testing_zher2k_batch.c testing_zherk_batch.c \
		 testing_zsymm_batch.c	 testing_zsyr2k_batch.c \
		 testing_zsyrk_batch.c testing_ztrmm_batch.c \
		 testing_ztrsm_batch.c 


AUX_COMMON  =  bblas_timer.c bblas_print_common.c

SOURCES       = $(BBLAS_SRC) $(TEST_SRC) $(OMP_SRC) $(MKL_SRC) $(CUDA_SRC) $(MAGMA_SRC)
SOURCES_Z      = $(SOURCES) $(BBLAS_COMMON) $(AUX_COMMON)
OBJECTS_Z     = $(SOURCES_Z:.c=.o) 

RM_FROM_NO_COMPLEX_SRC  =  zhemm_batch.c zher2k_batch.c \
		           zherk_batch.c 
RM_FROM_NO_COMPLEX_TEST =  testing_zhemm_batch.c \
	                   testing_zher2k_batch.c \
			   testing_zherk_batch.c \
		     	   bblas_omp/omp_zhemm_batch.c \
			   bblas_omp/omp_zher2k_batch.c bblas_omp/omp_zherk_batch.c

RM_FROM_NO_COMPLEX_OMP =  omp_zhemm_batch.c omp_zher2k_batch.c omp_zherk_batch.c

# single complex precision
SOURCES_C   = $(addprefix $(BBLAS_SRC_DIR)/, $(subst z,c,$(BBLAS_SRC_LIST)))
SOURCES_C  += $(addprefix $(BBLAS_TEST_DIR)/, $(subst z,c,$(TEST_SRC)))
SOURCES_C  += $(addprefix $(BBLAS_MAGMA_DIR)/, $(subst z,c,$(MAGMA_SRC_LIST)))
SOURCES_C  += $(addprefix $(BBLAS_CUDA_DIR)/, $(subst z,c,$(CUDA_SRC_LIST)))
SOURCES_C  += $(addprefix $(BBLAS_MKL_DIR)/,$(subst z,c,$(MKL_SRC_LIST)))
SOURCES_C  += $(addprefix $(BBLAS_OMP_DIR)/, $(subst z,c,$(OMP_SRC_LIST)))
SOURCES_C  += $(BBLAS_COMMON) $(AUX_COMMON)
OBJECTS_C   = $(SOURCES_C:.c=.o) 


# temporary variable to store sources 
#without routines specific to complex
SOURCES_REAL_SRC  = $(filter-out $(RM_FROM_NO_COMPLEX_SRC),$(BBLAS_SRC_LIST))
SOURCES_REAL_TEST = $(filter-out $(RM_FROM_NO_COMPLEX_TEST),$(TEST_SRC))

# double precision
SOURCES_D   = $(addprefix $(BBLAS_SRC_DIR)/, $(subst z,d,$(SOURCES_REAL_SRC)))
SOURCES_D  += $(addprefix $(BBLAS_TEST_DIR)/, $(subst z,d,$(SOURCES_REAL_TEST)))
SOURCES_D  += $(addprefix $(BBLAS_MAGMA_DIR)/, $(subst z,d,$(MAGMA_SRC_LIST)))
SOURCES_D  += $(addprefix $(BBLAS_CUDA_DIR)/, $(subst z,d,$(CUDA_SRC_LIST)))
SOURCES_D  += $(addprefix $(BBLAS_MKL_DIR)/,$(subst z,d,$(MKL_SRC_LIST)))
SOURCES_D  += $(addprefix $(BBLAS_OMP_DIR)/, $(subst z,d,$(filter-out $(RM_FROM_NO_COMPLEX_OMP), $(OMP_SRC_LIST))))
SOURCES_D  += $(BBLAS_COMMON) $(AUX_COMMON)
OBJECTS_D   = $(SOURCES_D:.c=.o) 

# single precision
SOURCES_S   = $(addprefix $(BBLAS_SRC_DIR)/, $(subst z,s,$(SOURCES_REAL_SRC)))
SOURCES_S  += $(addprefix $(BBLAS_TEST_DIR)/, $(subst z,s,$(SOURCES_REAL_TEST)))
SOURCES_S  += $(addprefix $(BBLAS_MAGMA_DIR)/, $(subst z,s,$(MAGMA_SRC_LIST)))
SOURCES_S  += $(addprefix $(BBLAS_CUDA_DIR)/, $(subst z,s,$(CUDA_SRC_LIST)))
SOURCES_S  += $(addprefix $(BBLAS_MKL_DIR)/,$(subst z,s,$(MKL_SRC_LIST)))
SOURCES_S  += $(addprefix $(BBLAS_OMP_DIR)/, $(subst z,s,$(filter-out $(RM_FROM_NO_COMPLEX_OMP), $(OMP_SRC_LIST))))
SOURCES_S  += $(BBLAS_COMMON) $(AUX_COMMON)
OBJECTS_S   = $(SOURCES_S:.c=.o) 

all: 
	$(MAKE) generate 
	$(MAKE) ztest stest dtest ctest 

.DEFAULT_GOAL := all

ztest: $(OBJECTS_Z) 
	$(CC) $(OBJECTS_Z) $(LDFLAGS)   -o $@  

stest: $(OBJECTS_S) 
	$(CC) $(OBJECTS_S) $(LDFLAGS)   -o $@  

dtest: $(OBJECTS_D) 
	$(CC) $(OBJECTS_D) $(LDFLAGS)   -o $@  

ctest: $(OBJECTS_C) 
	$(CC) $(OBJECTS_C) $(LDFLAGS)   -o $@  

.c.o:
	$(CC) $(CFLAGS) $(DEPS) $<   -o $@

clean:
	rm *.o $(BBLAS_OMP_DIR)/*.o $(BBLAS_SRC_DIR)/*.o  $(BBLAS_CUDA_DIR)/*.o $(BBLAS_MKL_DIR)/*.o $(BBLAS_MAGMA_DIR)/*.o
	rm ztest stest dtest ctest
