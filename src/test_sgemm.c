// Author: Samuel D. Relton
#include <stdio.h>
#include <cblas.h>
#include "sgemm_batched.h"

float A1[] = {
3, 1, 3,
1, 5, 9,
2, 6, 5
};

float B1[] = {
3, 1, 3,
1, 5, 9,
2, 6, 5
};

float C1[] = {
0, 0, 0,
0, 0, 0,
0, 0, 0
};

float A2[] = {
3, 1, 3,
1, 5, 9,
2, 6, 5
};

float B2[] = {
3, 1, 3,
1, 5, 9,
2, 6, 5
};

float C2[] = {
0, 0, 0,
0, 0, 0,
0, 0, 0
};

float ALPHA[] = {1.0, 1.0};
float BETA[] = {0.0, 0.0};

int LDA[] = {3, 3};
int LDB[] = {3, 3};
int LDC[] = {3, 3};

bblas_trans_t  TRANSA[] = {CblasNoTrans, CblasNoTrans};
bblas_trans_t  TRANSB[] = {CblasNoTrans, CblasNoTrans};

int M[] = {3, 3};
int N[] = {3, 3};
int K[] = {3, 3};

int BATCHCOUNT = 2;
bblas_batch_type_t BATCH_TYPE = BBLASFixed;

int INFO[] = {0, 0};

float *arrayA[] = {A1, A2};
float *arrayB[] = {B1, B2};
float *arrayC[] = {C1, C2};

int main()
{
	int i, j, k;

	sgemm_batched(TRANSA, TRANSB, M, N, K, ALPHA, arrayA,
		      LDA, arrayB, LDB, BETA, arrayC, LDC, BATCHCOUNT,
		      BATCH_TYPE, INFO);
	

	for (k = 0; k < BATCHCOUNT; k++) {
		putchar('\n');
		for (i=0; i<3; ++i) {
			for (j=0; j<3; ++j) {
				printf("%5.1f", arrayC[k][i*3+j]);

			}
			putchar('\n');
		}
	}

    return 0;
}
