/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */

double* my_solver(int N, double *A, double* B) {
	printf("NEOPT SOLVER\n");
	/* allocate memory for intermediary matrixes */
	double *AB = (double*) calloc(N * N, sizeof(double));
	double *ABAt = (double*) calloc(N * N, sizeof(double));
	double *BtBt = (double*) calloc(N * N, sizeof(double));
	/* result matrix */
	double *R = (double*) calloc(N * N, sizeof(double));

	int i,j,k;

	for (i = 0; i < N; i++) {
  		for (j = 0; j < N; j++) {
			/* calc AB */
			double *posij = AB + i * N + j;
      		*posij = 0.0;
			/* skip 0s operations of A */
			for (k = i; k < N; k++) {
				double *aik = A + i * N + k;
				double *bkj = B + k * N + j;
				(*posij) += (*aik) * (*bkj);
      		}
   		}
	}

	for (i = 0; i < N; i++) {
  		for (j = 0; j < N; j++) {
			/* calc BtBt */
			double *posij = BtBt + i * N + j;
      		*posij = 0.0;
      		for (k = 0; k < N; k++) {
				/* indexes are taken in reverse to simulate the transposition */
				double *aik = B + k * N + i;
				double *bkj = B + j * N + k;
				(*posij) += (*aik) * (*bkj);
      		}

			/* calc ABAt */
			posij = ABAt + i * N + j;
      		*posij = 0.0;
			/* skip 0s operations of A */
			for (k = j; k < N; k++) {
				double *aik = AB + i * N + k;
				/* indexes are taken in reverse to simulate the transposition */
				double *bkj = A + j * N + k;
				(*posij) += (*aik) * (*bkj);
      		}

			/* calc R coresponding to already calc value */
			*(R + i * N + j) = *posij + *(BtBt + i * N + j);
   		}
	}

	/* free memory */
	free(AB);
	free(ABAt);
	free(BtBt);
	return R;
}
