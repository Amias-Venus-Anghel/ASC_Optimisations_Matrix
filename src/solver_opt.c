/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	printf("OPT SOLVER\n");
	/* allocate memory for intermediary matrixes */
	double *AB = (double*) calloc(N * N, sizeof(double));
	double *ABAt = (double*) calloc(N * N, sizeof(double));
	double *BtBt = (double*) calloc(N * N, sizeof(double));
	/* result matrix */
	double *R = (double*) calloc(N * N, sizeof(double));

	register int i,j,k;

	/* calc AB */
	for (i = 0; i < N; i++) {
		/* skip 0s operations of A */
  		for (k = i; k < N; k ++) {
			register double *p = AB + i * N;
			register double aik = *(A + i * N + k);
			register double *bkj = B + k * N;

			for (j = 0; j < N; j++, p++) {
				*p += aik * (*bkj);
				bkj++;
      		}
   		}
	}

	for (i = 0; i < N; i++) {
  		for (j = 0; j < N; j++) {
			/* calc BtBt */
			register double suma = 0.0;
			register double *bkj = B + j * N;
      		for (k = 0; k < N; k++, bkj++) {
				/* indexes are taken in reverse to simulate the transposition */
				register double *aik = B + k * N + i;
				suma += (*aik) * (*bkj);
      		}
			*(BtBt + i * N + j) = suma;
		
			/* calc ABAt */
			register double suma2 = 0.0;
			/* skip 0s operations of A */
			for(k = j&(~3); k < N; k+=4) {
				register double *aik = AB + i * N + k;
				register double *bkj = A + j * N + k;
				/* calc 4 at once */
				suma2 += (*aik) * (*bkj);
				suma2 += (*(aik + 1)) * (*(bkj + 1));
				suma2 += (*(aik + 2)) * (*(bkj + 2));
				suma2 += (*(aik + 3)) * (*(bkj + 3));
			}
			*(ABAt + i * N + j) = suma2;
			/* calc R coresponding to already calc value */
			*(R + i * N + j) = suma + suma2;
   		}
	}
	
	/* free memory */
	free(AB);
	free(ABAt);
	free(BtBt);
	return R;	
}
