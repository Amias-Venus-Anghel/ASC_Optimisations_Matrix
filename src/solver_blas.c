/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include <cblas.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	printf("BLAS SOLVER\n");

	double *ABAt = (double*) calloc(N * N, sizeof(double));
	double *R = (double*) calloc(N * N, sizeof(double));

	/* copy B in ABAt */
	cblas_dcopy(N * N, B, 1, ABAt, 1);
	/* A * B */
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, 
		CblasNonUnit, N, N, 1, A, N, ABAt, N);
	/* AB * At */
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, 
		CblasNonUnit, N, N, 1, A, N, ABAt, N);
	/* Bt * Bt */
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, N, N, N, 1, B, N, B, N, 1, R, N);
	/* ABAt + BtBt */
	cblas_daxpy(N * N, 1.0, ABAt, 1, R, 1);
	
	/* free auxiliar memory */
	free (ABAt);

	return R;
}
