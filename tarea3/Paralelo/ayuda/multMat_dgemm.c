#include <stdio.h>

extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);

double *multMat_dgemm(int dimren,int dimcomp,int dimcol,double *mat1,double *mat2){

double *res;

char TRANSA, TRANSB;
int M, N, K;
double ALPHA, BETA;


res=malloc(sizeof(double)*dimren*dimcomp);




TRANSA = 'N';
TRANSB = 'N';
ALPHA = 1.0;
BETA = 0.0;

M=dimren;
K=dimcomp; 
N=dimcol;


dgemm_(&TRANSA, &TRANSB, &N, &M, &K, &ALPHA, mat2, &N, mat1, &K, &BETA, res, &N);



return res;

}


