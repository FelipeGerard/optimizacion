#include <stdio.h>
#include "matrizAleatoria.c"
#include "impMat.c"
#include "multMat.c"
extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);


void main(){
int dim1=3,dim2=3,dim3=3;

double *a,*b,*c,*d;

char TRANSA, TRANSB;
int M, N, K;
double ALPHA, BETA;



a=malloc(sizeof(double)*dim1*dim2);

b=malloc(sizeof(double)*dim2*dim3);

c=malloc(sizeof(double)*dim1*dim3);

matrizAleatoria(dim1,dim2,a);

impMat(dim1,dim2,a);


matrizAleatoria(dim2,dim3,b);

impMat(dim2,dim3,b);

multMat(dim1,dim2,dim3,a,b,c);


impMat(dim1,dim3,c);



TRANSA = 'N';
TRANSB = 'N';
ALPHA = 1.0;
BETA = 0.0;

M=dim1;
K=dim2; 
N=dim3;

d = malloc(sizeof(double)*dim1*dim3);

dgemm_(&TRANSA, &TRANSB, &N, &M, &K, &ALPHA, b, &N, a, &K, &BETA, d, &N);




impMat(dim1,dim3,d);



}


