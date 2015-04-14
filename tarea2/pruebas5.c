#include <stdio.h>
#include <lapacke.h>
/*extern 'C' {
#include <OpenBlas/include/cblas.h>
}*/
//extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);

void multMat(int nrow, int nshared, int ncol, double *A, double *B, double *C){
    for(int i=0; i<nrow; i++){
	for(int j=0; j<ncol; j++){
	    C[(i)*ncol + j] = 0.0;
	    for(int k=0; k<nshared; k++){
		C[(i)*ncol + j] += A[(i)*nshared + k]*B[(k)*ncol + j];
	    }
	}
    }
}

void printMat(int nrow, int ncol, double *A){
    for(int i=0; i<nrow; i++){
	for(int j=0; j<ncol; j++){
	    printf("[%d,%d] = %f\t", i, j, A[(j)*nrow + i]);
	}
	printf("\n");
    }   
}

void t(int nrow, int ncol, double *A, double *A_t){
    for(int i=0; i<ncol; i++){
	for(int j=0; j<nrow; j++){
	    A_t[(j)*ncol + i] = A[(i)*ncol + j];
	}
    }   
}

void main(){
    printf("Prueba en C:\n");
    int m=3, n=4, k=1;
    double *A, *B, *C, *D;
    A = malloc(sizeof(double)*m*k);
    B = malloc(sizeof(double)*k*n);
    C = malloc(sizeof(double)*m*n);
    D = malloc(sizeof(double)*m*n);
    for(int i=0; i<m; i++){
	for(int j=0; j<k; j++){
	    A[(j)*m + i] = (double)(i - j);
	    //	    C[i][j] = 10.0;
	}
    }
    for(int i=0; i<k; i++){
	for(int j=0; j<n; j++){
	    B[(j)*k + i] = (double)(i + j);
	    //	    C[i][j] = 10.0;
	}
    }
    for(int i=0; i<m; i++){
	for(int j=0; j<n; j++){
	   C[(j)*m + i] = 0.0; 
	   D[(j)*m + i] = 0.0; 
	    //	    C[i][j] = 10.0;
	}
    }
 
    //t(m, n, A, B);
    //B[n + 1] = 8;
    //    multMat(m, n, m, A, B, C);
    char ta = 'N';
    char tb = 'N';
    double alpha = 1.0;
    double beta = 0.0;
//    dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
    dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
    multMat(m, k, n, A, B, D);
    printf("A = \n");
    printMat(m,k,A);
    printf("\n");
    printf("B = \n");
    printMat(k,n,B);
    printf("\n");
    printf("D = \n");
    printMat(m,n,D);
    printf("\n");
    printf("C = \n");
    printMat(m,n,C);
    printf("\n");

    for(int p=0; p<m*n; p++){
	printf("%f\n", C[p]);
    }

}
