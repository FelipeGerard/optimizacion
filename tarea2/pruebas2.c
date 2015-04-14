#include <stdio.h>
#include <lapacke.h>
/*extern 'C' {
    #include <OpenBlas/include/cblas.h>
}*/
//extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
//extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);

void multMat(int nrow, int nshared, int ncol, double A[nrow][nshared], double B[nshared][ncol], double C[nrow][ncol]){
    for(int i=0; i<nrow; i++){
	for(int j=0; j<ncol; j++){
	    C[i][j] = 0.0;
	    for(int k=0; k<nshared; k++){
		C[i][j] = C[i][j] + A[i][k]*B[k][j];
	    }
	}
    }
}

void printMat(int nrow, int ncol, double A[nrow][ncol]){
    for(int i=0; i<nrow; i++){
	for(int j=0; j<ncol; j++){
	    printf("[%d,%d] = %f\t", i, j, A[i][j]);
	}
	printf("\n");
    }   
}

void t(int nrow, int ncol, double A[nrow][ncol], double A_t[ncol][nrow]){
    //A_t = malloc(sizeof(A));
    for(int i=0; i<ncol; i++){
	for(int j=0; j<nrow; j++){
	    A_t[i][j] = A[j][i];
	}
    }   
}

void main(){
    printf("Prueba en C:\n");
    int m=2, n=3;
    double A[m][n], B[n][m], C[m][m];
    for(int i=0; i<m; i++){
	for(int j=0; j<n; j++){
	    A[i][j] = (i)*n + j;
//	    C[i][j] = 10.0;
	}
    }
    t(m, n, A, B);
    B[2][1] = 8;
    multMat(m, n, m, A, B, C);
    printMat(m,n,A);
    printf("\n");
    printMat(n,m,B);
    printf("\n");
    printMat(m,m,C);
    printf("\n");
    for(int p=0; p<m*m; p++){
	printf("%f\n", *C[p]);
    }
}
