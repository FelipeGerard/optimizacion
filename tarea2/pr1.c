#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
extern void dgemm_(char *transaA, char *transaB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);

void printMat(int nrow, int ncol, double *A){
    for(int i=0; i<nrow; i++){
	for(int j=0; j<ncol; j++){
	    printf("[%d,%d] = %2.1f\t", i, j, A[j*nrow + i]);
	}
	printf("\n");
    }
}

void* matrizAleatoria(int nrow, int ncol,double *mat)
{
    int i,j;
    for(i=0;i<nrow;i++){
	for(j=0;j<ncol;j++){
	    mat[j*nrow+i]=rand()%10;
	}
    }


}

double* subMatrix(int nrow,int ncol, int num, int nproc, double *mat){

    double *matsub;
    int i,j;
    int i2, j2;
    printf("nrow=%d, ncol=%d, num=%d, nproc=%d", nrow, ncol, num, nproc);
    matsub=malloc(sizeof(double)*nrow*ncol/nproc);
    i2 = 0;
    for(i=0;i<nrow;i++){
	j2 = 0;
	for(j=ncol/nproc*num; j < ncol/nproc*(1+num); j++){
	    printf("i=%d, j=%d\ti2=%d, j2=%d, idx=%d, idx_sub=%d\tA=%2.1f\tAsub=%2.1f\n", i, j, i2, j2, j*nrow+i, j2*nrow+i2, mat[j*nrow + i], matsub[j2*nrow + i2]);
	    matsub[j2*nrow + i2] = mat[j*nrow + i];
	    j2++;
	}
	i2++;
    }
    return matsub;
}


void* zerosMat(int nrow,int ncol,double *mat){
    int i;
    for(i=0; i<nrow*ncol; i++){
	mat[i] = 0;
    }
}


int mod (int a, int b)
{
    if(b < 0) 
	return mod(-a, -b);   
    int ret = a % b;
    if(ret < 0)
	ret+=b;
    return ret;
}

void main(int argc, char *argv[]){

    int id,np, numpart;
    //    double *a,*b,*c;
    //    double *adist,*bdist,*cdist;
    sscanf(argv[1],"%d",&numpart);
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    printf("id = %d\n", id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    printf("==========================\n");

    int nrow=2, ncol=4;
    double *mat;
    mat = malloc(sizeof(double)*nrow*ncol);
    for(int p=0;p<nrow;p++){
	for(int q=0;q<ncol;q++){
	    //printf("%d\t%d", p, q);
	    mat[q*nrow+p]=q*nrow+p;
	}
    }
    double *submat;
    submat = malloc(sizeof(double)*2*2);
    submat = subMatrix(2, 4, numpart, 2, mat);
    printf("MAT = \n");
    printMat(2,4,mat);
    printf("SUBMAT = \n");
    printMat(2,2,submat);

    /*    submat1 = subMatrix(2, 4, 0, 2, mat);
	  submat2 = subMatrix(2, 4, 1, 2, mat);
	  printf("MAT = \n");
	  printMat(2,4,mat);
	  printf("SUBMAT 1 = \n");
	  printMat(2,2,submat1);
	  printf("SUBMAT 2 = \n");
	  printMat(2,2,submat2);
     */
    printf("==========================\n");
    free(mat);
    free(submat);
    MPI_Finalize();
}











