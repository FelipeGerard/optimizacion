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

    matsub=malloc(sizeof(double)*nrow/nproc*ncol/nproc);
    i2 = 0;
    printf("nrow=%d, ncol=%d, nproc=%d\n", nrow, ncol, nproc);
    for(i=nrow/nproc*num;i<nrow/nproc*(num+1);i++){
	j2 = 0;
	for(j=ncol/nproc*num; j < ncol/nproc*(1+num); j++){
	    printf("i=%d, j=%d\ti2=%d, j2=%d, idx=%d, idx_sub=%d\tA=%2.1f\tAsub=%2.1f\n", i, j, i2, j2, j*nrow+i, j2*(nrow/nproc)+i2, mat[j*nrow + i], matsub[j2*(nrow/nproc) + i2]);
	    matsub[j2*(nrow/nproc) + i2] = mat[j*nrow + i];
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
   

    int id,np;
    int dim1,dim2,dim3;
    double *a,*b,*c;
    double *adist,*bdist,*cdist;

    MPI_Init(&argc,&argv);

    sscanf(argv[1],"%d",&dim1);
    sscanf(argv[2], "%d",&dim2);
    sscanf(argv[3],"%d",&dim3);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if(id==0){
	printf("==========================\n");
	double *mat, *submat1, *submat2;
	mat = malloc(sizeof(double)*4*4);
	submat1 = malloc(sizeof(double)*2*4);
	submat2 = malloc(sizeof(double)*2*4);
	for(int p=0;p<4;p++){
	    for(int q=0;q<4;q++){
		mat[q*4+p]=q*4+p;
	    }
	}
	submat1 = subMatrix(4, 4, 1, 2, mat);
	submat2 = subMatrix(4, 4, 2, 2, mat);
	printMat(4,4,mat);
	printMat(2,4,submat1);
	printMat(2,4,submat2);

	printf("==========================\n");
    }

    if(id==0){

	a=malloc(sizeof(double)*dim1*dim2);
	b=malloc(sizeof(double)*dim2*dim3);
	matrizAleatoria(dim1,dim2,a);
	matrizAleatoria(dim2,dim3,b);
	printf("A = \n");
	printMat(dim1,dim2,a);
	printf("B = \n");
	printMat(dim1,dim2,b);
    }

    adist=malloc(sizeof(double)*(dim1/np)*dim2);
    bdist=malloc(sizeof(double)*(dim2/np)*dim3);


    MPI_Scatter(a,(dim1*dim2)/np,MPI_DOUBLE,adist,(dim1*dim2)/np,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Scatter(b,(dim2*dim3)/np,MPI_DOUBLE,bdist,(dim2*dim3)/np,MPI_DOUBLE,0,MPI_COMM_WORLD);


    //printMat(dim1/np,dim2,adist);

    //printMat(dim2/np,dim3,bdist);



    char TRANSA, TRANSB;
    int M, N, K;
    double ALPHA, BETA;

    TRANSA = 'N';
    TRANSB = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    M=dim1/np;
    K=dim2/np; 
    N=dim3;

    cdist=malloc(sizeof(double)*(dim1*dim3)/np);

    zerosMat(dim1/np,dim3,cdist);


    /*double *asub;


      asub=malloc(sizeof(double)*(dim1/np*dim2/np));

      asub=subMatrix(dim1,dim2,1,np,adist); //el 1 es i-jmodp

      printMat(dim1/np,dim2/np,asub);

      dgemm_(&TRANSA, &TRANSB, &N, &M, &K, &ALPHA, bdist, &N, asub, &K, &BETA, cdist, &N);

      printMat(dim1/np,dim3,cdist);*/


    int j;
    double *asub;

    asub=malloc(sizeof(double)*(dim1/np*dim2/np));

    for(j=0;j<np;j++){

	printf("\n---------------------------------\nID: %d\tj: %d\n", id, j);
	asub=subMatrix(dim1,dim2,mod(id-j,np),np,adist);
	dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, asub, &M, bdist, &K, &BETA, cdist, &M);

	printf("Ai = \n");
	printMat(dim1/np,dim2,adist);
	printf("Bi = \n");
	printMat(dim2/np,dim3,bdist);
	printf("Ci = \n");
	printMat(dim1/np,N,cdist);


	MPI_Send(bdist, (dim2*dim3)/np, MPI_DOUBLE, mod(id + 1,np), mod(id-j,np),MPI_COMM_WORLD);

	MPI_Recv(bdist,(dim2*dim3)/np, MPI_DOUBLE, mod(id - 1,np), mod(id-j-1,np), MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    }
    /*
       printf("\n==================================\nID: %d\n", id);
       printf("A = \n");
       printMat(M,K,adist);
       printf("B = \n");
       printMat(K,N,bdist);
       printf("C = \n");
       printMat(M,N,cdist);
     */



    if(id==0){
	free(a);
	free(b);
    }

    free(adist);
    free(bdist);
    free(asub);

    MPI_Finalize();
}











