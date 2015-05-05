#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "asignaCol.c"
#include "extraeCol.c"
#include "mod.c"
#include "prodPunto.c"
#include "zerosMat.c"
#include "constMat.c"
#include "impMat.c"
#include "impMatCSV.c"
#include "matrizId.c"
#include "multMat_dgemm.c"
#include "signo.c"

void *transp(int m, int n, double *mat, double *out){
    for(int i=0; i<m; i++){
	for(int j=0; j<n; j++){
	    out[j*m + i] = mat[i*n + j];
	}
    }
}

void *matrizPrueba(int dimren,int dimcol,double *mat){
    mat[0] = 0;
    mat[1] = 1;
    mat[2] = 2;
    mat[3] = 3;
    mat[4] = 4;
    mat[5] = 5;
    mat[6] = 6;
    mat[7] = 7;
    mat[8] = 8;
    mat[9] = 9;
    mat[10] = 10;
    mat[11] = 11;
    mat[12] = 12;
    mat[13] = 13;
    mat[14] = 14;
    mat[15] = 15;
}

void rotatePossession(int id, int np, int m, double *A_i, double *A_j, double *aux1, double *aux2, int verbose){
    if(verbose >= 1) printf("\nRotating possessions...\n");
    // Mandamos y recibimos Aj
    if(id == 0){
	if(verbose == 2) printf("\nSending Aj: %d --> %d", 0, 1); 
	MPI_Send(A_j, m, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	if(verbose == 2) printf("\nSuccess!\n");
    }else{
	if(verbose == 2) printf("\nSending Aj: %d --> %d", id, id-1); 
	MPI_Send(A_j, m, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD);
	if(verbose == 2) printf("\nSuccess!\n");
    }
    if(id == 1){
	if(verbose == 2) printf("\nReceiving Aj: %d --> %d", 0, 1); 
	MPI_Recv(aux1, m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(verbose == 2) printf("\nSuccess!\n");
	//impMat(m,1,A_j);
	//impMat(m,1,aux1);
    }else if(id < np-1){
	if(verbose == 2) printf("\nReceiving Aj: %d --> %d", id+1, id); 
	MPI_Recv(aux2, m, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(verbose == 2) printf("\nSuccess!\n");
    }
    // Esperamos a que se completen estas transferencias
    MPI_Barrier(MPI_COMM_WORLD);
    if(verbose == 2)    printf("\nCPU %d is through the first comm round...\n", id); 

    // Mandamos y recibimos Ai
    if(id == np-1){
	if(verbose == 2) printf("\nSending Ai: %d --> %d", np-1, np-1); 
	aux2 = A_i;
	if(verbose == 2) printf("\nReceiving Aj: %d --> %d", np-1, np-1); 
    }else if(id > 0){
	if(verbose == 2) printf("\nSending Ai: %d --> %d", id, id+1); 
	MPI_Send(A_i, m, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD);
    } 
    if(id > 1){
	if(verbose == 2) printf("\nReceiving Ai: %d --> %d", id-1, id); 
	MPI_Recv(aux1, m, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Asignamos Ai y Aj
    A_i = aux1;
    A_j = aux2;

    // Esperamos a que se completen todas las transferencias
    MPI_Barrier(MPI_COMM_WORLD);
}
void main(int argc,char **argv){
    int id,np,m,n; //dimensiones de la matriz a calcular su SVD
    double *Amat,*At,*V;
    int no_rot; //contador de no rotaciones
    char ban; //valor booleano
    int col_ort; //número de columnas ortogonales
    double prec; //precisión
    int sweeps;//contador de sweeps
    int maxsweeps;
    double *A_i, *A_j, *Adist, *Adt;
    double *V_i, *V_j, *Vdist, *Vdt;
    double precision;
    double prod,A_i_norma,A_j_norma,num,den,machineEps;

    int cha;
    cha = id;
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    machineEps = 1.e-16;
    no_rot = 0;
    ban = 'F';
    col_ort = n*(n-1)/2;
    precision = 1.e-8;
    sweeps = 0;
    maxsweeps = 0;


    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    MPI_Comm_size(MPI_COMM_WORLD, &np);


    if(id == 0){
	printf("m = %d\nn = %d\nnp = %d\n", m, n, np);
	Amat = malloc(sizeof(double)*m*n);
	At = malloc(sizeof(double)*n*m);
	V = malloc(sizeof(double)*n*n);
	//matrizAleatoria(m,n,Amat);
	zerosMat(m,n,Amat);//comentar esta línea para matrices aleatorias
	matrizPrueba(m,n,Amat);//comentar esta línea para matrices aleatorias
	transp(m,n,Amat,At);
	matrizId(n,n,V);
	printf("A = \n");
	impMat(m,n,Amat);
	printf("At = \n");
	impMat(m,n,At);
    }

    Adist = malloc(sizeof(double)*m*n/np);
    Vdist = malloc(sizeof(double)*n*n/np);
    Adt = malloc(sizeof(double)*m*n/np);
    Vdt = malloc(sizeof(double)*n*n/np);
    // En esta versión np debe ser igual a n/2
    MPI_Scatter(At, m*(n/np), MPI_DOUBLE, Adist, m*(n/np), MPI_DOUBLE, 0, MPI_COMM_WORLD); // n/np = 2 a fuerza
    MPI_Scatter(V, n*(n/np), MPI_DOUBLE, Vdist, n*(n/np), MPI_DOUBLE, 0, MPI_COMM_WORLD);


    double costeta,senteta;
    double ji,t;
    double cte;
    int i,j;

    printf("\n===============================================\n");
    printf("id = %d\n", id);

    transp(n/np,m,Adist,Adt);
    transp(n/np,n,Vdist,Vdt);
    A_i = extraeCol(m,n/np,0,Adt);
    A_j = extraeCol(m,n/np,1,Adt);
    V_i = extraeCol(n,n/np,0,Vdt);
    V_j = extraeCol(n,n/np,1,Vdt);

    printf("\nAdt = \n");
    impMat(m, n/np, Adt);
    printf("\nVdt = \n");
    impMat(n, n/np, Vdt);

    double *aux1, *aux2;
    aux1 = malloc(sizeof(double)*n);
    aux2 = malloc(sizeof(double)*n);


    //if(id == 0) printf("\n-------------- AFUERA\n");
    while(ban != 'T' && sweeps <= maxsweeps){
	no_rot = 0;

	//if(id == 0) printf("\n-------------- WHILE\n");
	for(int k = 0; k<2*(np-1); k++){
	    printf("\n-------------- FOR, sweep = %d, id = %d, k = %d\n", sweeps, id, k);
	    printf("Ai = ");
	    impMat(1,m,A_i);
	    printf("Aj = ");
	    impMat(1,m,A_j);
	    prod = prodPunto(m,A_i,A_j);

	    A_i_norma = sqrt(prodPunto(m,A_i,A_i));
	    A_j_norma = sqrt(prodPunto(m,A_j,A_j));

	    num = pow(A_j_norma,2)-pow(A_i_norma,2);
	    den = 2*prod;
	    if(fabs(den) <= machineEps*fabs(num)){
		costeta = 1;
		senteta = 0;
		no_rot += 1;
	    }
	    else{
		ji = num/den;
		t = signo(ji)/(fabs(ji)+sqrt(1+pow(ji,2)));	
		costeta = 1/sqrt(1+pow(t,2));
		senteta = costeta*t;
		cte = prod/(A_i_norma*A_j_norma);
		if(fabs(cte)>precision){
		    double *A_ij,*V_ij;
		    A_ij = constMat(m,2,A_i,A_j);
		    V_ij = constMat(n,2,V_i,V_j);
		    double *rotmatriz;
		    rotmatriz = malloc(sizeof(double)*2*2);
		    rotmatriz[0] = costeta;	
		    rotmatriz[1] = senteta;
		    rotmatriz[2] = -senteta;
		    rotmatriz[3] = costeta;
		    double *rot_A,*rot_V	;
		    rot_A = multMat_dgemm(m,2,2,A_ij,rotmatriz);
		    rot_V = multMat_dgemm(n,2,2,V_ij,rotmatriz);
		    A_i = extraeCol(m,2,0,rot_A);
		    A_j = extraeCol(m,2,1,rot_A);
		    V_i = extraeCol(n,2,0,rot_V);
		    V_j = extraeCol(n,2,1,rot_V);
		}
		else{
		    no_rot += np;
		}
		printf("Ai = ");
		impMat(1,m,A_i);
		printf("Aj = ");
		impMat(1,m,A_j);

	    } //fin else

	    rotatePossession(id, np, m, A_i, A_j, aux1, aux2, 0);
	    rotatePossession(id, np, n, V_i, V_j, aux1, aux2, 0);

	}//fin for


	sweeps += 1;

	if(no_rot == col_ort){
	    ban = 'T';
	}

    }//fin while

    printf("\nNúmero de rotaciones:%d\n",no_rot);
    printf("Número de sweeps:%d\n",sweeps);


    printf("V = \n");
    // impMat(n,n,V);

}
