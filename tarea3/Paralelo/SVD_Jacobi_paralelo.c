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

    m = atoi(argv[1]);
    n = atoi(argv[2]);
    machineEps = 1.e-16;
    no_rot = 0;
    ban = 'F';
    col_ort = n*(n-1)/2;
    precision = 1.e-8;
    sweeps = 0;
    maxsweeps = 10;


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
    
    printf("\nAdt = \n");
    impMat(m, n/np, Adt);
    printf("\nVdt = \n");
    impMat(n, n/np, Vdt);

    //if(id == 0) printf("\n-------------- AFUERA\n");
    while(ban != 'T' && sweeps <= maxsweeps){
	no_rot = 0;

	//if(id == 0) printf("\n-------------- WHILE\n");
	for(int k = 0; k<2*(np-1); k++){
	    printf("\n-------------- FOR, id = %d, k = %d\n", id, k);
	    printf("id = %d", id);
	    for(int p=0; p<m; p++){
		printf("\nAi[%d] = %f\t", p, A_i[p]);
	    }
	    printf("\n");
	    for(int p=0; p<m; p++){
		printf("\nAj[%d] = %f\t", p, A_j[p]);
	    }
	    prod = prodPunto(m,A_i,A_j);
	    //printf("\nprod = %f\n", prod);
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
		    double *V_i,*V_j;
		    V_i = extraeCol(n,n,k,Vdt);
		    V_j = extraeCol(n,n,j,Vdt);
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
		}
		else{
		    no_rot += np;
		}

		double *aux1, *aux2;
		aux1 = malloc(sizeof(double)*n);
		aux2 = malloc(sizeof(double)*n);
		// Mandamos y recibimos Aj
		if(id == 0){
		    MPI_Send(A_j, m, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		}else{
		    MPI_Send(A_j, m, MPI_DOUBLE, id-1, id, MPI_COMM_WORLD);
		    printf("\n************ MY NIGGA B-) ************\n");
		}
		if(id == 1){
		    MPI_Recv(aux1, m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}else{
		    MPI_Recv(aux2, m, MPI_DOUBLE, id+1, id+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// Mandamos y recibimos Ai
		if(id == np){
		    aux2 = A_i;
		}else if(id != 0){
		    MPI_Send(A_i, m, MPI_DOUBLE, id+1, id, MPI_COMM_WORLD);
		} 
		if(id != 0){
		    MPI_Recv(aux1, m, MPI_DOUBLE, id-1, id-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// Asignamos Ai y Aj
		A_i = aux1;
		A_j = aux2;


		////// Transferimos Adist
		// Transferimos A_j
		/*
		   if(id == 0){
		   MPI_Send(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		   }else if(id != np){
		   MPI_Send(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		   }
		   MPI_Recv(aux2, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);

		// Transferimos A_i
		if(id != 0){
		if(id != np){
		MPI_Send(A_i, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);
		}else{
		MPI_Send(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		}
		}
		MPI_Recv(aux1, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);




		if(id == 1){
		MPI_Recv(aux, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);
		}else if(id == np){
		MPI_Recv(aux, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);
		}else{
		MPI_Recv(aux, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);
		}


		MPI_Send(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		MPI_Send(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		MPI_Recv(aux, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);
		MPI_Send(A_i, m, MPI_DOUBLE, mod(id-1,np), mod(id+i,np), MPI_COMM_WORLD);
		MPI_Recv(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		A_i = aux;
		if(id != 0){
		MPI_Send(A_i, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		}
		if(id == 0){
		MPI_Send(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		MPI_Send(V_j, n, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		} else {
		MPI_Send(A_i, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		}
		if(id == 0){
		MPI_Recv(A_j, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		MPI_Recv(V_j, n, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		} else {
		MPI_Recv(A_i, m, MPI_DOUBLE, mod(id+1,np), mod(id-i,np), MPI_COMM_WORLD);
		}
		 */
	    } //fin else

	}//fin for


	sweeps += 1;

	if(no_rot == col_ort){
	    ban = 'T';
	}

    }//fin while

    printf("Número de rotaciones:%d\n",no_rot);
    printf("Número de sweeps:%d\n",sweeps);


    printf("V = \n");
    impMat(n,n,V);

}
