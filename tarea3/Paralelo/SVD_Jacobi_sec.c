#include <stdio.h>
#include <stdlib.h>
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

void *matrizPrueba(int dimren,int dimcol,double *mat){
    mat[0] = 1;
    mat[2] = 1;
    mat[4] = 1;
    mat[7] = 1;
    mat[8] = 1;
    mat[10] = 1;
    mat[12] = 1;
    mat[13] = 1;

}


void main(int argc,char **argv){
    int m,n; //dimensiones de la matriz a calcular su SVD
    double *Amat,*V;
    int no_rot; //contador de no rotaciones
    char ban; //valor booleano
    int col_ort; //número de columnas ortogonales
    double prec; //precisión
    int sweeps;//contador de sweeps
    int maxsweeps;
    double *A_i,*A_j;
    double *V_i,*V_j;
    double precision;


    m = atoi(argv[1]);

    n = atoi(argv[2]);
    Amat = malloc(sizeof(double)*m*n);
    //matrizAleatoria(m,n,Amat);
    zerosMat(m,n,Amat);//comentar esta línea para matrices aleatorias
    matrizPrueba(m,n,Amat);//comentar esta línea para matrices aleatorias

    printf("A = \n");
    impMat(m,n,Amat);

    V = malloc(sizeof(double)*n*n);
    matrizId(n,n,V);
    no_rot = 0;
    ban = 'F';
    col_ort = n*(n-1)/2;
    precision = 1.e-8;
    sweeps = 0;
    maxsweeps = 10;


    double prod,A_i_norma,A_j_norma,num,den,machineEps;
    machineEps = 1.e-16;
    double costeta,senteta;
    double ji,t;
    double cte;
    int i,j;



    while(ban != 'T' && sweeps <= maxsweeps){
	no_rot = 0;

	for(i = 0; i<(n-1); i++){
	    for(j = i+1;j<n;j++){
		A_i = extraeCol(m,n,i,Amat);
		A_j = extraeCol(m,n,j,Amat);

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
			double *V_i,*V_j;
			V_i = extraeCol(n,n,i,V);
			V_j = extraeCol(n,n,j,V);
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
			asignaCol(m,n,i,extraeCol(m,2,0,rot_A),Amat);	
			asignaCol(m,n,j,extraeCol(m,2,1,rot_A),Amat);
			asignaCol(n,n,i,extraeCol(n,2,0,rot_V),V);
			asignaCol(n,n,j,extraeCol(n,2,1,rot_V),V);
		    }
		    else{
			no_rot += 1;
		    }

		}


	    }//fin for

	}//fin for


	sweeps += 1;

	if(no_rot == col_ort){
	    ban = 'T';
	}

    }//fin while

    printf("Número de rotaciones:%d\n",no_rot);
    printf("Número de sweeps:%d\n",sweeps);


    printf("V = \n");
    impMatCSV(n,n,V);

}
