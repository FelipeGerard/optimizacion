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
#include "impMatTransp.c"
#include "matrizId.c"
#include "multMat_dgemm.c"
#include "signo.c"

void *transp(int m, int n, double *mat, double *out){
    for(int i=0; i<m; i++){
	for(int j=0; j<n; j++){
	    out[j*m + i] = mat[i*n + j];
	}
    }
    return NULL;
}

void matrizPrueba(int dimren,int dimcol,double *mat){
    int i;
    for(i = 0; i < dimren*dimcol; i++){
	mat[i] = i;
    }
}

void rotatePossession(int id, int np, int m, double *A_i, double *A_j, double *aux1, double *aux2, int verbose){
    int i = 0;

    if(verbose >= 1) printf("\nRotating possessions...\n");
    // Mandamos y recibimos Aj
    if(id == 0){
	if(verbose == 2) printf("\nSending Aj: %d --> %d\n", 0, 1); 
	MPI_Send(A_j, m, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	if(verbose == 2) printf("\nSuccess!\n");
    }else{
	if(verbose == 2) printf("\nSending Aj: %d --> %d\n", id, id-1); 
	MPI_Send(A_j, m, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD);
	if(verbose == 2) printf("\nSuccess!\n");
    }
    if(id == 1){
	if(verbose == 2) printf("\nReceiving Aj: %d --> %d\n", 0, 1); 
	MPI_Recv(aux1, m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(verbose == 2) printf("\nSuccess!\n");
    }else if(id < np-1){
	if(verbose == 2) printf("\nReceiving Aj: %d --> %d\n", id+1, id); 
	MPI_Recv(aux2, m, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(verbose == 2) printf("\nSuccess!\n");
    }
    // Esperamos a que se completen estas transferencias
    MPI_Barrier(MPI_COMM_WORLD);
    if(verbose == 2)    printf("\nCPU %d is through the first comm round...\n", id); 

    // Mandamos y recibimos Ai
    if(id == 0){
	for(i = 0; i<m; i++){
	    aux1[i] = A_i[i];
	}
    }else if(id == np-1){
	if(verbose == 2) printf("\nSending Ai: %d --> %d\n", np-1, np-1); 
	for(i = 0; i<m; i++){
	    aux2[i] = A_i[i];
	}
	if(verbose == 2) printf("\nReceiving Aj: %d --> %d\n", np-1, np-1); 
    }else{
	if(verbose == 2) printf("\nSending Ai: %d --> %d\n", id, id+1); 
	MPI_Send(A_i, m, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD);
    } 
    if(id > 1){
	if(verbose == 2) printf("\nReceiving Ai: %d --> %d\n", id-1, id); 
	MPI_Recv(aux1, m, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Asignamos Ai y Aj
    for(i = 0; i<m; i++){
	A_i[i] = aux1[i];
	A_j[i] = aux2[i];
    }

    // Esperamos a que se completen todas las transferencias
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc,char **argv){
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

    ///////////////////////////////
    ////// Inicializamos

    m = atoi(argv[1]);
    n = atoi(argv[2]);
    machineEps = 1.e-16;
    no_rot = 0;
    ban = 'F';
    col_ort = n*(n-1)/2;
    precision = 1.e-8;
    sweeps = 0;
    maxsweeps = atoi(argv[3]);


    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    MPI_Comm_size(MPI_COMM_WORLD, &np);


    if(id == 0){
	printf("\n================================================\n================================================");
	printf("\nm = %d\nn = %d\nnp = %d\nmaxsweeps = %d\n", m, n, np, maxsweeps);
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
	//printf("At = \n");
	//impMat(m,n,At);
	printf("V inicial = \n");
	impMat(n,n,V);
    }

    Adist = malloc(sizeof(double)*m*n/np);
    Vdist = malloc(sizeof(double)*n*n/np);
    Adt = malloc(sizeof(double)*m*n/np);
    Vdt = malloc(sizeof(double)*n*n/np);
    // Pasamos a cada nodo las columnas que le corresponden
    MPI_Scatter(At, m*(n/np), MPI_DOUBLE, Adist, m*(n/np), MPI_DOUBLE, 0, MPI_COMM_WORLD); // n/np = 2 a fuerza
    MPI_Scatter(V, n*(n/np), MPI_DOUBLE, Vdist, n*(n/np), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Indicadores de que nodo tiene que columna. En realidad no son necesarios
    double *coli, *colj, *colauxi, *colauxj;
    coli = malloc(sizeof(double)*1);
    colj = malloc(sizeof(double)*1);
    colauxi = malloc(sizeof(double)*1);
    colauxj = malloc(sizeof(double)*1);
    coli[0] = 2.0*id;
    colj[0] = 2.0*id + 1.0;

    double costeta,senteta;
    double ji,t;
    double cte;
    int i,j;

    transp(n/np,m,Adist,Adt);
    transp(n/np,n,Vdist,Vdt);
    A_i = extraeCol(m,n/np,0,Adt);
    A_j = extraeCol(m,n/np,1,Adt);
    V_i = extraeCol(n,n/np,0,Vdt);
    V_j = extraeCol(n,n/np,1,Vdt);

    //printf("\n===============================================\n");
    //printf("id = %d\n", id);
    //printf("\nAdt = \n");
    //impMat(m, n/np, Adt);
    //printf("\nVdt = \n");
    //impMat(n, n/np, Vdt);
    
    // Auxiliares para enviar y recibir 
    double *aux1, *aux2, *vaux1, *vaux2;
    aux1 = malloc(sizeof(double)*m);
    aux2 = malloc(sizeof(double)*m);
    vaux1 = malloc(sizeof(double)*n);
    vaux2 = malloc(sizeof(double)*n);


    //if(id == 0) printf("\n-------------- AFUERA\n");
    while(ban != 'T' && sweeps < maxsweeps){
	no_rot = 0;

	//if(id == 0) printf("\n-------------- WHILE\n");
	for(int k = 0; k<=2*(np-1); k++){

	    //printf("\n-------------- FOR, sweep = %d, id = %d, k = %d\n", sweeps, id, k);
	    //printf("s: %d, id: %d, k: %d, Ai = ", sweeps, id, k);
	    //impMat(1,m,A_i);
	    //printf("s: %d, id: %d, k: %d, Aj = ", sweeps, id, k);
	    //impMat(1,m,A_j);
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
	    } //fin else
	    // Intercambiamos columnas entre los nodos
	    rotatePossession(id, np, m, A_i, A_j, aux1, aux2, 0);
	    rotatePossession(id, np, n, V_i, V_j, vaux1, vaux2, 0);
	    rotatePossession(id, np, 1, coli, colj, vaux1, vaux2, 0);
	}//fin for


	//printf("\nid: %d salió del sweep %d", id, sweeps);
	sweeps += 1;

	if(no_rot == col_ort){
	    ban = 'T';
	}

    }//fin while

    // Esperamos a que todos los procesos terminen
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Juntamos las columnas de V
    double *V_ij, *VV;
    V_ij = malloc(sizeof(double)*2*n);
    if(id == 0){
	VV = malloc(sizeof(double)*n*n);
    }
    for(i = 0; i < n; i++){
	V_ij[i] = V_i[i];
    }
    for(i = 0; i < n; i++){
	V_ij[n + i] = V_j[i];
    }
    MPI_Gather(V_ij, 2*n, MPI_DOUBLE, VV, 2*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /////////////////////////
    /////// Imprimimos el resultado
    // Imprimimos V
    if(id == 0) {
	transp(n,n,VV,V);
	printf("\n================================================\n================================================");
	printf("\nNúmero de rotaciones:%d\n",no_rot);
	printf("Número de sweeps:%d\n",sweeps);
	printf("\n================================================");
	printf("\nV = \n");
	impMat(n,n,V);
    }

    if(id == 0){
	free(Amat);
	free(At);
	free(V);
    }
    free(A_i); free(A_j); free(Adist); free(Adt);
    free(V_i); free(V_j); free(Vdist); free(Vdt);
    free(aux1); free(aux2);
    free(vaux1); free(vaux2);
    //free(colauxi); free(colauxj);

    MPI_Finalize();
    return 0;
}
