#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void main(int argc, char *argv[]){

	double sum;
	int k,n,lim;

	//MPI:
	int tid,nthreads;
	double loop_min,loop_max,global_sum;

	sum = 0;

	sscanf(argv[1],"%d",&n);
	sscanf(argv[2],"%d",&lim);


	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &tid);

	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

	loop_min=1 + ((long)tid+0)*((long)lim)/((long)nthreads);
	loop_max=(((long)tid+1)*((long)lim)/(long)nthreads);

	for(k=loop_max;k>=loop_min;k--){
		sum+=1.0/pow(k,n);
	}

	char *cpu_name;
	cpu_name    = (char *)calloc(80,sizeof(char));
	gethostname(cpu_name,80);

	printf("machine ip = %s | id = %i | n = %ld, ..., %ld | local_sum = %f\n", cpu_name, tid, (long)loop_min, (long)loop_max, sum);

	MPI_Reduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if(tid == 0){
		printf("Suma global = %f", global_sum);
	}

	MPI_Finalize();

}
