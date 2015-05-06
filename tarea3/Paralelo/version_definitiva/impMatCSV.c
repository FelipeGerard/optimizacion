#include <stdio.h>


void impMatCSV(int dimren, int dimcol, double *mat){
    int i,j;
    for(i=0;i<dimren;i++){
	for(j=0;j<dimcol;j++){
		printf("%1.16f\t", mat[i*dimcol+j]);
		if(j < dimcol-1) printf(", ");
	}
	printf("\n");
    }
}
