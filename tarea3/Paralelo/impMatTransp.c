#include <stdio.h>

void impMatTransp(int dimren, int dimcol, double *mat){
    int i,j;
    for(j=0;j<dimcol;j++){
	for(i=0;i<dimren;i++){
	    printf("mat[%d] = %3.3f\t", i*dimcol+j,mat[i*dimcol+j]);
	}
	printf("\n");
    }
    printf("\n");
}
