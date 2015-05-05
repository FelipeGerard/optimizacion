#include <stdio.h>

void impMat(int dimren, int dimcol, double *mat){
    int i,j;
    for(i=0;i<dimren;i++){
	for(j=0;j<dimcol;j++){
	    printf("mat[%d]=%1.5f\t",(i*dimcol+j),mat[i*dimcol+j]);
	}
	printf("\n");
    }
    printf("\n");
}
