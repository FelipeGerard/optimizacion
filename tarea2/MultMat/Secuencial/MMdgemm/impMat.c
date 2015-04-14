#include <stdio.h>


void impMat(int dimren, int dimcol, double *mat){
int i,j;
for(i=0;i<dimren;i++){
for(j=0;j<dimcol;j++){
if((i*dimcol+j)%dimcol<dimcol-1){
printf("mat[%d]=%f\t",(i*dimcol+j),mat[i*dimcol+j]);
}
else{
printf("mat[%d]=%f\t\n",(i*dimcol+j),mat[i*dimcol+j]);
}

}
}
printf("\n");
}
