#include <stdio.h>
#include <stdlib.h>


void impMat(int dimren, int dimcol, double mat[dimren][dimcol]){
int i,j;
for(i=0;i<dimren;i++){
for(j=0;j<dimcol;j++){
if(j<dimcol-1){
printf("mat[%d][%d]=%f\t",i,j,mat[i][j]);
}
else{
printf("mat[%d][%d]=%f\n",i,j,mat[i][j]);
}
}

}
printf("\n");
}


void matrizAleatoria(int dimren, int dimcol,double mat[dimren][dimcol])
{
  int i,j;
  for(i=0;i<dimren;i++)
  {
    for(j=0;j<dimcol;j++)
    {
	mat[i][j]=rand()%5;
    }
  }
}


void multMat(int dimren, int dimcomp,int dimcol, double mat1[dimren][dimcomp],double mat2[dimcomp][dimcol],double matres[dimren][dimcomp] ){
int i,j,k;

for(i=0;i<dimren;i++){
for(k=0;k<dimcol;k++){
 matres[i][k]=0;
 for(j=0;j<dimcomp;j++){
  matres[i][k]=matres[i][k]+mat1[i][j]*mat2[j][k];
 }
}
}

}


void main(){
int dim1=3,dim2=3,dim3=3;
double a[dim1][dim2], b[dim2][dim3], c[dim1][dim3];
int i,j;


matrizAleatoria(dim1,dim2,a);

impMat(dim1,dim2,a);

matrizAleatoria(dim2,dim3,b);

impMat(dim2,dim3,b);

multMat(dim1,dim2,dim3,a,b,c);

impMat(dim1,dim3,c);


}

























