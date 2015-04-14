#include <stdio.h>
#include <stdlib.h>


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

void* matrizAleatoria(int dimren, int dimcol,double *mat)
{
  int i,j;
  for(i=0;i<dimren;i++)
  {
    for(j=0;j<dimcol;j++)
    {
	mat[i*dimcol+j]=rand()%5;
    }
  }


}


void multMat(int dimren, int dimcomp,int dimcol,double *mat1,double *mat2, double *matres){
int i,j,k;
for(i=0;i<dimren;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      matres[i*dimcol+j]=0;
      for(k=0;k<dimcomp;k++)
      {
        matres[i*dimcol+j]+=mat1[i*dimcomp+k]*mat2[k*dimcomp+j];
      }
    }
  }

}







void main(){
int dim1=3,dim2=3,dim3=3;

double *a,*b,*c;

a=malloc(sizeof(double)*dim1*dim2);

b=malloc(sizeof(double)*dim2*dim3);

c=malloc(sizeof(double)*dim1*dim3);

matrizAleatoria(dim1,dim2,a);

impMat(dim1,dim2,a);


matrizAleatoria(dim2,dim3,b);

impMat(dim2,dim3,b);

multMat(dim1,dim2,dim3,a,b,c);


impMat(dim1,dim3,c);




}




