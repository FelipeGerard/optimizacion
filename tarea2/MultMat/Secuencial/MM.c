#include <stdio.h>

void impMat(int dimren, int dimcol, int mat[dimren][dimcol]){
int i,j;
for(i=0;i<dimren;i++){
for(j=0;j<dimcol;j++){
if(j<dimcol-1){
printf("mat[%d][%d]=%d\t",i,j,mat[i][j]);
}
else{
printf("mat[%d][%d]=%d\n",i,j,mat[i][j]);
}
}

}
printf("\n");
}

void multMat(int dimren, int dimcomp,int dimcol, int mat1[dimren][dimcomp],int mat2[dimcomp][dimcol],int matres[dimren][dimcomp] ){
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
int a[dim1][dim2], b[dim2][dim3], c[dim1][dim3];
int i,j;


for(i=0;i<dim1;i++){
 for(j=0;j<dim2;j++){
  a[i][j]=2;
 }
}

for(i=0;i<dim2;i++){
 for(j=0;j<dim3;j++){
  b[i][j]=-5;
 }
}


impMat(dim1,dim2,a);

impMat(dim1,dim2,b);

multMat(dim1,dim2,dim3,a,b,c);

impMat(dim1,dim2,c);


}



