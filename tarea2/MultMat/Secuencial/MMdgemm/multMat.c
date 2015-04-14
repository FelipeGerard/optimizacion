

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
