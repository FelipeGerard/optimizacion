void *matrizId(int dimren,int dimcol,double *mat){

int i,j;
  for(i=0;i<dimren;i++)
  {
    for(j=0;j<dimcol;j++)
    {
	if(i==j){
	mat[i*dimcol+j]=1;
	}
	else{
	mat[i*dimcol+j]=0;
	}
    }
  }


}
