void *zerosMat(int dimren,int dimcol,double *mat){
    int i,j;
    for(i=0;i<dimren;i++){
	for(j=0;j<dimcol;j++){
	    mat[i*dimcol+j]=0;
	}
    }
    return NULL;
}
