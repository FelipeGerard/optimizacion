void *asignaCol(int dimren,int dimcol,int col, double *columna,double *mat){
    int i,j;

    for(i=0;i<dimren;i++){
	for(j=0;j<dimcol;j++){
	    if((i*dimcol+j)==(i*dimcol+col)){
		mat[i*dimcol+j]=columna[i];
	    }
	}
    }
    return NULL;
}
