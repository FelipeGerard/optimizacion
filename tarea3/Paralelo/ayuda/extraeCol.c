double *extraeCol(int dimren,int dimcol, int col,double *mat){
    int i,j;
    double *columna;

    columna=malloc(sizeof(double)*dimren);
    for(i=0;i<dimren;i++){
	for(j=0;j<dimcol;j++){
	    if((i*dimcol+j)==(i*dimcol+col)){
		columna[i]=mat[i*dimcol+j];
	    }
	}
    }
    return columna;
}


