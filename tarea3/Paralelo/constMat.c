double *constMat(int dimren,int dimcol, double *col1,double *col2){
    double *cbind;
    int i,j;
    cbind=malloc(sizeof(double)*dimren*2);

    for(i=0;i<dimren;i++){
	for(j=0;j<2;j++){
	    if(mod((i*dimcol+j),2)==0){
		cbind[i*dimcol+j]=col1[i];

	    }
	    else{
		cbind[i*dimcol+j]=col2[i];
	    }
	}
    }
    return cbind;
}


