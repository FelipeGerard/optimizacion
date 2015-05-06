double prodPunto(int dim,double *vec1,double *vec2){
    double res;
    int i;
    res=0.0;
    for(i=0;i<dim;i++){
	res=res+vec1[i]*vec2[i];
    }
    return res;
}
