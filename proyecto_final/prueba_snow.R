# ./configure --prefix /usr --enable-shared --enable-static
# make all install
# 
# sudo apt-get install r-cran-rmpi
# sudo R CMD INSTALL Rmpi --configure-args=--with-mpi=/usr/lib/openmpi
# sudo apt-get install r-cran-mpi

## MAC
# R CMD INSTALL Rmpi_0.6-3.tar.gz --configure-args=--with-mpi=`pwd`/OMPI
# install.packages("Rmpi", type="source",
#                  configure.args="--with-mpi=/Users/Felipe/algoritmos-de-gran-escala/paquetes/OMPI")
# install.packages('snow')
# R CMD INSTALL --configure-args="--with-Rmpi-include=`pwd`/OMPI/include --with-Rmpi-libpath=`pwd`/OMPI/lib --with-Rmpi-type=OPENMPI" Rmpi_0.6-3.tar.gz

#############################################
# brew install openmpi
# install.packages("Rmpi", type="source",
#############################################
                 
library(snow)
cl <- makeCluster(2, type = "MPI")
clusterEvalQ(cl, Sys.getenv("HOST"))

clusterApply(cl, 1:6, function(x) x+1)




stopCluster(cl)









