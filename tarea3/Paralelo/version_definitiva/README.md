
Compilar
-------------------
/ruta/a/mpi/bin/mpicc -std=c99 -o SVD_Jacobi_paralelo.o SVD_Jacobi_paralelo.c -lm -lblas


Correr
-------------------

* USO: `/ruta/a/mpi/bin/mpirun -H localhost -np num_procesos ./SVD_Jacobi_paralelo.o num_reng num_col max_sweeps`
* EJEMPLO: `/ruta/a/mpi/bin/mpirun -H localhost -np 2 ./SVD_Jacobi_paralelo.o 6 4 10`


NOTAS
-------------------

* Por ahora el método genera sus propias matrices, pero adentro se puede especificar alguna otra que se quiera.
