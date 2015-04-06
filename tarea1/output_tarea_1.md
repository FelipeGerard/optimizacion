------
Title: Tarea 1. Algoritmos de gran escala
Date: 6 de abril de 2015
Author: Andrés Villaseñor, Carlos González, Felipe Gerard
------

# 0. Setup

Dado que los tres integrantes del equipo tenemos Macs y que manejar los puertos de máquinas simuladas en Docker no es sencillo, optamos por utilizar máquinas de Amazon. Generamos tres instancias t2.micro, de las cuales una es el master y dos son esclavas. Configuramos esas máquinas y en ellas hicimos los ejercicios. Cabe mencionar que no son multicore, por lo que no hay ganancias si usamos varios procesos en una máquina. 

Utilizamos host0 como el master y host1 y host2 como esclavas. Por alguna razón no nos funcionó usar -hostfile, así que optamos por dar la lista de las máquinas a utilizar. Tendremos que investigar por qué no funcionó.

# 1. Hello world

Diversos ejercicios par ver que funcionara el cluster.

### hello-mpi.exe en el master (8 cores simulados)

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/hello_world$ mpirun -np 8 -H localhost -path ~/carpetaNodo/hello_world/ hello-mpi.exe
hello MPI user: from process = 6 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 0 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 4 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 3 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 2 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 1 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 7 on machine=ip-172-31-37-235, of NCPU=8 processes
hello MPI user: from process = 5 on machine=ip-172-31-37-235, of NCPU=8 processes
```

### hello-mpi.exe en paralelo en dos máquinas (8 cores simulados)

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/hello_world$ /opt/openmpi-1.8.4/bin/mpirun -np 8 -map-by core --host host1,host2 -path ~/carpetaNodo/hello_world/ hello-mpi.exe
hello MPI user: from process = 6 on machine=ip-172-31-37-236, of NCPU=8 processes
hello MPI user: from process = 2 on machine=ip-172-31-37-234, of NCPU=8 processes
hello MPI user: from process = 5 on machine=ip-172-31-37-236, of NCPU=8 processes
hello MPI user: from process = 3 on machine=ip-172-31-37-234, of NCPU=8 processes
hello MPI user: from process = 7 on machine=ip-172-31-37-236, of NCPU=8 processes
hello MPI user: from process = 1 on machine=ip-172-31-37-234, of NCPU=8 processes
hello MPI user: from process = 4 on machine=ip-172-31-37-236, of NCPU=8 processes
hello MPI user: from process = 0 on machine=ip-172-31-37-234, of NCPU=8 processes
mpi_user@ip-172-31-37-235:~/carpetaNFS/hello_world$
```

### hello-mpi-2.exe en el master (un solo core)

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/hello_world$ mpirun -H localhost hello-mpi-2.exe
0.000 tid=0 : hello MPI user: machine=ip-172-31-37-235 [NCPU=1]
```

### hello-mpi-2.exe en dos máquinas (8 cores simulados)

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/hello_world$ /opt/openmpi-1.8.4/bin/mpirun -np 8 -map-by core --host host1,host2 -path ~/carpetaNodo/hello_world/ hello-mpi-2.exe
0.000 tid=5 : hello MPI user: machine=ip-172-31-37-236 [NCPU=8]
0.000 tid=4 : hello MPI user: machine=ip-172-31-37-236 [NCPU=8]
0.000 tid=6 : hello MPI user: machine=ip-172-31-37-236 [NCPU=8]
0.000 tid=7 : hello MPI user: machine=ip-172-31-37-236 [NCPU=8]
0.000 tid=3 : hello MPI user: machine=ip-172-31-37-234 [NCPU=8]
0.000 tid=0 : hello MPI user: machine=ip-172-31-37-234 [NCPU=8]
0.000 tid=1 : hello MPI user: machine=ip-172-31-37-234 [NCPU=8]
0.000 tid=2 : hello MPI user: machine=ip-172-31-37-234 [NCPU=8]
```

# 2. Aproximación de serie $\sum_{n=1}^\infty \frac{1}{n^k}$

### Versión sin OMPI (hasta 1,000,000)

El número no es correcto. ¿Habrá algún *overflow*?

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time ./rzf.exe 2 1000000
archivo:./rzf.exe
exponente:2
exponente:2
limite:1000000

9.869598

real	0m0.013s
user	0m0.013s
sys	0m0.000s
```

### Versión original con OMPI en master (un core)

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time mpirun -H localhost rzfp_ej.exe 2 1000000000 | sort
1.644934

real	0m12.847s
user	0m12.588s
sys	0m0.111s
```

### Versión original con OMPI en paralelo (8 cores simulados en dos máquinas)

¡Algo anda mal! ¿Por qué hay un *overflow*? --> El problema es que los enteros no aguantan $10^9$ y dan la vuelta a los número negativos, así que tendremos que usar *long* y un poco más de formato.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ /opt/openmpi-1.8.4/bin/mpirun -np 8 -map-by core --host host1,host2 -path ~/carpetaNodo/rzf rzfp_ej.exe 2 1000000000
0.000000
0.000000
1.644934
0.000000
0.000000
inf
0.000000
0.000000
```

### Versión usando MPI_Reduce en master (un core)

Además de manejar los *long*, mejoramos el formato y agregamos MPI_Reduce para obtener la suma total sin importar el número de nodos ni el número máximo a sumar.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time mpirun -H localhost rzfp_modif.exe 2 1000000000 | sort -t "|" -k2,2
Suma global = 1.644934
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 1000000000 | local_sum = 1.644934

real	0m13.127s
user	0m12.848s
sys	0m0.130s
```

### Versión usando MPI_Reduce en master (4 cores simulados)

No ganamos nada con respecto a usar un solo core porque la máquina es de un solo core.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time mpirun -np 4 -H localhost rzfp_modif.exe 2 1000000000 | sort -t "|" -k2,2
Suma global = 1.644934
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 250000000 | local_sum = 1.644934
machine ip = ip-172-31-37-235 | id = 1 | n = 250000001, ..., 500000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 2 | n = 500000001, ..., 750000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 3 | n = 750000001, ..., 1000000000 | local_sum = 0.000000

real	0m13.129s
user	0m12.908s
sys	0m0.088s
```

### Versión usando MPI_Reduce en master (10 cores simulados)

Nuevamente no hay ganancias con respecto a usar un solo proceso.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time mpirun -np 10 -H localhost rzfp_modif.exe 2 1000000000 | sort -t "|" -k2,2
Suma global = 1.644934
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 100000000 | local_sum = 1.644934
machine ip = ip-172-31-37-235 | id = 1 | n = 100000001, ..., 200000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 2 | n = 200000001, ..., 300000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 3 | n = 300000001, ..., 400000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 4 | n = 400000001, ..., 500000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 5 | n = 500000001, ..., 600000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 6 | n = 600000001, ..., 700000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 7 | n = 700000001, ..., 800000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 8 | n = 800000001, ..., 900000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 9 | n = 900000001, ..., 1000000000 | local_sum = 0.000000

real	0m13.272s
user	0m13.060s
sys	0m0.098s
```

### Versión usando MPI_Reduce en paralelo (2 cores en 2 máquinas)

Hay muchísimas ganancias. Toma un poco más de la mitad que usando una máquina.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time /opt/openmpi-1.8.4/bin/mpirun --host host0,host1 -path /home/mpi_user/carpetaNodo/rzf rzfp_modif.exe 2 1000000000 | sort -t "|" -k2,2
Suma global = 1.644934
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 500000000 | local_sum = 1.644934
machine ip = ip-172-31-37-234 | id = 1 | n = 500000001, ..., 1000000000 | local_sum = 0.000000

real	0m7.177s
user	0m0.008s
sys	0m0.009s
```

### Versión usando MPI_Reduce en paralelo (4 cores simulados en 2 máquinas)

No hay ganancias sobre usar un proceso por máquina (y por lo tanto por core, en este caso)

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time /opt/openmpi-1.8.4/bin/mpirun -np 4 -map-by core --host host0,host1 -path /home/mpi_user/carpetaNodo/rzf rzfp_modif.exe 2 1000000000 | sort -t "|" -k2,2
Suma global = 1.644934
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 250000000 | local_sum = 1.644934
machine ip = ip-172-31-37-235 | id = 1 | n = 250000001, ..., 500000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 2 | n = 500000001, ..., 750000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 3 | n = 750000001, ..., 1000000000 | local_sum = 0.000000

real	0m7.305s
user	0m0.008s
```

### Versión usando MPI_Reduce en paralelo (10 cores simulados en 2 máquinas)

Nuevamente no hay ganancias porque las máquinas no son multicore.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time /opt/openmpi-1.8.4/bin/mpirun -np 10 -map-by core --host host0,host1 -path /home/mpi_user/carpetaNodo/rzf rzfp_modif.exe 2 1000000000 | sort -t "|" -k2,2
Suma global = 1.644934
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 100000000 | local_sum = 1.644934
machine ip = ip-172-31-37-235 | id = 1 | n = 100000001, ..., 200000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 2 | n = 200000001, ..., 300000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 3 | n = 300000001, ..., 400000000 | local_sum = 0.000000
machine ip = ip-172-31-37-235 | id = 4 | n = 400000001, ..., 500000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 5 | n = 500000001, ..., 600000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 6 | n = 600000001, ..., 700000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 7 | n = 700000001, ..., 800000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 8 | n = 800000001, ..., 900000000 | local_sum = 0.000000
machine ip = ip-172-31-37-234 | id = 9 | n = 900000001, ..., 1000000000 | local_sum = 0.000000

real	0m7.234s
user	0m0.005s
sys	0m0.012s
```

### Versión usando MPI_Reduce en master, pero hasta el término 1000 (un core)

Este ejercicio es para comparar el tiempo de una máquina vs. varias para pocos términos.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time mpirun -H localhost rzfp_modif.exe 2 1000 | sort
Suma global = 1.643935
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 1000 | local_sum = 1.643935

real	0m0.379s
user	0m0.054s
sys	0m0.143s
```

### Versión usando MPI_Reduce en paralelo, pero hasta el término 1000 (10 cores simulados en 2 máquinas)

Como vemos, si usamos pocos términos el *overhead* de usar paralelismo es mayor que las ganancias, por lo que no vale la pena. Otro objetivo de mostrar menos términos es que al usar muchos sólo la suma correspondiente a los primeros era significativa. Con pocos términos podemos apreciar que realmente MPI_Reduce está haciendo su trabajo correctamente.

```
mpi_user@ip-172-31-37-235:~/carpetaNFS/rzf$ time /opt/openmpi-1.8.4/bin/mpirun -np 10 -map-by core --host host0,host1 -path /home/mpi_user/carpetaNodo/rzf rzfp_modif.exe 2 1000 | sort -t "|" -k2,2
Suma global = 1.643935
machine ip = ip-172-31-37-235 | id = 0 | n = 1, ..., 100 | local_sum = 1.634984
machine ip = ip-172-31-37-235 | id = 1 | n = 101, ..., 200 | local_sum = 0.004963
machine ip = ip-172-31-37-235 | id = 2 | n = 201, ..., 300 | local_sum = 0.001660
machine ip = ip-172-31-37-235 | id = 3 | n = 301, ..., 400 | local_sum = 0.000831
machine ip = ip-172-31-37-235 | id = 4 | n = 401, ..., 500 | local_sum = 0.000499
machine ip = ip-172-31-37-234 | id = 5 | n = 501, ..., 600 | local_sum = 0.000333
machine ip = ip-172-31-37-234 | id = 6 | n = 601, ..., 700 | local_sum = 0.000238
machine ip = ip-172-31-37-234 | id = 7 | n = 701, ..., 800 | local_sum = 0.000178
machine ip = ip-172-31-37-234 | id = 8 | n = 801, ..., 900 | local_sum = 0.000139
machine ip = ip-172-31-37-234 | id = 9 | n = 901, ..., 1000 | local_sum = 0.000111

real	0m0.854s
user	0m0.006s
sys	0m0.012s
```







