Notas extra de instalación de OMPI
====================================

* IPs dentro de /etc/hosts de cada máquina (cambian cada que reiniciamos las máquinas de AWS)

```
54.69.144.172 host0
54.69.130.230 host1
54.69.144.222 host2
```

* Dentro de /etc/exports la parte `*(rw,sync)` no debe tener espacios.

* Según si es en master o en esclavos:

```
sudo apt-get install nfs-server
sudo apt-get install nfs-client

sudo mkdir /carpetaNFS
sudo mkdir /carpetaNodo
```

* Antes de sudo apt-get install build-essential, hay que hacer sudo apt-get update

* El punto es que la llave generada esté en .ssh/authorized_keys de cada máquina. Para ello juntamos todas las llaves
vía la carpeta compartida (simple cat de id_rsa.pub de cada máquina) y luego copiamos los permisos al .ssh/ de cada carpeta

* Hay que correr esto en la terminal para tener los comandos de OMPI sin dar rutas completas:

```
export PATH=/opt/openmpi-1.8.4/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi-1.8.4/lib
```

* Cuando corremos en varias máquinas hay que especificar qué mpirun hay que usar y el path dentro de los nodos

```
/opt/openmpi-1.8.4/bin/mpirun -np 8 -map-by core --host host1,host2 -path ~/carpetaNodo/hello_world/ hello-mpi.exe
```
