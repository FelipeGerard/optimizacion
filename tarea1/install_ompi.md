# IPs
54.69.144.172 host0
54.69.130.230 host1
54.69.144.222 host2

sudo apt-get install nfs-server
sudo apt-get install nfs-client

sudo mkdir /carpetaNFS
sudo mkdir /carpetaNodo

# Antes de sudo apt-get install build-essential, hay que hacer sudo apt-get update

# El punto es que la llave generada esté en .ssh/authorized_keys de cada máquina. Para ello juntamos todas las llaves
# vía la carpeta compartida y luego copiamos los permisos al .ssh/ de cada carpeta

export PATH=/opt/openmpi-1.8.4/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi-1.8.4/lib

/opt/openmpi-1.8.4/bin/mpirun -np 8 -map-by core --host host1,host2 -path ~/carpetaNodo/hello_world/ hello-mpi.exe
