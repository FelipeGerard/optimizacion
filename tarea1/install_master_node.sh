
#! /bin/bash

### Creamos la instancia de Docker del nodo master desde la compu
docker run -it --name master-node -v ~/optimizacion/node0:/home/mpi_user ubuntu:latest /bin/bash

### Dentro de Docker

# Creamos mpi_user y cambiamos a ese usuario
useradd -d /home/mpi_user -G sudo mpi_user
passwd mpi_user
# $: 1234
su mpi_user

# Editamos /etc/hosts
sudo vi /etc/hosts
# Agregamos
# 192.168.133.100	host0
# 192.168.133.101	host1

# Instalamos nfs-server
sudo apt-get install nfs-server

# Creamos la carpeta para NFS (como root)
exit
sudo mkdir /carpetaNFS

# Editamos /etc/exports (como root)
echo "/home/mpi_user/carpetaNFS *(rw,sync)" >> /etc/exports

# Reiniciamos el servicio (igual como root)
service nfs-kernel-server restart

# Instalamos openssh-server
sudo apt-get install openssh-server
