#! /bin/bash

### Creamos la instancia de Docker del nodo master desde la compu
docker run -it --name working-node -v ~/optimizacion/node1:/home/mpi_user ubuntu:latest /bin/bash

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

# Instalamos nfs-client
sudo apt-get install nfs-client

# Creamos la carpeta para NFS (en home de mpi_user) y 
sudo mkdir /carpetaNodo

# Montamos la carpeta de master en la del nodo
exit
sudo mount hostmaster:/carpetaNFS /carpetaNodo

# Reiniciamos el servicio (igual como root)
service nfs-kernel-server restart

# Instalamos openssh-server
sudo apt-get install openssh-server

# Loggear mpi_user con password
su mpi_user
ssh-keygen -t rsa
cd .ssh
cat id_rsa.pub >> authorized_keys
ssh-copy-id.
