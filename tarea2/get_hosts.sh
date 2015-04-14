#! /bin/bash

# Temp file with instance JSON
aws ec2 describe-instances > temp001

# Start cluster-1
instancias=$(< temp001 jq '.Reservations | .[] | .Instances | .[] | select(.SecurityGroups | .[].GroupName == "cluster-1") | .InstanceId' \
| awk 'BEGIN {r = ""} {r = r " " $0} END {print r}' | sed -e 's/"//g' -e 's/^ //')
echo $instancias
aws ec2 start-instances --instance-ids $(echo $instancias)
aws ec2 describe-instances > temp001

# Private IP address of master node
master_ip=$(< temp001 jq '.Reservations | .[] | .Instances | .[] | select(.SecurityGroups | .[].GroupName == "cluster-1") | select(.Tags | .[] | .Value == "master0") | .PrivateIpAddress' | sed 's/"//g')

# Public DNS of master node
master_pub_dns=$(< temp001 jq '.Reservations | .[] | .Instances | .[] | select(.SecurityGroups | .[].GroupName == "cluster-1") | select(.Tags | .[] | .Value == "master0") | .PublicDnsName' | sed 's/"//g')

# File of public DNSs of nodes (including master)
< temp001 jq '.Reservations | .[] | .Instances | .[] | select(.SecurityGroups | .[].GroupName == "cluster-1") | .PublicDnsName' | sed 's/"//g' \
> temp003

# New IPs and hosts for /etc/hosts
< temp001 \
| jq '.Reservations | .[] | .Instances | .[] | select(.SecurityGroups | .[].GroupName == "cluster-1") | .PrivateIpAddress' \
| sed 's/"//g' \
| awk -v master_ip="$master_ip" 'BEGIN { i = 1 } { if($1 == master_ip){print $1, "host0"} else {print $1, "host" i; i+=1}}' \
> temp002

# Delete old hosts from /etc/hosts
export SHELL=/bin/bash
parallel --nonall --slf temp003 "rm -f ./temp003"

# Add new IPs to /etc/hosts
parallel --nonall --slf temp003 --basefile temp002 "cat /etc/hosts | sed '/host[0-9]/d' > ./temp003; sudo sh -c 'cat ./temp002 ./temp003 > /etc/hosts'"

# Stop cluster-1
aws ec2 stop-instances --instance-ids $(echo $instancias)
