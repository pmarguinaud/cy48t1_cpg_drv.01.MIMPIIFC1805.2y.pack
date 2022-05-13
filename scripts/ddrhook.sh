#!/bin/bash

id=$1

./ddrhook.pl /scratch/work/marguina/tmp/arp.$id/0/drhook.txt /scratch/work/marguina/tmp/arp.$id/1/drhook.txt  > ddrhook01.$id.txt

./ddrhook.pl /scratch/work/marguina/tmp/arp.$id/0/drhook.txt /scratch/work/marguina/tmp/arp.$id/2/drhook.txt  > ddrhook02.$id.txt
./ddrhook.pl /scratch/work/marguina/tmp/arp.$id/2/drhook.txt /scratch/work/marguina/tmp/arp.$id/3/drhook.txt  > ddrhook23.$id.txt
./ddrhook.pl /scratch/work/marguina/tmp/arp.$id/3/drhook.txt /scratch/work/marguina/tmp/arp.$id/1/drhook.txt  > ddrhook31.$id.txt
./ddrhookN.pl /scratch/work/marguina/tmp/arp.$id/0/drhook.txt /scratch/work/marguina/tmp/arp.$id/2/drhook.txt  /scratch/work/marguina/tmp/arp.$id/1/drhook.txt > ddrhook021.$id.txt
