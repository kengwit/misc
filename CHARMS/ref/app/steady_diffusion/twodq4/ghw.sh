#! /bin/bash

x=0
i=0
while [ $i -lt 100 ] ; do
   echo $x $(expeval "ghw($x)")
   x=$(expeval "$x+0.01")
   i=$(($i+1))
done
