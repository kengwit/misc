#! /bin/bash

if [ $# -lt 1 ] ; then
   list="threedh8 threedt4"
else
   list=$*
fi

for dir in $(echo $list)
do
    cd $dir; 
    chmod +x *.sh
    ./test.sh
    cd ..
done
