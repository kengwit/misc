#! /bin/bash

if [ $# -lt 1 ] ; then
   list="oned twodq4 twodt3 twodt6 threedh8 threedt4 threedt10"
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
