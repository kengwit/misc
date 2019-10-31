#! /bin/bash

config_dir=../../../config
export PATH=$PATH:$config_dir
tools_dir=../../../tools
export PATH=$PATH:$tools_dir/bin/$(archdir)
if [ -x ./a.out  ] ; then
   a_out=./a.out
else
   a_out=./bin/$(archdir)/a.out
fi
echo Testing $a_out

clean.sh

for domain in Mandit Mandif
do
  $a_out $domain -nadapt 3 -pc_type jacobi -ksp_type cg -ksp_rtol 0.001
  hexvue.sh $domain 2
done


clean.sh
