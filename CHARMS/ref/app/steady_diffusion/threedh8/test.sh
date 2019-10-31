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

for domain in cutout ldomain 
do
  $a_out $domain
  hexvue.sh $domain 2
done

$a_out unit -nadapt 4
hexvue.sh unit 3
show_hier unit-hier3.egf

for domain in ringt ringf
do
  $a_out $domain -pc_type jacobi -ksp_type cg -nadapt 2
  hexvue.sh $domain 1
  show_hier ring-outline.egf $domain-hier1.egf
done

clean.sh
