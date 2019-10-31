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

for domain in sq one 
do
  $a_out $domain
  hexvue.sh $domain
done

$a_out adsq -nadapt 4
hexvue.sh adsq 3 


$a_out adone -nadapt 8
hexvue.sh adone 7


clean.sh
