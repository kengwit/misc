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

# 
for domain in ong slab cracker
do
  $a_out $domain -nadapt 2
  hexvue.sh $domain 1
done

clean.sh
