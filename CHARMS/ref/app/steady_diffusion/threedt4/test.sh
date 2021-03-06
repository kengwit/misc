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

for domain in cracker ong scube adscube
do
  $a_out $domain -nadapt 2
  hexvue.sh $domain 1
done

$a_out block1 -nadapt 1
hexvue.sh block1 0

clean.sh
