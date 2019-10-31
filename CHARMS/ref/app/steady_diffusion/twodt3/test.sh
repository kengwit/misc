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

rm -f *hvu
for domain in ldomain ldomain2 
do
  $a_out $domain
  hexvue -c $domain.hexvue
done

$a_out ucsd -nadapt 3
hexvue.sh ucsd 2

$a_out adldomain -nadapt 3
hexvue.sh adldomain 2 
$a_out charms -nadapt 5
hexvue.sh charms 4

clean.sh
