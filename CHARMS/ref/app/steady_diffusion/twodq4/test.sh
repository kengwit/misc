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
for domain in ldomain ldomain2 spiral spiral2 cool cool3
do
  $a_out $domain
  hexvue -c $domain.hexvue
done

$a_out oct -nadapt 3
hexvue.sh oct 2

$a_out spiralad -nadapt 3
hexvue.sh  spiralad 2

$a_out sqcrack -nadapt 6
hexvue.sh sqcrack 5
show_hier sqcrack-hier5.egf

$a_out nzebc -nadapt 3
hexvue.sh nzebc 2

clean.sh
