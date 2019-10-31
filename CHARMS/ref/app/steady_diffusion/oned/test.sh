#! /bin/bash

config_dir=../../../config
export PATH=$PATH:$config_dir
if [ -x ./a.out  ] ; then
   a_out=./a.out
else
   a_out=./bin/$(archdir)/a.out
fi
echo Testing $a_out

function cleanup {
  return
  domain=$1
  rm -f $domain*fpar $domain*.conn $domain*.fen $domain*.ebc
}
    
for nelems in 2 3 4
do
  echo Adaptive solution, $(($nelems+1)) nodes initially, uniform refinement
  grid.sh -nelems $nelems -len 100 -s -100 -k 1000 -ebci 100 -ebcf -133 -gradation 3 -adaptive
  prob=l$nelems
  echo Running $prob
  $a_out $prob | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup $prob
done

for nelems in 2 3
do
  echo Spike due to nonlinearly varying conductivity. Adaptive solution, $(($nelems+1)) nodes initially.
  grid.sh -nelems 3 -k "((x-4)+0.001)^2" -len 8 -prob spike$nelems -ebci 0 -ebcf 10000 -gradation "0.0001+0.2*(x-4)^2" -adaptive
  $a_out spike$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup spike$nelems
done

for nelems in 2 3
do
  echo Nonlinearly varying source term. Adaptive solution, $(($nelems+1)) nodes initially.
  grid.sh -nelems 3 -s "(6-x)^8" -len 8 -prob step$nelems -ebci -1000 -ebcf 3000 -gradation "0.01*(x-1.75)^4+0.0001" -adaptive
  $a_out step$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup step$nelems
done

for nelems in 2 4 8 
do
  echo Homogeneous EBC, positive source, $(($nelems+1)) nodes
  grid.sh -nelems $nelems -len 8 -s 1000 -k 10 -ebci 0 -ebcf 0
  $a_out l$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup l$nelems
done

for nelems in 8 
do
  echo Homogeneous EBC, negative source, $(($nelems+1)) nodes
  grid.sh -nelems $nelems -len 8 -s -1000 -k 10 -ebci 0 -ebcf 0
  $a_out l$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup l$nelems
done

for nelems in 17 
do
  echo Inhomogeneous EBC, positive source, $(($nelems+1)) nodes
  grid.sh -nelems $nelems -len 8 -s 100 -k 10 -ebci 100 -ebcf 0
  $a_out l$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup l$nelems
done

for nelems in 9 19 49
do
  echo Inhomogeneous EBC, positive source, $(($nelems+1)) nodes
  grid.sh -nelems $nelems -len 8 -s 100 -k 10 -ebci 0 -ebcf 100
  $a_out l$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup l$nelems
done

for nelems in 17 50 100
do
  echo Inhomogeneous EBC, no source, $(($nelems+1)) nodes
  grid.sh -nelems $nelems -len 8 -s 0 -k 10 -ebci 0 -ebcf 100
  $a_out l$nelems | xmgr -autoscale xy -noask -pipe 2>/dev/null
  cleanup l$nelems
done


