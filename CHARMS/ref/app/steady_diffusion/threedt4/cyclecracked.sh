#!/bin/bash

hsize=2

cat <<EOF > c.hexvue
! 
NEWFRAME WIDTH 320 HEIGHT 240 X 100 Y 100
edgeflag on
edgecolor grey50
hvuchattribs
hvuload cracked$hsize.hvu
!
HVUPICKSET 0 ! von Mises
! HVUSETDISP 0 0 0 sxsysz 0 0 3
! 
view iso
render wire
layer on   0
layer on   1
layer off 11
layer on  12
layer 0
!
fit
RENDER SHADE
render ldir -1 2 3
! fringe 1 18
VIEW REDRAW
! HVUQUIT
EOF

for sid in 1 2 3 4 5 6 7 8 9
do
  makecracked.sh $hsize -slist "$sid"
  a.out cracked$hsize
  hexvue -c c.hexvue
done


