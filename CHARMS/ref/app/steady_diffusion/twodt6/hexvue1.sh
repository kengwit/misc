#! /bin/bash

prob=$1; shift
if [ $# -ge 1 ] ; then
  step=$1; 
  hvufile=$prob-adapt$step.hvu
else
  step=
  hvufile=$prob.hvu
fi
hexvuefile=$prob.hexvue

cat <<EOF > $hexvuefile
! 
NEWFRAME WIDTH 480 HEIGHT 480 X 100 Y 100
edgeflag on
edgecolor black
shrink 1
hvuchattribs
hvuload $hvufile
!
HVUPICKSET 0 ! von Mises
HVUSETDISP 0 0 0 sxsysz 0 0 0
! 
view axes off
view bground black
view top
layer on   0
layer on   1
layer off 15
layer on  16
layer 0
!
fit
RENDER normal
render ldir -1 2 3
fringe -0.23 0.23
view zoom 0.99
view zoom 0.99
view zoom 0.99
VIEW REDRAW
hvupause 3
image $prob-$step.png
image $prob-$step.eps
! HVUQUIT
EOF


hexvue -c $hexvuefile 
rm -f $hexvuefile
