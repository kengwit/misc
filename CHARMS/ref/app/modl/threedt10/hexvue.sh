#! /bin/bash

prob=$1; shift
if [ $# -ge 1 ] ; then
  step=$1; shift
  hvufile=$prob-mode$step.hvu
else
  step=
  hvufile=$prob.hvu
fi
if [ $# -ge 1 ] ; then
  scale=$1; shift
else
  scale=1
fi

hexvuefile=$prob.hexvue


cat <<EOF > $hexvuefile
! 
NEWFRAME WIDTH 480 HEIGHT 480 X 100 Y 100
edgeflag on
edgecolor white
shrink 1
hvuchattribs
hvuload $hvufile
!
HVUPICKSET 3 ! z
HVUSETDISP 0 1 2 sxsysz $scale $scale $scale 
! 
view axes off
view bground grey50
view iso
layer on   0
layer on   1
layer off 15
layer on  16
layer 0
!
fit
RENDER SHADE
render ldir -1 2 3
! fringe 1 18
VIEW REDRAW
! HVUQUIT
EOF


hexvue -c $hexvuefile 
rm -f $hexvuefile
