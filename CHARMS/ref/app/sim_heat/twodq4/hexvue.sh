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
edgecolor white
shrink 1
hvuchattribs
hvuload $hvufile
!
HVUPICKSET 0 ! von Mises
HVUSETDISP 0 0 0 sxsysz 0 0 1
! 
view axes off
view bground grey50
view top
layer on   0
layer on   1
layer off 11
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
