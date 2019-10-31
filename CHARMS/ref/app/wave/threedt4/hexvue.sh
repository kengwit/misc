#! /bin/bash

prob=$1; shift
if [ $# -ge 1 ] ; then
  step=$1; 
  hvufile=$prob$step.xvu
else
  step=
  hvufile=$prob.xvu
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
HVUPICKSET 7 ! z
HVUSETDISP 0 1 2 sxsysz 1 1 1
! 
render wire
view axes off
view bground grey50
view nor -1 2 3
fit
render ldir 1 2 3
layer on   0
layer on   1
layer off 15
layer on  16
layer 0
!
RENDER shade dither
VIEW REDRAW
! HVUQUIT
EOF


hexvue -c $hexvuefile 
rm -f $hexvuefile
