#! /bin/bash

hvufile=$1; shift
hexvuefile=_sd.hexvue

cat <<EOF > $hexvuefile
! 
NEWFRAME WIDTH 480 HEIGHT 480 X 100 Y 100
edgeflag on
edgecolor white
shrink 1
hvuchattribs
hvuload $hvufile
!
HVUPICKSET 7 ! z
HVUSETDISP 0 1 2 sxsysz 10 10 10
! 
view axes off
view bground grey50
view load sd.view
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
