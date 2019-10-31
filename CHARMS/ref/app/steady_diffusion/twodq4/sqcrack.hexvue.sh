#! /bin/bash

prob=$1; shift
step=$1;

hexvuefile=prob
cat <<EOF >> hexvuefile=$1
! 
hvuload q.hvu
NEWFRAME WIDTH 320 HEIGHT 240 X 100 Y 100
! edgeflag on
! edgecolor white
! hvuchattribs
!
HVUPICKSET 0 ! von Mises
HVUSETDISP 0 0 0 sxsysz 0 0 10
! 
view load q.view
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