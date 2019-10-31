#! /bin/bash

prob=$1; shift
if [ $# -ge 1 ] ; then
  step=$1; 
  stem=$prob-adapt$step
else
  step=
  stem=$prob
fi
hvufile=$stem.hvu
hexvuefile=$prob.hexvue

cat <<EOF > $hexvuefile
! 
NEWFRAME WIDTH 480 HEIGHT 480 X 100 Y 100
edgeflag off
edgecolor white
shrink 1
hvuchattribs
hvuload $hvufile
!
HVUPICKSET 3 ! z
HVUSETDISP 0 1 2 sxsysz 1 1 1
! 
view axes off
view bground white
layer off 0
layer off 1
layer off 15
layer on  16
layer 1
layer on 1
!
render ldir -1 2 3
RENDER SHADE
hvupause 10
! x-cuts
view front
hvucutn 1 0 0
hvucutc -0.6 0 0
hvucrcut
hvucutc -0.59 0 0
hvucrcut
HVUPICKSET 0
fringe -0.156 0.135
fit
hvupause 10
image $stem-xcutx.png
HVUPICKSET 1
fringe -0.1 0.1
fit
hvupause 10
image $stem-xcuty.png
HVUPICKSET 2
fringe -0.413 0
fit
hvupause 10
image $stem-xcutz.png
HVUPICKSET 3
fringe 0 0.45
fit
hvupause 10
image $stem-xcuta.png
hvupause 10
hvucutrac ! remove all cuts
! y-cuts
view right
hvucutn 0 1 0
hvucutc 0 0.35 0
hvucrcut
hvucutc 0 0.34 0
hvucrcut
HVUPICKSET 0
fringe -0.156 0.135
fit
hvupause 10
image $stem-ycutx.png
HVUPICKSET 1
fringe -0.1 0.1
fit
hvupause 10
image $stem-ycuty.png
HVUPICKSET 2
fringe -0.413 0
fit
hvupause 10
image $stem-ycutz.png
HVUPICKSET 3
fringe 0 0.45
fit
hvupause 10
image $stem-ycuta.png
hvupause 10
hvucutrac ! remove all cuts
! z-cuts
view top
hvucutn 0 0 1
hvucutc 0 0 22.17
hvucrcut
hvucutc 0 0 22.18
hvucrcut
HVUPICKSET 0
fringe -0.156 0.135
fit
hvupause 10
image $stem-zcutx.png
HVUPICKSET 1
fringe -0.1 0.1
fit
hvupause 10
image $stem-zcuty.png
HVUPICKSET 2
fringe -0.413 0
fit
hvupause 10
image $stem-zcutz.png
HVUPICKSET 3
fringe 0 0.45
fit
hvupause 10
image $stem-zcuta.png
hvupause 10
HVUQUIT
EOF


hexvue -c $hexvuefile 
rm -f $hexvuefile
