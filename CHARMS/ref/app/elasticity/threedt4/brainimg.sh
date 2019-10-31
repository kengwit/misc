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
edgeflag on
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
layer on   0
layer on   1
layer off 15
layer on  16
layer 0
!
render ldir -1 2 3
RENDER SHADE
fringe 0 0.48
hvupause 20
!==============
! first w/o skull
! front
view front
fit
hvupause 20
image $stem-front.png
! iso
view iso
fit
hvupause 20
image $stem-iso.png
! right
view right
fit
hvupause 20
image $stem-right.png
! top
view top
fit
hvupause 20
image $stem-top.png
!==============
! now w/ skull
load skullsurf.egf
! front
view front
fit
hvupause 20
image $stem-frontws.png
! iso
view iso
fit
hvupause 20
image $stem-isows.png
! right
view right
fit
hvupause 20
image $stem-rightws.png
! top
view top
fit
hvupause 20
image $stem-topws.png
hvupause 20
HVUQUIT
EOF


hexvue -c $hexvuefile 
rm -f $hexvuefile
