#! /bin/bash

prob=$1; shift
if [ $# -ge 1 ] ; then
  step=$1; 
  hvufile=$prob-mode$step.hvu
else
  step=
  hvufile=$prob.hvu
fi
hexvuefile=$prob.hexvue
framename=$prob

maxscale=1000;

cx=152
cy=117
cz=0

nx=-1
ny=-1
nz=0
cat <<EOF > cut.view
212.935989 56.077068 109.842041
0.000000 0.000000 1.000000
0.933010 0.320440 0.163740
212.935989 56.077068 109.842041
431.375397 431.375397
EOF

nframes=10
ncycles=2
frame=0
while [ $frame -lt $(($ncycles*$nframes)) ]
do
    scale=$(expeval "sin(2*pi*$frame/$nframes)*$maxscale")
    echo Frame $frame, scale=$scale
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
HVUSETDISP 0 1 2 sxsysz $scale $scale $scale 
! 
view axes off
view bground white
layer off  0
layer on   1
layer off 15
layer off 16
layer 0
!
hvucutn $nx $ny $nz
hvucutc $cx $cy $cz
hvucutlr 1
hvucrcut
hvucutn 0 0 1
hvucutc $cx $cy 123
hvucutlr 1
hvucrcut
view load cut.view
fringe 0 0.03
RENDER SHADE
! render ldir -1 2 3
VIEW REDRAW
image $prob$frame.gif
HVUQUIT
EOF

    hexvue -c $hexvuefile 
#rm -f $hexvuefile
    frame=$(($frame+1))
done
