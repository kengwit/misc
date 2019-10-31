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
  maxfringe=$1; shift
else
  maxfringe=1
fi

hexvuefile=$prob.hexvue
framename=$prob

maxscale=0.51;

cx=162
cy=117
cz=0
cxd=162.1
cyd=117.1
czd=0

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
cat <<EOF > cut.view
211.122665 59.841419 107.086624
-0.732020 0.681254 0.006300
-0.820997 0.094430 0.563070
211.122665 59.841419 107.086624
404.241425 404.241425
EOF
cat <<EOF > cut.view
215.261627 63.738029 109.841454
-0.732020 0.681254 0.006300
-0.820997 0.094430 0.563070
207.367401 90.172089 148.866791
143.001648 143.001648
EOF


nframes=30
ncycles=1
frame=0
while [ $frame -lt $(($ncycles*$nframes)) ]
do
    scale=$(expeval "sin(2*pi*$frame/$nframes)*$maxscale")
    echo Frame $frame, scale=$scale

        if [ ${frame} -lt 10 ] ; then
                fnum=000${frame}
        elif [ ${frame} -lt 100 ] ; then
                fnum=00${frame}
        elif [ ${frame} -lt 1000 ] ; then
                fnum=0${frame}
        else
                fnum=${frame}
        fi

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
hvucutn $nx $ny $nz
hvucutc $cxd $cyd $czd
hvucutlr 1
hvucrcut
hvucutn 0 0 1
hvucutc $cx $cy 123
hvucutlr 1
hvucrcut
view load cut.view
!fringe 0 0.03
RENDER SHADE
render ldir -1 -2 3
VIEW REDRAW
image ${prob}-mode${step}-$fnum.gif
HVUQUIT
EOF

    hexvue -c $hexvuefile 
#rm -f $hexvuefile
    frame=$(($frame+1))
done
