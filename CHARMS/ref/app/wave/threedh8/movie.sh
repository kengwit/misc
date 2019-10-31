if [ $# -lt 1 ] ; then 
    echo 'Need arguments: prob hexvue_script [suffix]'
    exit 1
else
    prob=$1; shift
fi
if [ $# -lt 1 ] ; then 
    echo 'Need arguments: prob hexvue_script [suffix]'
    exit 1
else
    hexvue_script=$1; shift
fi
if [ $# -lt 1 ] ; then 
    suff=$$
else
    suff=$1; shift
fi
if [ $# -lt 1 ] ; then 
    first=0
else
    first=$1; shift
fi

mov=$prob-$suff.gif

create.frames.sh $prob -hexvue_script $hexvue_script -image_format gif -silent -first_frame $first
whirlgif -time 50 -loop 5 fram*gif > $mov
konqueror $mov