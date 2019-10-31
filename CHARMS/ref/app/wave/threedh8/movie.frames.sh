if [ $# -lt 1 ] ; then
    echo Expected argument
    exit 1
fi

f=$1

if [ ! -f $f.hexvue ] ; then 
    echo HexVue script $f.hexvue not found
    exit 1
fi

create.frames.sh $f -hexvue_script $f.hexvue -image_format gif -silent

