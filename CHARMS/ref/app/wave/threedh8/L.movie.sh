if [ $# -lt 1 ] ; then
    echo Expected argument
    exit 1
fi

f=$1
create.frames.sh $f -hexvue_script $f.hexvue -image_format png

