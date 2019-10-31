#! /bin/bash
#
#  defines
#
image_format=jpg
framename=frame
silent=f
first_frame=0
hexvue_script=hexvue.script
default_frame_script=_frame.sh
frame_script=$default_frame_script

#
#  usage function
#
function usage_and_exit {
echo "usage: $0 problem_prefix"
echo "  problem_prefix  "
echo "      The problem prefix is passed to the frame script as the first"
echo "      argument. It can be used to construct the name of a hexvue or for "
echo "      any other purpose; it is not used by $0 directly."
echo " One can also specify the following options:"
echo "  -image_format f"
echo "      One of jpg, gif, pict, tiff, eps (others might be also available;"
echo "      see the manual of the ImageMagick program \`import\'"
echo "  -frame_script fs"
echo "      This is the shell script that is run for each frame. The default"
echo "      is \"$frame_script\" (it is created by $0). The frame script is "
echo "      passed five arguments as described here: "
echo "      1. problem prefix,"
echo "      2. name of the hexvue script,"
echo "      3. name of a frame"
echo "      4. the frame number, and"
echo "      5. the frame image format."
echo "      The frame script is responsible (i) for creating the frame image, and "
echo "      (ii) for returning zero as its exit status when there should be "
echo "      additional frames, or non-zero when it created the last frame.  "
echo "      The default frame script looks like this (with some additional "
echo "      shell processing; see \"$frame_script\"):"
echo "          hexvue -l \${problem_prefix}\${frame}.hvu -c \${hexvue_script}"
echo "          mv junk.jpg \${framename}\${frame}.jpg"
echo "      As you can see, it runs hexvue, and it expects hexvue to write the"
echo "      file junk.jpg, which it then renames to \${framename}\${frame}.jpg"
echo "  -hexvue_script hs"
echo "      The HEXVUE script name is passed to the frame script as its"
echo "      second argument.  It is not used by $0 directly."
echo "  -frame_name f"
echo "      The frame name is the name of a JPEG file (or other graphics file)"
echo "      which constitutes an individual frame in the movie; the default name "
echo "      is \"$framename.N.jpg\", where \"N\" is the frame number.  "
echo "      Frame name is passed to the frame script as the third argument; it is not used 
by $0 directly."
echo "  -first_frame n"
echo "      The frames are numbered from zero if there is no previous frame present"
echo "      in the working directory; otherwise the "
echo "      frame number is set to be one more than the number of the"
echo "      last frame present.  The frame number is passed to the "
echo "      frame script as its fourth and last argument"
echo "  -silent"
echo "      Run silently"
echo "  -help"
echo "      Print this help text"
   exit 0
}

#
#  must have at least the problem prefix
#
if [ $# -lt 1 ] ; then
        usage_and_exit
fi
problem_prefix=$1; shift

#
#  get the additional arguments
#
while [ $# -ge 1 ] 
do
   a=$1; shift
   case $a in
      -frame_name)
         if [ $# -ge 1 ] ; then
            framename=$1; shift
         else
            usage_and_exit
         fi
         ;;
      -frame_script)
         if [ $# -ge 1 ] ; then
            frame_script=$1; shift
         else
            usage_and_exit
         fi
         ;;
      -first_frame)
         if [ $# -ge 1 ] ; then
            first_frame=$1; shift
         else
            usage_and_exit
         fi
         ;;
      -image_format)
         if [ $# -ge 1 ] ; then
            image_format=$1; shift
         else
            usage_and_exit
         fi
         ;;
      -hexvue_script)
         if [ $# -ge 1 ] ; then
            hexvue_script=$1; shift
         else
            usage_and_exit
         fi
         ;;
      -silent)
         silent=t
         ;;
      -help)
         usage_and_exit
         ;;
      *)
         usage_and_exit
         ;;
   esac
done


#
#  in case there are already some frames in here, try to figure out where to 
#  start with next frame numbers
#
frame=${first_frame}
while [ -f ${framename}${frame}.${image_format} ] 
do
        frame=$(($frame+1))
done
first_frame=$frame

#  the default frame script;  it is expected to create a frame
#  and to return either t when there is a next frame, or f when
#  the program is to stop.
#  The frame script is passed the following arguments:
#  1. problem_prefix
#  2. hexvue_script
#  3. frame name
#  4. frame number

cat <<EOF > $default_frame_script
if [ \$# -ne 5 ] ; then
        echo "usage: \$0 problem_prefix hexvue_script frame_name frame_number image_format"
        exit 1
fi
problem_prefix=\$1; shift
hexvue_script=\$1; shift
framename=\$1; shift
frame=\$1; shift
image_format=\$1; shift

if [ -f \${problem_prefix}\${frame}.hvu ] ;  then
   suffix=hvu
elif [  -f \${problem_prefix}\${frame}.xvu ] ;  then
   suffix=xvu
else
    echo No hexvue file found
    exit 1
fi
if [ -f \${problem_prefix}\${frame}.\${suffix} ] ;  then
        echo "post-processing \${problem_prefix}\${frame}..."
        hexvue -l \${problem_prefix}\${frame}.\${suffix}  -c \${hexvue_script}
        if [ \${frame} -lt 10 ] ; then
                fnum=000\${frame}
        elif [ \${frame} -lt 100 ] ; then
                fnum=00\${frame}
        elif [ \${frame} -lt 1000 ] ; then
                fnum=0\${frame}
        else
                fnum=\${frame}
        fi
        acounter=1
        mv junk.\${image_format} \${framename}\${fnum}.\${image_format}
        exit 0
else
        exit 1
fi
EOF
chmod +x $frame_script


if [ $silent != t ] ; then 
   echo "Problem prefix: $problem_prefix"
   echo "Frame name:     $framename"
   echo "First frame:    $first_frame"
   echo "HEXVUE script:  $hexvue_script"
   echo "Frame script:   $frame_script"
   echo "Continue? \c"
   read resp
   if [ x$resp = xn -o x$resp = no ] ; then
        exit 1
   fi
fi


frame=${first_frame}
next=0
while [ $next -eq 0 ] 
do
        $frame_script \
                ${problem_prefix} \
                ${hexvue_script} \
                ${framename} \
                $frame \
                $image_format 
        next=$?
        frame=$(($frame+1))
done
