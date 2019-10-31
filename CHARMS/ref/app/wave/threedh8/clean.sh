#! /bin/bash


testdata=false
err=true
pod=false
hist=false
interactive=false
aftertest=false

function usage_and_exit {
cat <<EOF
Usage: $(basename $0) [ options ]

Options:

-interactive|-i   ask for confirmation (default $interactive)

EOF
}

while [ $# -ge 1 ] 
do
   a=$1; shift
   case $a in
    -i|-interactive)
        interactive=true
        ;;
    -help)
        usage_and_exit
        ;;
    *)
        input=$a
        ;;
    esac
done

function ask {
    file=$1
    if [ $interactive = true ] ; then
       echo "Remove $file? [y|n] \c" >&2
       confirm=false
       while read resp
       do
           if [ "$resp" = "y" -o "$resp" = "yes" -o "$resp" = "YES" -o "$resp" = "Y" ] ; then
               confirm=true
               break
           fi
           if [ "$resp" = "n" -o "$resp" = "no" -o "$resp" = "NO" -o "$resp" = "N" ] ; then
               confirm=false
               break
           fi
       done
       echo $confirm
    else
       echo true
    fi
}

function del {
    for n in $*
    do
        echo Removing $n
        if [ -d $n ] ; then
           if [ $(ask $n) = true ] ; then
                rm -rf $n
           fi
        else
           if [ $(ask $n) = true ] ; then
                rm -f $n
           fi
        fi
    done
}

del core
del *.hvu
del *-hier*.egf
del *.vdt
del *.xmgr
del *.xvu
del tmp.*
del *~
