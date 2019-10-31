wfile=$1; shift
prob=$1

rm -f $prob.fen
rm -f $prob-gg*.conn

reading_Nodes=0
reading_Elements=0

while read s
do
    case "$s" in
    "Nodes")
        reading_Nodes=1
    ;;
    "Elements")
        reading_Elements=1
    ;;
    "&")
        reading_Nodes=0
        reading_Elements=0
    ;;
    *)
    if [ $reading_Nodes -ne 0 ] ; then
        echo "$s" >> $prob.fen
    elif [ $reading_Elements -ne 0 ] ; then 
        gg=$(echo "$s" | cut -d' ' -f22)
        echo $(echo "$s" | cut -d' ' -f2-21) >> $prob-gg$gg.conn
    else
        echo Ignoring "$s"
    fi
    ;;
    esac
done < $wfile
