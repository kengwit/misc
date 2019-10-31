#! /bin/bash

silent=t
nelems=15
len=1.0
ebci=0
ebcf=0
s=30000.0
k=100.0
prob_default=l
prob=
ref_fraction=0.5
adaptive=false
gradation=1

function usage_and_exit {
cat <<EOF
Usage: $(basename $0) [options]
Options:
  -silent         (be silent; default $silent)
  -nelems n       (number of elements; default $nelems)
  -len l          (length of the domain; default $len)
  -ebci v         (initial node EBC; default $ebci)
  -ebcf v         (final node EBC; default $ebcf)
  -s s            (heat source; default $s)
  -k k            (conductivity; default $k)
  -prob n         (problem name; default "$prob_default")
  -adaptive       (adaptive?; default "$adaptive")
  -ref_fraction r (refinement fraction; default "$ref_fraction")
  -gradation g    (gradation function; default "$len/$nelems")
EOF
exit 1
}

while [ $# -ge 1 ] 
do
   a=$1; shift
   case $a in
    -help)
        usage_and_exit
        ;;
    -nelems)
        if [ $# -ge 1 ] ; then
            nelems=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -len)
        if [ $# -ge 1 ] ; then
            len=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -ebci)
        if [ $# -ge 1 ] ; then
            ebci=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -ebcf)
        if [ $# -ge 1 ] ; then
            ebcf=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -s)
        if [ $# -ge 1 ] ; then
            s=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -k)
        if [ $# -ge 1 ] ; then
            k=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -prob)
        if [ $# -ge 1 ] ; then
            prob=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -ref_fraction)
        if [ $# -ge 1 ] ; then
            ref_fraction=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -gradation)
        if [ $# -ge 1 ] ; then
            gradation=$1; shift
        else
            usage_and_exit
        fi
        ;;
    -adaptive)
        adaptive=true
        ;;
    -silent)
        silent=true
        ;;
    *)
        ta="$args $a"
        args="$ta"
        ;;
    esac
done

prob_default=$prob_default$nelems
if [ x"$prob" = x ] ; then 
  domain=$prob_default
else
  domain=$prob
fi
fpar_file=$domain.fpar
fen_file=$domain.fen
ebc_file=$domain.ebc
gg_file=$domain-gg1.conn

if [ $silent != t ] ; then echo Writing $fpar_file; fi
cat <<EOF  > $fpar_file
! This is a sample input file for famuls.
!   
  obj algorithms
    obj steady_diffusion
      obj $domain ! instance
        ! the algorithm should work with submeshes of this mesh
        param gmesh = m1
        ! and here are the submeshes on which it should operate
        param gsubmeshes = sm1
        !
        obj gcell_groups
          obj gg1
            param implementation = default
            param material = generic_heat_diffusion
          end
        end
        !
        param les = default_les
        obj essential_bcs
          param ebc_file = $ebc_file
        end
        param print_sol = true
        obj adapt
          param adaptive = $adaptive
        end
      end
    end
    obj errest
      obj $domain
        obj gcell_groups
          obj gg1
            param implementation = default
          end
        end
      end
    end
    obj gradation
      obj $domain
        param mesh_size_func = $gradation
      end
    end
    obj refine
      obj $domain
        param ref_fraction = $ref_fraction
      end
    end
  end

  obj gmeshes
    obj m1
      param fen_file = $fen_file
      param gsubmeshes = sm1
    end
  end

  obj gsubmeshes
    obj sm1
      param gcell_groups = gg1
    end
  end

  obj gcell_groups
    obj gg1
      param type = line_l2
      param conn_file = $gg_file
    end
  end

  obj materials
    obj generic_heat_diffusion
      param type = diffusion_iso
      param conductivity = $k
      param internal_heat_generation_density = $s
    end
  end

  obj les
    obj default_les
      param fles_solver = gep
    end
  end
EOF


if [ $silent != t ] ; then echo Writing $fen_file; fi
rm -f $fen_file
n=1;
while [ $n -le $(($nelems+1)) ] ;
do
  echo $n $(expeval "($n-1)*$len/$nelems") 0 0 >> $fen_file;
  n=$(($n+1)); 
done


if [ $silent != t ] ; then echo Writing $ebc_file; fi
rm -f $ebc_file
echo "function \"(1-x/$len)*$ebci+(x/$len)*$ebcf\"" >> $ebc_file
echo "box 0 0 0 0 0 0 inflate 0.001" >> $ebc_file
echo "box $len 0 0 $len 0 0 inflate 0.001" >> $ebc_file


if [ $silent != t ] ; then echo Writing $gg_file; fi
rm -f $gg_file
n=1;
while [ $n -le $nelems ] ;
do
  echo $n $(($n+1)) >> $gg_file
  n=$(($n+1)); 
done

if [ $silent != t ] ; then
  echo Exact solution for homogeneous boundary conditions:
  echo "   Max phi = " $(expeval "($s/$k)/2*($len/2)^2")
  echo "Run: a.out $domain"
fi

