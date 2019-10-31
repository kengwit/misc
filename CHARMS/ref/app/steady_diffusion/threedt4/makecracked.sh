#!/bin/bash

function usage_and_exit {
cat <<EOF
Usage: $(basename $0) hsize
EOF
exit 0
}

if [ $# -lt 1 ] ; then
    usage_and_exit
fi
hsize=$1

prob=cracked
slist=""

args=""
while [ $# -ge 1 ] 
do
   a=$1; shift
   case $a in
    -help)
        usage_and_exit
        ;;
    -slist)
        if [ $# -ge 1 ] ; then
            slist=$1; shift
        else
            usage_and_exit
        fi
        ;;
    *)
        ta="$args $a"
        args="$ta"
        ;;
    esac
done

rm -f $prob.cmd
cat <<EOF > $prob.cmd
load "$prob.egf"
selc inside
new name "mesh_view" x 474 y 0 width 300 height 300
render wire
fit
prune_elongated 0
prune_small 0
! vltol 0.0001
! variable_vltol 
! use_layer_as_surf_tag
presplit "$hsize"
! subdiv_split
smesh
crease_angle 60
corner_angle 45
bess
! no print requested
model tshow
new name "show_view" x 800 y 0 width 300 height 300 model tshow wire 
view bground black
view top
view preserve_vup
show te
show ts_orient
view iso
fit
! interactive
make_tr 1
vdt $prob$hsize.vin
! quit
EOF
egfbess -c $prob.cmd

vdt.sh -gen_interior_pts -i $prob$hsize -dump_data > $prob$hsize.vdt

if [ x"$slist" = x ] ; then
  cmd="vdt2famuls.sh $prob$hsize.vdt $prob$hsize"
else
  cmd="vdt2famuls.sh $prob$hsize.vdt $prob$hsize -slist \"$slist\""
fi

echo
echo Running $cmd
eval $cmd


cat <<EOF > $prob$hsize.fpar
  obj algorithms
    obj steady_diffusion
      obj $prob$hsize ! instance
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
          param ebc_file = $prob$hsize.ebc
        end
        obj adapt
          param adaptive = false
        end
      end
    end
    obj hexvue
      obj hexvue
        obj gcell_groups
          obj gg1
            param implementation = default
          end
        end
        param hexvue_file = $prob$hsize.hvu
      end
    end
    obj gradation
      obj $prob$hsize
!        param mesh_size_func = 0.005+200*(x*y)^2
        param mesh_size_func = 0.0005+0.2*(x^2+y^2)
      end
    end
    obj refine
      obj $prob$hsize
        param true_hierarchical = false
        param ref_fraction = 0.25
      end
    end
  end

  obj gmeshes
    obj m1
      param fen_file = $prob$hsize.fen
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
      param type = solid_t4
      param conn_file = $prob$hsize-gg1.conn
    end
  end

  obj materials
    obj generic_heat_diffusion
      param conductivity = 100
      param internal_heat_generation_density = 10000
    end
  end

  obj les
    obj default_les
      param fles_solver = spetsc ! smeschach ! 
    end
  end
EOF


