load "mand12.egf"
selc inside
new name "mesh_view" x 474 y 0 width 300 height 300
render wire
fit
prune_elongated 0
prune_small 0
! vltol 0.0001
! variable_vltol 
! use_layer_as_surf_tag
! presplit 1
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
