a.out brain5000bl -pc_type jacobi -ksp_type cg -ksp_rtol 1e-2 -ksp_monitor -whier 1 -nadapt 3 -wdispmag 1	
a.out brain5000blh1 -pc_type jacobi -ksp_type cg -ksp_rtol 1e-2 -ksp_monitor -whier 1 -nadapt 3 -wdispmag 1	
brainimg.sh brain5000blh1; brainimg.sh brain5000bl; brainimg.sh brain5000blh1 2; brainimg.sh brain5000bl 2
