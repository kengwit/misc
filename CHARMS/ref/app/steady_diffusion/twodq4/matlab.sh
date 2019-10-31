cat <<EOF > init.m
addpath ../../../tools/matlab
famuls_init
EOF
echo "Running matlab: execute init.m as the first thing"
matlab -nodesktop -nosplash
