#!/bin/bash
#
cond=4
declare -a array=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
for k in $(seq 0 1)
do
    echo "Step $k of 1"
for i in $(seq 0 15)
do
    taskset -c $i ./exe_map_fma_native aux_$i.dat $cond $i>> time_out_$cond.7_native_fma.dat &
    array[$i]=$!
done
for kk in ${array[@]}
do
    wait $kk
done
done
#
