##!/bin/bash
for ((i=1;i<=1;i++))
do
  bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_05e(${i})" -logfile "runall_Cluster05e${i}.out"
  sleep 30
done
