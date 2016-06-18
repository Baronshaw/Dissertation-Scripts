##!/bin/bash
for ((i=144;i<=144;i++))
do
  bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_08b(${i})" -logfile "runall_Cluster_08b_${i}.out"
  sleep 30
done
