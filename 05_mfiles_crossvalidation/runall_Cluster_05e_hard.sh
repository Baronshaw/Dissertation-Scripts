##!/bin/bash
for ((i=1;i<=10;i++))
do
  bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_05e_hard(${i})" -logfile "runall_Cluster05e_hard${i}.out"
  sleep 30
done
