##!/bin/bash
for ((i=1;i<=10;i++))
do
  bsub -x -q week -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2013a/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_12d(${i})" -logfile "runall_Cluster12d${i}.out"
  sleep 30
done
