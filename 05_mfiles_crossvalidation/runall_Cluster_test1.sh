##!/bin/bash
for ((i=2;i<=4;i++))
do 
  bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012a/matlab -nodesktop -nosplash -singleCompThread -r "runall_test1(${i})" -logfile "runall_Cluster_test1${i}.out"
  sleep 30
done
