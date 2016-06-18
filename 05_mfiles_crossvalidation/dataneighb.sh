##!/bin/bash
for ((i=1999;i<=2010;i++))
do
  bsub -x -q week -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_05f(${i})" -logfile "dataneighb${i}.out"
  sleep 30
done
