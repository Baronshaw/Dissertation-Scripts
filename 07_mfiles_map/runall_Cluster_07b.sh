##!/bin/bash
for ((i=1;i<=12;i++))
do
  bsub -x -q hour -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "BMEmaps(${i})" -logfile "runall_Cluster07b${i}.out"
  sleep 30
done
