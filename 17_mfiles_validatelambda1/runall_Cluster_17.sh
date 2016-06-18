##!/bin/bash
for ((i=1;i<=10;i++))
do
  bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "generatelambda(${i})" -logfile "generatelambda${i}.out"
  sleep 60
done
