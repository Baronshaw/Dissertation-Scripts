##!/bin/bash
bsub -x -q hour -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_08a" -logfile runall_Cluster_08a.out"

