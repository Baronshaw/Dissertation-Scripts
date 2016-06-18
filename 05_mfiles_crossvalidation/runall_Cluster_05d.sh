##!/bin/bash
bsub -x -q day -n 3 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_05d" -logfile "runall_Cluster05d.out"
