##!/bin/bash
bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_07" -logfile "runall_Cluster07.out"
