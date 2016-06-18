##!/bin/bash
bsub -x -q week -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_02b" -logfile "runall_Cluster_02b.out"
 
