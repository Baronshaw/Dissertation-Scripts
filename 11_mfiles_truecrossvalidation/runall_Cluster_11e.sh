##!/bin/bash
bsub -x -q hour -n 3 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_11e" -logfile "runall_Cluster11e.out"
 
