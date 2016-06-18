##!/bin/bash
bsub -x -q day -n 3 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_11f" -logfile "runall_Cluster11f.out"
 
