##!/bin/bash
bsub -x -q week /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_02a" -logfile "runall_Cluster_02a.out"
 
