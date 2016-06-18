##!/bin/bash
bsub -q day matlab -nodesktop -nosplash -singleCompThread -r "runall_Cluster_01" -logfile "runall_Cluster_01.out"
