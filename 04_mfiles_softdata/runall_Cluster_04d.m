function [] = runall_Cluster_04d()
% this function will run all the data in the 04_mfiles_softdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../04_mfiles_softdata

% gather all the soft data by year
gatherSoftData();

end