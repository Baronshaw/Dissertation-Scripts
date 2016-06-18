function [] = runall_Cluster_08a()
% this function will run all the data in the 08_mfiles_QAQC folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../08_mfiles_QAQC

% the mock locations needed for QA analysis
mocklocations();
todeliver();

end