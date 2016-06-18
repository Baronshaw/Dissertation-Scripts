function [] = runall_Cluster_02b()
% this function will run all the data in the 02_mfiles_offset folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../02_mfiles_offset

% calculate the final mean trend for ALL the obs data
calcMeanTrendfinal();

end