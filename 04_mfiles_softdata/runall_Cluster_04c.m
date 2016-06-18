function [] = runall_Cluster_04c(SOFTIDX)
% this function will run all the data in the 04_mfiles_softdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../04_mfiles_softdata

% calculate the final long mean trend for ALL the soft data 
calcMeanTrendfinal_soft(SOFTIDX);

end