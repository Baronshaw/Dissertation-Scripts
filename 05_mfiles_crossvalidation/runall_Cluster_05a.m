function [] = runall_Cluster_05a(FOLDIDX)
% this function will run all the data in the 05_mfiles_softdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../05_mfiles_crossvalidation

% calculate the spatial 10 fold cross validation 
getXvalidation_10fold(0,1,1,FOLDIDX);

end