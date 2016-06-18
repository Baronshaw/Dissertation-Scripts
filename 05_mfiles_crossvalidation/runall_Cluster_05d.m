function [] = runall_Cluster_05d()
% this function will run all the data in the 05_mfiles_softdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../05_mfiles_crossvalidation

% calculate the spatial 10 fold cross validation statistics
getXvalStatistics_10fold(0,1,1);
getXvalStatistics_10fold(0,0,1);
getXvalStatistics_10fold(1,0,1);

% show the results of the spatial 10 fold cross valdiation statistics
for i = 1:5
    displayXvalStatistics_10fold(i);
end

% % Xval functions/testing different scenarios
% Xvalforcedisolation(); % needs editing 2/17/2014
% displayXvalforcedisolation(); % needs editing 2/17/2014

end