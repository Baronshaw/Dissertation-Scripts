function [] = runall_Cluster_02a()
% this function will run all the data in the 02_mfiles_offset folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../02_mfiles_offset

% explore mean trend parameters
getMeanTrend();

% plot mean trend parameters
plotMeanTrend();

% final mean trend parameters
getplotMeanTrendfinal();

end