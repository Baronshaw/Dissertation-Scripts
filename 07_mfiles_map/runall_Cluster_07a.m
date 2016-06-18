function [] = runall_Cluster_07a()
% this function will run all the data in the 07_mfiles_maps folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../07_mfiles_map

for i = 13:6:154
    WHIMS_est(i,2,1);
    WHIMS_est(i,2,0);
end

% % calculate and display time series and randomly selected monitoring stations
% BMEtimeseries(100,0,1,1);
% BMEtimeseries(100,0,0,1);
% BMEtimeseries(100,1,0,1);

end