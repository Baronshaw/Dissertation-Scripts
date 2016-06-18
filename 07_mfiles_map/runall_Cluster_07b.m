function [] = runall_Cluster_07b()
% this function will run all the data in the 07_mfiles_maps folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../07_mfiles_map

temp = datenum(1999,1,1):datenum(2010,12,31);
temp = datevec(temp);
idx = temp(:,2) == 1 | temp(:,2) == 7;
daydisp = unique(temp(idx,1:2),'rows');
 
% display maps created in BMEmaps
for i = 1:size(daydisp,1)
    %displayBMEmaps([daydisp(i,1:2) 1],0,1,1);
    %displayBMEmaps([daydisp(i,1:2) 1],0,0,1);
    displayBMEmaps([daydisp(i,1:2) 1],1,0,1);
end

% % calculate and display time series and randomly selected monitoring stations
% BMEtimeseries(100,0,1,1);
% BMEtimeseries(100,0,0,1);
% BMEtimeseries(100,1,0,1);

end