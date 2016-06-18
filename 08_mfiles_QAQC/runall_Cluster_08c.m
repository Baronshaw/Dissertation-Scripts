function [] = runall_Cluster_08c()
% this function will run all the data in the 08_mfiles_QAQC folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../08_mfiles_QAQC

vals = datevec(datenum(1998,1:160,1));
tryruns = find(vals(:,2)==1|vals(:,2)==7);

% % QA time series for 1998
% QAplots(10:15,0,'overlap');
% close all;
% QAplots(10:15,1,'overlap');
% close all;

% % QA time series for the summer
% QAplots(102:104,0,'sum06');
% close all;
% QAplots(102:104,1,'sum06');
% close all;

% QA time series for the entire time period for a few locations
for i = 13:6:156
    QAplots(i:i+5,0,'whole',0);
    close all;
    QAplots(i:i+5,1,'whole',0);
    close all;
end

% % set of QA plots of maps
% for i = 3:length(tryruns)-1
%     for j = 0:2
%         disp(vals(tryruns(i),:));
%         QAplots2(tryruns(i),j);
%         close all;
%     end
% end

end