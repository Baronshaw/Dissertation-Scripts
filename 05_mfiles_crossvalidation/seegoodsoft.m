function [] = seegoodsoft()
% see where the soft data was actually beneficial / not beneficial

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
ch = pd;
cMS = unique(ch(:,1:2),'rows');
temp = datevec(ch(:,3));
ckall = cell(12,1);
for i = 1:12
    ckall{i} = cell(length(cMS),1);
    for j = 1:length(cMS)
        disp([i j]);
        idxtemp = temp(:,1) == 2001 & temp(:,2) == i; 
        idx = cMS(j,1) == ch(:,1) & cMS(j,2) == ch(:,2) & idxtemp;
        ckall{i}{j} = ch(idx,:);
    end
end

% load all the LOOCV results for hard only 
for i = 1:12  
    load(sprintf('../matfiles/Xval_LOOCV_mon%d__nosoft_long_gauss.mat',i));
    zkall_soft{i,1} = zk_madd;
    zhall{i,1} = zh_Xval;
    allhards{i,1} = allhard;
    blah = 5;
end

% load all the LOOCV results for hard + soft

% get all the places where the soft decreased the absolute error and look
% at the % decrease in realtive error

% in the best areas, look at obs, hard est, hard neighb, lambda1s,
% lambda2s, hard and soft est

% is there a consistent trend with regions that have the best performance?
% Or times of the year that have the best performance?

% look as some S-curves. Do I see a clear pattern is the S-curves that led
% lambda1/lambda2 to be helpful?

% now do the same for the worst performance

% problem: are these too specific to be informative/see trends?










end