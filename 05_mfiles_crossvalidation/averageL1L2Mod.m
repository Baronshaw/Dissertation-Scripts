function [] = averageL1L2Mod()
% this function will average lambda1, lambda2, and the modeled values for
% each space/time location. This "averaged" value can then be used to
% calculate more traditional statistics to use from comparison

% load each soft data day
n = 3; 
deltaT = 365; 
numbins = 10; 
minpnts = 150; 
negval = 0;
temp = datevec(datenum(2001,1,1):datenum(2001,12,31));
dayWHI = temp(:,1:3);
dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
sumLam1 = cell(length(dayWHIdisps),1);
sumMod = cell(length(dayWHIdisps),1);
sumLam2 = cell(length(dayWHIdisps),1);
allLocs = cell(length(dayWHIdisps),1);
MESCurve = cell(length(dayWHIdisps),1);
NMESCurve = cell(length(dayWHIdisps),1);
SESCurve = cell(length(dayWHIdisps),1);
RSCurve = cell(length(dayWHIdisps),1);
R2SCurve = cell(length(dayWHIdisps),1);
MAESCurve = cell(length(dayWHIdisps),1);
NMAESCurve = cell(length(dayWHIdisps),1); 
RMSESCurve = cell(length(dayWHIdisps),1);
NRMSESCurve = cell(length(dayWHIdisps),1);
L2SCurve = cell(length(dayWHIdisps),1);
for i = 1:length(dayWHIdisps)
    disp(i);
    load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
        dayWHIdisps(i),deltaT,n,minpnts,numbins,negval));
    
    % average modeled
    midpnts = arrayfun(@(j) (perctile_data(j,2:end)-perctile_data(j,1:end-1))./2 ...
        + perctile_data(j,1:end-1), 1:length(perctile_data),'UniformOutput',false);
    midpnts = cell2mat(midpnts');
    sumMod{i,1} = mean(midpnts,2);
    
    % average lambda1
    sumLam1{i,1} = mean(mean_Obs,2);
    
    % average lambda2
    sumLam2{i,1} = size(var_Obs,2).*sum(sqrt(var_Obs),2);
    
    % locations
    allLocs{i,1} = [CTMlocs repmat(datenum(dayWHI(i,:)),size(CTMlocs,1),1)];
    
    % all statistics per S-curve
    % is this different from bias(sumMod-sumLam1) versus bias(sumMode(bin)-sumLam1(bin))?
    MESCurvetemp = arrayfun(@(j) (1./length(midpnts(j,1:end))).*sum(midpnts(j,1:end)-mean_Obs(j,1:end)), ...
        1:length(midpnts),'UniformOutput',false);
    MESCurve{i,1} = cell2mat(MESCurvetemp');
    NMESCurvetemp = arrayfun(@(j) sum(midpnts(j,1:end)-mean_Obs(j,1:end))/sum(mean_Obs(j,1:end)), ...
        1:length(midpnts),'UniformOutput',false);
    NMESCurve{i,1} = cell2mat(NMESCurvetemp');
    SESCurvetemp = arrayfun(@(j) sqrt((1/9).*sum((midpnts(j,1:end)-mean_Obs(j,1:end))-MESCurve{i,1}(j)).^2), ...
        1:length(midpnts),'UniformOutput',false); % mind the hardcoded 1/9!
    SESCurve{i,1} = cell2mat(SESCurvetemp');
    RSCurvetemp = arrayfun(@(j) corr(mean_Obs(j,1:end)',midpnts(j,1:end)','type','Pearson'), ...
        1:length(midpnts),'UniformOutput',false);    
    RSCurve{i,1} = cell2mat(RSCurvetemp');
    R2SCurve{i,1} = RSCurve{i,1}.^2;
    MAESCurvetemp = arrayfun(@(j) (1./length(midpnts(j,1:end))).*sum(abs(midpnts(j,1:end)-mean_Obs(j,1:end))), ...
        1:length(midpnts),'UniformOutput',false);
    MAESCurve{i,1} = cell2mat(MAESCurvetemp');
    NMAESCurvetemp = arrayfun(@(j) sum(abs(midpnts(j,1:end)-mean_Obs(j,1:end)))./sum(mean_Obs(j,1:end)), ...
        1:length(midpnts),'UniformOutput',false);
    NMAESCurve{i,1} = cell2mat(NMAESCurvetemp');
    RMSESCurvetemp = arrayfun(@(j) sqrt((1./length(midpnts(j,1:end))).*sum(midpnts(j,1:end)-mean_Obs(j,1:end)).^2), ...
        1:length(midpnts),'UniformOutput',false);
    RMSESCurve{i,1} = cell2mat(RMSESCurvetemp');
    NRMSESCurvetemp = arrayfun(@(j) sqrt(sum(midpnts(j,1:end)-mean_Obs(j,1:end)).^2)./sqrt(sum(mean_Obs(j,1:end))), ...
        1:length(midpnts),'UniformOutput',false);
    NRMSESCurve{i,1} = cell2mat(NRMSESCurvetemp');
    L2SCurvetemp = arrayfun(@(j) (1./length(midpnts(j,1:end))).*sum(sqrt(var_Obs(j,1:end))), ...
        1:length(midpnts),'UniformOutput',false);
    L2SCurve{i,1} = cell2mat(L2SCurvetemp');
    

end
sumLam1 = cell2mat(sumLam1);
sumMod = cell2mat(sumMod);
sumLam2 = cell2mat(sumLam2);
allLocs = cell2mat(allLocs);
MESCurve = cell2mat(MESCurve);
NMESCurve = cell2mat(NMESCurve);
SESCurve = cell2mat(SESCurve);
RSCurve = cell2mat(RSCurve);
R2SCurve = cell2mat(R2SCurve);
MAESCurve = cell2mat(MAESCurve);
NMAESCurve = cell2mat(NMAESCurve); 
RMSESCurve = cell2mat(RMSESCurve);
NRMSESCurve = cell2mat(NRMSESCurve);
L2SCurve = cell2mat(L2SCurve);

save('../matfiles/avgL1L2Mod.mat','sumLam1','sumMod','sumLam2','allLocs', ...
    'MESCurve','NMESCurve','SESCurve','RSCurve','R2SCurve','MAESCurve', ...
    'NMAESCurve','RMSESCurve','NRMSESCurve','L2SCurve');

end