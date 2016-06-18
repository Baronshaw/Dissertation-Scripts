function [] = correctedTraditionalVSSCurve()
% this function will go through the Scurves and for each grid calculate the
% traditional model performance statistics and will calculate the SCurve
% model performance statistics. I can then go through each day and create
% maps for 2001 and then for selected locations, create time series for
% each of the model perfomrnace statistics

% loop through each day of soft data from 2001
n = 3; exstr = '';
deltaT = 365; 
deltaT = 90; exstr = '_T90';
numbins = 10; 
minpnts = 150; 
negval = 0;
temp = datevec(datenum(2001,1,1):datenum(2001,12,31));
dayWHI = temp(:,1:3);
dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
allLocs = cell(length(dayWHIdisps),1);
TME = cell(length(dayWHIdisps),1); TNME = cell(length(dayWHIdisps),1);
TSE = cell(length(dayWHIdisps),1); TR = cell(length(dayWHIdisps),1);
TR2 = cell(length(dayWHIdisps),1); TMAE = cell(length(dayWHIdisps),1);
TNMAE = cell(length(dayWHIdisps),1); TRMSE = cell(length(dayWHIdisps),1);
TNRMSE = cell(length(dayWHIdisps),1);
SME = cell(length(dayWHIdisps),1); SNME = cell(length(dayWHIdisps),1);
SSE = cell(length(dayWHIdisps),1); SR = cell(length(dayWHIdisps),1);
SR2 = cell(length(dayWHIdisps),1); SMAE = cell(length(dayWHIdisps),1);
SNMAE = cell(length(dayWHIdisps),1); SRMSE = cell(length(dayWHIdisps),1);
SNRMSE = cell(length(dayWHIdisps),1);
for i = 1:length(dayWHIdisps)
    disp(i);
    load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
        dayWHIdisps(i),deltaT,n,minpnts,numbins,negval));
    len = size(valMod,1);
    
    % all space/time locations
    allLocs{i,1} = [CTMlocs repmat(datenum(dayWHI(i,:)),size(CTMlocs,1),1)];
    
    % traditional statistics
    TME{i,1} = cell2mat( arrayfun(@(j) (1./sum(~isnan(valObs(j,:)))).*nansum(valMod(j,:)-valObs(j,:)), ...
        1:len,'UniformOutput',false) )';
    TNME{i,1} = cell2mat( arrayfun(@(j) nansum(valMod(j,:)-valObs(j,:))/nansum(valObs(j,:)), ...
        1:len,'UniformOutput',false) )';
    TSE{i,1} = cell2mat( arrayfun(@(j) sqrt((1./(sum(~isnan(valObs(j,:)))-1)).*nansum((valMod(j,:)-valObs(j,:))-TME{i,1}(j)).^2), ...
        1:len,'UniformOutput',false) )';     
    TR{i,1} = cell2mat( arrayfun(@(j) corr(valObs(j,~isnan(valObs(j,:)))',valMod(j,~isnan(valObs(j,:)))','type','Pearson'), ...
        1:len,'UniformOutput',false) )';    
    TR2{i,1} = TR{i,1}.^2; 
    TMAE{i,1} = cell2mat( arrayfun(@(j) (1./sum(~isnan(valObs(j,:)))).*nansum(abs(valMod(j,:)-valObs(j,:))), ...
        1:len,'UniformOutput',false) )';
    TNMAE{i,1} = cell2mat( arrayfun(@(j) nansum(abs(valMod(j,:)-valObs(j,:)))./nansum(valObs(j,:)), ...
        1:len,'UniformOutput',false) )'; 
    TRMSE{i,1} = cell2mat( arrayfun(@(j) sqrt((1./sum(~isnan(valObs(j,:)))).*nansum(valMod(j,:)-valObs(j,:)).^2), ...
        1:len,'UniformOutput',false) )';
    TNRMSE{i,1} = cell2mat( arrayfun(@(j) sqrt(nansum(valMod(j,:)-valObs(j,:)).^2)./sqrt(nansum(valObs(j,:))), ...
        1:len,'UniformOutput',false) )';

    % average modeled
    midpnts = arrayfun(@(j) (perctile_data(j,2:end)-perctile_data(j,1:end-1))./2 ...
        + perctile_data(j,1:end-1), 1:length(perctile_data),'UniformOutput',false);
    midpnts = cell2mat(midpnts');
    len2 = size(midpnts,2);
    
    % S-curve statistics
    SME{i,1} = cell2mat( arrayfun(@(j) (1./len2).*sum(midpnts(j,:)-mean_Obs(j,:)), ...
        1:len,'UniformOutput',false) )';
    SNME{i,1} = cell2mat( arrayfun(@(j) sum(midpnts(j,:)-mean_Obs(j,:))/sum(mean_Obs(j,:)), ...
        1:len,'UniformOutput',false) )';
    SSE{i,1} = cell2mat( arrayfun(@(j) nanmean(var_Obs(j,:)), ...
        1:len,'UniformOutput',false) )'; 
    SR{i,1} = cell2mat( arrayfun(@(j) corr(mean_Obs(j,~isnan(mean_Obs(j,:)))',midpnts(j,~isnan(mean_Obs(j,:)))','type','Pearson'), ...
        1:len,'UniformOutput',false) )';    
    SR2{i,1} = SR{i,1}.^2;
    SMAE{i,1} = cell2mat( arrayfun(@(j) (1./len2).*sum(abs(midpnts(j,:)-mean_Obs(j,:))), ...
        1:len,'UniformOutput',false) )';
    SNMAE{i,1} = cell2mat( arrayfun(@(j) sum(abs(midpnts(j,:)-mean_Obs(j,:)))./sum(mean_Obs(j,:)), ...
        1:len,'UniformOutput',false) )';
    SRMSE{i,1} = cell2mat( arrayfun(@(j) sqrt((1./len2).*sum(midpnts(j,:)-mean_Obs(j,:)).^2), ...
        1:len,'UniformOutput',false) )';
    SNRMSE{i,1} = cell2mat( arrayfun(@(j) sqrt(sum(midpnts(j,:)-mean_Obs(j,:)).^2)./sqrt(sum(mean_Obs(j,:))), ...
        1:len,'UniformOutput',false) )';
    
end

allLocs = cell2mat(allLocs);
TME = cell2mat(TME); 
TNME = cell2mat(TNME);
TSE = cell2mat(TSE); 
TR = cell2mat(TR);
TR2 = cell2mat(TR2); 
TMAE = cell2mat(TMAE);
TNMAE = cell2mat(TNMAE); 
TRMSE = cell2mat(TRMSE);
TNRMSE = cell2mat(TNRMSE);
SME = cell2mat(SME); 
SNME = cell2mat(SNME);
SSE = cell2mat(SSE); 
SR = cell2mat(SR);
SR2 = cell2mat(SR2); 
SMAE = cell2mat(SMAE);
SNMAE = cell2mat(SNMAE); 
SRMSE = cell2mat(SRMSE);
SNRMSE = cell2mat(SNRMSE);

% save results
save(sprintf('../matfiles/correctTradVSSCurve%s.mat',exstr),'allLocs', ...
    'TME','TNME','TSE','TR','TR2','TMAE','TNMAE','TRMSE','TNRMSE', ...
    'SME','SNME','SSE','SR','SR2','SMAE','SNMAE','SRMSE','SNRMSE');

end