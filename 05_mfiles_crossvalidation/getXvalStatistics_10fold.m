function [] = getXvalStatistics_10fold(soft,constant,gauss)
% this function will calculate the cross validation statistics for the 10
% fold cross validation

if nargin < 1, soft = 0; end % soft data or not
if nargin < 2, constant = 1; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 0, constr = '_long'; else constr = '_constant'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

for i = 1:10
    
    % loading results
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',i,softstr,constr,gaussstr));
    zkall{i,1} = zk_madd;
    subZall{i,1} = zh_Xval;
    ckxy{i,1} = ck(:,1:2);
    ckall{i,1} = ck(:,3);
    vkall{i,1} = vk;

    % all the statistics by fold
    idx = ~isnan(zk_madd);
    peari{i,1} = (corr(zk_madd(idx),zh_Xval(idx),'type','Pearson')).^2;
    speari{i,1} = (corr(zk_madd(idx),zh_Xval(idx),'type','Spearman')).^2;
    MEi{i,1} =nanmean((zk_madd(idx)-zh_Xval(idx)));
    MAEi{i,1} = nanmean(abs(zk_madd(idx)-zh_Xval(idx)));
    MSEi{i,1} = nanmean((zk_madd(idx)-zh_Xval(idx)).^2);
    RMSEi{i,1} = sqrt(nanmean((zk_madd(idx)-zh_Xval(idx)).^2));
    MRi{i,1} = nanmean(sqrt(vk(idx)));
    SEi{i,1} = sqrt( sum((mean(zk_madd(idx))-zk_madd(idx)).^2)./(sum(idx)-1) ) ./ sqrt(sum(idx));
    raui{i,1} = corr(zk_madd(idx),zh_Xval(idx),'type','Pearson');
    MSi{i,1} = nanmean( ((zk_madd(idx)-zh_Xval(idx)) - nanmean(zk_madd(idx)-zh_Xval(idx)))./nanvar(zk_madd(idx)-zh_Xval(idx)) );
    RMSSi{i,1} = sqrt(nanmean(( ((zk_madd(idx)-zh_Xval(idx)) - nanmean(zk_madd(idx)-zh_Xval(idx)))./nanvar(zk_madd(idx)-zh_Xval(idx)) ).^2));
    varesti{i,1} = nanvar(zk_madd(idx));
    varvali{i,1} = nanvar(zh_Xval(idx));
    fold_number{i,1} = i*ones(sum(idx),1);
    
end

zkall = cell2mat(zkall);
subZall = cell2mat(subZall);
ckall = cell2mat(ckall);
vkall = cell2mat(vkall);
fold_number = cell2mat(fold_number);
ckxy = cell2mat(ckxy);
[yr mo da] = datevec(ckall);
idx = ~isnan(zkall);

% all the statistics
pear = (corr(zkall(idx),subZall(idx),'type','Pearson')).^2;
spear = (corr(zkall(idx),subZall(idx),'type','Spearman')).^2;
ME =nanmean((zkall(idx)-subZall(idx)));
MAE = nanmean(abs(zkall(idx)-subZall(idx)));
MSE = nanmean((zkall(idx)-subZall(idx)).^2);
RMSE = sqrt(nanmean((zkall(idx)-subZall(idx)).^2));
MR = nanmean(sqrt(vkall(idx)));
SE = sqrt( sum((mean(zkall(idx))-zkall(idx)).^2)./(sum(idx)-1) ) ./ sqrt(sum(idx));
rau = corr(zkall(idx),subZall(idx),'type','Pearson');
MS = nanmean( ((zkall(idx)-subZall(idx)) - nanmean(zkall(idx)-subZall(idx)))./nanvar(zkall(idx)-subZall(idx)) );
RMSS = sqrt(nanmean(( ((zkall(idx)-subZall(idx)) - nanmean(zkall(idx)-subZall(idx)))./nanvar(zkall(idx)-subZall(idx)) ).^2));
varest = nanvar(zkall(idx));
varval = nanvar(subZall(idx));

% by year
years = NaN*ones(length(subZall),1);
for i = 1999:2010
    
    idx = i == yr;
    years(idx) = i;
    
    % all the statistics
    pearyr{i-1998,1} = (corr(zkall(idx),subZall(idx),'type','Pearson')).^2;
    spearyr{i-1998,1} = (corr(zkall(idx),subZall(idx),'type','Spearman')).^2;
    MEyr{i-1998,1} =nanmean((zkall(idx)-subZall(idx)));
    MAEyr{i-1998,1} = nanmean(abs(zkall(idx)-subZall(idx)));
    MSEyr{i-1998,1} = nanmean((zkall(idx)-subZall(idx)).^2);
    RMSEyr{i-1998,1} = sqrt(nanmean((zkall(idx)-subZall(idx)).^2));
    MRyr{i-1998,1} = nanmean(sqrt(vkall(idx)));
    SEyr{i-1998,1} = sqrt( sum((mean(zkall(idx))-zkall(idx)).^2)./(sum(idx)-1) ) ./ sqrt(sum(idx));
    rauyr{i-1998,1} = corr(zkall(idx),subZall(idx),'type','Pearson');
    MSyr{i-1998,1} = nanmean( ((zkall(idx)-subZall(idx)) - nanmean(zkall(idx)-subZall(idx)))./nanvar(zkall(idx)-subZall(idx)) );
    RMSSyr{i-1998,1} = sqrt(nanmean(( ((zkall(idx)-subZall(idx)) - nanmean(zkall(idx)-subZall(idx)))./nanvar(zkall(idx)-subZall(idx)) ).^2));
    varestyr{i-1998,1} = nanvar(zkall(idx));
    varvalyr{i-1998,1} = nanvar(subZall(idx));
    
end

% saving results
spatial_coordinate_x = ckxy(:,1);
spatial_coordinate_y = ckxy(:,2);
[yr mo da] = datevec(ckall);
time_t = yr*10000 + mo*100 + da;
observed_pm2p5_daily = subZall;
BME_mean = zkall;
BME_variance = vkall;

if constant == 1
    save('../matfiles/Xval_information_constant','spatial_coordinate_x',...
        'spatial_coordinate_y','time_t','observed_pm2p5_daily','years', ...
        'BME_mean','BME_variance','fold_number');
else
    save('../matfiles/Xval_information_long','spatial_coordinate_x',...
        'spatial_coordinate_y','time_t','observed_pm2p5_daily','years', ...
        'BME_mean','BME_variance','fold_number');
end

% save results
alltogether = {ME;MAE;MSE;pear;spear;MR;SE;rau;MS;RMSS;varval;varest};
save(sprintf('../matfiles/Xval_10fold%s%s%s_results.mat',softstr,constr,gaussstr),...
    'ME','MAE','MSE','pear','spear','RMSE','MR','SE','rau','MS','RMSS','varval','varest','alltogether',...
    'MEi','MAEi','MSEi','peari','speari','RMSEi','MRi','SEi','raui','MSi','RMSSi','varvali','varesti', ...
    'MEyr','MAEyr','MSEyr','pearyr','spearyr','RMSEyr','MRyr','SEyr','rauyr','MSyr','RMSSyr','varvalyr','varestyr');

end