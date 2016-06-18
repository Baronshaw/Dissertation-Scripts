function [] = runall_Cluster_11e()
% this function will perform all the cross-validation results for each of
% the folds

% load BME
cd ../BMELIB2.0b
startup();
cd ../11_mfiles_truecrossvalidation

%%% calculate results %%%

% gathering all the data
for i = 1:10   
    % loading results
    load(sprintf('../matfiles/Xval_true10fold_fold%d_hardonly.mat',i));
    zkall{i,1} = zk_madd;
    smoothingParam = [900000 300000 100 50];
    load(sprintf('../matfiles/meanTrend_true10fold_%d_%d_%d_%d_%d.mat',i,smoothingParam));
    zhall{i,1} = zh_Xval;
    ckall{i,1} = ch_Xval;
    vkall{i,1} = vk;  
    foldall{i,1} = i*ones(length(zkall{i,1}),1);
end
zkall = cell2mat(zkall);
idx = ~isnan(zkall); zkall = zkall(idx);
zhall = cell2mat(zhall); zhall = zhall(idx);
ckall = cell2mat(ckall); ckall = ckall(idx,:);
vkall = cell2mat(vkall); vkall = vkall(idx);
foldall = cell2mat(foldall); foldall = foldall(idx);

% overall daily
n_sites = size(unique(ckall(:,1:2),'rows'),1);
n_observations = size(ckall,1);
errors = zkall - zhall;
RMSE = sqrt( mean(errors.^2) );
MAE = mean(abs(errors));
ME = mean(errors);
r2 = (corr(zkall,zhall,'type','Pearson')).^2;
MS = mean(errors./sqrt(vkall));
RMSS = std(errors./sqrt(vkall));
MR = mean(sqrt(vkall));
std_est = std(zkall);
std_obs = std(zhall);
QAr2 = ( ( std_est.^2 + std_obs.^2 - (RMSE.^2-ME.^2) )./(2*std_est*std_obs) )^2;
X = unique(ckall(:,1:2),'rows')'; Y = unique(ckall(:,1:2),'rows')';
D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
temp = sort(D,2);
dists_info = mean(temp(:,2))/1000;

aberr = abs(errors./zhall);

% by station daily
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = unisID(i,1) == ckall(:,1) & unisID(i,2) == ckall(:,2);
    n_sites_sID(i) = size(unique(ckall(idx,1:2),'rows'),1);
    n_observations_sID(i) = size(ckall(idx,:),1);
    errors = zkall(idx) - zhall(idx);
    RMSE_sID(i) = sqrt( mean(errors.^2) );
    MAE_sID(i) = mean(abs(errors));
    ME_sID(i) = mean(errors);
    r2_sID(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
    MS_sID(i) = mean(errors./sqrt(vkall(idx)));
    RMSS_sID(i) = std(errors./sqrt(vkall(idx)));
    MR_sID(i) = mean(sqrt(vkall(idx)));
    std_est_sID(i) = std(zkall(idx));
    std_obs_sID(i) = std(zhall(idx));
    QAr2_sID(i) = ( ( std_est_sID(i).^2 + std_obs_sID(i).^2 - (RMSE_sID(i).^2-ME_sID(i).^2) )./(2*std_est_sID(i)*std_obs_sID(i)) )^2;
    X = unique(ckall(idx,1:2),'rows')'; Y = unique(ckall(idx,1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    if numel(temp)>1, dists_info_sID(i) = mean(temp(:,2))/1000; end
end

% by fold daily
unifold = unique(foldall);
for i = 1:length(unifold)
    
    idx = unifold(i) == foldall;
    
    n_sites_fold(i) = size(unique(ckall(idx,1:2),'rows'),1);
    n_observations_fold(i) = size(ckall(idx,:),1);
    errors = zkall(idx) - zhall(idx);
    RMSE_fold(i) = sqrt( mean(errors.^2) );
    MAE_fold(i) = mean(abs(errors));
    ME_fold(i) = mean(errors);
    r2_fold(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
    MS_fold(i) = mean(errors./sqrt(vkall(idx)));
    RMSS_fold(i) = std(errors./sqrt(vkall(idx)));
    MR_fold(i) = mean(sqrt(vkall(idx)));
    std_est_fold(i) = std(zkall(idx));
    std_obs_fold(i) = std(zhall(idx));
    QAr2_fold(i) = ( ( std_est_fold(i).^2 + std_obs_fold(i).^2 - (RMSE_fold(i).^2-ME_fold(i).^2) )./(2*std_est_fold(i)*std_obs_fold(i)) )^2;
    X = unique(ckall(idx,1:2),'rows')'; Y = unique(ckall(idx,1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    dists_info_fold(i) = mean(temp(:,2))/1000;
    
end

% by month daily
temp = datevec(ckall(:,3));
mos = temp(:,2);
unimos = unique(mos);
for i = 1:length(unimos)
    
    idx = unimos(i) == mos;
    
    n_sites_mos(i) = size(unique(ckall(idx,1:2),'rows'),1);
    n_observations_mos(i) = size(ckall(idx,:),1);
    errors = zkall(idx) - zhall(idx);
    RMSE_mos(i) = sqrt( mean(errors.^2) );
    MAE_mos(i) = mean(abs(errors));
    ME_mos(i) = mean(errors);
    r2_mos(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
    MS_mos(i) = mean(errors./sqrt(vkall(idx)));
    RMSS_mos(i) = std(errors./sqrt(vkall(idx)));
    MR_mos(i) = mean(sqrt(vkall(idx)));
    std_est_mos(i) = std(zkall(idx));
    std_obs_mos(i) = std(zhall(idx));
    QAr2_mos(i) = ( ( std_est_mos(i).^2 + std_obs_mos(i).^2 - (RMSE_mos(i).^2-ME_mos(i).^2) )./(2*std_est_mos(i)*std_obs_mos(i)) )^2;    
    X = unique(ckall(idx,1:2),'rows')'; Y = unique(ckall(idx,1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    dists_info_mos(i) = mean(temp(:,2))/1000;
    
end

% by season daily
clear idx;
idx{1} = mos == 12 | mos == 1 | mos == 2;
idx{2} = mos == 3 | mos == 4 | mos == 5;
idx{3} = mos == 6 | mos == 7 | mos == 8;
idx{4} = mos == 9 | mos == 10 | mos == 11;
for i = 1:length(idx)
    
    n_sites_sea(i) = size(unique(ckall(idx{i},1:2),'rows'),1);
    n_observations_sea(i) = size(ckall(idx{i},:),1);
    errors = zkall(idx{i}) - zhall(idx{i});
    RMSE_sea(i) = sqrt( mean(errors.^2) );
    MAE_sea(i) = mean(abs(errors));
    ME_sea(i) = mean(errors);
    r2_sea(i) = (corr(zkall(idx{i}),zhall(idx{i}),'type','Pearson')).^2;
    MS_sea(i) = mean(errors./sqrt(vkall(idx{i})));
    RMSS_sea(i) = std(errors./sqrt(vkall(idx{i})));
    MR_sea(i) = mean(sqrt(vkall(idx{i})));
    std_est_sea(i) = std(zkall(idx{i}));
    std_obs_sea(i) = std(zhall(idx{i}));
    QAr2_sea(i) = ( ( std_est_sea(i).^2 + std_obs_sea(i).^2 - (RMSE_sea(i).^2-ME_sea(i).^2) )./(2*std_est_sea(i)*std_obs_sea(i)) )^2;    
    X = unique(ckall(idx{i},1:2),'rows')'; Y = unique(ckall(idx{i},1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    dists_info_sea(i) = mean(temp(:,2))/1000;
    
end

% by year daily
temp = datevec(ckall(:,3));
yrs = temp(:,1);
uniyr = unique(yrs);
for i = 1:length(uniyr)
    
    idx = uniyr(i) == yrs;
    
    n_sites_year(i) = size(unique(ckall(idx,1:2),'rows'),1);
    n_observations_year(i) = size(ckall(idx,:),1);
    errors = zkall(idx) - zhall(idx);
    RMSE_year(i) = sqrt( mean(errors.^2) );
    MAE_year(i) = mean(abs(errors));
    ME_year(i) = mean(errors);
    r2_year(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
    MS_year(i) = mean(errors./sqrt(vkall(idx)));
    RMSS_year(i) = std(errors./sqrt(vkall(idx)));
    MR_year(i) = mean(sqrt(vkall(idx)));
    std_est_year(i) = std(zkall(idx));
    std_obs_year(i) = std(zhall(idx));
    QAr2_year(i) = ( ( std_est_year(i).^2 + std_obs_year(i).^2 - (RMSE_year(i).^2-ME_year(i).^2) )./(2*std_est_year(i)*std_obs_year(i)) )^2;    
    X = unique(ckall(idx,1:2),'rows')'; Y = unique(ckall(idx,1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    dists_info_year(i) = mean(temp(:,2))/1000;
        
end

% by EPA region daily
% go through original station locations
for i = 1:length(uniyr)
    load(sprintf('../datafiles/Observed_PM2p5/MasterDaily_PM2p5_%d.mat',uniyr(i)));
    alllon{i,1} = longitude;
    alllat{i,1} = latitude;
    allID{i,1} = location;
end
alllon = cell2mat(alllon);
alllat = cell2mat(alllat);
allID = cell2mat(allID);
[uniloc uniidx] = unique([alllon alllat],'rows');
unilon = uniloc(:,1);
unilat = uniloc(:,2);
uniID = allID(uniidx);

% convert them to coordinate system
load('../05_mfiles_crossvalidation/projexample.mat');
load('../09_mfiles_projections/Projections.mat');
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load('../09_mfiles_projections/Ellipsoids.mat');
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
cd ../09_mfiles_projections
projuniID = ell2lambertcc([unilon unilat],'whiproj2001'); 
cd ../11_mfiles_truecrossvalidation
% match estimations locations with station id
sIDall = NaN*ones(size(ckall,1),1);
for i = 1:length(projuniID)
    idx = projuniID(i,1) == ckall(:,1) & projuniID(i,2) == ckall(:,2);
    sIDall(idx) = uniID(i);
end
stateall = floor(sIDall./10^7);
% match id's with EPA region
% from https://aqs.epa.gov/aqsweb/codes/data/StateCountyCodes.csv
M = csvread('../05_mfiles_crossvalidation/StateCountyCodes_mod.csv');
temp = unique(M,'rows');
stateIDs = temp(:,1);
EPAregionsIDs = temp(:,2);
EPAall = NaN*ones(size(ckall,1),1);
for i = 1:length(stateIDs)
    idx = stateIDs(i) == stateall;
    EPAall(idx) = EPAregionsIDs(i);
end

% by state daily
unistate = unique(stateall);
for i = 1:length(unistate)
    
    idx = unistate(i) == stateall;
    
    n_sites_state(i) = size(unique(ckall(idx,1:2),'rows'),1);
    n_observations_state(i) = size(ckall(idx,:),1);
    errors = zkall(idx) - zhall(idx);
    RMSE_state(i) = sqrt( mean(errors.^2) );
    MAE_state(i) = mean(abs(errors));
    ME_state(i) = mean(errors);
    r2_state(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
    MS_state(i) = mean(errors./sqrt(vkall(idx)));
    RMSS_state(i) = std(errors./sqrt(vkall(idx)));
    MR_state(i) = mean(sqrt(vkall(idx)));
    std_est_state(i) = std(zkall(idx));
    std_obs_state(i) = std(zhall(idx));
    QAr2_state(i) = ( ( std_est_state(i).^2 + std_obs_state(i).^2 - (RMSE_state(i).^2-ME_state(i).^2) )./(2*std_est_state(i)*std_obs_state(i)) )^2;    
    X = unique(ckall(idx,1:2),'rows')'; Y = unique(ckall(idx,1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    if numel(temp) > 1, dists_info_state(i) = mean(temp(:,2))/1000; end
    
end

% actual EPA calculation
uniEPA = unique(EPAall);
for i = 1:length(uniEPA)-1
    
    idx = uniEPA(i) == EPAall;
    
    n_sites_EPA(i) = size(unique(ckall(idx,1:2),'rows'),1);
    n_observations_EPA(i) = size(ckall(idx,:),1);
    errors = zkall(idx) - zhall(idx);
    RMSE_EPA(i) = sqrt( mean(errors.^2) );
    MAE_EPA(i) = mean(abs(errors));
    ME_EPA(i) = mean(errors);
    r2_EPA(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
    MS_EPA(i) = mean(errors./sqrt(vkall(idx)));
    RMSS_EPA(i) = std(errors./sqrt(vkall(idx)));
    MR_EPA(i) = mean(sqrt(vkall(idx)));
    std_est_EPA(i) = std(zkall(idx));
    std_obs_EPA(i) = std(zhall(idx));
    QAr2_EPA(i) = ( ( std_est_EPA(i).^2 + std_obs_EPA(i).^2 - (RMSE_EPA(i).^2-ME_EPA(i).^2) )./(2*std_est_EPA(i)*std_obs_EPA(i)) )^2;
    X = unique(ckall(idx,1:2),'rows')'; Y = unique(ckall(idx,1:2),'rows')';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    temp = sort(D,2);
    dists_info_EPA(i) = mean(temp(:,2))/1000;
    
end

% averging data to yearly
[uniyrsID uniidx] = unique([ckall(:,1:2) yrs],'rows');
ckyrs = uniyrsID;
zkyrs = NaN*ones(length(uniyrsID),1);
zhyrs = NaN*ones(length(uniyrsID),1);
vkyrs = NaN*ones(length(uniyrsID),1);
foldyrs = foldall(uniidx);
EPAyrs = EPAall(uniidx);
for i = 1:size(uniyrsID,1);
    if mod(i,100)==0, disp([i length(uniyrsID)]); end
    idx = uniyrsID(i,1) == ckall(:,1) & uniyrsID(i,2) == ckall(:,2) & uniyrsID(i,3) == yrs;
    zkyrs(i) = mean(zkall(idx));
    zhyrs(i) = mean(zhall(idx));  
    vkyrs(i) = mean(vkall(idx))/sum(idx);
end

% overall yearly
n_sitesyrs = size(unique(uniyrsID(:,1:2),'rows'),1);
n_observationsyrs = size(uniyrsID,1);
errors = zkyrs - zhyrs;
RMSEyrs = sqrt( mean(errors.^2) );
MAEyrs = mean(abs(errors));
MEyrs = mean(errors);
r2yrs = (corr(zkyrs,zhyrs,'type','Pearson')).^2;
MSyrs = mean(errors./sqrt(vkyrs));
RMSSyrs = std(errors./sqrt(vkyrs));
MRyrs = mean(sqrt(vkyrs));
std_estyrs = std(zkyrs);
std_obsyrs = std(zhyrs);
QAr2yrs = ( ( std_estyrs.^2 + std_obsyrs.^2 - (RMSEyrs.^2-MEyrs.^2) )./(2*std_estyrs*std_obsyrs) )^2;

% by fold yearly
uniyrsfold = unique(foldyrs);
for i = 1:length(uniyrsfold)
    
    idx = uniyrsfold(i) == foldyrs;
    
    n_sites_foldyrs(i) = size(unique(uniyrsID(idx,1:2),'rows'),1);
    n_observations_foldyrs(i) = size(uniyrsID(idx,:),1);
    errors = zkyrs(idx) - zhyrs(idx);
    RMSE_foldyrs(i) = sqrt( mean(errors.^2) );
    MAE_foldyrs(i) = mean(abs(errors));
    ME_foldyrs(i) = mean(errors);
    r2_foldyrs(i) = (corr(zkyrs(idx),zhyrs(idx),'type','Pearson')).^2;
    MS_foldyrs(i) = mean(errors./sqrt(vkyrs(idx)));
    RMSS_foldyrs(i) = std(errors./sqrt(vkyrs(idx)));
    MR_foldyrs(i) = mean(sqrt(vkyrs(idx)));
    std_est_foldyrs(i) = std(zkyrs(idx));
    std_obs_foldyrs(i) = std(zhyrs(idx));
    QAr2_foldyrs(i) = ( ( std_est_foldyrs(i).^2 + std_obs_foldyrs(i).^2 - ...
        (RMSE_foldyrs(i).^2-ME_foldyrs(i).^2) )./(2*std_est_foldyrs(i)*std_obs_foldyrs(i)) )^2;
    
end

% by year yearly
uniyrsyr = unique(uniyrsID(:,3));
for i = 1:length(uniyrsyr)
    
    idx = uniyrsyr(i) == uniyrsID(:,3);
    
    n_sites_yearyrs(i) = size(unique(uniyrsID(idx,1:2),'rows'),1);
    n_observations_yearyrs(i) = size(uniyrsID(idx,:),1);
    errors = zkyrs(idx) - zhyrs(idx);
    RMSE_yearyrs(i) = sqrt( mean(errors.^2) );
    MAE_yearyrs(i) = mean(abs(errors));
    ME_yearyrs(i) = mean(errors);
    r2_yearyrs(i) = (corr(zkyrs(idx),zhyrs(idx),'type','Pearson')).^2;
    MS_yearyrs(i) = mean(errors./sqrt(vkyrs(idx)));
    RMSS_yearyrs(i) = std(errors./sqrt(vkyrs(idx)));
    MR_yearyrs(i) = mean(sqrt(vkyrs(idx)));
    std_est_yearyrs(i) = std(zkyrs(idx));
    std_obs_yearyrs(i) = std(zhyrs(idx));
    QAr2_yearyrs(i) = ( ( std_est_yearyrs(i).^2 + std_obs_yearyrs(i).^2 - ...
        (RMSE_yearyrs(i).^2-ME_yearyrs(i).^2) )./(2*std_est_yearyrs(i)*std_obs_yearyrs(i)) )^2;
    
end

% by EPA region yearly
uniyrsEPA = unique(EPAyrs);
for i = 1:length(uniyrsEPA)
    
    idx = uniyrsEPA(i) == EPAyrs;
    
    n_sites_EPAyrs(i) = size(unique(uniyrsID(idx,1:2),'rows'),1);
    n_observations_EPAyrs(i) = size(uniyrsID(idx,:),1);
    errors = zkyrs(idx) - zhyrs(idx);
    RMSE_EPAyrs(i) = sqrt( mean(errors.^2) );
    MAE_EPAyrs(i) = mean(abs(errors));
    ME_EPAyrs(i) = mean(errors);
    r2_EPAyrs(i) = (corr(zkyrs(idx),zhyrs(idx),'type','Pearson')).^2;
    MS_EPAyrs(i) = mean(errors./sqrt(vkyrs(idx)));
    RMSS_EPAyrs(i) = std(errors./sqrt(vkyrs(idx)));
    MR_EPAyrs(i) = mean(sqrt(vkyrs(idx)));
    std_est_EPAyrs(i) = std(zkyrs(idx));
    std_obs_EPAyrs(i) = std(zhyrs(idx));
    QAr2_EPAyrs(i) = ( ( std_est_EPAyrs(i).^2 + std_obs_EPAyrs(i).^2 - ...
        (RMSE_EPAyrs(i).^2-ME_EPAyrs(i).^2) )./(2*std_est_EPAyrs(i)*std_obs_EPAyrs(i)) )^2;
    
end

% save results
save(sprintf('Xval_true10fold_results_hardonly.mat'),...
    'n_sites','n_observations','RMSE','MAE','ME','r2','MS','RMSS','MR','std_est',...
    'std_obs','QAr2','n_sites_fold','n_observations_fold','RMSE_fold',...
    'MAE_fold','ME_fold','r2_fold','MS_fold','RMSS_fold','MR_fold',...
    'std_est_fold','std_obs_fold','QAr2_fold','n_sites_year','n_observations_year',...
    'RMSE_year','MAE_year','ME_year','r2_year','MS_year','RMSS_year','MR_year',...
    'std_est_year','std_obs_year','QAr2_year','n_sites_EPA','n_observations_EPA',...
    'RMSE_EPA','MAE_EPA','ME_EPA','r2_EPA','MS_EPA','RMSS_EPA','MR_EPA',...
    'std_est_EPA','std_obs_EPA','QAr2_EPA','n_sitesyrs','n_observationsyrs',...
    'RMSEyrs','MAEyrs','MEyrs','r2yrs','MSyrs','RMSSyrs','MRyrs','std_estyrs','std_obsyrs','QAr2yrs',...
    'n_sites_foldyrs','n_observations_foldyrs','RMSE_foldyrs','MAE_foldyrs',...
    'ME_foldyrs','r2_foldyrs','MS_foldyrs','RMSS_foldyrs','MR_foldyrs','std_est_foldyrs','std_obs_foldyrs',...
    'QAr2_foldyrs','n_sites_yearyrs','n_observations_yearyrs','RMSE_yearyrs',...
    'MAE_yearyrs','ME_yearyrs','r2_yearyrs','MS_yearyrs','RMSS_yearyrs','MR_yearyrs','std_est_yearyrs','std_obs_yearyrs',...
    'QAr2_yearyrs','n_sites_EPAyrs','n_observations_EPAyrs','RMSE_EPAyrs','MAE_EPAyrs',...
    'ME_EPAyrs','r2_EPAyrs','MS_EPAyrs','RMSS_EPAyrs','MR_EPAyrs','std_est_EPAyrs','std_obs_EPAyrs','QAr2_EPAyrs',...
    'dists_info','dists_info_fold','dists_info_year','dists_info_EPA',...
    'n_sites_mos','n_observations_mos','RMSE_mos','MAE_mos','ME_mos','r2_mos','MS_mos',...
    'RMSS_mos','MR_mos','std_est_mos','std_obs_mos','QAr2_mos','dists_info_mos',...
    'n_sites_sea','n_observations_sea','RMSE_sea','MAE_sea','ME_sea','r2_sea','MS_sea',...
    'RMSS_sea','MR_sea','std_est_sea','std_obs_sea','QAr2_sea','dists_info_sea',...
    'n_sites_state','n_observations_state','RMSE_state','MAE_state','ME_state','r2_state','MS_state',...
    'RMSS_state','MR_state','std_est_state','std_obs_state','QAr2_state','dists_info_state',...
    'n_sites_sID','n_observations_sID','RMSE_sID','MAE_sID','ME_sID','r2_sID',...
    'MS_sID','RMSS_sID','MR_sID','std_est_sID','std_obs_sID','QAr2_sID',...
    'ckall','zhall','zkall','vkall','foldall','EPAall',...
    'ckyrs','zhyrs','zkyrs','vkyrs','foldyrs','EPAyrs'); 

%%% display results %%%

% load data
load(sprintf('Xval_true10fold_results_hardonly.mat'));

% header string
strval = 'n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2';

% table for all data daily
n = [n_sites n_observations RMSE MAE ME r2 MS RMSS MR std_est std_obs QAr2];
outid = fopen(sprintf('Xval_true10fold_alldaily_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_alldaily_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by fold daily
n = [n_sites_fold' n_observations_fold' RMSE_fold' MAE_fold' ME_fold' r2_fold' ...
    MS_fold' RMSS_fold' MR_fold' std_est_fold' std_obs_fold' QAr2_fold'];
outid = fopen(sprintf('Xval_true10fold_folddaily_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_folddaily_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by month daily
n = [n_sites_mos' n_observations_mos' RMSE_mos' MAE_mos' ME_mos' r2_mos' ...
    MS_mos' RMSS_mos' MR_mos' std_est_mos' std_obs_mos' QAr2_mos'];
outid = fopen(sprintf('Xval_true10fold_monthdaily_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_monthdaily_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by season daily
n = [n_sites_sea' n_observations_sea' RMSE_sea' MAE_sea' ME_sea' r2_sea' ...
    MS_sea' RMSS_sea' MR_sea' std_est_sea' std_obs_sea' QAr2_sea'];
outid = fopen(sprintf('Xval_true10fold_seasondaily_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_seasondaily_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by year daily
n = [n_sites_year' n_observations_year' RMSE_year' MAE_year' ME_year' r2_year' ...
    MS_year' RMSS_year' MR_year' std_est_year' std_obs_year' QAr2_year'];
outid = fopen('Xval_true10fold_yeardaily_hardonly.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_yeardaily_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by EPA region daily
n = [n_sites_EPA' n_observations_EPA' RMSE_EPA' MAE_EPA' ME_EPA' r2_EPA' MS_EPA' ...
    RMSS_EPA' MR_EPA' std_est_EPA' std_obs_EPA' QAr2_EPA'];
outid = fopen('Xval_true10fold_regiondaily_hardonly.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_regiondaily_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table for all data yearly
n = [n_sitesyrs n_observationsyrs RMSEyrs MAEyrs MEyrs r2yrs MSyrs RMSSyrs MRyrs std_estyrs std_obsyrs ...
    QAr2yrs];
outid = fopen(sprintf('Xval_true10fold_allyearly_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_allyearly_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by fold yearly
n = [n_sites_foldyrs' n_observations_foldyrs' RMSE_foldyrs' MAE_foldyrs' ME_foldyrs' r2_foldyrs' ...
    MS_foldyrs' RMSS_foldyrs' MR_foldyrs' std_est_foldyrs' std_obs_foldyrs' QAr2_foldyrs'];
outid = fopen(sprintf('Xval_true10fold_foldyearly_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_foldyearly_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by year yearly
n = [n_sites_yearyrs' n_observations_yearyrs' RMSE_yearyrs' MAE_yearyrs' ...
    ME_yearyrs' r2_yearyrs' MS_yearyrs' RMSS_yearyrs' MR_yearyrs' std_est_yearyrs' std_obs_yearyrs' QAr2_yearyrs'];
outid = fopen(sprintf('Xval_true10fold_yearyearly_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_yearyearly_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by EPA region yearly
n = [n_sites_EPAyrs' n_observations_EPAyrs' RMSE_EPAyrs' MAE_EPAyrs' ME_EPAyrs' ...
    r2_EPAyrs' MS_EPAyrs' RMSS_EPAyrs' MR_EPAyrs' std_est_EPAyrs' std_obs_EPAyrs' QAr2_EPAyrs'];
outid = fopen(sprintf('Xval_true10fold_regionyearly_hardonly.csv'),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xval_true10fold_regionyearly_hardonly.csv'), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

end