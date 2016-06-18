function [] = getXvalStatistics_LOOCV_test1(soft,constant,gauss)
% this function will calculate the cross validation statistics for the 
% LOOCV cross validation

if nargin < 1, soft = 0; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% gathering ck
if ~exist('../matfiles/Xval_LOOCV_ckall.mat')
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
    save('../matfiles/Xval_LOOCV_ckall.mat','ckall');
else
    load('../matfiles/Xval_LOOCV_ckall.mat');
end

% formatting ck
len = size(ckall{1},1);
ckalltemp = cell(size(ckall,1)*len,1);
n = 1;
for i = 1:12
    for j = 1:len
        ckalltemp{n,1} = ckall{i}{j};
        n = n + 1;
    end
end
ckall = cell2mat(ckalltemp);

% gathering all the data
for i = 1:12  
    load(sprintf('../matfiles/Xval_LOOCV_mon%d_%s%s%s_test1.mat',i,softstr,constr,gaussstr));
    zkall{i,1} = cell2mat(zk_madd);
    zhall{i,1} = cell2mat(zh_Xval);
    vkall{i,1} = cell2mat(vk);
end
zkall = cell2mat(zkall);
idx = ~isnan(zkall); zkall = zkall(idx);
zhall = cell2mat(zhall); zhall = zhall(idx);
ckall = ckall(idx,:);
vkall = cell2mat(vkall); vkall = vkall(idx);

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
        
    % by year/station daily
    uniIDyr = unique(ckall(idx,1:2),'rows');
    zksub = zkall(idx); zhsub = zhall(idx);
    for j = 1:n_sites_year(i)
        idxIDyr = uniIDyr(j,1) == ckall(idx,1) & uniIDyr(j,2) == ckall(idx,2);
        errors = zksub(idxIDyr) - zhsub(idxIDyr);
        RMSE_IDyr{i}(j) = mean(errors.^2); 
    end
    
end

% by EPA region daily
% go through original station locations
cd Observed_PM2p5
for i = 1:length(uniyr)
    load(sprintf('MasterDaily_PM2p5_%d.mat',uniyr(i)));
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
cd ..

% convert them to coordinate system
load projexample.mat;
cd 09_mfiles_projections
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load Ellipsoids.mat
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
projuniID = ell2lambertcc([unilon unilat],'whiproj2001'); 
cd ..
% match estimations locations with station id
sIDall = NaN*ones(size(ckall,1),1);
for i = 1:length(projuniID)
    idx = projuniID(i,1) == ckall(:,1) & projuniID(i,2) == ckall(:,2);
    sIDall(idx) = uniID(i);
end
stateall = floor(sIDall./10^7);
% match id's with EPA region
% from https://aqs.epa.gov/aqsweb/codes/data/StateCountyCodes.csv
M = csvread('StateCountyCodes_mod.csv');
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
save(sprintf('Xval_LOOCV%s%s%s_results_test1.mat',softstr,constr,gaussstr),...
    'n_sites','n_observations','RMSE','MAE','ME','r2','MS','RMSS','MR','std_est',...
    'std_obs','QAr2','n_sites_year','n_observations_year',...
    'RMSE_year','MAE_year','ME_year','r2_year','MS_year','RMSS_year','MR_year',...
    'std_est_year','std_obs_year','QAr2_year','n_sites_EPA','n_observations_EPA',...
    'RMSE_EPA','MAE_EPA','ME_EPA','r2_EPA','MS_EPA','RMSS_EPA','MR_EPA',...
    'std_est_EPA','std_obs_EPA','QAr2_EPA','n_sitesyrs','n_observationsyrs',...
    'RMSEyrs','MAEyrs','MEyrs','r2yrs','MSyrs','RMSSyrs','MRyrs','std_estyrs','std_obsyrs','QAr2yrs',...
    'n_sites_yearyrs','n_observations_yearyrs','RMSE_yearyrs',...
    'MAE_yearyrs','ME_yearyrs','r2_yearyrs','MS_yearyrs','RMSS_yearyrs','MR_yearyrs','std_est_yearyrs','std_obs_yearyrs',...
    'QAr2_yearyrs','n_sites_EPAyrs','n_observations_EPAyrs','RMSE_EPAyrs','MAE_EPAyrs',...
    'ME_EPAyrs','r2_EPAyrs','MS_EPAyrs','RMSS_EPAyrs','MR_EPAyrs','std_est_EPAyrs','std_obs_EPAyrs','QAr2_EPAyrs',...
    'dists_info','dists_info_year','dists_info_EPA',...
    'n_sites_mos','n_observations_mos','RMSE_mos','MAE_mos','ME_mos','r2_mos','MS_mos',...
    'RMSS_mos','MR_mos','std_est_mos','std_obs_mos','QAr2_mos','dists_info_mos',...
    'n_sites_sea','n_observations_sea','RMSE_sea','MAE_sea','ME_sea','r2_sea','MS_sea',...
    'RMSS_sea','MR_sea','std_est_sea','std_obs_sea','QAr2_sea','dists_info_sea',...
    'n_sites_state','n_observations_state','RMSE_state','MAE_state','ME_state','r2_state','MS_state',...
    'RMSS_state','MR_state','std_est_state','std_obs_state','QAr2_state','dists_info_state',...
    'n_sites_sID','n_observations_sID','RMSE_sID','MAE_sID','ME_sID','r2_sID',...
    'MS_sID','RMSS_sID','MR_sID','std_est_sID','std_obs_sID','QAr2_sID',...
    'RMSE_IDyr',...
    'ckall','zhall','zkall','vkall','EPAall',...
    'ckyrs','zhyrs','zkyrs','vkyrs','EPAyrs'); 

end