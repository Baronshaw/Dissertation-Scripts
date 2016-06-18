function [] = getXvalStatistics_LOOCV(soft,constant,gauss,forceddist)
% this function will calculate the cross validation statistics for the 
% LOOCV cross validation

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, forceddist = 0; end % distance in m

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
temp2005 = datevec(pd(:,3)); temp2005 = temp2005(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckall = cell(size(cMS,1),1);
for i = 1:size(cMS,1)
    disp(i);
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckall{i,1} = ch(idx & temp2001==2001,:);
    ckall{i,1} = ch(idx & temp2005==2005,:);
end

% loading results
load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_L2RMSE.mat', ...
    softstr,constr,gaussstr,floor(forceddist./1000)));
zkall = zk_madd;
zhall = zh_Xval;
vkall = vk;     

zkall = cell2mat(zkall);
idx = ~isnan(zkall); zkall = zkall(idx);
zhall = cell2mat(zhall); zhall = zhall(idx);
ckall = cell2mat(ckall); ckall = ckall(idx,:);
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

allstats = [n_observations n_sites RMSE MAE ME r2 MS RMSS MR std_est std_obs QAr2 ];

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
    
end

% save results
save(sprintf('Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_results_timezone.mat', ...
    softstr,constr,gaussstr,floor(forceddist./1000)),...
    'n_sites','n_observations','RMSE','MAE','ME','r2','MS','RMSS','MR','std_est',...
    'std_obs','QAr2',...
    'n_sites_year','n_observations_year',...
    'RMSE_year','MAE_year','ME_year','r2_year','MS_year','RMSS_year','MR_year',...
    'std_est_year','std_obs_year','QAr2_year'); 

end