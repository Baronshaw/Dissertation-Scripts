function [] = exploreHighError()
% this function will investigate the high errors with the soft data and see
% if it is/isn't coming from the soft data itself
 
constant = 0; 
gauss = 1; 

if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% gathering all the data
for i = 1:10   
    
    % loading results
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',i,'_nosoft',constr,gaussstr));
    zkallh{i,1} = zk_madd;
    zhallh{i,1} = zh_Xval;
    ckallh{i,1} = ck;
    vkallh{i,1} = vk;  
    foldallh{i,1} = i*ones(length(zkallh{i,1}),1);
    
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',i,'_soft',constr,gaussstr));
    zkalls{i,1} = zk_madd;
    zhalls{i,1} = zh_Xval;
    ckalls{i,1} = ck;
    vkalls{i,1} = vk;  
    foldalls{i,1} = i*ones(length(zkalls{i,1}),1);
    
end

zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx);
zhallh = cell2mat(zhallh); zhallh = zhallh(idx);
ckallh = cell2mat(ckallh); ckallh = ckallh(idx,:);
vkallh = cell2mat(vkallh); vkallh = vkallh(idx);
foldallh = cell2mat(foldallh); foldallh = foldallh(idx);

zkalls = zkalls(idx);
zhalls = cell2mat(zhalls); zhalls = zhalls(idx);
ckalls = cell2mat(ckalls); ckalls = ckalls(idx,:);
vkalls = cell2mat(vkalls); vkalls = vkalls(idx);
foldalls = cell2mat(foldalls); foldalls = foldalls(idx);

% looping through each soft data year
temp = datevec(ckalls(:,3));
yrs = temp(:,1);
uniyr = unique(yrs);
soft_years = [2001:2002 2005:2007];
allAEs = abs(zkalls - zhalls);

% look at a histogram of errors
for i = 1:length(soft_years)
    idx = soft_years(i) == yrs;
    figure; hold on;
    hist(allAEs(idx),100);
    higherr(i) = prctile(allAEs(idx),99); % 99%tile of errors
    title(sprintf('abs errors for %d',soft_years(i)));
end

% fill in high errors with the hard data results 
zkall_update = zkalls;
vkall_update = vkalls;
for i = 1:length(soft_years)
    idx = soft_years(i) == yrs & allAEs >= higherr(i);
    zkall_update(idx) = zkallh(idx);
    vkall_update(idx) = vkallh(idx);
end
blah = 5;

% recalculate all cross-validation statistics by year and compare the two
% by year daily
temp = datevec(ckalls(:,3));
yrs = temp(:,1);
uniyr = unique(yrs);
for i = 1:length(uniyr)
    
    idx = uniyr(i) == yrs;
    
    n_sites_year_update(i) = size(unique(ckalls(idx,1:2),'rows'),1);
    n_observations_year_update(i) = size(ckalls(idx,:),1);
    errors = zkall_update(idx) - zhalls(idx);
    RMSE_year_update(i) = sqrt( mean(errors.^2) );
    MAE_year_update(i) = mean(abs(errors));
    ME_year_update(i) = mean(errors);
    r2_year_update(i) = (corr(zkall_update(idx),zhalls(idx),'type','Pearson')).^2;
    MS_year_update(i) = mean(errors./sqrt(vkall_update(idx)));
    RMSS_year_update(i) = std(errors./sqrt(vkall_update(idx)));
    MR_year_update(i) = mean(sqrt(vkall_update(idx)));
    std_est_year_update(i) = std(zkall_update(idx));
    std_obs_year_update(i) = std(zhalls(idx));
    QAr2_year_update(i) = ( ( std_est_year_update(i).^2 + std_obs_year_update(i).^2 - ...
        (RMSE_year_update(i).^2-ME_year_update(i).^2) )./(2*std_est_year_update(i)*std_obs_year_update(i)) )^2;    
    
end

% hopefully there should be very little difference, supporting my inutition
% that it is not a few large errors that are contributing to poor
% statistics

% load previous results
load(sprintf('Xval_10fold%s%s%s_results.mat','_nosoft',constr,gaussstr));
n_sites_yearh = n_sites_year; n_observations_yearh = n_observations_year; 
RMSE_yearh = RMSE_year; MAE_yearh = MAE_year; ME_yearh = ME_year;
r2_yearh = r2_year; MS_yearh = MS_year; RMSS_yearh = RMSS_year;
MR_yearh = MR_year; std_est_yearh = std_est_year; std_obs_yearh = std_obs_year;
QAr2_yearh = QAr2_year;
load(sprintf('Xval_10fold%s%s%s_results.mat','_soft',constr,gaussstr));
n_sites_years = n_sites_year; n_observations_years = n_observations_year; 
RMSE_years = RMSE_year; MAE_years = MAE_year; ME_years = ME_year;
r2_years = r2_year; MS_years = MS_year; RMSS_years = RMSS_year;
MR_years = MR_year; std_est_years = std_est_year; std_obs_years = std_obs_year;
QAr2_years = QAr2_year;

% header string
strval = 'RMSE,,,MAE,,,ME,,,r2,,,MS,,,RMSS,,,MR,,,std est,,';

% put all the results in a table
n = [RMSE_yearh' RMSE_years' RMSE_year_update' ...
    MAE_yearh' MAE_years' MAE_year_update' ...
    ME_yearh' ME_years' ME_year_update' ...
    r2_yearh' r2_years' r2_year_update' ...
    MS_yearh' MS_years' MS_year_update' ...
    RMSS_yearh' RMSS_years' RMSS_year_update' ...
    MR_yearh' MR_years' MR_year_update' ...
    std_est_yearh' std_est_years' std_est_year_update'];
outid = fopen('Xval_10fold_yeardaily_update.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('Xval_10fold_yeardaily_update.csv', ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

end