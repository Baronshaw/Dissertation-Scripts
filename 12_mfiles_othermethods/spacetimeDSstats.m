function [] = spacetimeDSstats(additive,multiplicative)
% this function will take the results from the Xval results from the
% space/time downscaler and calculate statistics and put them into table 
% form

if nargin < 1, additive = 'ind'; end 
if nargin < 2, multiplicative = 'ind'; end

% load data
for i = 2001:2002
    load(sprintf('../matfiles/prepCTMandObs_%d.mat',i));
    Modall{i-2000,1} = Mod; Obsall{i-2000,1} = Obs;
    coordObsall{i-2000,1} = coordObs; cht{i-2000,1} = yrmodaObs;
    load(sprintf('../matfiles/prepCTM_%d.mat',i));
    distCTMvall{i-2000,1} = distCTMv; dailyCTMvall{i-2000,1} = dailyCTMv;
    yrmodaCTMvall{i-2000,1} = yrmodaCTMv;
end
zhm_paired = cell2mat(Modall); zho_paired = cell2mat(Obsall);
ch_paired = cell2mat(coordObsall); cht_paired = cell2mat(cht);
chm_all = cell2mat(distCTMvall); zhm_all = cell2mat(dailyCTMvall);
chtm_all = cell2mat(yrmodaCTMvall);
yrall = floor(cht_paired./10000);
moall = floor((cht_paired - yrall*10000)/100);
daall = cht_paired  - yrall*10000 - moall*100;
cht_paired = datenum(yrall,moall,daall);
chtm_all = datenum(chtm_all(:,1),chtm_all(:,2),chtm_all(:,3));
unidates = unique(cht_paired);

% load Xval data
load(sprintf('matdata/LOOCE_DS_add_%s_muli_%s.mat',additive,multiplicative));

% set up data: location day actual modeled estimated
ckall_Xval = cell2mat(ck);
zhX = cell2mat(zhX); zk = cell2mat(zk); vk = cell2mat(vk);
cktall_Xval = ckall_Xval(:,3);
ckall_Xval = ckall_Xval(:,1:2);
[lia lib] = ismember([ch_paired cht_paired],[ckall_Xval cktall_Xval],'rows');
zh_Xval = zhX(lib);
zk_Xval = zk(lib);
vk_Xval = vk(lib);

% calculate statistics 1) by day and 2) overall
n_sites = size(unique(ckall_Xval(:,1:2),'rows'),1);
n_observations = size(ckall_Xval,1);
errors = zk_Xval - zh_Xval;
RMSE = sqrt( mean(errors.^2) );
MAE = mean(abs(errors));
ME = mean(errors);
r2 = (corr(zk_Xval,zh_Xval,'type','Pearson')).^2;
MS = mean(errors./sqrt(vk_Xval));
RMSS = std(errors./sqrt(vk_Xval));
MR = mean(sqrt(vk_Xval));
std_est = std(zk_Xval);
std_obs = std(zh_Xval);
QAr2 = ( ( std_est.^2 + std_obs.^2 - (RMSE.^2-ME.^2) )./(2*std_est*std_obs) )^2;

% you can calculate the statistics for each day, but this makes more sense
% for the static DS only
for i = 1:length(unidates)
    idx = cktall_Xval == unidates(i);
    n_sites_d(i,1) = size(unique(ckall_Xval(idx,1:2),'rows'),1);
    n_observations_d(i,1) = size(ckall_Xval(idx,:),1);
    errors = zk_Xval(idx) - zh_Xval(idx);
    RMSE_d(i,1) = sqrt( mean(errors.^2) );
    MAE_d(i,1) = mean(abs(errors));
    ME_d(i,1) = mean(errors);
    r2_d(i,1) = (corr(zk_Xval(idx),zh_Xval(idx),'type','Pearson')).^2;
    MS_d(i,1) = mean(errors./sqrt(vk_Xval(idx)));
    RMSS_d(i,1) = std(errors./sqrt(vk_Xval(idx)));
    MR_d(i,1) = mean(sqrt(vk_Xval(idx)));
    std_est_d(i,1) = std(zk_Xval(idx));
    std_obs_d(i,1) = std(zh_Xval(idx));
    QAr2_d(i,1) = ( ( std_est_d(i,1).^2 + std_obs_d(i,1).^2 - (RMSE_d(i,1).^2-ME_d(i,1).^2) )./(2*std_est_d(i,1)*std_obs_d(i,1)) )^2;
end

% save Xval results
save(sprintf('matdata/LOOCE_DS_table_add_%s_muli_%s.mat',additive,multiplicative), ...
    'ckall_Xval','cktall_Xval','zh_Xval','zk_Xval','vk_Xval', ...
    'n_sites','n_observations', 'RMSE','MAE','ME','r2','MS','RMSS','MR','std_est','std_obs','QAr2', ... 
    'unidates','n_sites_d','n_observations_d','RMSE_d','MAE_d','ME_d','r2_d', ...
    'MS_d','RMSS_d','MR_d','std_est_d','std_obs_d','QAr2_d');

% put results in a table form
strval = 'n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2'; % header string
n = [n_sites n_observations RMSE MAE ME r2 MS RMSS MR std_est std_obs QAr2];
outid = fopen(sprintf('tables/Xval_DS_overall_add_%s_muli_%s.csv',additive,multiplicative),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/Xval_DS_overall_add_%s_muli_%s.csv',additive,multiplicative),n,'delimiter',',','precision',6,'-append','roffset',1)

strval = 'year,month,day,n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2'; % header string
temp = datevec(unidates);
n = [temp(:,1) temp(:,2) temp(:,3) n_sites_d n_observations_d RMSE_d MAE_d ME_d r2_d MS_d RMSS_d MR_d std_est_d std_obs_d QAr2_d];
outid = fopen(sprintf('tables/Xval_DS_daily_add_%s_muli_%s.csv',additive,multiplicative),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('tables/Xval_DS_daily_add_%s_muli_%s.csv',additive,multiplicative),n,'delimiter',',','precision',6,'-append','roffset',1)

end