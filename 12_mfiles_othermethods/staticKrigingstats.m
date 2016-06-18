function [] = staticKrigingstats()
% this function will take the results of the RAMPXval and calculcate
% statistics and create a table

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckall_Xval = cell(size(cMS,1),1);
for i = 1:size(cMS,1)
    disp(i);
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckall_Xval{i,1} = ch(idx&temp2001==2001|idx&temp2001==2002,:);
end

% load Xval data
softstr = '_soft'; 
constr = '_long'; 
gaussstr = '_gauss'; 
forceddist = 0;
load(sprintf('matdata/Xvalforcediso_LOOCV_kriging_foriso%dkm_timezone.mat',floor(forceddist./1000)));

% set up data: location day actual modeled estimated
zk_Xval = zk_madd;
vk_Xval = vk;     
zk_Xval = cell2mat(zk_Xval);
idx = ~isnan(zk_Xval); zk_Xval = zk_Xval(idx);
zh_Xval = cell2mat(zh_Xval); zh_Xval = zh_Xval(idx);
ckall_Xval = cell2mat(ckall_Xval); ckall_Xval = ckall_Xval(idx,:);
vk_Xval = cell2mat(vk_Xval); vk_Xval = vk_Xval(idx);
unidates = unique(ckall_Xval(:,3));

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

for i = 1:length(unidates)
    disp(i);
    idx = ckall_Xval(:,3) == unidates(i);
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
save('matdata/LOOCE_kriging_table.mat','ckall_Xval','zh_Xval','zk_Xval','vk_Xval', ...
    'n_sites','n_observations', 'RMSE','MAE','ME','r2','MS','RMSS','MR','std_est','std_obs','QAr2', ... 
    'unidates','n_sites_d','n_observations_d','RMSE_d','MAE_d','ME_d','r2_d', ...
    'MS_d','RMSS_d','MR_d','std_est_d','std_obs_d','QAr2_d');

% put results in a table form
strval = 'n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2'; % header string
n = [n_sites n_observations RMSE MAE ME r2 MS RMSS MR std_est std_obs QAr2];
outid = fopen('tables/Xval_kriging_overall.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/Xval_kriging_overall.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

strval = 'year,month,day,n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2'; % header string
temp = datevec(unidates);
n = [temp(:,1) temp(:,2) temp(:,3) n_sites_d n_observations_d RMSE_d MAE_d ME_d r2_d MS_d RMSS_d MR_d std_est_d std_obs_d QAr2_d];
outid = fopen('tables/Xval_kriging_daily.csv','w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite('tables/Xval_kriging_daily.csv',n,'delimiter',',','precision',6,'-append','roffset',1)

end