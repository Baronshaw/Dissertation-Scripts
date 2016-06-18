function [] = displayXvalforcedisolationStatistics_10fold_test1(soft,constant,gauss)
% this function will display the results of the 10 fold cross validation
% statistics

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% load data
load(sprintf('Xvalforcediso_10fold%s%s%s_results_test1.mat',softstr,constr,gaussstr));

% header string
strval = 'n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2';

% table for all data daily
n = [n_sites n_observations RMSE MAE ME r2 MS RMSS MR std_est std_obs QAr2];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_alldaily_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_alldaily_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by fold daily
n = [n_sites_fold' n_observations_fold' RMSE_fold' MAE_fold' ME_fold' r2_fold' ...
    MS_fold' RMSS_fold' MR_fold' std_est_fold' std_obs_fold' QAr2_fold'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_folddaily_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_folddaily_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by month daily
n = [n_sites_mos' n_observations_mos' RMSE_mos' MAE_mos' ME_mos' r2_mos' ...
    MS_mos' RMSS_mos' MR_mos' std_est_mos' std_obs_mos' QAr2_mos'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_monthdaily_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_monthdaily_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by season daily
n = [n_sites_sea' n_observations_sea' RMSE_sea' MAE_sea' ME_sea' r2_sea' ...
    MS_sea' RMSS_sea' MR_sea' std_est_sea' std_obs_sea' QAr2_sea'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_seasondaily_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_seasondaily_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by year daily
n = [n_sites_year' n_observations_year' RMSE_year' MAE_year' ME_year' r2_year' ...
    MS_year' RMSS_year' MR_year' std_est_year' std_obs_year' QAr2_year'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_yeardaily_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_yeardaily_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by EPA region daily
n = [n_sites_EPA' n_observations_EPA' RMSE_EPA' MAE_EPA' ME_EPA' r2_EPA' MS_EPA' ...
    RMSS_EPA' MR_EPA' std_est_EPA' std_obs_EPA' QAr2_EPA'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_regiondaily_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_regiondaily_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table for all data yearly
n = [n_sitesyrs n_observationsyrs RMSEyrs MAEyrs MEyrs r2yrs MSyrs RMSSyrs MRyrs std_estyrs std_obsyrs ...
    QAr2yrs];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_allyearly_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_allyearly_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by fold yearly
n = [n_sites_foldyrs' n_observations_foldyrs' RMSE_foldyrs' MAE_foldyrs' ME_foldyrs' r2_foldyrs' ...
    MS_foldyrs' RMSS_foldyrs' MR_foldyrs' std_est_foldyrs' std_obs_foldyrs' QAr2_foldyrs'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_foldyearly_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_foldyearly_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by year yearly
n = [n_sites_yearyrs' n_observations_yearyrs' RMSE_yearyrs' MAE_yearyrs' ...
    ME_yearyrs' r2_yearyrs' MS_yearyrs' RMSS_yearyrs' MR_yearyrs' std_est_yearyrs' std_obs_yearyrs' QAr2_yearyrs'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_yearyearly_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_yearyearly_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by EPA region yearly
n = [n_sites_EPAyrs' n_observations_EPAyrs' RMSE_EPAyrs' MAE_EPAyrs' ME_EPAyrs' ...
    r2_EPAyrs' MS_EPAyrs' RMSS_EPAyrs' MR_EPAyrs' std_est_EPAyrs' std_obs_EPAyrs' QAr2_EPAyrs'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_regionyearly_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_regionyearly_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

%%% adding median distance between stations

% header string
strval = 'n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2,med dist';

% table for all data daily
n = [n_sites n_observations RMSE MAE ME r2 MS RMSS MR std_est std_obs QAr2 dists_info];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_alldailydists_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_alldailydists_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by fold daily
n = [n_sites_fold' n_observations_fold' RMSE_fold' MAE_fold' ME_fold' r2_fold' ...
    MS_fold' RMSS_fold' MR_fold' std_est_fold' std_obs_fold' QAr2_fold' dists_info_fold'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_folddailydists_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_folddailydists_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by year daily
n = [n_sites_year' n_observations_year' RMSE_year' MAE_year' ME_year' r2_year' ...
    MS_year' RMSS_year' MR_year' std_est_year' std_obs_year' QAr2_year' dists_info_year'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_yeardailydists_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_yeardailydists_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% table by EPA region daily
n = [n_sites_EPA' n_observations_EPA' RMSE_EPA' MAE_EPA' ME_EPA' r2_EPA' MS_EPA' ...
    RMSS_EPA' MR_EPA' std_est_EPA' std_obs_EPA' QAr2_EPA' dists_info_EPA'];
outid = fopen(sprintf('Xvalforcediso_10fold%s%s%s_regiondailydists_test1.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_10fold%s%s%s_regiondailydists_test1.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

end