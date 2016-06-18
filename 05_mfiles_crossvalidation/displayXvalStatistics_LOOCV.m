function [] = displayXvalStatistics_LOOCV(soft,constant,gauss)
% this function will display the results of the 10 fold cross validation
% statistics

if nargin < 1, soft = 0; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

forceddist = 0:100000:1000000;
yearz = 1999:2010;
forceddist = 0;
yearz = 2001;

% load data
for i = 1:length(forceddist)
    
    load(sprintf('Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_results_timezone.mat', ...
        softstr,constr,gaussstr,floor(forceddist(i)./1000)));
    n_sitesz(i,1) = n_sites;
    n_observationsz(i,1) = n_observations;
    RMSEz(i,1) = RMSE;
    MAEz(i,1) = MAE;
    MEz(i,1) = ME;
    r2z(i,1) = r2;
    MSz(i,1) = MS;
    RMSSz(i,1) = RMSS;
    MRz(i,1) = MR;
    std_estz(i,1) = std_est;
    std_obsz(i,1) = std_obs;
    QAr2z(i,1) = QAr2;
    
    n_sites_yearz{i,1} = n_sites_year';
    n_observations_yearz{i,1} = n_observations_year';
    RMSE_yearz{i,1} = RMSE_year';
    MAE_yearz{i,1} = MAE_year';
    ME_yearz{i,1} = ME_year';
    r2_yearz{i,1} = r2_year';
    MS_yearz{i,1} = MS_year';
    RMSS_yearz{i,1} = RMSS_year';
    MR_yearz{i,1} = MR_year';
    std_est_yearz{i,1} = std_est_year';
    std_obs_yearz{i,1} = std_obs_year';
    QAr2_yearz{i,1} = QAr2_year';
    
end

% header string
strval = 'dist,n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2';

% table for all data daily
n = [floor(forceddist/1000)' n_sitesz n_observationsz RMSEz MAEz MEz r2z MSz RMSSz MRz std_estz std_obsz QAr2z];
outid = fopen(sprintf('Xvalforcediso_LOOCV%s%s%s_foriso_alldaily.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_LOOCV%s%s%s_foriso_alldaily.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

% header string
strval = 'year,dist,n sites,n obs,RMSE,MAE,ME,r2,MS,RMSS,MR,std est,std obs,QA r2';

% table by year daily
temp = repmat(floor(forceddist/1000),length(n_sites_year),1);
tempyr = repmat(yearz',length(forceddist),1);
n = [tempyr temp(:) cell2mat(n_sites_yearz) cell2mat(n_observations_yearz) ...
    cell2mat(RMSE_yearz) cell2mat(MAE_yearz) cell2mat(ME_yearz) cell2mat(r2_yearz) ...
    cell2mat(MS_yearz) cell2mat(RMSS_yearz) cell2mat(MR_yearz) cell2mat(std_est_yearz) ...
    cell2mat(std_obs_yearz) cell2mat(QAr2_yearz)];
outid = fopen(sprintf('Xvalforcediso_LOOCV%s%s%s_foriso_yeardaily.csv',softstr,constr,gaussstr),'w+');
fprintf(outid,'%s',strval);
fclose(outid);
dlmwrite(sprintf('Xvalforcediso_LOOCV%s%s%s_foriso_yeardaily.csv',softstr,constr,gaussstr), ...
    n,'delimiter',',','precision',6,'-append','roffset',1)

end