function [] = staticDS_estimate(est_what)
% this function will create estimates at given space/time locations for the
% static DS methods. This will create estimates as well as values for said
% parameters

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "staticDS_estimate" -logfile "runall_Cluster12.out"

if nargin < 1, est_what = 2; end % 1=est,2=beta0,3=beta1
if est_what == 1 
    est_what_str = 'est';
elseif est_what == 2
    est_what_str = 'beta0';
elseif est_what == 3;
    est_what_str = 'beta1';
end
    
% loading BME function
cd ../BMELIB2.0b
startup
cd ../12_mfiles_othermethods

% load paired modeled and observed data (years 2001 and 2002 for now)
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

% modifying dates
yrall = floor(cht_paired./10000);
moall = floor((cht_paired - yrall*10000)/100);
daall = cht_paired  - yrall*10000 - moall*100;
cht_paired = datenum(yrall,moall,daall);
chtm_all = datenum(chtm_all(:,1),chtm_all(:,2),chtm_all(:,3));
unidates = unique(cht_paired);
len = length(unidates);

% once a month
unidates =  datenum([repmat(2001,24,1) [1:24]' repmat(1,24,1)]);
len = length(unidates);

% load parameter files
load('matdata/DS_results_v1.mat');

% kriging parameters
cs = []; zs = []; vs = [];
nhmax = 70; nsmax = 0;
order = NaN;
options = BMEoptions;
options(3) = 150000;

%%% estimates
% looping through each day
zk_mtr = cell(len,1);
vk = cell(len,1);
zk = cell(len,1);
ckall = cell(len,1);
%matlabpool open 12
for i = 1:len
    disp([i len]);
    % getting data from specific day
    idx = cht_paired == unidates(i);
    chsub = ch_paired(idx,:); chtsub = cht_paired(idx); zhosub = zho_paired(idx); zhmsub = zhm_paired(idx);
    idx = chtm_all == unidates(i);
    chmallsub = chm_all(idx,:); chtmallsub = chtm_all(idx); zhmallsub = zhm_all(idx);

    % removing mean trend from the data
    zhosub_mtr = zhosub - (beta0c{i} + beta1c{i}.*zhmsub);  
    
    % intialize parameters
    zk_mtr{i} = NaN*ones(length(chsub),1); 
    vk{i} = NaN*ones(length(chsub),1); 
    zk{i} = NaN*ones(length(chsub),1);
    ckall{i} = chmallsub;  

    % covariance parameters
    ch = chsub;
    zh = zhosub_mtr;
    A11 = betac{i,1}(1); A12 = betac{i,1}(2); A22 = betac{i,1}(3); 
    phi0 = betac{i,1}(4); phi1 = betac{i,1}(5);  tau2 = tau2c{i,1}; 
    zhmp = prctile(zhmsub,0:10:100);
    problems = NaN*ones(1,10);
    for k = 1:10
        problems(k) = (zhmp(k+1)-zhmp(k))./2 + zhmp(k); 
    end
    xB1 = problems(1); xB2 = problems(2); 
    covmodel = {'nuggetC','exponentialC_DS','exponentialC_DS'};
    covparam = {[tau2] [(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1]};
    dmax = [2000000];

    % find modeled values of the observed data
    temp1 = arrayfun(@(x) find(abs(x-zhmp)==min(abs(x-zhmp))),zhmsub,'UniformOutput',false);
    temp2 = find( cell2mat(cellfun(@(x) length(x)>1,temp1,'UniformOutput',false)) == 1);
    for k = 1:length(temp2)
        temp1{temp2(k)} = temp1{temp2(k)}(1);
    end
    zhm = zhmp(cell2mat(temp1));
    xB1M = repmat(zhm,length(zhm),1); xB2M = repmat(zhm',1,length(zhm));
    covparamdd_DS = {[tau2] [(A11-A12.*xB1M).*(A11-A12.*xB2M)] [(A12-A22.*xB1M).*(A12-A22.*xB2M)]};

    % find modeled values of the estimation locations
    if ismember(ckall{i},ch,'rows') == 1 % estimation location in the data 
        [ia ib] = ismember(ckall{i},ch,'rows');
        ck_mod = zhmsub(ib);
    elseif isequal(ckall{i},chmallsub) % estimation locations are modeled locations
        ck_mod = zhmallsub;
    else % new estimation location
        cMS1 = ckall{i}(:,1:2); cMS2 = chmallsub;
        DMS = sqrt(bsxfun(@plus,dot(cMS1,cMS1,2),dot(cMS2,cMS2,2)')-2*(cMS1*cMS2'));
        ck_mod = zhmallsub(DMS==min(DMS),:);
    end
    temp1 = arrayfun(@(x) find(abs(x-zhmp)==min(abs(x-zhmp))),ck_mod,'UniformOutput',false);
    temp2 = find( cell2mat(cellfun(@(x) length(x)>1,temp1,'UniformOutput',false)) == 1);
    for k = 1:length(temp2)
        temp1{temp2(k)} = temp1{temp2(k)}(1);
    end
    zkm = zhmp(cell2mat(temp1));
    
    % defining covariance parameters
    if est_what == 1        
        covparamkd_DS = {[tau2] [(A11-A12.*zkm')*(A11-A12.*zhm)] [(A12-A22.*zkm')*(A12-A22.*zhm)]};
        covparamkk_DS = {[tau2] [(A11-A12.*zkm')*(A11-A12.*zkm)] [(A12-A22.*zkm')*(A12-A22.*zkm)]};
    elseif est_what == 2
        covmodelkd_DS = {'exponentialC_DS','exponentialC_DS'};
        covparamkd = {[(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1]};
        covparamkd_DS = {[A11*(A11-A12.*zhm)] [A12*(A12-A22.*zhm)]};
        covmodelkk_DS = {'exponentialC','exponentialC'};
        covparamkk_DS = {[A11.^2 3.*phi0] [A12.^2 3*phi1]};
    elseif est_what == 3
        covmodelkd_DS = {'exponentialC_DS','exponentialC_DS'};
        covparamkd = {[(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1]};
        covparamkd_DS = {[A12*(A11-A12.*zhm)] [A22*(A12-A22.*zhm)]};
        covmodelkk_DS = {'exponentialC','exponentialC'};
        covparamkk_DS = {[A12.^2 3.*phi0] [A22.^2 3*phi1]};
    end

    % calculate an estimate
    % options(1) = 1;   
    if est_what == 1
        tic
        [zk_mtr{i},vk{i}]=krigingME_DS(ckall{i}(:,1:2),ch,cs,zh,zs,vs,covmodel,covparam, ...
            covparamdd_DS,covparamkd_DS,covparamkk_DS,nhmax,nsmax,dmax,order,options);
        toc
        zk{i} = zk_mtr{i} + (beta0c{i}+beta1c{i}.*ck_mod); % add mean trend
    elseif est_what == 2
        tic
        [zk_mtr{i},vk{i}]=krigingME_DSB(ckall{i}(:,1:2),ch,cs,zh,zs,vs,covmodel,covparam, ...
            covparamdd_DS,covmodelkd_DS,covparamkd,covparamkd_DS,covmodelkk_DS,covparamkk_DS, ...
            nhmax,nsmax,dmax,order,options);
        toc
        zk{i} = zk_mtr{i} + (beta0c{i}); % add mean trend
    elseif est_what == 3
        tic
        [zk_mtr{i},vk{i}]=krigingME_DSB(ckall{i}(:,1:2),ch,cs,zh,zs,vs,covmodel,covparam, ...
            covparamdd_DS,covmodelkd_DS,covparamkd,covparamkd_DS,covmodelkk_DS,covparamkk_DS, ...
            nhmax,nsmax,dmax,order,options);
        toc
        zk{i} = zk_mtr{i} + (beta1c{i}); % add mean trend
    end

end
%matlabpool close

% saving LOOCE
save(sprintf('matdata/estimation_DS_%s.mat',est_what_str),'zk_mtr','zk','vk','ckall','unidates');

end