function [] = spacetimeDS_estimate(est_what,additive,multiplicative)
% this function will create estimates at given space/time locations for the
% static DS methods. This will create estimates as well as values for said
% parameters

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "spacetimeDS_estimate(1,'ind','ind')" -logfile "runall_Cluster12.out"

if nargin < 1, est_what = 3; end % 1=est,2=beta0,3=beta1
if nargin < 2, additive = 'dyn'; end 
if nargin < 3, multiplicative = 'ind'; end
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
    load(sprintf('../matfiles/prepCTM_%d.mat',i));
    distCTMvall{i-2000,1} = distCTMv; dailyCTMvall{i-2000,1} = dailyCTMv;
    yrmodaCTMvall{i-2000,1} = yrmodaCTMv;
end
chm_all = cell2mat(distCTMvall); zhm_all = cell2mat(dailyCTMvall);
chtm_all = cell2mat(yrmodaCTMvall);
chtm_all = datenum(chtm_all(:,1),chtm_all(:,2),chtm_all(:,3));

% once a month
unidates =  datenum([repmat(2001,24,1) [1:24]' repmat(1,24,1)]);
len = length(unidates);

% load parameter files
load(sprintf('matdata/STDS_results_add_%s_muli_%s.mat',additive,multiplicative));
problemz = problems;

% calculate all recursive betas for dyn add bias
if strcmp(additive,'dyn')
    for i = 1:length(tME)
        if i == 1, 
            all_beta0t(i,1) = beta0t;
            all_beta1t(i,1) = beta1t;
        else
            all_beta0t(i,1) = rho0*all_beta0t(i-1,1) + eta0t; 
            all_beta1t(i,1) = rho1*all_beta1t(i-1,1) + eta1t;
        end
    end
else
    all_beta0t = NaN*ones(length(tME),1);
    all_beta1t = NaN*ones(length(tME),1);
end
    
% kriging parameters
cs = []; zs = []; vs = [];
nhmax = 70; nsmax = 0;
order = NaN;
options = BMEoptions;
options(3) = 150000;
zh = zho_mtr;

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
    idxh = cht == unidates(i);
    idxm = chtm_all == unidates(i);
    chmallsub = chm_all(idxm,:); chtmallsub = chtm_all(idxm); zhmallsub = zhm_all(idxm);

    % intialize parameters
    zk_mtr{i} = NaN*ones(length(chmallsub),1); 
    vk{i} = NaN*ones(length(chmallsub),1); 
    zk{i} = NaN*ones(length(chmallsub),1);
    ckall{i} = [chmallsub chtmallsub];  

    % covariance parameters
    idxT = unique(cht) == unidates(i);
    zhmp = problemz{idxT};
    problems = NaN*ones(1,10);
    for k = 1:10
        problems(k) = (zhmp(k+1)-zhmp(k))./2 + zhmp(k); 
    end
    xB1 = problems(1); xB2 = problems(2); 
    covmodel = {'nuggetC/nuggetC','exponentialC_DS/nuggetC','exponentialC_DS/nuggetC'};
    covparam = {[tau2 tau2] [(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0 tau2] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1 tau2]};
    dmax3 = (3.*phi0+3.*phi1)/(tau2+tau2);
    dmax = [2000000 20 dmax3];

    % find modeled values of the observed data
    temp1 = arrayfun(@(x) find(abs(x-zhmp)==min(abs(x-zhmp))),zhm(idxh),'UniformOutput',false);
    temp2 = find( cell2mat(cellfun(@(x) length(x)>1,temp1,'UniformOutput',false)) == 1);
    for k = 1:length(temp2)
        temp1{temp2(k)} = temp1{temp2(k)}(1);
    end
    zhmM = zhmp(cell2mat(temp1));
    xB1M = repmat(zhmM,length(zhmM),1); xB2M = repmat(zhmM',1,length(zhmM));
    covparamdd_DS = {[tau2] [(A11-A12.*xB1M).*(A11-A12.*xB2M)] [(A12-A22.*xB1M).*(A12-A22.*xB2M)]};

    % find modeled values of the estimation locations
    if ismember(ckall{i}(:,1:2),ch,'rows') == 1 % estimation location in the data 
        [ia ib] = ismember(ckall{i},ch,'rows');
        ck_mod = zhm_paired(ib);
    elseif isequal(ckall{i}(:,1:2),chmallsub) % estimation locations are modeled locations
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
        covparamkd_DS = {[tau2] [(A11-A12.*zkm')*(A11-A12.*zhmM)] [(A12-A22.*zkm')*(A12-A22.*zhmM)]};
        covparamkk_DS = {[tau2] [(A11-A12.*zkm')*(A11-A12.*zkm)] [(A12-A22.*zkm')*(A12-A22.*zkm)]};
    elseif est_what == 2
        covmodelkd_DS = {'exponentialC_DS/nuggetC','exponentialC_DS/nuggetC'};
        covparamkd = {[(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0 tau2] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1 tau2]};
        covparamkd_DS = {[A11*(A11-A12.*zhmM) tau2] [A12*(A12-A22.*zhmM) tau2]};
        covmodelkk_DS = {'exponentialC/nuggetC','exponentialC/nuggetC'};
        covparamkk_DS = {[A11.^2 3.*phi0 tau2] [A12.^2 3*phi1 tau2]};        
    elseif est_what == 3
        covmodelkd_DS = {'exponentialC_DS/nuggetC','exponentialC_DS/nuggetC'};
        covparamkd = {[(A11-A12.*xB1).*(A11-A12.*xB2) 3.*phi0 tau2] [(A12-A22.*xB1).*(A12-A22.*xB2) 3*phi1 tau2]};
        covparamkd_DS = {[A12*(A11-A12.*zhmM) tau2] [A22*(A12-A22.*zhmM) tau2]};
        covmodelkk_DS = {'exponentialC/nuggetC','exponentialC/nuggetC'};
        covparamkk_DS = {[A12.^2 3.*phi0 tau2] [A22.^2 3*phi1 tau2]};  
    end

    % calculate an estimate
    options(1) = 1;
    if est_what == 1
        tic
        [zk_mtr{i},vk{i}]=krigingME_DS(ckall{i},[ch cht],cs,zh,zs,vs,covmodel,covparam, ...
            covparamdd_DS,covparamkd_DS,covparamkk_DS,nhmax,nsmax,dmax,order,options);
        toc
        zk{i} = zk_mtr{i} + (ck_mod); % add mean trend
    elseif est_what == 2       
        tic
        [zk_mtr{i},vk{i}]=krigingME_DSB(ckall{i},[ch cht],cs,zh,zs,vs,covmodel,covparam, ...
            covparamdd_DS,covmodelkd_DS,covparamkd,covparamkd_DS,covmodelkk_DS,covparamkk_DS, ...
            nhmax,nsmax,dmax,order,options);
        toc
        if strcmp(additive,'ind')
            zk{i} = zk_mtr{i} + (beta0t); % add mean trend
        else strcmp(additive,'dyn')
            idxdyn = unidates(i) == tME;
            zk{i} = zk_mtr{i} + (rho0.*all_beta0t(idxdyn) + eta0t); % add mean trend
        end
    elseif est_what == 3        
        tic
        [zk_mtr{i},vk{i}]=krigingME_DSB(ckall{i},[ch cht],cs,zh,zs,vs,covmodel,covparam, ...
            covparamdd_DS,covmodelkd_DS,covparamkd,covparamkd_DS,covmodelkk_DS,covparamkk_DS, ...
            nhmax,nsmax,dmax,order,options);
        toc
        if strcmp(additive,'ind')
            zk{i} = zk_mtr{i} + (beta1t); % add mean trend
        else strcmp(additive,'dyn')
            idxdyn = unidates(i) == tME;
            zk{i} = zk_mtr{i} + (rho1.*all_beta1t(idxdyn) + eta1t); % add mean trend
        end
    end

end
%matlabpool close

% saving LOOCE
save(sprintf('matdata/estimation_stDS_%s_add_%s_muli_%s.mat',est_what_str,additive,multiplicative), ...
    'zk_mtr','zk','vk','ckall','unidates');

end