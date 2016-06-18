function [] = spacetimeDownScaler(additive,multiplicative)
% this function will estimate the parameters of the static downscaler
% method and attempt to impliment some sort of cross validation

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "spacetimeDownScaler" -logfile "runall_Cluster12.out"

if nargin < 1, additive = 'dyn'; end 
if nargin < 2, multiplicative = 'ind'; end

cd ../BMELIB2.0b
startup
cd ../12_mfiles_othermethods

% load paired modeled and observed data (years 2001 and 2002 for now)
for i = 2001:2002
    load(sprintf('../matfiles/prepCTMandObs_%d.mat',i));
    Modall{i-2000,1} = Mod; Obsall{i-2000,1} = Obs;
    coordObsall{i-2000,1} = coordObs; cht{i-2000,1} = yrmodaObs;
end
zhm = cell2mat(Modall); zho = cell2mat(Obsall);
ch = cell2mat(coordObsall); cht = cell2mat(cht);

% modifying dates
yrall = floor(cht./10000);
moall = floor((cht - yrall*10000)/100);
daall = cht  - yrall*10000 - moall*100;
cht = datenum(yrall,moall,daall);
unidates = unique(cht);
len = length(unique(cht));

% convert from s/t vector to s/t grid
[ZHO,cMS,tME,nanratio]=valstv2stg([ch cht],zho);
[ZHM,cMS,tME,nanratio]=valstv2stg([ch cht],zhm);

% estimation of parameters: mean trend
if strcmp(additive,'ind') & strcmp(multiplicative,'ind')
    
    % estimate beta0t and beta1t
    p = polyfit(zhm,zho,1);
    beta0t = p(2);
    beta1t = p(1); 
    
    % substract mean trend from data
    zho_mt = beta0t + beta1t*zhm;
    zho_mtr = zho - zho_mt;
    T = length(tME);

elseif strcmp(additive,'dyn') & strcmp(multiplicative,'ind')
    
    % estimate beta0t, beta1t, rho0, rho1, eta0t, eta1t
    beta0 = [0 0 0 0 0 0];
    LB = [-Inf -Inf -1 -1 -Inf -Inf];
    UB = [Inf Inf 1 1 Inf Inf];
    % running objective function
    options = optimset('Display','iter');
    [beta,fval,exitflag,output] = fminsearchbnd(@(beta)myErrFunMST(beta,ZHO,ZHM),beta0,LB,UB,options);
    beta0t = beta(1); beta1t = beta(2);
    rho0 = beta(3); rho1 = beta(4);
    eta0t = beta(5); eta1t = beta(6);
    
    % subtracting mean trend from data
    [cm T] = size(ZHO);
    for i = 1:T
        if i == 1, 
            all_beta0t(i,1) = beta0t;
            all_beta1t(i,1) = beta1t;
        else
            all_beta0t(i,1) = rho0*all_beta0t(i-1,1) + eta0t; 
            all_beta1t(i,1) = rho1*all_beta1t(i-1,1) + eta1t;
        end
    end
    ZHO_mt = rho0.*repmat(all_beta0t',cm,1) + eta0t + rho1.*repmat(all_beta1t',cm,1).*ZHM + eta1t.*ZHM;
    ZHO_mtr = ZHO - ZHO_mt;
    [ch_mtr,zho_mtr]=valstg2stv(ZHO_mtr,cMS,tME);
    [dummy,zhm]=valstg2stv(ZHM,cMS,tME);
    [dummy,zho_mt]=valstg2stv(ZHO_mt,cMS,tME);
    [dummy,zho]=valstg2stv(ZHO,cMS,tME);
    idx = isnan(zho_mtr);
    ch = ch_mtr(~idx,1:2); 
    cht = ch_mtr(~idx,3);
    zho_mtr = zho_mtr(~idx);    
    zhm = zhm(~idx);
    zho_mt = zho_mt(~idx);
    zho = zho(~idx);
       
elseif strcmp(additive,'ind') & strcmp(multiplicative,'dyn')
else
end

% covariance
% estimate of tau2
temp = zho_mt - zho;
a = prctile(temp,25); b = prctile(temp,75);
tau2 = var(temp(temp>=a&temp<=b));

% prepping objective function
rLags = cell(T,1);
Crtests = cell(T,1);
problems = cell(T,1);
lenp = 10;
for i = 1:T
    disp(i);
    idx = cht == tME(i);
    if sum(idx) > 0 
        [rLags{i,1},Crtests{i,1},problems{i,1}] = prepErrFunST(ch(idx,:),cht(idx),zho_mtr(idx),zhm(idx),lenp);
        problems{i,1}(end+1) = tau2;
    end
end
save('matdata/tempsave.mat','rLags','Crtests','problems');
load('matdata/tempsave.mat');
beta0 = [3 3 0 50000 50000];
LB = [0 0 -Inf 0 0];
UB = [Inf Inf Inf Inf Inf];

% running objective function
tic
options = optimset('Display','iter','TolFun',1e-2);
[beta,fval,exitflag,output] = fminsearchbnd(@(beta)myErrFunCST(beta,rLags,Crtests,problems),beta0,LB,UB,options);
toc
A11 = beta(1);
A12 = beta(2);
A22 = beta(3);
phi0 = beta(4);
phi1 = beta(5);

% save results
if strcmp(additive,'ind') & strcmp(multiplicative,'ind')
    
    save(sprintf('matdata/STDS_results_add_%s_muli_%s.mat',additive,multiplicative), ...
        'ZHO','ZHM','cMS','tME','beta0t','beta1t', ...
        'tau2','A11','A12','A22','phi0','phi1','rLags','Crtests','problems', ...
        'zho','zhm','ch','cht','zho_mt','zho_mtr');

elseif strcmp(additive,'dyn') & strcmp(multiplicative,'ind')
    
    save(sprintf('matdata/STDS_results_add_%s_muli_%s.mat',additive,multiplicative), ...
        'ZHO','ZHM','cMS','tME','beta0t','beta1t','rho0','rho1','eta0t','eta1t', ...
        'tau2','A11','A12','A22','phi0','phi1','rLags','Crtests','problems', ...
        'zho','zhm','ch','cht','zho_mt','zho_mtr');
    
elseif strcmp(additive,'ind') & strcmp(multiplicative,'dyn')
else
end


end