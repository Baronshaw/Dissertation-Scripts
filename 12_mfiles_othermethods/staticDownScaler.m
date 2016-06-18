function [] = staticDownScaler()
% this function will estimate the parameters of the static downscaler
% method and attempt to impliment some sort of cross validation

% parameters to estimate: beta0, beta1, A11, A12, A22, tau2, phi0, phi1

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "staticDownScaler" -logfile "runall_Cluster12.out"

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

% initialize variables
beta0c = cell(len,1);
beta1c = cell(len,1);
tau2c = cell(len,1);
zhmpc = cell(len,1);
betac = cell(len,1);
beta0cg = cell(len,1);
beta1cg = cell(len,1);
rLags = cell(len,1);
Crtests = cell(len,1);
problems = cell(len,1);

% loop through every day to estimate variables
for i = 1:len
    disp(i);
    % getting data from specific day
    idx = cht == unidates(i);
    chsub = ch(idx,:); chtsub = cht(idx); zhosub = zho(idx); zhmsub = zhm(idx);
    
    % estimate beta0 and beta1
    p = polyfit(zhmsub,zhosub,1);
    beta0c{i} = p(2);
    beta1c{i} = p(1); 
    
    % removing mean trend from the data
    zhosub_mtr = zhosub - (beta0c{i} + beta1c{i}.*zhmsub);

    % estimate of tau2
    temp = (beta0c{i}+beta1c{i}.*zhmsub)-zhosub;
    a = prctile(temp,25); b = prctile(temp,75);
    tau2c{i} = var(temp(temp>=a&temp<=b));
    
    % percentile of modeled data
    zhmpc{i} = prctile(zhmsub,0:10:100);
    zhmp = prctile(zhmsub,0:10:100);
    lenp = length(zhmp) - 1;
    
    % prepping objective function
    [rLags{i,1},Crtests{i,1},problems{i,1}] = prepErrFun(chsub,chtsub,zhosub_mtr,zhmsub,zhmp,lenp);
    problems{i,1}(end+1) = tau2c{i};
    beta0 = [3 3 0 50000 50000];
    LB = [0 0 -Inf 0 0];
    UB = [Inf Inf Inf Inf Inf];
    
    % running objective function
    tic
    options = optimset('Display','iter');
    [betac{i,1},fval,exitflag,output] = fminsearchbnd(@(beta)myErrFunC(beta,rLags{i,1},Crtests{i,1},problems{i,1}),beta0,LB,UB,options);
    toc
    
    % GLS estimate of beta0 and beta1
    try
        [beta0cg{i,1},beta1cg{i,1}] = getGLSB(betac{i,1},problems{i,1},chsub,zhosub,zhmsub,zhmp);
    catch
        beta0cg{i,1} = NaN; betacg{i,1} = NaN;
    end

end

% save results
save('matdata/DS_results_v1.mat','zhmpc','cht','beta0c','beta1c','tau2c','betac','beta0cg','beta1cg','problems','rLags','Crtests');

end