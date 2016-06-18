function [] = seeNeib_one(include_mean)
% this function will get the soft data neighborhood of all data locations

% bsub -x -q week -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "seeNeib_one(1)" -logfile "seeNeib_one1.out"
% bsub -x -q week -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "seeNeib_one(0)" -logfile "seeNeib_one0.out"

if nargin < 1, include_mean = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../05_mfiles_crossvalidation

soft = 1; % soft data or not
constant = 0;  % constant offset or not
gauss = 1; % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');

% load covariance information
if constant == 0
    load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
elseif constant == 1
    load('../matfiles/covmod_constant_joint exponential exponential_joint.mat');
end
covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};
covparam = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};
dmax = [2000000 365 f.alp*f.ar1/f.at1 + (1-f.alp)*f.ar2/f.at2];

% gathering all the soft data location  
CTMyears = [2001 2002 2005:2007]; 
for i = 1:length(CTMyears) 
    disp(CTMyears(i));
    load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',CTMyears(i)));
    cssA{i,1} = css;
    limiA{i,1} = limi;

    lambda1A{i,1} = lambda1;
    lambda2A{i,1} = lambda2;

    load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 300000 100 50],CTMyears(i)));
    [lia lib] = ismember(pI,cssA{i,1},'rows');
    lib(lib==0)=[];
    if include_mean == 1
        temp = mI(lib);
        temp(isnan(temp)) = 0;
        lambda1A_mrmv{i,1} = lambda1A{i,1} - temp;
    else
        lambda1A_mrmv{i,1} = lambda1A{i,1};
    end

end  

css = cell2mat(cssA);
limi = cell2mat(limiA);

lambda1 = cell2mat(lambda1A_mrmv); 
lambda2 = cell2mat(lambda2A); 
    
% other BME parameters
softpdftype = 1; 
nhmax = 7;
nsmax = 3;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

% perform BME calculations
len = size(ch,1); 
alllambda1 = cell(len,1); 
alllambda2 = cell(len,1); 
allcs = cell(len,1);

matlabpool open 12
parfor i = 1:size(ch,1)
    disp(i);
    
    tic
    [allcs{i,1},alllambda1{i,1},dh,sumnslocal,index]=neighbours2a(ch(i,:),css,lambda1,nsmax,dmax);
    alllambda2{i,1}=lambda2(index);
    toc

end
matlabpool close

% saving results
save(sprintf('../matfiles/GetSoftNeib_%d.mat',include_mean), ...
    'ch','nhmax','nsmax','alllambda1','alllambda2','allcs');

end