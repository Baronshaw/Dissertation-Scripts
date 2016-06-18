function [] = runall_Cluster_11d1(FOLDIDX)
% this function will do the estimation for 2005 for the true Xval

if nargin < 1, FOLDIDX = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../11_mfiles_truecrossvalidation

% load all data
smoothingParam = [900000 300000 100 50];
load(sprintf('../matfiles/meanTrend_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam));
zh_mtr_nonXval = zh_nonXval - mI_nonXval;

% load soft data
smoothingParam = [2*900000 300000 100 50];
load(sprintf('../matfiles/meanTrendSoft_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam));
CTMyears = [2001 2002 2005:2007]; 
for i = 1:length(CTMyears) 
    disp(CTMyears(i));
    load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',CTMyears(i)));
    cssA{i,1} = css;
    limiA{i,1} = limi;

    lambda1A{i,1} = lambda1;
    lambda2A{i,1} = lambda2;

    smoothingParam = [2*900000 300000 100 50];
    load(sprintf('../matfiles/meanTrendSoft_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam));
    [lia lib] = ismember(round(pI),round(cssA{i,1}),'rows');
    lib(lib==0)=[];
    lambda1A_mtr{i,1} = lambda1A{i,1}-mIsoft_nonXval(lib);
end 
css = cell2mat(cssA);
lambda1 = cell2mat(lambda1A_mtr); 
lambda2 = cell2mat(lambda2A);

% load covariance model
load(sprintf('../matfiles/covmod_true10fold_%d.mat',FOLDIDX));
covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};
covparam = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};
dmax = [2000000 365 f.alp*f.ar1/f.at1 + (1-f.alp)*f.ar2/f.at2];

% other BME parameters
softpdftype = 1; 
nhmax = 7;
nsmax = 3;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

% added 1/7/2015 to calculate hard only
css = [];
lambda1 = [];
lambda2 = [];
nsmax = 0;

% calculate estimation at cross validation locations
[zk,vk,temp1,temp2,temp1a,temp1b,allhard,alllambda1,alllambda2,allch,allcs] ...
        = krigingME2_correct_parallel(ch_Xval,ch_nonXval,css,...
        zh_mtr_nonXval,lambda1,lambda2,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    
% calculate and add back mean trend
% mean trend calculation
cd ../10_mfiles_newmeantrend
tic
smoothingParam = [900000 300000 100 50];
[mI_Xval]=expKernelSmooth_stv_parallel(ch_nonXval,zh_nonXval,smoothingParam,ch_Xval);
toc
cd ../11_mfiles_truecrossvalidation 
zk_madd = zk + mI_Xval;

% updated 1/7/2015 to calculate hard only
save(sprintf('../matfiles/Xval_true10fold_fold%d_hardonly.mat',FOLDIDX), ...
        'ch_Xval','nhmax','nsmax','zk','vk','temp1','temp2','temp1a','temp1b', ...
        'allhard','alllambda1','alllambda2','allch','allcs','mI_Xval','zk_madd');

end