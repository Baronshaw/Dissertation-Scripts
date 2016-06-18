function [] = staticKriging_estimate()
% this function will create estimates at given space/time locations for the
% static Kriging method. 

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "staticKriging_estimate" -logfile "runall_Cluster12k.out"
    
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

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));

% added later
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
zh_mrmv = zh - mI;

% load covariance information
load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};
covparam = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};
dmax = [2000000 0 99999999999999]; % spatial only

% gathering all the soft data locations
CTMyears = [2001 2002];  
for i = 1:length(CTMyears) 
    disp(CTMyears(i));
    load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',CTMyears(i)));
    cssA{i,1} = css;
    limiA{i,1} = limi;
    lambda1A{i,1} = lambda1;
    lambda2A{i,1} = lambda2;

    load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 300000 100 50],CTMyears(i)));
    pI = round(pI); cssA{i} = round(cssA{i}); % dealing with computing artifacts
    [lia lib] = ismember(pI,cssA{i,1},'rows');
    lib(lib==0)=[];
    mI(isnan(mI)) = 0; % added 9/28/2014
    lambda1A_mrmv{i,1} = lambda1A{i,1}-mI(lib);        
end  
    
css = cell2mat(cssA);
limi = cell2mat(limiA);
lambda1 = cell2mat(lambda1A_mrmv); 
lambda2 = cell2mat(lambda2A); 

% other BME parameters
softpdftype = 1; 
nhmax = 7;
nsmax = 0;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

%%% estimates
% looping through each day
zk_mtr = cell(len,1);
vk = cell(len,1);
temp1 = cell(len,1);
temp2 = cell(len,1);
temp1a = cell(len,1);
temp1b = cell(len,1);
allhard = cell(len,1);
alllambda1 = cell(len,1);
alllambda2 = cell(len,1);
allch = cell(len,1);
allcs = cell(len,1);
zk = cell(len,1);
ckall = cell(len,1);
matlabpool open 12
for i = 1:len
    disp([i len]);
    % getting data from specific day
    idx = chtm_all == unidates(i);
    chmallsub = chm_all(idx,:); chtmallsub = chtm_all(idx);

    % intialize parameters
    ckall{i} = [chmallsub chtmallsub];  

    [zk_mtr{i,1},vk{i,1},temp1{i,1},temp2{i,1},temp1a{i,1},temp1b{i,1}, ...
            allhard{i,1},alllambda1{i,1},alllambda2{i,1},allch{i,1},allcs{i,1}] ...
            =krigingME2_correct_parallel(ckall{i},ch,css,...
            zh_mrmv,lambda1,lambda2,covmodel,covparam,nhmax,nsmax,dmax,order,options);

end
matlabpool close

% add mean trend back to cross-validation points
% reload to not confuse mI hard and mI soft
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50)); 
pd = pd; zd = zd;
matlabpool open 12
parfor i = 1:length(ckall)
    disp(i);
    cd ../10_mfiles_newmeantrend
    [mI]=expKernelSmooth_stv(pd,zd,[900000,300000,100,50],ckall{i});
    cd ../12_mfiles_othermethods
    zk{i} = zk_mtr{i} + mI;
end
matlabpool close

% saving results
save(sprintf('matdata/estimation_kriging_foriso%dkm_timezone.mat',0), ...
    'nhmax','nsmax','temp1','temp2','temp1a','temp1b', ...
    'zk_mtr','zk','vk','ckall','unidates','mI');

end