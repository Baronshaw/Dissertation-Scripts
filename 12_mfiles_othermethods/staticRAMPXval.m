function [] = staticRAMPXval()
% this function will create a static (spatial only) version of the RAMP
% method and perform a LOOCV

% only use the following if parallel computing is needed
% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "staticRAMPXval" -logfile "runall_Cluster12.out"

% loading BME function
cd ../BMELIB2.0b
startup
cd ../12_mfiles_othermethods

soft = 1; 
constant = 0; 
gauss = 1;
softstr = '_soft'; 
constr = '_long'; 
gaussstr = '_gauss'; 
forceddist = 0;

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
nsmax = 3;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

% perform BME calculations
len = size(cMS,1);
zk = cell(len,1); vk = cell(len,1); temp1 = cell(len,1); temp2 = cell(len,1);
temp1a = cell(len,1); temp1b = cell(len,1); allhard = cell(len,1);
alllambda1 = cell(len,1); alllambda2 = cell(len,1); allch = cell(len,1);
allcs = cell(len,1);

matlabpool open 12
for i = 1:size(cMS,1)
    disp(i);
    idxdist = sqrt( (cMS(i,1)-ch(:,1)).^2 + (cMS(i,2)-ch(:,2)).^2 ) > forceddist;
    idxsoft = sqrt( (cMS(i,1)-css(:,1)).^2 + (cMS(i,2)-css(:,2)).^2 ) < 1000000; % hard coded

    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ck = ch(idx&temp2001==2001|idx&temp2001==2002,:); % added 9/14/2014

    ch_nonXval = ch(~idx & idxdist,:);
    zh_mrmv_nonXval = zh_mrmv(~idx & idxdist,:);
    tic
    if ~isempty(ck)
        L2RMSE = lambda2(idxsoft); % added 11/5/2014
        [zk{i,1},vk{i,1},temp1{i,1},temp2{i,1},temp1a{i,1},temp1b{i,1}, ...
            allhard{i,1},alllambda1{i,1},alllambda2{i,1},allch{i,1},allcs{i,1}] ...
            =krigingME2_correct_parallel(ck,ch_nonXval,css(idxsoft,:),...
            zh_mrmv_nonXval,lambda1(idxsoft),L2RMSE,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    end
    toc
end
matlabpool close

% temporary save
save(sprintf('matdata/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone.mat',softstr,constr,gaussstr,floor(forceddist./1000)), ...
        'ck','nhmax','nsmax','zk','vk','temp1','temp2','temp1a','temp1b', ...
        'allhard','alllambda1','alllambda2','allch','allcs','forceddist');

% add mean trend back to cross-validation points
% reload to not confuse mI hard and mI soft
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50)); 
for i = 1:size(cMS,1)
    disp(i);
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ck = ch(idx&temp2001==2001|idx&temp2001==2002,:); % added 9/14/2014
    
    [liA locB] = ismember(ck,ch,'rows');
    locB(locB==0) = []; % so I can use locB 
    if ~isempty(locB)
        if gauss == 0
            zk_madd = moments(:,1) + mI(locB);
            mI_Xval = mI(locB);
            zh_Xval = zh(locB);
        elseif gauss == 1 & constant == 0
            zk_madd{i,1} = zk{i,1} + mI(locB);
            mI_Xval{i,1} = mI(locB); % offset
            zh_Xval{i,1} = zh(locB);
        elseif gauss == 1 & constant == 1
            zk_madd{i,1} = zk + mO;
            mI_Xval{i,1} = mO*ones(length(ck),1); % offset
            zh_Xval{i,1} = zh(locB);
        end
    end
end

% saving data
save(sprintf('matdata/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone.mat',softstr,constr,gaussstr,floor(forceddist./1000)), ...
    'ck','nhmax','nsmax','zk','vk','temp1','temp2','temp1a','temp1b', ...
    'allhard','alllambda1','alllambda2','allch','allcs','zk_madd','zh_Xval','mI_Xval','forceddist');

end