function [] = BMEmaps(daydisp,soft,constant,gauss)
% this will calculate the posterior BME mean and variance for maps

if nargin < 1, daydisp = [2001 7 1]; end
if nargin < 2, soft = 1; end % soft data or not
if nargin < 3, constant = 0; end % constant offset or not
if nargin < 4, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../07_mfiles_map

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;

% subtract mean trend from all the data
if constant == 0
    zh_mrmv = zh - mI;
else
    mO = mean(zh);
    zh_mrmv = zh - mO;
end

% load covariance information
if constant == 0
    load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
elseif constant == 1
    load('../matfiles/covmod_constant_joint exponential exponential_joint.mat');
end
covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};
covparam = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};
dmax = [2000000 365 f.alp*f.ar1/f.at1 + (1-f.alp)*f.ar2/f.at2];

% estimation locations
load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',2001));
cMSCTM = unique(css(:,1:2),'rows');
ck = [cMSCTM repmat(datenum(daydisp),length(cMSCTM),1)]; 

% gathering all the soft data locations
if soft == 1
    
    CTMyears = [2001 2002]; % excluding 2005:2007 bc not working out
    for i = 1:length(CTMyears) 
        disp(CTMyears(i));
        load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',CTMyears(i)));
        cssA{i,1} = css;
        limiA{i,1} = limi;
        
        if gauss == 0
            probdensA{i,1} = probdens;
        else
            lambda1A{i,1} = lambda1;
            lambda2A{i,1} = lambda2;
        end
        
        load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 300000 100 50],CTMyears(i)));
        [lia lib] = ismember(round(pI),round(cssA{i,1}),'rows');
        lib(lib==0)=[];
        
        if gauss == 0 & constant == 0
            probdensA_mrmv{i,1} = NaN*ones(size(probdensA{i,1}));
            for j = 1:size(lib,1)                
                probdensA_mrmv{i,1}(j,:) = probdensA{i,1}(j,:)-mI(lib(j)); 
            end
        elseif gauss == 1 & constant == 0
            mI(isnan(mI)) = 0; % added 7/9/2015
            lambda1A_mrmv{i,1} = lambda1A{i,1}-mI(lib);
        else
            disp('something strange happened');
        end
        
    end  
    
    css = cell2mat(cssA);
    limi = cell2mat(limiA);
    
    if gauss == 0, probdens = cell2mat(probdensA_mrmv);
    else lambda1 = cell2mat(lambda1A_mrmv); lambda2 = cell2mat(lambda2A); end
    
else
    
    css = [];
    lambda1 = [];
    lambda2 = [];
    
end

% other BME parameters
softpdftype = 1; 
nhmax = 7;
if soft == 1, nsmax = 3; else nsmax = 0; end
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

% perform BME calculations
if gauss == 0

    len = length(css);
    nl = 13*ones(len,1);
    % perfrom BMEprobaMoments
    [tempa tempb] = unique(ch,'rows');
    [moments,info]=BMEprobaMoments2_parallel(ck,ch(tempb,:),css,zh_mrmv(tempb),softpdftype, ...
        nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    
elseif gauss == 1
    
    [zk,vk,temp1,temp2,temp1a,temp1b]=krigingME2_correct_parallel(ck,ch,css,...
        zh_mrmv,lambda1,lambda2,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    
end
save('temp.mat','zk','vk','temp1','temp2','temp1a','temp1b');
% add mean trend back to cross-validation points
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50)); 
cd ../10_mfiles_newmeantrend
smoothingParam = [2*900000 300000 100 50]; 
[mI]=expKernelSmooth_stv_parallel(pd,zd,smoothingParam,ck);
toc
cd ../07_mfiles_map
if gauss == 0
    zk_madd = moments(:,1) + mI;
elseif gauss == 1 & constant == 0
    zk_madd = zk + mI;
elseif gauss == 1 & constant == 1
    zk_madd = zk + mO;
end

% saving results
if gauss == 1
    save(sprintf('../matfiles/BMEmaps_%d_%0.2d_%0.2d%s%s%s.mat',daydisp(1),daydisp(2),daydisp(3), ...
        softstr,constr,gaussstr),'ck','nhmax','nsmax','zk','vk','zk_madd','mI');
else
    save(sprintf('../matfiles/BMEmaps_%d_%0.2d_%0.2d%s%s%s.mat',daydisp(1),daydisp(2),daydisp(3), ...
        softstr,constr,gaussstr),'ck','nhmax','nsmax','moments','info','zk_madd','mI');
end

end