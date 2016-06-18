function [] = getXvalforcedisolation_LOOCV(soft,constant,gauss,forceddist)
% this will calculate the spatial LOO cross validation for all the 
% observed data with an increaseing spatial radius

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, forceddist = 0; end % units = meters

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% loading BME function
cd ../BMELIB2.0b
startup
cd ../05_mfiles_crossvalidation

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));

% added later
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
% temp2005 = datevec(pd(:,3)); temp2005 = temp2005(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');

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

% gathering all the soft data locations
if soft == 1
    
    CTMyears = [2001 2002 2005]; % [2001 2002 2005:2007]; 
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
        pI = round(pI); cssA{i} = round(cssA{i}); % dealing with computing artifacts
        [lia lib] = ismember(pI,cssA{i,1},'rows');
        lib(lib==0)=[];
        
        if gauss == 0 & constant == 0
            probdensA_mrmv{i,1} = NaN*ones(size(probdensA{i,1}));
            for j = 1:size(lib,1)                
                probdensA_mrmv{i,1}(j,:) = probdensA{i,1}(j,:)-mI(lib(j)); 
            end
        elseif gauss == 1 & constant == 0
            mI(isnan(mI)) = 0; % added 9/28/2014
            lambda1A_mrmv{i,1} = lambda1A{i,1}-mI(lib);
        else
            disp('something strange happened');
        end
        
    end  
    
    css = cell2mat(cssA);
    limi = cell2mat(limiA);
    
    if gauss == 0, probdens = cell2mat(probdensA_mrmv);
    else
        lambda1 = cell2mat(lambda1A_mrmv); 
        lambda2 = cell2mat(lambda2A); 
    end
    
else
    
    css = [];
    lambda1 = [];
    lambda2 = [];
    
end

% other BME parameters
softpdftype = 1; 
nhmax = 7;
nsmax = 3;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

% perform BME calculations
if gauss == 0

    len = length(css);
    nl = 13*ones(len,1);
    % perfrom BMEprobaMoments
    [tempa tempb] = unique(ch_nonXval,'rows');
    [moments,info]=BMEprobaMoments2_parallel(ck,ch_nonXval(tempb,:),css,zh_mrmv_nonXval(tempb),softpdftype, ...
        nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    
elseif gauss == 1
    len = size(cMS,1);
    ckXval = cell(len,1); zk = cell(len,1); vk = cell(len,1); temp1 = cell(len,1); temp2 = cell(len,1);
    temp1a = cell(len,1); temp1b = cell(len,1); allhard = cell(len,1);
    alllambda1 = cell(len,1); alllambda2 = cell(len,1); allch = cell(len,1);
    allcs = cell(len,1);
    matlabpool open 12
    for i = 1:size(cMS,1)
        disp(i);
        idxdist = sqrt( (cMS(i,1)-ch(:,1)).^2 + (cMS(i,2)-ch(:,2)).^2 ) > forceddist;
        if soft == 1
            idxsoft = sqrt( (cMS(i,1)-css(:,1)).^2 + (cMS(i,2)-css(:,2)).^2 ) < 1000000; % hard coded
        end
        % above line added for compuational time
        
        idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
        ck = ch(idx,:); ckXval{i,1} = ck; % added 9/14/2014
        ch_nonXval = ch(~idx & idxdist,:);
        zh_mrmv_nonXval = zh_mrmv(~idx & idxdist,:);
        tic
        if soft == 1
            if ~isempty(ck)
                [zk{i,1},vk{i,1},temp1{i,1},temp2{i,1},temp1a{i,1},temp1b{i,1}, ...
                    allhard{i,1},alllambda1{i,1},alllambda2{i,1},allch{i,1},allcs{i,1}] ...
                    =krigingME2_correct_parallel(ck,ch_nonXval,css(idxsoft,:),...
                    zh_mrmv_nonXval,lambda1(idxsoft),lambda2(idxsoft),covmodel,covparam,nhmax,nsmax,dmax,order,options);
            end
        else
            nsmax = 0;
            [zk{i,1},vk{i,1},temp1{i,1},temp2{i,1},temp1a{i,1},temp1b{i,1}, ...
                allhard{i,1},alllambda1{i,1},alllambda2{i,1},allch{i,1},allcs{i,1}] ...
                =krigingME2_correct_parallel(ck,ch_nonXval,css,...
                zh_mrmv_nonXval,lambda1,lambda2,covmodel,covparam,nhmax,nsmax,dmax,order,options);
        end
        toc
            
    end
    matlabpool close
    
end

% temporary save
save(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_20150408.mat',softstr,constr,gaussstr,floor(forceddist./1000)), ...
        'ckXval','nhmax','nsmax','zk','vk','temp1','temp2','temp1a','temp1b', ...
        'allhard','alllambda1','alllambda2','allch','allcs','forceddist');

% add mean trend back to cross-validation points
% reload to not confuse mI hard and mI soft
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50)); 
for i = 1:size(cMS,1)
    
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ck = ch(idx,:); % added 9/14/2014

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
if gauss == 1
    save(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_20150408.mat',softstr,constr,gaussstr,floor(forceddist./1000)), ...
        'ckXval','nhmax','nsmax','zk','vk','temp1','temp2','temp1a','temp1b', ...
        'allhard','alllambda1','alllambda2','allch','allcs','zk_madd','zh_Xval','mI_Xval','forceddist');
elseif gauss == 0
    save(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',fold,softstr,constr,gaussstr), ...
        'ck','nhmax','nsmax','moments','info','zk_madd','mI_Xval','zh_Xval');
end

end