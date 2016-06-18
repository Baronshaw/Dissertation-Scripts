function [] = WHIMS_est(monthsdate,testrun,soft,constant,gauss)
% this is the main function for the WHIMS run 
tic
if nargin < 1, monthsdate = 96; end 
if nargin < 2, testrun = 1; end
if nargin < 3, soft = 1; end % soft data or not
if nargin < 4, constant = 0; end % constant offset or not
if nargin < 5, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

if testrun == 1 
    teststr = 'test_'; 
    teststr2 = 'QAQC';
elseif testrun == 0
    teststr = ''; 
    teststr2 = 'est';
elseif testrun == 2
    teststr = 'BMEest'; 
    teststr2 = 'est';
end

% sh test2.sh
% ##!/bin/bash
% for ((i=10;i<=160;i++))
% do
%   bsub -q week matlab -nodesktop -nosplash -singleCompThread -r "analysis(${i})$
% done

% begining variables
inputdates1 = datevec(datenum(1998,monthsdate,1));
inputdates2 = datevec(datenum(1998,monthsdate+1,1)-1);
startdate = inputdates1(:,1:3);
enddate = inputdates2(:,1:3);

daydisp1 = startdate(1)*10000 + startdate(2)*100 + startdate(3);
daydisp2 = enddate(1)*10000 + enddate(2)*100 + enddate(3);

if testrun == 1
    % load mock locations
    load('../matfiles_QAQC/forQA.mat');
    partt = datenum(startdate(1),startdate(2),startdate(3)):datenum(enddate(1),enddate(2),enddate(3));
    allpartt = repmat(partt,length(simux),1);
    ck = [ repmat([simux simuy],length(partt),1) allpartt(:) ];
elseif testrun == 0
    % load participant locations
    load('../matfiles_est/partdata.mat');
    partt = datenum(startdate(1),startdate(2),startdate(3)):datenum(enddate(1),enddate(2),enddate(3));
    allpartt = repmat(partt,length(partx),1);
    ck = [ repmat([partx party],length(partt),1) allpartt(:) ];
end

% loading mean trend removed data
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

% gathering all the soft data locations
if soft == 1
    
    CTMyears = [2001 2002 2005:2007]; 
    for i = 1:length(CTMyears) 
        disp(CTMyears(i));
        load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',CTMyears(i)));
        
        % for the sake of efficiency, only pick nearest temporal locations
        idx_sclose = css(:,3) >= datenum(startdate(1),startdate(2),startdate(3)) - 365 & ...
            css(:,3) <= datenum(enddate(1),enddate(2),enddate(3)) + 365;
        
        cssA{i,1} = css(idx_sclose,:);
        limiA{i,1} = limi(idx_sclose,:);
        
        if gauss == 0
            probdensA{i,1} = probdens(idx_sclose,:);
        else
            lambda1A{i,1} = lambda1(idx_sclose,:);
            lambda2A{i,1} = lambda2(idx_sclose,:);
        end
        
        load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 300000 100 50],CTMyears(i)));
        [lia lib] = ismember(pI,cssA{i,1},'rows');
        lib(lib==0)=[];
        
        if gauss == 0 & constant == 0
            probdensA_mrmv{i,1} = NaN*ones(size(probdensA{i,1}));
            for j = 1:size(lib,1)                
                probdensA_mrmv{i,1}(j,:) = probdensA{i,1}(j,:)-mI(lib(j)); 
            end
        elseif gauss == 1 & constant == 0
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

% for the sake of efficiency, only pick nearest temporal locations
idx_hclose = ch(:,3) >= datenum(startdate(1),startdate(2),startdate(3)) - 365 & ...
    ch(:,3) <= datenum(enddate(1),enddate(2),enddate(3)) + 365;

% other BME parameters
softpdftype = 1; 
nhmax = 7;
nsmax = 3;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

if testrun == 2 % BME background estimation
    % finding the bounds of US and which day has the most hard data
    load('../09_mfiles_projections/USAcontiguous.mat');
    cd ../09_mfiles_projections
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../08_mfiles_QAQC
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    xg=linspace(lax(1),lax(2),148);
    yg=linspace(lax(3),lax(4),112);
    [xg yg]=meshgrid(xg,yg); 
    testdays = datenum(startdate(1),startdate(2),12:18);
    bestday = cell2mat( arrayfun(@(x) sum(ch(:,3)==x),testdays,'UniformOutput',false) );
    ckt = testdays(bestday==max(bestday)); ckt = ckt(1); % if repeates
    ck = [xg(:) yg(:) repmat(ckt,length(xg(:)),1)];
end

% perform BME calculations
if gauss == 0

    len = length(css);
    nl = 13*ones(len,1);
    [moments,info]=BMEprobaMoments2_parallel(ck,ch(idx_hclose,:),css,zh_mrmv(idx_hclose),softpdftype, ...
        nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    
elseif gauss == 1
    
    [zk,vk,temp1,temp2,temp1a,temp1b,allhard,alllambda1,alllambda2,allch,allcs]=...
            krigingME2_correct_parallel(ck,ch(idx_hclose,:),css,...
            zh_mrmv(idx_hclose),lambda1,lambda2,covmodel,covparam,nhmax,nsmax,dmax,order,options);
        
end 

% caclulate mean trend of estimation locations
smoothingParam = [900000 300000 100 50];
cd ../10_mfiles_newmeantrend
[mItest]=expKernelSmooth_stv_parallel(ch(idx_hclose,:),zh(idx_hclose),smoothingParam,ck);
cd ../06_mfiles_estimation
save(sprintf('../matfiles_%s/%smI_%d_%d.mat',teststr2,teststr,daydisp1,daydisp2),'mItest','ck');

% save results
if gauss == 1
    save(sprintf('../matfiles_%s/%sWHIMS_%d_%d%s%s%s.mat', ...
        teststr2,teststr,daydisp1,daydisp2,softstr,constr,gaussstr), ...
        'ck','nhmax','nsmax','zk','vk','temp1','temp2','temp1a','temp1b',...
        'allhard','alllambda1','alllambda2','allch','allcs');
elseif gauss == 0
    save(sprintf('../matfiles_%s/%sWHIMS_%d_%d%s%s%s.mat', ...
        teststr2,teststr,daydisp1,daydisp2,softstr,constr,gaussstr), ...
        'ck','nhmax','nsmax','moments','info');
end

end