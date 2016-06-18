function [] = SeeEffectOfLambda1(soft,constant,gauss,yearz)
% this function will go through ever ck and see the closest lambda1. This
% will determine if lambda1 SHOULD have helped or hurt the estimation of ck

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, yearz = 1999; end

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;

% loading soft data
CTMyears = [2001 2002 2005:2007]; 
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
    [lia lib] = ismember(pI,cssA{i,1},'rows');
    lib(lib==0)=[];

end  
css = cell2mat(cssA);
limi = cell2mat(limiA);
if gauss == 0, probdens = cell2mat(probdensA);
else lambda1 = cell2mat(lambda1A); lambda2 = cell2mat(lambda2A); end

% load covariance information
if constant == 0
    load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
elseif constant == 1
    load('../matfiles/covmod_constant_joint exponential exponential_joint.mat');
end
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

% loop through and get all the values in the neighborhoods
zs = lambda1;
vs = lambda2;
   
idxh = ch(:,3) <= datenum(yearz,12,31) & ch(:,3) >= datenum(yearz,1,1);
idxs = css(:,3) < datenum(yearz+1,12,31) & css(:,3) > datenum(yearz-1,1,1);

ck = ch(idxh,:);
zslocal = cell(sum(idxh),1);
vslocal = cell(sum(idxh),1);
vs = vs(idxs);

if sum(idxs) > 0
    parfor i = 1:length(ck)
        ck0 = ck(i,:);
        tic
        [cslocal,zslocal{i},dh,sumnslocal,index]=neighbours2a(ck0,css(idxs,:),zs(idxs),nsmax,dmax);
        vslocal{i}=vs(index);
        toc
    end
end

save(sprintf('dataneighb_%d.mat',yearz),'ck','ch','zh','zslocal','vslocal');

end