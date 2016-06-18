function [] = Xvalforcedisolation(years,soft,constant,gauss,percent,forceddist)
% this function will randomly select a certain percentage of monitoring
% stations (as defined by 'percent') and will isolate that monitor by
% removing all stations within a distance of 'forceddist'

if nargin < 1, years = 2001; end % year of time series
if nargin < 2, soft = 0; end % soft data or not
if nargin < 3, constant = 0; end % constant offset or not
if nargin < 4, gauss = 1; end % gaussian soft data or not
if nargin < 5, percent = 5; end
if nargin < 6, forceddist = 500; end % distance in kilometers

if soft == 0, softstr = '_nosoft'; else softstr = ''; end
if constant == 1, constr = '_constant'; else constr = ''; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = ''; end

% the 'long' offset is the final offset chosen
if constant == 1
    load('covariance_test/expcov__constantoffset_allyears.mat');
else
    load(sprintf('covariance_test/expcov__%s_allyears.mat','long'));
end

% years to look at time series
years = 2001;
load('distprepsoft_new.mat');

% gathering the hard data
zh = Obsg{1}(:); 
cht = repmat(tMEO{1},length(Obsg{1}),1); cht = cht(:);
chty = floor(cht./10000); chtm = floor( (cht-(chty*10000))./100 );
chtd = cht - chty*10000 - chtm*100;
cht = datenum(chty,chtm,chtd);
chs = repmat(cMSObs{1},length(tMEO{1}),1);
ch = [ chs cht ];
idx = isnan(zh);
zh(idx) = []; ch(idx,:) = [];
smoothingParam = [1400000 300000 100 50];
if ~exist('test_hardoffset_new.mat')
    cd mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv(ch,zh,smoothingParam,ch);
    toc
    cd ..
    save('test_hardoffset_new.mat','mI');
else
    load('test_hardoffset_new.mat');
end
if constant == 0
    zh = zh - mI;
else
    zh = zh - mO;
end

% calculate offset at CTM centriods
load('distprepsoft_new.mat');
zd = Obsg{1}(:);
pdt = repmat(tMEO{1},length(cMSObs{1}),1); pdt = pdt(:);
pdty = floor(pdt./10000); pdtm = floor( (pdt-(pdty*10000))./100 );
pdtd = pdt - pdty*10000 - pdtm*100;
pdt = datenum(pdty,pdtm,pdtd);
pd = [ repmat(cMSObs{1},length(tMEO{1}),1) pdt ];
idx = isnan(zd);
zd(idx) = []; pd(idx,:) = [];
pIt = repmat(tMECTM{1},length(cMSCTM{1}),1); pIt = pIt(:);
pIty = floor(pIt./10000); pItm = floor( (pIt-(pIty*10000))./100 );
pItd = pIt - pIty*10000 - pItm*100;
pIt = datenum(pIty,pItm,pItd);
pI = [ repmat(cMSCTM{1},length(tMECTM{1}),1) pIt ];
clear Modg Obsg % making space
smoothingParam = [1400000 300000 100 50];

if ~exist('test_softoffset_new2.mat')
    cd mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv(pd,zd,smoothingParam,pI);
    toc
    cd ..
    save('test_softoffset_new2.mat','mI','pI');
else
    load('test_softoffset_new2.mat');
end

% gathering all the soft data
if soft == 1
    css = pI;
    lambda1 = cell(length(tMECTM{1}),1);
    lambda2 = cell(length(tMECTM{1}),1);
    for i = 1:length(tMECTM{1})
        disp(i);
        load(sprintf('matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d_new.mat',tMECTM{1}(i),365,3,150,10,0));
        lambda1{i} = meanGivMod(:,1);
        lambda2{i} = varGivMod(:,1);
    end
    lambda1 = cell2mat(lambda1) - mI;
    lambda2 = cell2mat(lambda2);
    lambda2(lambda2<0) = 3;
elseif soft == 0
    css = [];
    lambda1 = [];
    lambda2 = [];
end

% load covariance information
if constant == 0
    load('covariance_test/covmodjoint_r_long_joint exponential exponential.mat');
elseif constant == 1
    load('covariance_test/covmodjoint_constantoffset_joint exponential exponential.mat');
end
covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};
covparam = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};
dmax = [1000000 365 f.alp*f.ar1/f.at1 + (1-f.alp)*f.ar2/f.at2];

% other BME parameters
softpdftype = 1; 
nhmax = 7;
nsmax = 3;
order = NaN;
options = BMEoptions;
options(1) = 1;
options(3) = 150000;

%%% Xvaldiation part

% picking the cross-validation locations
rand('seed',0);
len = floor((percent/100)*length(cMSObs{1})); % number of monitors for Xval
randidx = randperm(length(cMSObs{1}));
ckXval = cMSObs{1}(randidx(1:len),:); % Xval locations
[r c] = size(ckXval);

for i = 1:r
    disp([i r]);
    idxdist = sqrt( (ckXval(i,1)-cMSObs{1}(:,1)).^2 + (ckXval(i,2)-cMSObs{1}(:,2)).^2 );
    idxdist = idxdist <= forceddist*1000;

    % removing nearby stations
    [locA locB] = ismember(ch(:,1:2),cMSObs{1}(idxdist,:),'rows');
    chtemp = ch(~locA,:);
    zhtemp = zh(~locA);
    [locA locB] = ismember(ch(:,1:2),ckXval(i,:),'rows');
    subZ = zh(locA);
    ck{i} = ch(locA,:);

    if gauss == 0
        % finding probdens of the soft data given mean and variance
        len = length(lambda1);
        limi = arrayfun(@linspace,lambda1-1.96.*sqrt(lambda2),lambda1+1.96.*sqrt(lambda2),13*ones(len,1),...
            'uni',0);
        probdens = arrayfun(@truncNorm,limi,lambda1,lambda2,-mI,'uni',0);
        limi = cell2mat(limi);
        probdens = cell2mat(probdens); probdens = probdens(:,2:end);
        nl = 13*ones(len,1);

        % perfrom BMEprobaMoments
        [moments{i},info{i}]=BMEprobaMoments2(ck{i},chtemp,css,zhtemp,softpdftype, ...
            nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    elseif gauss == 1	
        [zk{i},vk{i},temp1{i},temp2{i}]=krigingME2(ck{i},chtemp,css,zhtemp,lambda1,lambda2, ...
            covmodel,covparam,nhmax,nsmax,dmax,order,options);
    end

    % add mean trend back
    load('test_hardoffset_new.mat'); % load again to not confuse mI from softdata and mI from harddata
    [liA locB] = ismember(ck{i},ch,'rows');
    locB(locB==0) = []; % so I can use locB 
    if gauss == 0
        zkresoff{i} = moments{i}(:,1) + mI(locB);
    elseif gauss == 1 & constant == 0
        zkresoff{i} = zk{i} + mI(locB);
    elseif gauss == 1 & constant == 1
        zkresoff{i} = zk{i} + mO;
    end
    mItimeSeries{i} = mI(locB); % offset
    subZresoff{i} = subZ + mI(locB);

end

if gauss == 1
    save(sprintf('Xval_results/XvalforcedisolationTESTING_per%dfor%d_%d%s%s%s.mat',percent,forceddist,years,softstr,constr,gaussstr), ...
        'ck','nhmax','nsmax','zk','vk','temp1','temp2','zkresoff','subZresoff','mItimeSeries');
elseif gauss == 0
    save(sprintf('Xval_results/XvalforcedisolationTESTING_per%dfor%d_%d%s%s%s.mat',percent,forceddist,years,softstr,constr,gaussstr), ...
    'ck','nhmax','nsmax','moments','info','zkresoff','subZresoff','mItimeSeries');
end

end