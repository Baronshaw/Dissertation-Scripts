function [overall,years,regions,yearregions,stations,yearstations, ...
    nears,yearnears,folds,yearfolds,uniyr,uniR,uniMS,unifold] = stat_MP_bme(soft,forceddist,loo)
% this function will calculate the cross validation statistics 

% parameters
if nargin < 1, soft = 1; end 
if nargin < 2, forceddist = 0; end
if nargin < 3, loo = 0; end % loo == 0 = 'LOO'; loo == 1 = '10fold'; loo == 2 = 'true10fold';

if soft == 0 
    softstr = '_nosoft'; 
elseif soft == 1 
    softstr = '_soft'; 
elseif soft == 2
    softstr = '';
elseif soft == 3
    softstr = '';
elseif soft == 4
    softstr = '_add_ind_muli_ind';
elseif soft == 5
    softstr = '_add_dyn_muli_ind';
end

constr = '_long';
gausstr = '_gauss';
 
% loading results
if loo == 0
    if soft < 2 % krig, RAMP
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_20150408.mat', ...
            softstr,constr,gausstr,floor(forceddist./1000)));
        zk = cell2mat(zk_madd);
        zh = cell2mat(zh_Xval);
        vk = cell2mat(vk);
        ck = cell2mat(ckXval);
    elseif soft == 2 % CAMP
    elseif soft == 3 % staticDS
        load(sprintf('../12_mfiles_othermethods/matdata/LOOCE_DSforcedisolation_%dkm%s.mat',floor(forceddist./1000),softstr));
        temp1 = cell2mat(cellfun(@length,zk,'UniformOutput',false));
        temp2 = cell2mat(arrayfun(@(x,y) repmat(x,y,1),unidates,temp1,'UniformOutput',false));
        zk = cell2mat(zk);
        zh = cell2mat(zhall); zh = zh(:,1);
        vk = cell2mat(vk);
        ck = cell2mat(ckall);
        ck = [ck temp2];
    elseif soft > 3 % staticDS
        load(sprintf('../12_mfiles_othermethods/matdata/LOOCE_DSforcedisolation_%dkm%s.mat',floor(forceddist./1000),softstr));
        zk = cell2mat(zk);
        zh = cell2mat(zhX);
        vk = cell2mat(vk);
        ck = cell2mat(ck);
    end
elseif loo == 1
    for i = 1:10
        load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s_timezone_20150408.mat', ...
            i,softstr,constr,gausstr));
        zks{i,1} = zk_madd;
        zh{i,1} = zh_Xval;
        vks{i,1} = vk;
        cks{i,1} = ck;
        fold{i,1} = i*ones(length(ck),1);
    end
    zk = cell2mat(zks); zh = cell2mat(zh); vk = cell2mat(vks); ck = cell2mat(cks); fold = cell2mat(fold);
elseif loo == 2
    if soft == 0, softstr = '_hardonly'; elseif soft == 1, softstr = ''; end
    for i = 1:10
        smoothingParam = [900000 300000 100 50];
        load(sprintf('../matfiles/meanTrend_true10fold_%d_%d_%d_%d_%d.mat',i,smoothingParam));
        zh{i,1} = zh_Xval;
        load(sprintf('../matfiles/Xval_true10fold_fold%d%s.mat',i,softstr));
        zks{i,1} = zk_madd;
        vks{i,1} = vk;
        ck{i,1} = ch_Xval;
        fold{i,1} = i*ones(length(ck{i,1}),1);
    end
    zk = cell2mat(zks); zh = cell2mat(zh); vk = cell2mat(vks); ck = cell2mat(ck); fold = cell2mat(fold);
end

% info about region, fold
[PMs,EPAall] = get_info(ck);

% overall
[num,mObs,mMod,mBias,nBias,nmBias,fBias,mErr,nErr,nmErr,fErr,R,R2,sBias,msBias, ...
    rmsBias,nrmsBias,mDsBias,m2DmsBias,s2DmsBias,beta1,vObs,vMod] = ...
    calcstats(zh,ck,zk,vk);
overall = { num,mObs,mMod,mBias,nBias,nmBias,fBias,mErr,nErr,nmErr,fErr,R,R2,sBias,msBias, ...
    rmsBias,nrmsBias,mDsBias,m2DmsBias,s2DmsBias,beta1,vObs,vMod };

% years
dummy = datevec(ck(:,3));
uniyr = unique(dummy(:,1));
numY=NaN*ones(length(uniyr),1);mObsY=NaN*ones(length(uniyr),1);
mModY=NaN*ones(length(uniyr),1);mBiasY=NaN*ones(length(uniyr),1);
nBiasY=NaN*ones(length(uniyr),1);nmBiasY=NaN*ones(length(uniyr),1);
fBiasY=NaN*ones(length(uniyr),1);mErrY=NaN*ones(length(uniyr),1);
nErrY=NaN*ones(length(uniyr),1);nmErrY=NaN*ones(length(uniyr),1);
fErrY=NaN*ones(length(uniyr),1);RY=NaN*ones(length(uniyr),1);
R2Y=NaN*ones(length(uniyr),1);sBiasY=NaN*ones(length(uniyr),1);
msBiasY=NaN*ones(length(uniyr),1);rmsBiasY=NaN*ones(length(uniyr),1);
nrmsBiasY=NaN*ones(length(uniyr),1);mDsBiasY=NaN*ones(length(uniyr),1);
m2DmsBiasY=NaN*ones(length(uniyr),1);s2DmsBiasY=NaN*ones(length(uniyr),1);
beta1Y=NaN*ones(length(uniyr),1);vObsY=NaN*ones(length(uniyr),1);
vModY=NaN*ones(length(uniyr),1);
for i = 1:length(uniyr)
    idx = uniyr(i) == dummy(:,1);
    [numY(i,1),mObsY(i,1),mModY(i,1),mBiasY(i,1),nBiasY(i,1),nmBiasY(i,1), ...
        fBiasY(i,1),mErrY(i,1),nErrY(i,1),nmErrY(i,1),fErrY(i,1),RY(i,1), ...
        R2Y(i,1),sBiasY(i,1),msBiasY(i,1),rmsBiasY(i,1),nrmsBiasY(i,1), ...
        mDsBiasY(i,1),m2DmsBiasY(i,1),s2DmsBiasY(i,1),beta1Y(i,1),vObsY(i,1),vModY(i,1)] = ...
        calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
end
years = { numY,mObsY,mModY,mBiasY,nBiasY,nmBiasY,fBiasY,mErrY,nErrY,nmErrY,fErrY,RY,R2Y,sBiasY,msBiasY, ...
    rmsBiasY,nrmsBiasY,mDsBiasY,m2DmsBiasY,s2DmsBiasY,beta1Y,vObsY,vModY };

% region (overall)
uniR = unique(EPAall);
numR=NaN*ones(length(uniR),1);mObsR=NaN*ones(length(uniR),1);
mModR=NaN*ones(length(uniR),1);mBiasR=NaN*ones(length(uniR),1);
nBiasR=NaN*ones(length(uniR),1);nmBiasR=NaN*ones(length(uniR),1);
fBiasR=NaN*ones(length(uniR),1);mErrR=NaN*ones(length(uniR),1);
nErrR=NaN*ones(length(uniR),1);nmErrR=NaN*ones(length(uniR),1);
fErrR=NaN*ones(length(uniR),1);RR=NaN*ones(length(uniR),1);
R2R=NaN*ones(length(uniR),1);sBiasR=NaN*ones(length(uniR),1);
msBiasR=NaN*ones(length(uniR),1);rmsBiasR=NaN*ones(length(uniR),1);
nrmsBiasR=NaN*ones(length(uniR),1);mDsBiasR=NaN*ones(length(uniR),1);
m2DmsBiasR=NaN*ones(length(uniR),1);s2DmsBiasR=NaN*ones(length(uniR),1);
beta1R=NaN*ones(length(uniR),1);vObsR=NaN*ones(length(uniR),1);
vModR=NaN*ones(length(uniR),1);
for i = 1:length(uniR)
    idx = uniR(i) == EPAall;
    [numR(i,1),mObsR(i,1),mModR(i,1),mBiasR(i,1),nBiasR(i,1),nmBiasR(i,1), ...
        fBiasR(i,1),mErrR(i,1),nErrR(i,1),nmErrR(i,1),fErrR(i,1),RR(i,1), ...
        R2R(i,1),sBiasR(i,1),msBiasR(i,1),rmsBiasR(i,1),nrmsBiasR(i,1), ...
        mDsBiasR(i,1),m2DmsBiasR(i,1),s2DmsBiasR(i,1),beta1R(i,1),vObsR(i,1),vModR(i,1)] = ...
        calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
end
regions = { numR,mObsR,mModR,mBiasR,nBiasR,nmBiasR,fBiasR,mErrR,nErrR,nmErrR,fErrR,RR,R2R,sBiasR,msBiasR, ...
    rmsBiasR,nrmsBiasR,mDsBiasR,m2DmsBiasR,s2DmsBiasR,beta1R,vObsR,vModR };

% region (by year)
numYR=NaN*ones(length(uniyr),length(uniR));mObsYR=NaN*ones(length(uniyr),length(uniR));
mModYR=NaN*ones(length(uniyr),length(uniR));mBiasYR=NaN*ones(length(uniyr),length(uniR));
nBiasYR=NaN*ones(length(uniyr),length(uniR));nmBiasYR=NaN*ones(length(uniyr),length(uniR));
fBiasYR=NaN*ones(length(uniyr),length(uniR));mErrYR=NaN*ones(length(uniyr),length(uniR));
nErrYR=NaN*ones(length(uniyr),length(uniR));nmErrYR=NaN*ones(length(uniyr),length(uniR));
fErrYR=NaN*ones(length(uniyr),length(uniR));RYR=NaN*ones(length(uniyr),length(uniR));
R2YR=NaN*ones(length(uniyr),length(uniR));sBiasYR=NaN*ones(length(uniyr),length(uniR));
msBiasYR=NaN*ones(length(uniyr),length(uniR));rmsBiasYR=NaN*ones(length(uniyr),length(uniR));
nrmsBiasYR=NaN*ones(length(uniyr),length(uniR));mDsBiasYR=NaN*ones(length(uniyr),length(uniR));
m2DmsBiasYR=NaN*ones(length(uniyr),length(uniR));s2DmsBiasYR=NaN*ones(length(uniyr),length(uniR));
beta1YR=NaN*ones(length(uniyr),length(uniR));vObsYR=NaN*ones(length(uniyr),length(uniR));
vModYR=NaN*ones(length(uniyr),length(uniR));
for i = 1:length(uniyr)
    for j = 1:length(uniR)
        idx = uniyr(i) == dummy(:,1) & uniR(j) == EPAall;
        [numYR(i,j),mObsYR(i,j),mModYR(i,j),mBiasYR(i,j),nBiasYR(i,j),nmBiasYR(i,j), ...
            fBiasYR(i,j),mErrYR(i,j),nErrYR(i,j),nmErrYR(i,j),fErrYR(i,j),RYR(i,j), ...
            R2YR(i,j),sBiasYR(i,j),msBiasYR(i,j),rmsBiasYR(i,j),nrmsBiasYR(i,j), ...
            mDsBiasYR(i,j),m2DmsBiasYR(i,j),s2DmsBiasYR(i,j),beta1YR(i,j),vObsYR(i,j),vModYR(i,j)] = ...
            calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
    end
end
yearregions = { numYR,mObsYR,mModYR,mBiasYR,nBiasYR,nmBiasYR,fBiasYR,mErrYR,nErrYR,nmErrYR,fErrYR,RYR,R2YR,sBiasYR,msBiasYR, ...
    rmsBiasYR,nrmsBiasYR,mDsBiasYR,m2DmsBiasYR,s2DmsBiasYR,beta1YR,vObsYR,vModYR };

% by station (overall)
uniMS = unique(ck(:,1:2),'rows');
numMS=NaN*ones(length(uniMS),1);mObsMS=NaN*ones(length(uniMS),1);
mModMS=NaN*ones(length(uniMS),1);mBiasMS=NaN*ones(length(uniMS),1);
nBiasMS=NaN*ones(length(uniMS),1);nmBiasMS=NaN*ones(length(uniMS),1);
fBiasMS=NaN*ones(length(uniMS),1);mErrMS=NaN*ones(length(uniMS),1);
nErrMS=NaN*ones(length(uniMS),1);nmErrMS=NaN*ones(length(uniMS),1);
fErrMS=NaN*ones(length(uniMS),1);RMS=NaN*ones(length(uniMS),1);
R2MS=NaN*ones(length(uniMS),1);sBiasMS=NaN*ones(length(uniMS),1);
msBiasMS=NaN*ones(length(uniMS),1);rmsBiasMS=NaN*ones(length(uniMS),1);
nrmsBiasMS=NaN*ones(length(uniMS),1);mDsBiasMS=NaN*ones(length(uniMS),1);
m2DmsBiasMS=NaN*ones(length(uniMS),1);s2DmsBiasMS=NaN*ones(length(uniMS),1);
beta1MS=NaN*ones(length(uniMS),1);vObsMS=NaN*ones(length(uniMS),1);
vModMS=NaN*ones(length(uniMS),1);
[temp uniidx] = ismember(ck(:,1:2),uniMS,'rows');
for i = 1:length(uniMS)
    idx = uniMS(i,1) == ck(:,1) & uniMS(i,2) == ck(:,2);
    [numMS(i,1),mObsMS(i,1),mModMS(i,1),mBiasMS(i,1),nBiasMS(i,1),nmBiasMS(i,1), ...
        fBiasMS(i,1),mErrMS(i,1),nErrMS(i,1),nmErrMS(i,1),fErrMS(i,1),RMS(i,1), ...
        R2MS(i,1),sBiasMS(i,1),msBiasMS(i,1),rmsBiasMS(i,1),nrmsBiasMS(i,1), ...
        mDsBiasMS(i,1),m2DmsBiasMS(i,1),s2DmsBiasMS(i,1),beta1MS(i,1),vObsMS(i,1),vModMS(i,1)] = ...
        calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
end
stations = { numMS,mObsMS,mModMS,mBiasMS,nBiasMS,nmBiasMS,fBiasMS,mErrMS,nErrMS,nmErrMS,fErrMS,RMS,R2MS,sBiasMS,msBiasMS, ...
    rmsBiasMS,nrmsBiasMS,mDsBiasMS,m2DmsBiasMS,s2DmsBiasMS,beta1MS,vObsMS,vModMS };

% by station (by year)
numYMS=NaN*ones(length(uniyr),length(uniMS));mObsYMS=NaN*ones(length(uniyr),length(uniMS));
mModYMS=NaN*ones(length(uniyr),length(uniMS));mBiasYMS=NaN*ones(length(uniyr),length(uniMS));
nBiasYMS=NaN*ones(length(uniyr),length(uniMS));nmBiasYMS=NaN*ones(length(uniyr),length(uniMS));
fBiasYMS=NaN*ones(length(uniyr),length(uniMS));mErrYMS=NaN*ones(length(uniyr),length(uniMS));
nErrYMS=NaN*ones(length(uniyr),length(uniMS));nmErrYMS=NaN*ones(length(uniyr),length(uniMS));
fErrYMS=NaN*ones(length(uniyr),length(uniMS));RYMS=NaN*ones(length(uniyr),length(uniMS));
R2YMS=NaN*ones(length(uniyr),length(uniMS));sBiasYMS=NaN*ones(length(uniyr),length(uniMS));
msBiasYMS=NaN*ones(length(uniyr),length(uniMS));rmsBiasYMS=NaN*ones(length(uniyr),length(uniMS));
nrmsBiasYMS=NaN*ones(length(uniyr),length(uniMS));mDsBiasYMS=NaN*ones(length(uniyr),length(uniMS));
m2DmsBiasYMS=NaN*ones(length(uniyr),length(uniMS));s2DmsBiasYMS=NaN*ones(length(uniyr),length(uniMS));
beta1YMS=NaN*ones(length(uniyr),length(uniMS));vObsYMS=NaN*ones(length(uniyr),length(uniMS));
vModYMS=NaN*ones(length(uniyr),length(uniMS));
for i = 1:length(uniyr)
    for j = 1:length(uniMS)
        idx = uniyr(i) == dummy(:,1) & uniMS(j,1) == ck(:,1) & uniMS(j,2) == ck(:,2);
        [numYMS(i,j),mObsYMS(i,j),mModYMS(i,j),mBiasYMS(i,j),nBiasYMS(i,j),nmBiasYMS(i,j), ...
            fBiasYMS(i,j),mErrYMS(i,j),nErrYMS(i,j),nmErrYMS(i,j),fErrYMS(i,j),RYMS(i,j), ...
            R2YMS(i,j),sBiasYMS(i,j),msBiasYMS(i,j),rmsBiasYMS(i,j),nrmsBiasYMS(i,j), ...
            mDsBiasYMS(i,j),m2DmsBiasYMS(i,j),s2DmsBiasYMS(i,j),beta1YMS(i,j),vObsYMS(i,j),vModYMS(i,j)] = ...
            calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
    end
end
yearstations = { numYMS,mObsYMS,mModYMS,mBiasYMS,nBiasYMS,nmBiasYMS,fBiasYMS,mErrYMS,nErrYMS,nmErrYMS,fErrYMS,RYMS,R2YMS,sBiasYMS,msBiasYMS, ...
    rmsBiasYMS,nrmsBiasYMS,mDsBiasYMS,m2DmsBiasYMS,s2DmsBiasYMS,beta1YMS,vObsYMS,vModYMS };

% distance to nearest neighbor (overall)
X = uniMS'; Y = uniMS';
D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
sorted = sort(D,2);
nearest = sorted(:,2);
prcnear = prctile(nearest,0:10:100);
nearestall = nearest(uniidx);
numN=NaN*ones(10,1);mObsN=NaN*ones(10,1);
mModN=NaN*ones(10,1);mBiasN=NaN*ones(10,1);
nBiasN=NaN*ones(10,1);nmBiasN=NaN*ones(10,1);
fBiasN=NaN*ones(10,1);mErrN=NaN*ones(10,1);
nErrN=NaN*ones(10,1);nmErrN=NaN*ones(10,1);
fErrN=NaN*ones(10,1);RN=NaN*ones(10,1);
R2N=NaN*ones(10,1);sBiasN=NaN*ones(10,1);
msBiasN=NaN*ones(10,1);rmsBiasN=NaN*ones(10,1);
nrmsBiasN=NaN*ones(10,1);mDsBiasN=NaN*ones(10,1);
m2DmsBiasN=NaN*ones(10,1);s2DmsBiasN=NaN*ones(10,1);
beta1N=NaN*ones(10,1);vObsN=NaN*ones(10,1);
vModN=NaN*ones(10,1);
for i = 1:10
    idx = nearestall >= prcnear(i) & nearestall < prcnear(i+1);
    [numN(i,1),mObsN(i,1),mModN(i,1),mBiasN(i,1),nBiasN(i,1),nmBiasN(i,1), ...
        fBiasN(i,1),mErrN(i,1),nErrN(i,1),nmErrN(i,1),fErrN(i,1),RN(i,1), ...
        R2N(i,1),sBiasN(i,1),msBiasN(i,1),rmsBiasN(i,1),nrmsBiasN(i,1), ...
        mDsBiasN(i,1),m2DmsBiasN(i,1),s2DmsBiasN(i,1),beta1N(i,1),vObsN(i,1),vModN(i,1)] = ...
        calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
end
nears = { numN,mObsN,mModN,mBiasN,nBiasN,nmBiasN,fBiasN,mErrN,nErrN,nmErrN,fErrN,RN,R2N,sBiasN,msBiasN, ...
    rmsBiasN,nrmsBiasN,mDsBiasN,m2DmsBiasN,s2DmsBiasN,beta1N,vObsN,vModN };

% distance to nearest neighbor (by year)
numYN=NaN*ones(length(uniyr),10);mObsYN=NaN*ones(length(uniyr),10);
mModYN=NaN*ones(length(uniyr),10);mBiasYN=NaN*ones(length(uniyr),10);
nBiasYN=NaN*ones(length(uniyr),10);nmBiasYN=NaN*ones(length(uniyr),10);
fBiasYN=NaN*ones(length(uniyr),10);mErrYN=NaN*ones(length(uniyr),10);
nErrYN=NaN*ones(length(uniyr),10);nmErrYN=NaN*ones(length(uniyr),10);
fErrYN=NaN*ones(length(uniyr),10);RYN=NaN*ones(length(uniyr),10);
R2YN=NaN*ones(length(uniyr),10);sBiasYN=NaN*ones(length(uniyr),10);
msBiasYN=NaN*ones(length(uniyr),10);rmsBiasYN=NaN*ones(length(uniyr),10);
nrmsBiasYN=NaN*ones(length(uniyr),10);mDsBiasYN=NaN*ones(length(uniyr),10);
m2DmsBiasYN=NaN*ones(length(uniyr),10);s2DmsBiasYN=NaN*ones(length(uniyr),10);
beta1YN=NaN*ones(length(uniyr),10);vObsYN=NaN*ones(length(uniyr),10);
vModYN=NaN*ones(length(uniyr),10);
for i = 1:length(uniyr)
    for j = 1:10
        idx = uniyr(i) == dummy(:,1) & nearestall >= prcnear(j) & nearestall < prcnear(j+1);
        [numYN(i,j),mObsYN(i,j),mModYN(i,j),mBiasYN(i,j),nBiasYN(i,j),nmBiasYN(i,j), ...
            fBiasYN(i,j),mErrYN(i,j),nErrYN(i,j),nmErrYN(i,j),fErrYN(i,j),RYN(i,j), ...
            R2YN(i,j),sBiasYN(i,j),msBiasYN(i,j),rmsBiasYN(i,j),nrmsBiasYN(i,j), ...
            mDsBiasYN(i,j),m2DmsBiasYN(i,j),s2DmsBiasYN(i,j),beta1YN(i,j),vObsYN(i,j),vModYN(i,j)] = ...
            calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
    end
end
yearnears = { numYN,mObsYN,mModYN,mBiasYN,nBiasYN,nmBiasYN,fBiasYN,mErrYN,nErrYN,nmErrYN,fErrYN,RYN,R2YN,sBiasYN,msBiasYN, ...
    rmsBiasYN,nrmsBiasYN,mDsBiasYN,m2DmsBiasYN,s2DmsBiasYN,beta1YN,vObsYN,vModYN };

% by fold (overall)
if loo == 0 
    folds = NaN;
    unifold = NaN;
elseif loo == 1 | loo == 2
    unifold = unique(fold);
    numF=NaN*ones(length(unifold),1);mObsF=NaN*ones(length(unifold),1);
    mModF=NaN*ones(length(unifold),1);mBiasF=NaN*ones(length(unifold),1);
    nBiasF=NaN*ones(length(unifold),1);nmBiasF=NaN*ones(length(unifold),1);
    fBiasF=NaN*ones(length(unifold),1);mErrF=NaN*ones(length(unifold),1);
    nErrF=NaN*ones(length(unifold),1);nmErrF=NaN*ones(length(unifold),1);
    fErrF=NaN*ones(length(unifold),1);RF=NaN*ones(length(unifold),1);
    R2F=NaN*ones(length(unifold),1);sBiasF=NaN*ones(length(unifold),1);
    msBiasF=NaN*ones(length(unifold),1);rmsBiasF=NaN*ones(length(unifold),1);
    nrmsBiasF=NaN*ones(length(unifold),1);mDsBiasF=NaN*ones(length(unifold),1);
    m2DmsBiasF=NaN*ones(length(unifold),1);s2DmsBiasF=NaN*ones(length(unifold),1);
    beta1F=NaN*ones(length(unifold),1);vObsF=NaN*ones(length(unifold),1);
    vModF=NaN*ones(length(unifold),1);
    for i = 1:length(unifold)
        idx = unifold(i) == fold;
        [numF(i,1),mObsF(i,1),mModF(i,1),mBiasF(i,1),nBiasF(i,1),nmBiasF(i,1), ...
            fBiasF(i,1),mErrF(i,1),nErrF(i,1),nmErrF(i,1),fErrF(i,1),RF(i,1), ...
            R2F(i,1),sBiasF(i,1),msBiasF(i,1),rmsBiasF(i,1),nrmsBiasF(i,1), ...
            mDsBiasF(i,1),m2DmsBiasF(i,1),s2DmsBiasF(i,1),beta1F(i,1),vObsF(i,1),vModF(i,1)] = ...
            calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
    end
    folds = { numF,mObsF,mModF,mBiasF,nBiasF,nmBiasF,fBiasF,mErrF,nErrF,nmErrF,fErrF,RF,R2F,sBiasF,msBiasF, ...
        rmsBiasF,nrmsBiasF,mDsBiasF,m2DmsBiasF,s2DmsBiasF,beta1F,vObsF,vModF };
end

% by fold (by year)
if loo == 0
    yearfolds = NaN;
elseif loo == 1 | loo == 2
    numYF=NaN*ones(length(uniyr),length(unifold));mObsYF=NaN*ones(length(uniyr),length(unifold));
    mModYF=NaN*ones(length(uniyr),length(unifold));mBiasYF=NaN*ones(length(uniyr),length(unifold));
    nBiasYF=NaN*ones(length(uniyr),length(unifold));nmBiasYF=NaN*ones(length(uniyr),length(unifold));
    fBiasYF=NaN*ones(length(uniyr),length(unifold));mErrYF=NaN*ones(length(uniyr),length(unifold));
    nErrYF=NaN*ones(length(uniyr),length(unifold));nmErrYF=NaN*ones(length(uniyr),length(unifold));
    fErrYF=NaN*ones(length(uniyr),length(unifold));RYF=NaN*ones(length(uniyr),length(unifold));
    R2YF=NaN*ones(length(uniyr),length(unifold));sBiasYF=NaN*ones(length(uniyr),length(unifold));
    msBiasYF=NaN*ones(length(uniyr),length(unifold));rmsBiasYF=NaN*ones(length(uniyr),length(unifold));
    nrmsBiasYF=NaN*ones(length(uniyr),length(unifold));mDsBiasYF=NaN*ones(length(uniyr),length(unifold));
    m2DmsBiasYF=NaN*ones(length(uniyr),length(unifold));s2DmsBiasYF=NaN*ones(length(uniyr),length(unifold));
    beta1YF=NaN*ones(length(uniyr),length(unifold));vObsYF=NaN*ones(length(uniyr),length(unifold));
    vModYF=NaN*ones(length(uniyr),length(unifold));
    for i = 1:length(uniyr)
        for j = 1:length(unifold)
            idx = uniyr(i) == dummy(:,1) & unifold(j) == fold;
            [numYF(i,j),mObsYF(i,j),mModYF(i,j),mBiasYF(i,j),nBiasYF(i,j),nmBiasYF(i,j), ...
                fBiasYF(i,j),mErrYF(i,j),nErrYF(i,j),nmErrYF(i,j),fErrYF(i,j),RYF(i,j), ...
                R2YF(i,j),sBiasYF(i,j),msBiasYF(i,j),rmsBiasYF(i,j),nrmsBiasYF(i,j), ...
                mDsBiasYF(i,j),m2DmsBiasYF(i,j),s2DmsBiasYF(i,j),beta1YF(i,j),vObsYF(i,j),vModYF(i,j)] = ...
                calcstats(zh(idx),ck(idx),zk(idx),vk(idx));
        end
    end
    yearfolds = { numYF,mObsYF,mModYF,mBiasYF,nBiasYF,nmBiasYF,fBiasYF,mErrYF,nErrYF,nmErrYF,fErrYF,RYF,R2YF,sBiasYF,msBiasYF, ...
        rmsBiasYF,nrmsBiasYF,mDsBiasYF,m2DmsBiasYF,s2DmsBiasYF,beta1YF,vObsYF,vModYF };
end

end