function [] = stat_MP_cmaq(yrz)
% this function will calculate measures of CMAQ model performance

% Traditional Measures of Model Performance Include:
% from http://www.epa.gov/cleanairactbenefits/jan10/812_AQM_mpe_memo.pdf
% 1) number of paired modeled and obs
% 2) mean obs value
% 3) mean modeled value
% 4) mean bias
% 5) normalized bias
% 6) normalized mean bias
% 7) fractional bias
% 8) mean error
% 9) normalized error
% 10) normalized mean error
% 11) fractional error 
% 12) correlation
% proposed new measures
% 13) correlation squared
% 14) standard bias
% 15) mean squared bias
% 16) root mean squared bias
% 17) normalized root mean squared bias
% 18) mean bias/standard bias
% 19) mean bias squared/mean squared bias
% 20) variance of bias/mean squared bias
% other measures
% 21) beta1
% 22) variance of obs
% 23) variance of mod

% seperate by: 1) percentile, 2) season, 3) region, 4) network, 5) land

if nargin < 1, yrz = 2001; end

% load CMAQ paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yrz));

% load data information 
load('matfiles/info.mat')

% overall
num = size(Obs,1); 
mObs = mean(Obs); 
mMod = mean(Mod); 
mBias = mean(Mod-Obs);
idxff = Obs > 0; % adding a little fudge factor
nBias = 100.*(1./num).*sum((Mod(idxff)-Obs(idxff))./Obs(idxff)); 
nmBias = sum(Mod-Obs)./sum(Obs);
fBias = 100.*(1./num).*sum((Mod-Obs)./(0.5.*(Mod+Obs))); 
mErr = mean(abs(Mod-Obs)); 
nErr = 100.*(1./num).*sum(abs(Mod(idxff)-Obs(idxff))./Obs(idxff));
nmErr = sum(abs(Mod-Obs))./sum(Obs); 
fErr = 100.*(1./num).*sum(abs(Mod-Obs)./(0.5.*(Mod+Obs))); 
R = corr(Mod,Obs); 
R2 = R.^2;
sBias = std(Mod-Obs);
msBias = mean((Mod-Obs).^2);
rmsBias = sqrt(mean((Mod-Obs).^2));
nrmsBias = sqrt(mean((Mod-Obs).^2))./mean(Obs);
mDsBias = mBias./sBias;
m2DmsBias = (mean(Mod-Obs)).^2./mean((Mod-Obs).^2);
s2DmsBias = var(Mod-Obs)./mean((Mod-Obs).^2);
dummy = polyfit(Obs,Mod,1);
beta1 = dummy(1);
vObs = var(Obs);
vMod = var(Mod);

% by percentile
uniPMP = unique(PMd);
for i = 1:length(unique(PMdec))
    idx = PMdec == i;
    numP(i) = size(Obs(idx),1); 
    mObsP(i) = mean(Obs(idx)); 
    mModP(i) = mean(Mod(idx)); 
    mBiasP(i) = mean(Mod(idx)-Obs(idx));
    idxff = Obs > 0; % adding a little fudge factor
    nBiasP(i) = 100.*(1./numP(i)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
    nmBiasP(i) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
    fBiasP(i) = 100.*(1./numP(i)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErrP(i) = mean(abs(Mod(idx)-Obs(idx))); 
    nErrP(i) = 100.*(1./numP(i)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
    nmErrP(i) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
    fErrP(i) = 100.*(1./numP(i)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    RP(i) = corr(Mod(idx),Obs(idx));
    R2P(i) = RP(i).^2;
    sBiasP(i) = std(Mod(idx)-Obs(idx));
    msBiasP(i) = mean((Mod(idx)-Obs(idx)).^2);
    rmsBiasP(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
    nrmsBiasP(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
    mDsBiasP(i) = mBiasP(i)./sBiasP(i);
    m2DmsBiasP(i) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
    s2DmsBiasP(i) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1P(i) = dummy(1);
    vObsP(i) = var(Obs(idx));
    vModP(i) = var(Mod(idx));
end

% by season
% season: 1 = winter, 2 = spring, 3 = summer, 4 = fall
uniPMS = unique(PMs);
for i = 1:length(unique(PMs))
    idx = PMs == i;
    numS(i) = size(Obs(idx),1); 
    mObsS(i) = mean(Obs(idx)); 
    mModS(i) = mean(Mod(idx)); 
    mBiasS(i) = mean(Mod(idx)-Obs(idx));
    idxff = Obs > 0; % adding a little fudge factor
    nBiasS(i) = 100.*(1./numS(i)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); % fix obs = Inf
    nmBiasS(i) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
    fBiasS(i) = 100.*(1./numS(i)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErrS(i) = mean(abs(Mod(idx)-Obs(idx))); 
    nErrS(i) = 100.*(1./numS(i)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
    nmErrS(i) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
    fErrS(i) = 100.*(1./numS(i)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    RS(i) = corr(Mod(idx),Obs(idx)); 
    R2S(i) = RS(i).^2;
    sBiasS(i) = std(Mod(idx)-Obs(idx));
    msBiasS(i) = mean((Mod(idx)-Obs(idx)).^2);
    rmsBiasS(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
    nrmsBiasS(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
    mDsBiasS(i) = mBiasS(i)./sBiasS(i);
    m2DmsBiasS(i) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
    s2DmsBiasS(i) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1S(i) = dummy(1);
    vObsS(i) = var(Obs(idx));
    vModS(i) = var(Mod(idx));
end

% by rural/urban/suburban
% rural = 2, suburban = 3, urban = 5
for i = [2 3 5]
    idx = PMl == i;
    numL(i) = size(Obs(idx),1); 
    mObsL(i) = mean(Obs(idx)); 
    mModL(i) = mean(Mod(idx)); 
    mBiasL(i) = mean(Mod(idx)-Obs(idx));
    idxff = Obs > 0; % adding a little fudge factor
    nBiasL(i) = 100.*(1./numL(i)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
    nmBiasL(i) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
    fBiasL(i) = 100.*(1./numL(i)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErrL(i) = mean(abs(Mod(idx)-Obs(idx))); 
    nErrL(i) = 100.*(1./numL(i)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
    nmErrL(i) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
    fErrL(i) = 100.*(1./numL(i)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    RL(i) = corr(Mod(idx),Obs(idx)); 
    R2L(i) = RL(i).^2;
    sBiasL(i) = std(Mod(idx)-Obs(idx));
    msBiasL(i) = mean((Mod(idx)-Obs(idx)).^2);
    rmsBiasL(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
    nrmsBiasL(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
    mDsBiasL(i) = mBiasL(i)./sBiasL(i);
    m2DmsBiasL(i) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
    s2DmsBiasL(i) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1L(i) = dummy(1);
    vObsL(i) = var(Obs(idx));
    vModL(i) = var(Mod(idx));
end

% by network
% 'Rank' 1-4,6 = AQS, 'Rank' 5 = IMPROVE
for i = 1:2
    if i == 1, idx = PMr ~= 5; else idx = PMr == 5; end;
    numR(i) = size(Obs(idx),1); 
    mObsR(i) = mean(Obs(idx)); 
    mModR(i) = mean(Mod(idx)); 
    mBiasR(i) = mean(Mod(idx)-Obs(idx));
    idxff = Obs > 0; % adding a little fudge factor
    nBiasR(i) = 100.*(1./numR(i)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
    nmBiasR(i) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
    fBiasR(i) = 100.*(1./numR(i)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErrR(i) = mean(abs(Mod(idx)-Obs(idx))); 
    nErrR(i) = 100.*(1./numR(i)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
    nmErrR(i) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
    fErrR(i) = 100.*(1./numR(i)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    RR(i) = corr(Mod(idx),Obs(idx));
    R2R(i) = RR(i).^2;
    sBiasR(i) = std(Mod(idx)-Obs(idx));
    msBiasR(i) = mean((Mod(idx)-Obs(idx)).^2);
    rmsBiasR(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
    nrmsBiasR(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
    mDsBiasR(i) = mBiasR(i)./sBiasR(i);
    m2DmsBiasR(i) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
    s2DmsBiasR(i) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1R(i) = dummy(1);
    vObsR(i) = var(Obs(idx));
    vModR(i) = var(Mod(idx));
end

% by EPA region
uniPME = unique(EPAall);
for i = 1:length(unique(EPAall))
    idx = EPAall == i;
    numE(i) = size(Obs(idx),1); 
    mObsE(i) = mean(Obs(idx)); 
    mModE(i) = mean(Mod(idx)); 
    mBiasE(i) = mean(Mod(idx)-Obs(idx));
    idxff = Obs > 0; % adding a little fudge factor
    nBiasE(i) = 100.*(1./numE(i)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
    nmBiasE(i) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
    fBiasE(i) = 100.*(1./numE(i)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErrE(i) = mean(abs(Mod(idx)-Obs(idx))); 
    nErrE(i) = 100.*(1./numE(i)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
    nmErrE(i) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
    fErrE(i) = 100.*(1./numE(i)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    RE(i) = corr(Mod(idx),Obs(idx));  
    R2E(i) = RE(i)^2;
    sBiasE(i) = std(Mod(idx)-Obs(idx));
    msBiasE(i) = mean((Mod(idx)-Obs(idx)).^2);
    rmsBiasE(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
    nrmsBiasE(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
    mDsBiasE(i) = mBiasE(i)./sBiasE(i);
    m2DmsBiasE(i) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
    s2DmsBiasE(i) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1E(i) = dummy(1);
    vObsE(i) = var(Obs(idx));
    vModE(i) = var(Mod(idx));
end

% by station
uniPMM = unique(coordObs,'rows');
for i = 1:length(uniPMM)
    idx = coordObs(:,1) == uniPMM(i,1) & coordObs(:,2) == uniPMM(i,2);
    numM(i) = size(Obs(idx),1); 
    mObsM(i) = mean(Obs(idx)); 
    mModM(i) = mean(Mod(idx)); 
    mBiasM(i) = mean(Mod(idx)-Obs(idx));
    idxff = Obs > 0; % adding a little fudge factor
    nBiasM(i) = 100.*(1./numM(i)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
    nmBiasM(i) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
    fBiasM(i) = 100.*(1./numM(i)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErrM(i) = mean(abs(Mod(idx)-Obs(idx))); 
    nErrM(i) = 100.*(1./numM(i)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
    nmErrM(i) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
    fErrM(i) = 100.*(1./numM(i)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    RM(i) = corr(Mod(idx),Obs(idx)); 
    R2M(i) = RM(i).^2;
    sBiasM(i) = std(Mod(idx)-Obs(idx));
    msBiasM(i) = mean((Mod(idx)-Obs(idx)).^2);
    rmsBiasM(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
    nrmsBiasM(i) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
    mDsBiasM(i) = mBiasM(i)./sBiasM(i);
    m2DmsBiasM(i) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
    s2DmsBiasM(i) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1M(i) = dummy(1);
    vObsM(i) = var(Obs(idx));
    vModM(i) = var(Mod(idx));
end

% by day
dayz = datenum(yr,mo,da);
unidayz = unique(dayz);
for i = 1:length(unidayz)
    idx = dayz == unidayz(i);
    numD(i) = sum(idx);
    mObsD(i) = mean(Obs(idx));
    mModD(i) = mean(Mod(idx));
    BiasD(i) = mModD(i) - mObsD(i);
    ErrD(i) = abs(mModD(i) - mObsD(i));
end

% save results
save('matfiles/traditional_performance.mat', ...
    'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
    'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
    'uniPMP','numP','mObsP','mModP','mBiasP','nBiasP','nmBiasP','fBiasP','mErrP','nErrP','nmErrP','fErrP','RP', ...   
    'R2P','sBiasP','msBiasP','rmsBiasP','nrmsBiasP','mDsBiasP','m2DmsBiasP','s2DmsBiasP','beta1P','vObsP','vModP', ...
    'uniPMS' ,'numS','mObsS' ,'mModS','mBiasS','nBiasS','nmBiasS','fBiasS' ,'mErrS','nErrS','nmErrS','fErrS','RS', ...  
    'R2S','sBiasS','msBiasS','rmsBiasS','nrmsBiasS','mDsBiasS','m2DmsBiasS','s2DmsBiasS','beta1S','vObsS','vModS', ...
    'numL','mObsL' ,'mModL','mBiasL','nBiasL','nmBiasL','fBiasL','mErrL' ,'nErrL','nmErrL','fErrL','RL', ...   
    'R2L','sBiasL','msBiasL','rmsBiasL','nrmsBiasL','mDsBiasL','m2DmsBiasL','s2DmsBiasL','beta1L','vObsL','vModL', ...
    'numR','mObsR','mModR' ,'mBiasR','nBiasR' ,'nmBiasR','fBiasR' ,'mErrR' ,'nErrR','nmErrR','fErrR','RR', ...    
    'R2R','sBiasR','msBiasR','rmsBiasR','nrmsBiasR','mDsBiasR','m2DmsBiasR','s2DmsBiasR','beta1R','vObsR','vModR', ...
    'uniPME','numE','mObsE','mModE','mBiasE','nBiasE','nmBiasE','fBiasE' ,'mErrE','nErrE','nmErrE','fErrE','RE', ...   
    'R2E','sBiasE','msBiasE','rmsBiasE','nrmsBiasE','mDsBiasE','m2DmsBiasE','s2DmsBiasE','beta1E','vObsE','vModE', ...
    'uniPMM','numM','mObsM','mModM','mBiasM','nBiasM','nmBiasM','fBiasM' ,'mErrM','nErrM','nmErrM','fErrM','RM', ...
    'R2M','sBiasM','msBiasM','rmsBiasM','nrmsBiasM','mDsBiasM','m2DmsBiasM','s2DmsBiasM','beta1M','vObsM','vModM', ...
    'unidayz','numD','mObsD','mModD','BiasD','ErrD'); 

end