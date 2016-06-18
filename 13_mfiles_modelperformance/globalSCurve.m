function [] = globalSCurve(yrz,numbins,modplots,negval)
% this function will create a global SCurve for 2001 (de Nazelle method)

% parameters subject to change
if nargin < 1, yrz = 2001; end
if nargin < 2, numbins = 10; end % number of bins for each plot
if nargin < 3, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 4, negval = 0; end % 0 = there are no negative predicted values

% load CMAQ paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yrz));

% load CMAQ modeled data
load(sprintf('../matfiles/prepCTM_%d.mat',yrz))

% perciles for each grid 
perctile_data = prctile(Mod,linspace(0,100,numbins+1))';

% find the mean and variance of each bin
for j = 1:numbins
    disp(j);
    idx = Mod >= perctile_data(j) & Mod < perctile_data(j+1);    
    
    num(j,1) = nansum(idx);
    mObs(j,1) = nanmean(Obs(idx)); 
    mMod(j,1) = nanmean(Mod(idx)); 
    mBias(j,1) = nanmean(Mod(idx)-Obs(idx));
    idxff = Obs(idx) <= 0; Obs2 = Obs(idx); Mod2  = Mod(idx); Obs2(idxff) = NaN; Mod2(idxff) = NaN; % adding fudge factor
    nBias(j,1) = 100.*(1./num(j,1)).*nansum((Mod2-Obs2)./Obs2); 
    nmBias(j,1) = nansum(Mod(idx)-Obs(idx))./nansum(Obs(idx));
    fBias(j,1) = 100.*(1./num(j,1)).*nansum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    mErr(j,1) = nanmean(abs(Mod(idx)-Obs(idx))); 
    nErr(j,1) = 100.*(1./num(j,1)).*nansum(abs(Mod2-Obs2)./Obs2); 
    nmErr(j,1) = nansum(abs(Mod(idx)-Obs(idx)))./nansum(Obs(idx)); 
    fErr(j,1) = 100.*(1./num(j,1)).*nansum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
    R(j,1) = corr(Mod(idx),Obs(idx));    
    R2(j,1) = R(j,1).^2;        
    sBias(j,1) = nanstd(Mod(idx)-Obs(idx));
    msBias(j,1) = nanmean((Mod(idx)-Obs(idx)).^2);
    rmsBias(j,1) = sqrt(nanmean((Mod(idx)-Obs(idx)).^2));
    nrmsBias(j,1) = sqrt(nanmean((Mod(idx)-Obs(idx)).^2))./nanmean(Obs(idx));
    mDsBias(j,1) = mBias(j,1)./sBias(j,1);
    m2DmsBias(j,1) = (mBias(j,1)).^2./msBias(j,1); 
    s2DmsBias(j,1) = (sBias(j,1)).^2./msBias(j,1);
    dummy = polyfit(Obs(idx),Mod(idx),1);
    beta1(j,1) = dummy(1);
    vObs(j,1) = nanvar(Obs(idx));
    vMod(j,1) = nanvar(Mod(idx));
       
end 

% fixed modeled values
numGivMod = interp1(mMod,num,modplots,'linear','extrap');
mObsGivMod = interp1(mMod,mObs,modplots,'linear','extrap');
mModGivMod = interp1(mMod,mMod,modplots,'linear','extrap');
mBiasGivMod = interp1(mMod,mBias,modplots,'linear','extrap');
nBiasGivMod = interp1(mMod,nBias,modplots,'linear','extrap');
nmBiasGivMod = interp1(mMod,nmBias,modplots,'linear','extrap');
fBiasGivMod = interp1(mMod,fBias,modplots,'linear','extrap');
mErrGivMod = interp1(mMod,mErr,modplots,'linear','extrap');
nErrGivMod = interp1(mMod,nErr,modplots,'linear','extrap');
nmErrGivMod = interp1(mMod,nmErr,modplots,'linear','extrap'); 
fErrGivMod = interp1(mMod,fErr,modplots,'linear','extrap');
RGivMod = interp1(mMod,R,modplots,'linear','extrap');
R2GivMod = interp1(mMod,R2,modplots,'linear','extrap'); 
sBiasGivMod = interp1(mMod,sBias,modplots,'linear','extrap');
msBiasGivMod = interp1(mMod,msBias,modplots,'linear','extrap');
rmsBiasGivMod = interp1(mMod,rmsBias,modplots,'linear','extrap'); 
nrmsBiasGivMod = interp1(mMod,nrmsBias,modplots,'linear','extrap'); 
mDsBiasGivMod = interp1(mMod,mDsBias,modplots,'linear','extrap');
m2DmsBiasGivMod = interp1(mMod,m2DmsBias,modplots,'linear','extrap');
s2DmsBiasGivMod = interp1(mMod,s2DmsBias,modplots,'linear','extrap');
beta1GivMod = interp1(mMod,beta1,modplots,'linear','extrap');
vObsGivMod = interp1(mMod,vObs,modplots,'linear','extrap');
vModGivMod = interp1(mMod,vMod,modplots,'linear','extrap');

% if values are negative, make them zero, but show S-curves
if negval == 0
    mMod(mMod<0) = 0;
    mModGivMod(mModGivMod<0) = 0;
end

% save results
save('matfiles/globalSCurve.mat','dailyCTMv','distCTMv','yrmodaCTMv', ...
    'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
    'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
    'numGivMod','mObsGivMod','mModGivMod','mBiasGivMod','nBiasGivMod','nmBiasGivMod', ...
    'fBiasGivMod','mErrGivMod','nErrGivMod','nmErrGivMod','fErrGivMod','RGivMod', ... 
    'R2GivMod','sBiasGivMod','msBiasGivMod','rmsBiasGivMod','nrmsBiasGivMod','mDsBiasGivMod', ...
    'm2DmsBiasGivMod','s2DmsBiasGivMod','beta1GivMod','vObsGivMod','vModGivMod', ...
    'yrz','numbins','modplots','negval');

end