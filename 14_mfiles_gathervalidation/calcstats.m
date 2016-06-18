function [num,mObs,mMod,mBias,nBias,nmBias,fBias, ...
    mErr,nErr,nmErr,fErr,R,R2,sBias,msBias, ...
    rmsBias,nrmsBias,mDsBias,m2DmsBias,s2DmsBias,beta1,vObs,vMod] = ...
    calcstats(zh,ck,zk,vk)
% calculcate all the crossvalidation statistics

% initialize parameters
if nargin < 1, zh = rand(100,1); end
if nargin < 2, ck = rand(100,3); end
if nargin < 3, zk = rand(100,1); end
if nargin < 4, vk = rand(100,1); end

idx = ~isnan(zh) & ~isnan(zk) & zk < 640;
Obs = zh(idx); Mod = zk(idx);
num = size(Obs,1); 
mObs = mean(Obs); 
mMod = mean(Mod); 
mBias = mean(Mod-Obs);
idxff = Obs > 0.000001; % adding a little fudge factor
nBias = 100.*(1./num).*sum((Mod(idxff)-Obs(idxff))./Obs(idxff)); 
nmBias = sum(Mod-Obs)./sum(Obs);
fBias = 100.*(1./num).*sum((Mod-Obs)./(0.5.*(Mod+Obs))); 
mErr = mean(abs(Mod-Obs)); 
nErr = 100.*(1./num).*sum(abs(Mod(idxff)-Obs(idxff))./Obs(idxff));
nmErr = sum(abs(Mod-Obs))./sum(Obs); 
fErr = 100.*(1./num).*sum(abs(Mod-Obs)./(0.5.*(Mod+Obs))); 
if isempty(Mod), R = NaN; else R = corr(Mod,Obs); end
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

end