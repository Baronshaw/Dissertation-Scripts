function [] = stat_MP_cmaq_grid(yrz)
% this function estimates the measures of traditional model performance on
% a grid

if nargin < 1, yrz = 2001; end

% loop through each day of the year
yrNday = datevec(datenum(yrz,1,1):datenum(yrz,12,31));
yrNday = yrNday(:,1).*10^4 + yrNday(:,2).*10^2 + yrNday(:,3);

% variables 
load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(1)));
r = length(CTMlocs); c = length(yrNday);
num = NaN*ones(r,c); mObs = NaN*ones(r,c); mMod = NaN*ones(r,c);
mBias = NaN*ones(r,c); nBias = NaN*ones(r,c); nmBias = NaN*ones(r,c);
fBias = NaN*ones(r,c); mErr = NaN*ones(r,c); nErr = NaN*ones(r,c);
nmErr = NaN*ones(r,c); fErr = NaN*ones(r,c); R = NaN*ones(r,c);
R2 = NaN*ones(r,c); sBias = NaN*ones(r,c); msBias = NaN*ones(r,c);
rmsBias = NaN*ones(r,c); nrmsBias = NaN*ones(r,c); mDsBias = NaN*ones(r,c);
m2DmsBias = NaN*ones(r,c); s2DmsBias = NaN*ones(r,c); beta1 = NaN*ones(r,c);
vObs = NaN*ones(r,c); vMod = NaN*ones(r,c);

% load soft data files
for i = 1:length(yrNday)
    disp(i);
    load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));
    Mod = valMod;
    Obs = valObs;
    idx = ~isnan(Mod);
    
    num(:,i) = sum(idx,2);
    mObs(:,i) = nanmean(Obs,2); 
    mMod(:,i) = nanmean(Mod,2); 
    mBias(:,i) = nanmean(Mod-Obs,2);
    idxff = Obs <= 0; Obs2 = Obs; Mod2  = Mod; Obs2(idxff) = NaN; Mod2(idxff) = NaN; % adding fudge factor
    nBias(:,i) = 100.*(1./num(:,i)).*nansum((Mod2-Obs2)./Obs2,2); 
    nmBias(:,i) = nansum(Mod-Obs,2)./nansum(Obs,2);
    fBias(:,i) = 100.*(1./num(:,i)).*nansum((Mod-Obs)./(0.5.*(Mod+Obs)),2); 
    mErr(:,i) = nanmean(abs(Mod-Obs),2); 
    nErr(:,i) = 100.*(1./num(:,i)).*nansum(abs(Mod2-Obs2)./Obs2,2); 
    nmErr(:,i) = nansum(abs(Mod-Obs),2)./nansum(Obs,2); 
    fErr(:,i) = 100.*(1./num(:,i)).*nansum(abs(Mod-Obs)./(0.5.*(Mod+Obs)),2); 

    idxj = ~isnan(Mod);
    for j = 1:r
        R(j,i) = corr(Mod(j,idxj(j,:))',Obs(j,idxj(j,:))');
    end
    
    R2(:,i) = R(:,i).^2;
        
    sBias(:,i) = nanstd(Mod-Obs,[],2);
    msBias(:,i) = nanmean((Mod-Obs).^2,2);
    rmsBias(:,i) = sqrt(nanmean((Mod-Obs).^2,2));
    nrmsBias(:,i) = sqrt(nanmean((Mod-Obs).^2,2))./nanmean(Obs,2);
    mDsBias(:,i) = mBias(:,i)./sBias(:,i);
    m2DmsBias(:,i) = (mBias(:,i)).^2./msBias(:,i); 
    s2DmsBias(:,i) = (sBias(:,i)).^2./msBias(:,i);
    dummy = arrayfun(@(x) polyfit(Obs(x,~isnan(Obs(x,:)))',Mod(x,~isnan(Obs(x,:)))',1), 1:r, 'UniformOutput', false)';
    dummy2 = cell2mat(dummy);
    beta1(:,i) = dummy2(:,1);
    vObs(:,i) = nanvar(Obs,[],2);
    vMod(:,i) = nanvar(Mod,[],2);
    
end

% save results
save('matfiles/traditional_performance_grid.mat', ...
    'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
    'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
    'CTMlocs','yrNday'); 

end