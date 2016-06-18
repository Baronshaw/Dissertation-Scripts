function [] = stat_MP_cmaqII(yrz)
% this function will calculate measures of CMAQ model performance by
% combining two characteristics

% seperate by: 1) percentile, 2) season, 3) network/land
% network/land = IMPROVE, AQS rural, AQS suburban, AQS urban
%   IMPROVE: PMr == 5, AQS rural: PMr ~= 5 & PMl == 2, AQS suburban: PMr ~=5 &
%       PMl == 3, AQS urban: PMr ~= 5 & PMl == 5

if nargin < 1, yrz = 2001; end

% load CMAQ paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yrz));

% load data information 
load('matfiles/info.mat')

% percentile/season
uniPMP = unique(PMd);
uniPMS = unique(PMs);
for i = 1:length(unique(PMdec))    
    for j = 1:length(unique(PMs))
        idx = PMdec == i & PMs == j;
        numPS(i,j) = size(Obs(idx),1); 
        mObsPS(i,j) = mean(Obs(idx)); 
        mModPS(i,j) = mean(Mod(idx)); 
        mBiasPS(i,j) = mean(Mod(idx)-Obs(idx));
        idxff = Obs > 0; % adding a little fudge factor
        nBiasPS(i,j) = 100.*(1./numPS(i,j)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
        nmBiasPS(i,j) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
        fBiasPS(i,j) = 100.*(1./numPS(i,j)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        mErrPS(i,j) = mean(abs(Mod(idx)-Obs(idx))); 
        nErrPS(i,j) = 100.*(1./numPS(i,j)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
        nmErrPS(i,j) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
        fErrPS(i,j) = 100.*(1./numPS(i,j)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        RPS(i,j) = corr(Mod(idx),Obs(idx));
        R2PS(i,j) = RPS(i,j).^2;
        sBiasPS(i,j) = std(Mod(idx)-Obs(idx));
        msBiasPS(i,j) = mean((Mod(idx)-Obs(idx)).^2);
        rmsBiasPS(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
        nrmsBiasPS(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
        mDsBiasPS(i,j) = mBiasPS(i,j)./sBiasPS(i,j);
        m2DmsBiasPS(i,j) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
        s2DmsBiasPS(i,j) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
        dummy = polyfit(Obs(idx),Mod(idx),1);
        beta1PS(i,j) = dummy(1);
        vObsPS(i,j) = var(Obs(idx));
        vModPS(i,j) = var(Mod(idx));
    end
end

% percentile/network/land
uniPMP = unique(PMd);
for i = 1:length(unique(PMdec))
    for j = 1:4
        if j == 1, idx = PMdec == i & PMr == 5;
        elseif j == 2, idx = PMdec == i & PMr ~= 5 & PMl == 2;
        elseif j == 3, idx = PMdec == i & PMr ~=5 & PMl == 3;
        elseif j == 4, idx = PMdec == i & PMr ~= 5 & PMl == 5; end
        numPN(i,j) = size(Obs(idx),1); 
        mObsPN(i,j) = mean(Obs(idx)); 
        mModPN(i,j) = mean(Mod(idx)); 
        mBiasPN(i,j) = mean(Mod(idx)-Obs(idx));
        idxff = Obs > 0; % adding a little fudge factor
        nBiasPN(i,j) = 100.*(1./numPN(i,j)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
        nmBiasPN(i,j) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
        fBiasPN(i,j) = 100.*(1./numPN(i,j)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        mErrPN(i,j) = mean(abs(Mod(idx)-Obs(idx))); 
        nErrPN(i,j) = 100.*(1./numPN(i,j)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
        nmErrPN(i,j) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
        fErrPN(i,j) = 100.*(1./numPN(i,j)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        RPN(i,j) = corr(Mod(idx),Obs(idx));
        R2PN(i,j) = RPN(i,j).^2;
        sBiasPN(i,j) = std(Mod(idx)-Obs(idx));
        msBiasPN(i,j) = mean((Mod(idx)-Obs(idx)).^2);
        rmsBiasPN(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
        nrmsBiasPN(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
        mDsBiasPN(i,j) = mBiasPN(i,j)./sBiasPN(i,j);
        m2DmsBiasPN(i,j) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
        s2DmsBiasPN(i,j) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
        dummy = polyfit(Obs(idx),Mod(idx),1);
        beta1PN(i,j) = dummy(1);
        vObsPN(i,j) = var(Obs(idx));
        vModPN(i,j) = var(Mod(idx));
    end
end

% season/network/land
uniPMS = unique(PMs);
for i = 1:length(unique(PMs))
    for j = 1:4
        if j == 1, idx = PMs == i & PMr == 5;
        elseif j == 2, idx = PMdec == i & PMr ~= 5 & PMl == 2;
        elseif j == 3, idx = PMdec == i & PMr ~=5 & PMl == 3;
        elseif j == 4, idx = PMdec == i & PMr ~= 5 & PMl == 5; end
        numSN(i,j) = size(Obs(idx),1); 
        mObsSN(i,j) = mean(Obs(idx)); 
        mModSN(i,j) = mean(Mod(idx)); 
        mBiasSN(i,j) = mean(Mod(idx)-Obs(idx));
        idxff = Obs > 0; % adding a little fudge factor
        nBiasSN(i,j) = 100.*(1./numSN(i,j)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
        nmBiasSN(i,j) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
        fBiasSN(i,j) = 100.*(1./numSN(i,j)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        mErrSN(i,j) = mean(abs(Mod(idx)-Obs(idx))); 
        nErrSN(i,j) = 100.*(1./numSN(i,j)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
        nmErrSN(i,j) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
        fErrSN(i,j) = 100.*(1./numSN(i,j)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        RSN(i,j) = corr(Mod(idx),Obs(idx));
        R2SN(i,j) = RSN(i,j).^2;
        sBiasSN(i,j) = std(Mod(idx)-Obs(idx));
        msBiasSN(i,j) = mean((Mod(idx)-Obs(idx)).^2);
        rmsBiasSN(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
        nrmsBiasSN(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
        mDsBiasSN(i,j) = mBiasSN(i,j)./sBiasSN(i,j);
        m2DmsBiasSN(i,j) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
        s2DmsBiasSN(i,j) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
        dummy = polyfit(Obs(idx),Mod(idx),1);
        beta1SN(i,j) = dummy(1);
        vObsSN(i,j) = var(Obs(idx));
        vModSN(i,j) = var(Mod(idx));
    end
end

% by station/season
uniPMM = unique(coordObs,'rows');
for i = 1:length(uniPMM)
    for j = 1:4
        idx = coordObs(:,1) == uniPMM(i,1) & coordObs(:,2) == uniPMM(i,2) & PMs == j;
        numMS(i,j) = size(Obs(idx),1); 
        mObsMS(i,j) = mean(Obs(idx)); 
        mModMS(i,j) = mean(Mod(idx)); 
        mBiasMS(i,j) = mean(Mod(idx)-Obs(idx));
        idxff = Obs > 0; % adding a little fudge factor
        nBiasMS(i,j) = 100.*(1./numMS(i,j)).*sum((Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff)); 
        nmBiasMS(i,j) = sum(Mod(idx)-Obs(idx))./sum(Obs(idx));
        fBiasMS(i,j) = 100.*(1./numMS(i,j)).*sum((Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        mErrMS(i,j) = mean(abs(Mod(idx)-Obs(idx))); 
        nErrMS(i,j) = 100.*(1./numMS(i,j)).*sum(abs(Mod(idx&idxff)-Obs(idx&idxff))./Obs(idx&idxff));
        nmErrMS(i,j) = sum(abs(Mod(idx)-Obs(idx)))./sum(Obs(idx)); 
        fErrMS(i,j) = 100.*(1./numMS(i,j)).*sum(abs(Mod(idx)-Obs(idx))./(0.5.*(Mod(idx)+Obs(idx)))); 
        if sum(idx)>0
            RMS(i,j) = corr(Mod(idx),Obs(idx)); 
        else
            RMS(i,j) = NaN;
        end
        R2MS(i,j) = RMS(i,j).^2;
        sBiasMS(i,j) = std(Mod(idx)-Obs(idx));
        msBiasMS(i,j) = mean((Mod(idx)-Obs(idx)).^2);
        rmsBiasMS(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2));
        nrmsBiasMS(i,j) = sqrt(mean((Mod(idx)-Obs(idx)).^2))./mean(Obs(idx));
        mDsBiasMS(i,j) = mBiasMS(i,j)./sBiasMS(i,j);
        m2DmsBiasMS(i,j) = (mean(Mod(idx)-Obs(idx))).^2./mean((Mod(idx)-Obs(idx)).^2);
        s2DmsBiasMS(i,j) = var(Mod(idx)-Obs(idx))./mean((Mod(idx)-Obs(idx)).^2);
        dummy = polyfit(Obs(idx),Mod(idx),1);
        if sum(idx)>0
            beta1MS(i,j) = dummy(1);
        else
            betaMS(i,j) = NaN;
        end
        vObsMS(i,j) = var(Obs(idx));
        vModMS(i,j) = var(Mod(idx));
    end
end

% save results
save('matfiles/traditional_performanceII.mat', ...
    'numPS','mObsPS','mModPS','mBiasPS','nBiasPS' ,'nmBiasPS','fBiasPS', ...
    'mErrPS','nErrPS','nmErrPS','fErrPS','RPS','R2PS','sBiasPS', ...
    'msBiasPS','rmsBiasPS','nrmsBiasPS','mDsBiasPS','m2DmsBiasPS', ...
    's2DmsBiasPS','beta1PS','vObsPS','vModPS', ...
    'numPN','mObsPN','mModPN','mBiasPN','nBiasPN' ,'nmBiasPN','fBiasPN', ...
    'mErrPN','nErrPN','nmErrPN','fErrPN','RPN','R2PN','sBiasPN', ...
    'msBiasPN','rmsBiasPN','nrmsBiasPN','mDsBiasPN','m2DmsBiasPN', ...
    's2DmsBiasPN','beta1PN','vObsPN','vModPN', ...
    'numSN','mObsSN','mModSN','mBiasSN','nBiasSN' ,'nmBiasSN','fBiasSN', ...
    'mErrSN','nErrSN','nmErrSN','fErrSN','RSN','R2SN','sBiasSN', ...
    'msBiasSN','rmsBiasSN','nrmsBiasSN','mDsBiasSN','m2DmsBiasSN', ...
    's2DmsBiasSN','beta1SN','vObsSN','vModSN', ...
    'numMS','mObsMS','mModMS','mBiasMS','nBiasMS' ,'nmBiasMS','fBiasMS', ...
    'mErrMS','nErrMS','nmErrMS','fErrMS','RMS','R2MS','sBiasMS', ...
    'msBiasMS','rmsBiasMS','nrmsBiasMS','mDsBiasMS','m2DmsBiasMS', ...
    's2DmsBiasMS','beta1MS','vObsMS','vModMS'); 

end