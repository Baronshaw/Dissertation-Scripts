function [] = stat_MP_cmaq_gridbin(yrz)
% this function estimates the measures of traditional model performance on
% a grid within each bin of an SCurve

if nargin < 1, yrz = 2001; end

% loop through each day of the year
yrNday = datevec(datenum(yrz,1,1):datenum(yrz,12,31));
yrNday = yrNday(:,1).*10^4 + yrNday(:,2).*10^2 + yrNday(:,3);
len = length(yrNday);

% variables 
load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(1)));
r = length(CTMlocs); c = length(yrNday);
[rp cp] = size(perctile_data);
num = cell(len,1); mObs = cell(len,1); mMod = cell(len,1);
mBias = cell(len,1); nBias = cell(len,1); nmBias = cell(len,1);
fBias = cell(len,1); mErr = cell(len,1); nErr = cell(len,1);
nmErr = cell(len,1); fErr = cell(len,1); R = cell(len,1);
R2 = cell(len,1); sBias = cell(len,1); msBias = cell(len,1);
rmsBias = cell(len,1); nrmsBias = cell(len,1); mDsBias = cell(len,1);
m2DmsBias = cell(len,1); s2DmsBias = cell(len,1); beta1 = cell(len,1);
vObs = cell(len,1); vMod = cell(len,1);

% load soft data files
for i = 1:len
    tic
    disp(i);
    temp = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));
    Modtemp = temp.valMod;
    [rm cm] = size(Modtemp);
    
    idx = arrayfun(@(x,y) Modtemp>repmat(perctile_data(:,x),1,cm) & ...
        Modtemp<repmat(perctile_data(:,y),1,cm),1:cp-1,2:cp,'UniformOutput',false); 
    
    num{i} = NaN*ones(r,cp-1); mObs{i} = NaN*ones(r,cp-1); mMod{i} = NaN*ones(r,cp-1);
    mBias{i} = NaN*ones(r,cp-1); nBias{i} = NaN*ones(r,cp-1); nmBias{i} = NaN*ones(r,cp-1);
    fBias{i} = NaN*ones(r,cp-1); mErr{i} = NaN*ones(r,cp-1); nErr{i} = NaN*ones(r,cp-1);
    nmErr{i} = NaN*ones(r,cp-1); fErr{i} = NaN*ones(r,cp-1); R{i} = NaN*ones(r,cp-1);
    R2{i} = NaN*ones(r,cp-1); sBias{i} = NaN*ones(r,cp-1); msBias{i} = NaN*ones(r,cp-1);
    rmsBias{i} = NaN*ones(r,cp-1); nrmsBias{i} = NaN*ones(r,cp-1); mDsBias{i} = NaN*ones(r,cp-1);
    m2DmsBias{i} = NaN*ones(r,cp-1); s2DmsBias{i} = NaN*ones(r,cp-1); beta1{i} = NaN*ones(r,cp-1);
    vObs{i} = NaN*ones(r,cp-1); vMod{i} = NaN*ones(r,cp-1);

    for j = 1:cp-1
        disp([i j]);
        Mod = temp.valMod; Mod(~idx{j}) = NaN;
        Obs = temp.valObs; Obs(~idx{j}) = NaN;

        num{i}(:,j) = sum(idx{j},2);
        mObs{i}(:,j) = nanmean(Obs,2); 
        mMod{i}(:,j) = nanmean(Mod,2); 
        mBias{i}(:,j) = nanmean(Mod-Obs,2);
        idxff = Obs <= 0; Obs2 = Obs; Mod2  = Mod; Obs2(idxff) = NaN; Mod2(idxff) = NaN; % adding fudge factor
        nBias{i}(:,j) = 100.*(1./num{i}(:,j)).*nansum((Mod2-Obs2)./Obs2,2); 
        nmBias{i}(:,j) = nansum(Mod-Obs,2)./nansum(Obs,2);
        fBias{i}(:,j) = 100.*(1./num{i}(:,j)).*nansum((Mod-Obs)./(0.5.*(Mod+Obs)),2); 
        mErr{i}(:,j) = nanmean(abs(Mod-Obs),2); 
        nErr{i}(:,j) = 100.*(1./num{i}(:,j)).*nansum(abs(Mod2-Obs2)./Obs2,2); 
        nmErr{i}(:,j) = nansum(abs(Mod-Obs),2)./nansum(Obs,2); 
        fErr{i}(:,j) = 100.*(1./num{i}(:,j)).*nansum(abs(Mod-Obs)./(0.5.*(Mod+Obs)),2); 

        idxk = ~isnan(Mod);
        for k = 1:r
            if sum(idxk(k,:)) > 0
                R{i}(k,j) = corr(Mod(k,idxk(k,:))',Obs(k,idxk(k,:))');
            end
        end
        
        R2{i}(:,j) = R{i}(:,j).^2;

        sBias{i}(:,j) = nanstd(Mod-Obs,[],2);
        msBias{i}(:,j) = nanmean((Mod-Obs).^2,2);
        rmsBias{i}(:,j) = sqrt(nanmean((Mod-Obs).^2,2));
        nrmsBias{i}(:,j) = sqrt(nanmean((Mod-Obs).^2,2))./nanmean(Obs,2);
        mDsBias{i}(:,j) = mBias{i}(:,j)./sBias{i}(:,j);
        m2DmsBias{i}(:,j) = (mBias{i}(:,j)).^2./msBias{i}(:,j); 
        s2DmsBias{i}(:,j) = (sBias{i}(:,j)).^2./msBias{i}(:,j);
        dummy = arrayfun(@(x) polyfit(Obs(x,~isnan(Obs(x,:)))',Mod(x,~isnan(Obs(x,:)))',1), 1:r, 'UniformOutput', false)';
        dummy2 = cell2mat(dummy);
        beta1{i}(:,j) = dummy2(:,1);
        vObs{i}(:,j) = nanvar(Obs,[],2);
        vMod{i}(:,j) = var(Mod,[],2);
    
    end
    toc
end

% save results
save('matfiles/traditional_performance_gridbin.mat', ...
    'num','mObs','mMod','mBias','nBias','nmBias','fBias','mErr','nErr','nmErr','fErr','R', ... 
    'R2','sBias','msBias','rmsBias','nrmsBias','mDsBias','m2DmsBias','s2DmsBias','beta1','vObs','vMod', ...
    'CTMlocs','yrNday'); 

end