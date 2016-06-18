function [] = runCAMP(yr)
% this function will impliment the CAMP method

if nargin < 1, yr = 2001; end

modplots = 0:5:50;

% gather all modeled and obs for 2001
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yr));

% 10 percentile bins
prctbins = prctile(Mod,0:10:100);

% mean/variance in each bin
for i = 1:length(prctbins)-1
    mean_Mod(i,1) = mean(Mod(Mod>=prctbins(i)&Mod<prctbins(i+1)));
    mean_Obs(i,1) = mean(Obs(Mod>=prctbins(i)&Mod<prctbins(i+1)));
    var_Obs(i,1) = var(Obs(Mod>=prctbins(i)&Mod<prctbins(i+1)));
end

% load all modeled data for 2001
load(sprintf('../matfiles/prepCTM_%d.mat',yr));

% calculate lambda1 and lambda2 for all modeled values
meanGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
[r c] = size(meanGivMod);
for i = 1:c
    if i == 1
        meanGivMod(:,i) = interp1(mean_Mod,mean_Obs,dailyCTMv,'linear','extrap');
        varGivMod(:,i) = interp1(mean_Mod,var_Obs,dailyCTMv,'linear','extrap');
    else
        meanGivMod(:,i) = interp1(mean_Mod,mean_Obs,repmat(modplots(i-1),r,1),'linear','extrap');
        varGivMod(:,i) = interp1(mean_Mod,var_Obs,repmat(modplots(i-1),r,1),'linear','extrap');
    end
end

% if values are negative, make them zero, but show S-curves
meanGivMod(meanGivMod<0) = 0;

% save results
save(sprintf('CAMPmethod_%d.mat',yr),'modplots','meanGivMod','varGivMod', ...
    'mean_Mod','mean_Obs','var_Obs','dailyCTMv','distCTMv','yrmodaCTMv');

end