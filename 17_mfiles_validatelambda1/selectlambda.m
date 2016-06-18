function [] = selectlambda(yr2select)
% this function will gather the resulting lambda1 and lambda2 for Paper One
% for fixed modeld values and use it as the "selected" lambda1 and lambda2
% for the simulation

% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "selectlambda" -logfile "selectlambda.out"

if nargin < 1, yr2select = 2002; end

% gathering all lambda1s and lambda2s for fixed modeled values
n = 3; 
deltaT = 365; 
numbins = 10; 
minpnts = 150; 
negval = 0;

temp = datevec(datenum(yr2select,1,1):datenum(yr2select,12,31));
dayWHI = temp(:,1:3);
dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
[r c] = size(dayWHI);
css = cell(r,1);
meanGivMod_all = cell(r,1);
varGivMod_all = cell(r,1);
modplots_all = cell(r,1);
dailyCTMg_all = cell(r,1);

for i = 1:r
    disp([yr2select i]);
    load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat',dayWHIdisps(i),deltaT,n,minpnts,numbins,negval));
    css{i,1} = [cMSCTMv repmat(datenum(dayWHI(i,1),dayWHI(i,2),dayWHI(i,3)),length(dailyCTMg),1)];
    meanGivMod_all{i,1} = meanGivMod(:,2:end);
    varGivMod_all{i,1} = abs(varGivMod(:,2:end)); % no there's no zeros
    modplots_all{i,1} = repmat(modplots,length(dailyCTMg),1);
    dailyCTMg_all{i,1} = dailyCTMg;
end

css = cell2mat(css);
meanGivMod_all = cell2mat(meanGivMod_all);
varGivMod_all = cell2mat(varGivMod_all);
modplots_all = cell2mat(modplots_all);
dailyCTMg_all = cell2mat(dailyCTMg_all);

meanGivCMAQ_all = NaN*ones(length(css),1);
varGivCMAQ_all = NaN*ones(length(css),1);

% get the lambdas given CMAQ modeled value from above
matlabpool open 12
parfor i = 1:length(css)
    if mod(i,1000)==0, disp(i); end
    meanGivCMAQ_all(i) = interp1(modplots_all(i,:),meanGivMod_all(i,:),dailyCTMg_all(i),'linear','extrap');
    varGivCMAQ_all(i) = interp1(modplots_all(i,:),varGivMod_all(i,:),dailyCTMg_all(i),'linear','extrap');
end
matlabpool close

save(sprintf('PM2p5_meanGivMod_yr%d.mat',yr2select), ...
    'css','meanGivMod_all','varGivMod_all','modplots_all', ...
    'meanGivCMAQ_all','varGivCMAQ_all');

end