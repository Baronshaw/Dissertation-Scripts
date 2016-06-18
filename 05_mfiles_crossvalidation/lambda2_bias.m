function [] = lambda2_bias()
% this function will take all the neighborhoods of an estimation location
% and determine if a smaller lambda2 implies a lower bias

% load LOO estimation information
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd-mI;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckall = cell(size(cMS,1),1);
for i = 1:size(cMS,1)
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckall{i,1} = ch(idx,:);
end
softstr = '_soft'; constr = '_long'; gaussstr = '_gauss'; 
load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm.mat',softstr,constr,gaussstr,0));
zkall = zk; % mean trend removed
zhall = zh; % mean trend removed
vkall = vk;     
zkall = cell2mat(zkall);
idx = ~isnan(zkall); zkall = zkall(idx);
zhall = zhall(idx);
ckall = cell2mat(ckall); ckall = ckall(idx,:);
vkall = cell2mat(vkall); vkall = vkall(idx);

% load neighborhood information 
load('../matfiles/GetSoftNeib_1.mat');
allcs = allcs(idx);
alllambda1 = alllambda1(idx); ix=cellfun(@isempty,alllambda1); alllambda1(ix)={NaN*ones(3,1)}; 
alllambda2 = alllambda2(idx); ix=cellfun(@isempty,alllambda2); alllambda2(ix)={NaN*ones(3,1)}; 

% limit everything to 2001 and 2002
dates = datevec(ckall(:,3));
dates = dates(:,1);
idx = dates == 2001 | dates == 2002;

alllambda1 = alllambda1(idx); alllambda2 = alllambda2(idx); allcs = allcs(idx);
zhall = zhall(idx); zkall = zkall(idx);

matlam1 = cell2mat(alllambda1);
matlam2 = cell2mat(alllambda2);

% group lambda1 values by percentile 
lam1per = prctile(matlam1,0:100);
lam2per = cell(length(lam1per)-1,1);

% within each percentile group of lambda1's, compare the corresponding``
% lambda2's and corresponding bias
for i = 1:length(lam1per)-1
    
    idx = lam1per(i) <= matlam1 & lam1per(i+1) > matlam1;
    lam2per{i} = matlam2(idx);
    
    hardidx = ceil(find(idx==1)./3); % '3' because of soft data neighb
    
    % for each 'lam2per' in group, see the corresponding bias, keep in mind
    % that one estimate can include up to 3 lambda1 percentiles 
    biasper = abs(zhall(hardidx) - zkall(hardidx));
    
    % plot lam2per versus bias, hopefully there is a trend
    figure; hold on;
    plot(lam2per{i},biasper,'b.');
    title(sprintf('lambda1 versus bias for lam1 = (%f to %f)',lam1per(i),lam1per(i+1)));
    xlabel('lambda2');
    ylabel('bias');
    
    
end

end