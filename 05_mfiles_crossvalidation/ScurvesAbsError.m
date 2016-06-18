function [] = ScurvesAbsError(soft,constant,gauss,forceddist)
% this function will look at the relationship between absolute error and
% bin # and bin width in Scurves

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, forceddist = 0; end % distance in m

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% loading hard data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckall = cell(size(cMS,1),1);
for i = 1:size(cMS,1)
    disp(i);
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckall{i,1} = ch(idx & temp2001==2001,:);
end

% loading soft data
load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_n6T90.mat', ...
    softstr,constr,gaussstr,floor(forceddist./1000)));
zkall = zk_madd;
zhall = zh_Xval;
vkall = vk;     

zkall = cell2mat(zkall);
idx = ~isnan(zkall); zkall = zkall(idx);
zhall = cell2mat(zhall); zhall = zhall(idx);
ckall = cell2mat(ckall); ckall = ckall(idx,:);
vkall = cell2mat(vkall); vkall = vkall(idx);

% maybe check out how many locations didn't use soft data at all
blah = cell(length(allcs),1);
for i = 1:length(allcs)
    blah{i} = cell2mat(allcs{i});
end
testing = cellfun(@isempty,blah,'UniformOutput',false);
testing2 = cell2mat(testing);
idx = find(testing2==1);
figure; hold on;
plot(cMS(:,1),cMS(:,2),'bo');
plot(cMS(idx,1),cMS(idx,2),'ro');

% is this soft data missingness consistent across different soft data
% parameters?

% are the resulting estimates the same for over a 1/3 of values when
% compared with hard data?

% if yes, why is there no soft data used for all these estimation
% locations?

% calculate absolute bias to find when absolute bias is highest
absbias = abs(zkalls-zhall) - abs(zkallh-zkall); % <0 good, >0 bad

% retain the cell structure of estimates to more easily load soft data

% loop through several (all?) instances of absolute bias
% or loop through each day with soft data (that might be faster)
% for each of the 3 soft data, check the location/day

% find the days in a cell structure for the soft data pertaining to that
% day

% load the soft data for that day

% check to see for each space/time estimation location, for each three soft
% data points the: 1) bin number, bin width, average bin width (?)

% plot to see if there is some sort of trend between absolute error and bin
% number and bin width








end