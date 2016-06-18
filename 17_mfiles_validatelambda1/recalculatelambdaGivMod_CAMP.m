function [] = recalculatelambdaGivMod_CAMP()
% this function will perform a recalculation of CAMP and find lambda1 and
% lambda2 for mixed modeled values and for the modeled values within the
% CMAQ grid itself

% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "recalculatelambdaGivMod_CAMP" -logfile "recalculatelambdaGivMod_CAMP.out"

% loading BME function
cd ../BMELIB2.0b
startup
cd ../17_mfiles_validatelambda1

modplots = 0:5:50;

% load data
load('recalculatelambdaGivPer_CAMP.mat');

% load all modeled values for 2001
load(sprintf('../matfiles/prepCTM_%d.mat',2001)); 

% matching CMAQ locations to recalculated locations
if ~exist('recalculate_dailyvCTMorder.mat')
    cssCMAQ = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];
    [aidx bidx] = ismember(css,cssCMAQ,'rows');
    dailyCTMvorder = dailyCTMv(bidx);
    disp(sprintf('proportion of aidx = %0.2f',sum(aidx)/length(aidx)));
    save('recalculate_dailyvCTMorder.mat','dailyCTMvorder')
else
    load('recalculate_dailyvCTMorder.mat')
end

% final variables
lambda1_recalcModAll = NaN*ones(length(css),length(modplots)+1);
lambda2_recalcModAll = NaN*ones(length(css),length(modplots)+1);
lambda1_recalcsimuModAll = NaN*ones(length(css),length(modplots)+1);
lambda2_recalcsimuModAll = NaN*ones(length(css),length(modplots)+1);

% loop through all modeled values locations
matlabpool open 12
parfor i = 1:length(css)
    if mod(i,1000)==0, disp(i); end
    lambda1_recalcModAll(i,:) = interp1(meanMod_recalc(i,:),lambda1_recalc(i,:),[modplots dailyCTMvorder(i)],'linear','extrap');
    lambda2_recalcModAll(i,:) = interp1(meanMod_recalc(i,:),lambda2_recalc(i,:),[modplots dailyCTMvorder(i)],'linear','extrap');
    lambda1_recalcsimuModAll(i,:) = interp1(meanMod_recalc(i,:),lambda1_recalcsimu(i,:),[modplots dailyCTMvorder(i)],'linear','extrap');
    lambda2_recalcsimuModAll(i,:) = interp1(meanMod_recalc(i,:),lambda2_recalcsimu(i,:),[modplots dailyCTMvorder(i)],'linear','extrap');
end
matlabpool close

lambda1_recalcMod = lambda1_recalcModAll(:,1:end-1);
lambda2_recalcMod = lambda2_recalcModAll(:,1:end-1);
lambda1_recalcsimuMod = lambda1_recalcsimuModAll(:,1:end-1);
lambda2_recalcsimuMod = lambda2_recalcsimuModAll(:,1:end-1);

lambda1_recalcModGrid = lambda1_recalcModAll(:,end);
lambda2_recalcModGrid = lambda2_recalcModAll(:,end);
lambda1_recalcsimuModGrid = lambda1_recalcsimuModAll(:,end);
lambda2_recalcsimuModGrid = lambda2_recalcsimuModAll(:,end);

% save results
save('recalculatelambdaGivMod_CAMP.mat','css','modplots','dailyCTMvorder', ...
    'lambda1_recalcMod','lambda2_recalcMod','lambda1_recalcModGrid','lambda2_recalcModGrid', ...
    'lambda1_recalcsimuMod','lambda2_recalcsimuMod','lambda1_recalcsimuModGrid','lambda2_recalcsimuModGrid');

end
