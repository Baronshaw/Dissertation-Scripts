function [] = getXvalStatistics_LOOCV_summary(soft,constant,gauss,forceddist)
% this function will calculate the cross validation statistics for the 
% LOOCV cross validation across all soft data methods attempted

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, forceddist = 0; end % distance in m

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckallC = cell(size(cMS,1),1);
ckallhC = cell(size(cMS,1),1);
for i = 1:size(cMS,1)
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckallC{i,1} = ch(idx & temp2001==2001,:);
    ckallhC{i,1} = ch(idx,:);
end

softtry = { '' ; '' ; '_n6' ; '_n12' ; '_n6T90' ; '_lam2b' ; '_1p5lam2' ; ...
    '_2lam2' ; '_LinReg' ; '_n6delT365' ; '_n12delT365' };

for i = 1:length(softtry)
    disp(i);
    % loading results
    if i == 1
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm.mat', ...
            '_nosoft',constr,gaussstr,floor(forceddist./1000)));
    else
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s.mat', ...
            softstr,constr,gaussstr,floor(forceddist./1000),softtry{i}));
    end
    zkall = zk_madd;
    zhall = zh_Xval;
    vkall = vk;     

    zkall = cell2mat(zkall);
    zhall = cell2mat(zhall);
    vkall = cell2mat(vkall);
    
    if i == 1
        ckalls = cell2mat(ckallC); ckallh = cell2mat(ckallhC); 
        [lia lib] = ismember(ckalls,ckallh,'rows');
        zkall = zkall(lib);
        zhall = zhall(lib);
        vkall = vkall(lib);
    end
    idx = ~isnan(zkall); zkall = zkall(idx);
    ckalls = cell2mat(ckallC); ckalls = ckalls(idx,:);
    zhall = zhall(idx);
    vkall = vkall(idx);
    
    % overall 
    n_sites = size(unique(ckalls(:,1:2),'rows'),1);
    n_observations = size(ckalls,1);
    errors = zkall - zhall;
    RMSE(i) = sqrt( mean(errors.^2) );
    MAE(i) = mean(abs(errors));
    ME(i) = mean(errors);
    r2(i) = (corr(zkall,zhall,'type','Pearson')).^2;
    MS(i) = mean(errors./sqrt(vkall));
    RMSS(i) = std(errors./sqrt(vkall));
    MR(i) = mean(sqrt(vkall));
    std_est(i) = std(zkall);
    std_obs(i) = std(zhall);
    QAr2(i) = ( ( std_est(i).^2 + std_obs(i).^2 - (RMSE(i).^2-ME(i).^2) )./(2*std_est(i)*std_obs(i)) )^2;

end
allstats = [RMSE' MAE' ME' r2' MS' RMSS' MR' std_est' std_obs' QAr2'];

% save results
save(sprintf('Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_results_summary.mat', ...
    softstr,constr,gaussstr,floor(forceddist./1000)),...
    'n_sites','n_observations','RMSE','MAE','ME','r2','MS','RMSS','MR','std_est',...
    'std_obs','QAr2','softtry','allstats'); 

end