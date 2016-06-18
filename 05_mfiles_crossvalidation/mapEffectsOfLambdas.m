function [] = mapEffectsOfLambdas()
% this function will create plots/histograms/maps that show the results of
% SeeEffectOfLambda1.m. Hopefully I'll be se able to isolate to see
% when/where/how many times lambda1 SHOULD have been helpful and if it was
% actually helpful.

constant = 0; 
gauss = 1; 

if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;

% loading all hard and soft data results
for i = 1:10   
    
    % loading results
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',i,'_nosoft',constr,gaussstr));
    zkallh{i,1} = zk_madd;
    zhallh{i,1} = zh_Xval;
    ckallh{i,1} = ck;
    vkallh{i,1} = vk;  
    foldallh{i,1} = i*ones(length(zkallh{i,1}),1);
    
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',i,'_soft',constr,gaussstr));
    zkalls{i,1} = zk_madd;
    zhalls{i,1} = zh_Xval;
    ckalls{i,1} = ck;
    vkalls{i,1} = vk;  
    foldalls{i,1} = i*ones(length(zkalls{i,1}),1);
    
end

zkallh = cell2mat(zkallh); 
zkalls = cell2mat(zkalls);
zhallh = cell2mat(zhallh); 
ckallh = cell2mat(ckallh); 
vkallh = cell2mat(vkallh); 
foldallh = cell2mat(foldallh); 
zhalls = cell2mat(zhalls); 
ckalls = cell2mat(ckalls); 
vkalls = cell2mat(vkalls); 
foldalls = cell2mat(foldalls); 

% loop through each year with soft data
soft_years = [2001:2002 2005 2006:2007];
nsmax = 3;
for i = 1:length(soft_years)
    
    % load lambda neighborhoods
    load(sprintf('dataneighb_%d.mat',soft_years(i)));
    zh_lambda = NaN*ones(length(ck),1);
    zkh_lambda = NaN*ones(length(ck),1);
    lambda_help = NaN*ones(length(ck),nsmax);
    absdiff_help = NaN*ones(length(ck),nsmax);
    reldiff_help = NaN*ones(length(ck),nsmax);
    lambda_help2 = NaN*ones(length(ck),nsmax);
    zks_lambda = NaN*ones(length(ck),1);
    actuallambda_help = NaN*ones(length(ck),1);
    varlambda = NaN*ones(length(ck),1);
    
    [lia lib] = ismember(ck,ckallh,'rows');
    lib(lib==0) = [];
    % go through each location
    for j = 1:length(ck)
        
        if mod(j,1000) == 0, disp(j); end

        % get obs for that location
        zh_lambda(j) = zhallh(lib(j));

        % get hard zk for that location
        zkh_lambda(j) = zkallh(lib(j));

        % get lambda1 for that location

        % should lambda1 have helped? yes/no, i.e. is the bias between lambda1 and 
        % obs < the bias between hard zk and obs, also by how much? Calculate
        % difference and relative difference.
        if ~isempty(zslocal{j})
            lambda_help(j,:) = abs(zslocal{j}-zh_lambda(j))<abs(zkh_lambda(j)-zh_lambda(j));
            absdiff_help(j,:) = abs(zslocal{j}-zh_lambda(j));
            reldiff_help(j,:) = 100*(zslocal{j}-zh_lambda(j))./zh_lambda(j);
            if zkh_lambda(j)-zh_lambda(j) > 0,
                lambda_help2(j,:) = zslocal{j} <= zkh_lambda(j);
            elseif zkh_lambda(j)-zh_lambda(j) < 0,
                lambda_help2(j,:) = zslocal{j} >= zkh_lambda(j);
            end
        end

        % was my prediction correct? why or why not? This is where we check out
        % lambda2.
        zks_lambda(j) = zkalls(lib(j));
        actuallambda_help(j) = abs(zks_lambda(j)-zh_lambda(j))<abs(zkh_lambda(j)-zh_lambda(j));

        % also can check a running theory that all the lambda1's (and
        % lambda2's) in a neighborhood at all equal, meaning if one lambda1 is
        % bad, all lambda1's will be bad
        varlambda(j) = var(zslocal{j});
    
    end
    
    % save results
    save(sprintf('help_lambda_%d.mat',soft_years(i)),'zh_lambda','zkh_lambda', ...
        'lambda_help','absdiff_help','reldiff_help','lambda_help2','zks_lambda', ...
        'actuallambda_help','varlambda');
   
end

end