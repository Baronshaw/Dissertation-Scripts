function [] = quantSoftIncDec()
% this function will quanitify the typical bias increase/reduction that
% comes from including the soft data by year

% some potential outcomes include 1) see if certain % increases or decreases
% are typically in a given geographical location or time of the year 2) see
% the typically increase/decrease due to soft data (maybe only a 5%
% decrease in bias is typical, we'll see)

constant = 0; 
gauss = 1; 

if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

soft_years = [2001:2002 2005:2007]; 

% gathering all the data
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

zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx);
zhallh = cell2mat(zhallh); zhallh = zhallh(idx);
ckallh = cell2mat(ckallh); ckallh = ckallh(idx,:);
vkallh = cell2mat(vkallh); vkallh = vkallh(idx);
foldallh = cell2mat(foldallh); foldallh = foldallh(idx);

zkalls = zkalls(idx);
zhalls = cell2mat(zhalls); zhalls = zhalls(idx);
ckalls = cell2mat(ckalls); ckalls = ckalls(idx,:);
vkalls = cell2mat(vkalls); vkalls = vkalls(idx);
foldalls = cell2mat(foldalls); foldalls = foldalls(idx);

temp = datevec(ckalls(:,3));
yrs = temp(:,1);

% looping through each soft year
for i = 1:length(soft_years)
    
    idx = soft_years(i) == yrs;
    
    % looking at absoluate error, count # of times it increased/decreased
    allAEs = abs(zkalls(idx) - zhalls(idx));
    allAEh = abs(zkallh(idx) - zhallh(idx));    
    idxbetter = allAEs<allAEh;
  
    % difference in absolute bias
    alldiff = allAEs - allAEh;
    figure; hold on;
    hist(alldiff(alldiff>-60),100);
    
    % what's the best metric to show: the difference in abs bias, rel diff?
    allperdiff = 100*(allAEs-allAEh)./allAEs;
    idxint = allperdiff > -900;
    figure; hold on;
    hist(allperdiff(idxint),100);
    
    % loop through each station and count when soft was better
    uniallincdiff = [];
    uniallperincdiff = [];
    unialldecdiff = [];
    uniallperdecdiff = [];
    cksub = ckalls(idx,:);
    unisub = unique(cksub(:,1:2),'rows');
    for j = 1:size(unisub)
        uniidx = unisub(j,1) == cksub(:,1) & unisub(j,2) == cksub(:,2);
        if sum(uniidx&idxbetter)>0
            uniallincdiff(j) = nanmean(allAEs(uniidx&idxbetter) - allAEh(uniidx&idxbetter));
            uniallperincdiff(j) = nanmean(100*(allAEs(uniidx&idxbetter)-allAEh(uniidx&idxbetter)) ...
                ./allAEs(uniidx&idxbetter));
        end
        if sum(uniidx&~idxbetter)>0
            unialldecdiff(j) = nanmean(allAEs(uniidx&~idxbetter) - allAEh(uniidx&~idxbetter));        
            uniallperdecdiff(j) = nanmean(100*(allAEs(uniidx&~idxbetter)-allAEh(uniidx&~idxbetter)) ...
                ./allAEs(uniidx&~idxbetter));
        end
    end
    
    % map showing increase diff across the US
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    colorplot(unisub,uniallincdiff','hot',Property,Value);
    cax = [prctile(uniallincdiff,5) prctile(uniallincdiff,95)];
    caxis(cax);
    colorbar;
    title(sprintf('diff for ea loc in %d where soft data reduced abs bias',soft_years(i)));
    
    % map showing decrease diff across the US
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    colorplot(unisub,unialldecdiff','hot',Property,Value);
    cax = [prctile(unialldecdiff,5) prctile(unialldecdiff,95)];
    caxis(cax);
    colorbar;
    title(sprintf('diff for ea loc in %d where soft data increased abs bias',soft_years(i)));
    
    % map showing increase % across the US
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    colorplot(unisub,uniallperincdiff','hot',Property,Value);
    cax = [prctile(uniallperincdiff,5) prctile(uniallperincdiff,95)];
    caxis(cax);
    colorbar;
    title(sprintf('%% for ea loc in %d where soft data reduced abs bias',soft_years(i)));
    
    % map showing decrease % across the US
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    colorplot(unisub,uniallperdecdiff','hot',Property,Value);
    cax = [prctile(uniallperdecdiff,5) prctile(uniallperdecdiff,95)];
    caxis(cax);
    colorbar;
    title(sprintf('%% for ea loc in %d where soft data increased abs bias',soft_years(i)));
    
    % create the same plots for each day
    uniallincdiffday = [];
    uniallperincdiffday = [];
    unialldecdiffday = [];
    uniallperdecdiffday = [];
    uniday = unique(cksub(:,3)) - datenum(soft_years(i)-1,12,31);
    cksubtemp = cksub - datenum(soft_years(i)-1,12,31);
    unicountday = [];
    for j = 1:size(uniday,1)
        uniidx = uniday(j) == cksubtemp(:,3);
        if sum(uniidx&idxbetter)>0
            uniallincdiffday(j) = nanmean(allAEs(uniidx&idxbetter) - allAEh(uniidx&idxbetter));
            uniallperincdiffday(j) = nanmean(100*(allAEs(uniidx&idxbetter)-allAEh(uniidx&idxbetter)) ...
                ./allAEs(uniidx&idxbetter));
        end
        if sum(uniidx&~idxbetter)>0
            unialldecdiffday(j) = nanmean(allAEs(uniidx&~idxbetter) - allAEh(uniidx&~idxbetter));        
            uniallperdecdiffday(j) = nanmean(100*(allAEs(uniidx&~idxbetter)-allAEh(uniidx&~idxbetter)) ...
                ./allAEs(uniidx&~idxbetter));
        end
    end
    
    % time series of diff increase
    figure; hold on;
    plot(uniday,uniallincdiffday,'b--');
    title(sprintf('diff of ea day in %d where soft data reduced abs bias',soft_years(i)));
    
    % time series of diff decrease
    figure; hold on;
    plot(uniday,unialldecdiffday,'b--');
    title(sprintf('diff of ea day in %d where soft data increased abs bias',soft_years(i)));
    
    % time series of % increase
    figure; hold on;
    plot(uniday,uniallperincdiffday,'b--');
    title(sprintf('%% of ea day in %d where soft data reduced abs bias',soft_years(i)));
    
    % time series of % decrease
    figure; hold on;
    plot(uniday,uniallperdecdiffday,'b--');
    title(sprintf('%% of ea day in %d where soft data increased abs bias',soft_years(i)));
           
end

end