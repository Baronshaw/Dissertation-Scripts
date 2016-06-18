function [] = exploreEffectsOfLambdas()
% this function will explore the effects of the soft data on the estimate,
% specifically if there is a geographical/temporal trend to performance

constant = 0; 
gauss = 1; 

if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

soft_years = [2001:2002 2005 2006:2007]; 

% gathering all the data
for i = 1:10   
    
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s.mat',i,'_soft',constr,gaussstr));
    zkallsall{i,1} = zk_madd;
    zhallsall{i,1} = zh_Xval;
    ckallsall{i,1} = ck;
    vkalls{i,1} = vk;  
    foldalls{i,1} = i*ones(length(zkallsall{i,1}),1);
    
end

% zkallsall = cell2mat(zkallsall); idx = ~isnan(zkallsall);
% zkallsall = zkallsall(idx);
% zhallsall = cell2mat(zhallsall); zhallsall = zhallsall(idx);
% ckallsall = cell2mat(ckallsall); ckallsall = ckallsall(idx,:);
% vkalls = cell2mat(vkalls); vkalls = vkalls(idx);
% foldalls = cell2mat(foldalls); foldalls = foldalls(idx);

zkallsall = cell2mat(zkallsall); 
zhallsall = cell2mat(zhallsall); 
ckallsall = cell2mat(ckallsall); 
vkalls = cell2mat(vkalls); 
foldalls = cell2mat(foldalls); 

temp = datevec(ckallsall(:,3));
yrs = temp(:,1);

% looping through each soft year
for i = 1:length(soft_years)
    
    zkalls = zkallsall(yrs==soft_years(i));
    zhalls = zhallsall(yrs==soft_years(i));
    ckalls = ckallsall(yrs==soft_years(i),:);
    % load results
    load(sprintf('help_lambda_%d.mat',soft_years(i)));
    load(sprintf('dataneighb_%d.mat',soft_years(i)));
    actuallambda_help2 = ~actuallambda_help;
    test1 = nansum(lambda_help2,2); 
    % I = improved, E = estimate, W = worsened, S = soft
    
    idx = [];
    idx{1} = actuallambda_help == 1 & test1 == 0; % soft improved, but shouldn't have   
    idx{2} = actuallambda_help == 1 & test1 == 3; % soft improved, and should have    
    idx{3} = actuallambda_help == 0 & test1 == 0; % soft worsened, and should have       
    idx{4} = actuallambda_help == 0 & test1 == 3; % soft worsened, but shouldn't have   

    for j = 1:length(idx)
    
        % looking at overall absoluate error
        allAE = abs(zkalls(idx{j}) - zhalls(idx{j}));
        figure; hold on;
        hist(allAE,100);
        print(gcf,'-painters','-dpng','-r600',sprintf('overallabsbias_yr%d_cat%d.png',soft_years(i),j));

        % looking at overall relative error
        allperdiff = 100*(zkalls(idx{j}) - zhalls(idx{j}))./zkalls(idx{j});
        figure; hold on;
        hist(allperdiff,100);
        print(gcf,'-painters','-dpng','-r600',sprintf('overallrelbias_yr%d_cat%d.png',soft_years(i),j));

        % loop through each station and count when soft was better
        cksub = ckalls(idx{j},:);
        zkallssub = zkalls(idx{j});
        zhallssub = zhalls(idx{j});
        unisub = unique(cksub(:,1:2),'rows');
        stationAE = NaN*ones(length(unisub),1);
        stationperdiff = NaN*ones(length(unisub),1);
        for k = 1:size(unisub)
            uniidx = unisub(k,1) == cksub(:,1) & unisub(k,2) == cksub(:,2);
            stationAE(k) = nanmean( abs(zkallssub(uniidx) - zhallssub(uniidx)) );
            stationperdiff(k) = nanmean( 100*(zkallssub(uniidx) - zhallssub(uniidx))./zkallssub(uniidx) );
        end

        % map showing abs bias across the US
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
        for k = 1:length(X)
            cd ../09_mfiles_projections
            states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
            cd ../05_mfiles_crossvalidation
            plot(states(:,1),states(:,2),'k-');
        end
        % colorplot
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',5,[0 0 0]};
        colorplot(unisub,stationAE,'hot',Property,Value);
        cax = [prctile(stationAE,5) prctile(stationAE,95)];
        caxis(cax);
        colorbar;
        title(sprintf('abs bias for ea loc in %d for category %d',soft_years(i),j)); 
        print(gcf,'-painters','-dpng','-r600',sprintf('overallabsbias_yr%d_cat%d.png',soft_years(i),j));

        % map showing relative bias across the US
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
        for k = 1:length(X)
            cd ../09_mfiles_projections
            states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
            cd ../05_mfiles_crossvalidation
            plot(states(:,1),states(:,2),'k-');
        end
        % colorplot
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',5,[0 0 0]};
        colorplot(unisub,stationperdiff,'hot',Property,Value);
        cax = [prctile(stationperdiff,5) prctile(stationperdiff,95)];
        caxis(cax);
        colorbar;
        title(sprintf('rel diff for ea loc in %d for category %d',soft_years(i),j));
        print(gcf,'-painters','-dpng','-r600',sprintf('overallrelbias_yr%d_cat%d.png',soft_years(i),j));

        % create the same plots for each day
        cksub = ckalls(idx{j},:);
        zkallssub = zkalls(idx{j});
        zhallssub = zhalls(idx{j});
        uniday = unique(cksub(:,3)) - datenum(soft_years(i)-1,12,31);
        cksubtemp = cksub - datenum(soft_years(i)-1,12,31);
        dayAE = NaN*ones(length(uniday),1);
        dayperdiff = NaN*ones(length(uniday),1);
        for k = 1:size(uniday,1)
            uniidx = uniday(k) == cksubtemp(:,3);
            dayAE(k) = nanmean( abs(zkallssub(uniidx) - zhallssub(uniidx)) );
            dayperdiff(k) = nanmean( 100*(zkallssub(uniidx) - zhallssub(uniidx))./zkallssub(uniidx) );
        end

        % time series of abs bias
        figure; hold on;
        plot(uniday,dayAE,'b--');
        title(sprintf('abs bias ea day in %d for category %d',soft_years(i),j));
        print(gcf,'-painters','-dpng','-r600',sprintf('TSabsbias_yr%d_cat%d.png',soft_years(i),j));

        % time series of rel bias
        figure; hold on;
        plot(uniday,dayperdiff,'b--');
        title(sprintf('rel bias ea day in %d for category %d',soft_years(i),j));
        print(gcf,'-painters','-dpng','-r600',sprintf('TSrelbias_yr%d_cat%d.png',soft_years(i),j));
    
    end
    
    close all
           
end

end