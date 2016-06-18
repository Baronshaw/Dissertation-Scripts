function [] = displayMeanTrendSoft(FOLDIDX,daydisp)
% this function will create maps for a given fold on certain days to
% compare the difference between the previous mean trend and the true 10
% fold cross validation for the soft data

if nargin < 1, FOLDIDX = 1; end
if nargin < 2, daydisp = [7 1]; end % day in each year to display

% load previous mean trend
smoothingParam = [900000 300000 100 50];
CTMyears = [2001:2002 2005:2007];
for i = 1:length(CTMyears)
    load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',smoothingParam,CTMyears(i)));
    mI_prev{i,1} = mI;
    pI_prev{i,1} = pI;
    zd_prev{i,1} = zd;
    pd_prev{i,1} = pd;
end
mI_prev = cell2mat(mI_prev); pI_prev = cell2mat(pI_prev);
zd_prev = cell2mat(zd_prev); pd_prev = cell2mat(pd_prev);
temp_prev = datevec(pI_prev(:,3));

% load new mean trend
smoothingParam = [2*900000 300000 100 50];
load(sprintf('../matfiles/meanTrendSoft_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam));

% getting unique days to display
temp = datevec(pI(:,3));
idx = temp(:,2) == daydisp(1) & temp(:,3) == daydisp(2);
toplot = unique(temp(idx,1:3),'rows');

% loop through the first of each month and plot both figures
for i = 1:size(toplot,1)
    
    % plot the previous mean trend first
    idx = temp_prev(:,1) == toplot(i,1) & temp_prev(:,2) == toplot(i,2) & temp_prev(:,3) == toplot(i,3);
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../11_mfiles_truecrossvalidation    
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(mI_prev,5) prctile(mI_prev,95)];
    plotField(pI_prev(idx,1:2),mI_prev(idx),[ -3000000 3000000 -2000000 1500000 ],plotax);   
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../11_mfiles_truecrossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    caxis(cax);
    colorbar;
    title(sprintf('Previous soft mean trend for %d %d %d',toplot(i,:)));
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('PrevMeanTrendSoft_%d_%d_%d.pdf',toplot(i,:)));
    print(gcf,'-painters','-dpng','-r600',sprintf('PrevMeanTrendSoft_%d_%d_%d.png',toplot(i,:)));
    
    % plot the new mean trend first
    idx = temp(:,1) == toplot(i,1) & temp(:,2) == toplot(i,2) & temp(:,3) == toplot(i,3);
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../11_mfiles_truecrossvalidation    
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(mI_prev,5) prctile(mI_prev,95)];
    plotField(pI(idx,1:2),mIsoft_nonXval(idx),[ -3000000 3000000 -2000000 1500000 ],plotax);   
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../11_mfiles_truecrossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    caxis(cax);
    colorbar;
    title(sprintf('New soft mean trend for fold %d %d %d %d',FOLDIDX,toplot(i,:)));
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('NewMeanTrendSoft_fold%d_%d_%d_%d.pdf',FOLDIDX,toplot(i,:)));
    print(gcf,'-painters','-dpng','-r600',sprintf('NewMeanTrendSoft_fold%d_%d_%d_%d.png',FOLDIDX,toplot(i,:)));
        
end

end