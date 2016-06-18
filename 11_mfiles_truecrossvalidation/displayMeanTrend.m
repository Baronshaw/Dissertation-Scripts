function [] = displayMeanTrend(FOLDIDX,daydisp)
% this function will create maps for a given fold on certain days to
% compare the difference between the previous mean trend and the true 10
% fold cross validation

if nargin < 1, FOLDIDX = 1; end
if nargin < 2, daydisp = [7 1]; end % day in each year to display

% load previous mean trend
smoothingParam = [900000 300000 100 50];
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',smoothingParam));
mI_prev = mI; pI_prev = pI; zd_prev = zd; pd_prev = pd;
temp_prev = datevec(pI_prev(:,3));

% load new mean trend
smoothingParam = [900000 300000 100 50];
load(sprintf('../matfiles/meanTrend_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam));

% getting unique days to display
temp = datevec(ch_nonXval(:,3));
idx = temp(:,2) == daydisp(1) & temp(:,3) == daydisp(2);
toplot = unique(temp(idx,1:3),'rows');

% loop through the first of each month and plot both figures
for i = 1:size(toplot,1)
    
    % plot the previous mean trend first
    idx = temp_prev(:,1) == toplot(i,1) & temp_prev(:,2) == toplot(i,2) & temp_prev(:,3) == toplot(i,3);
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../11_mfiles_truecrossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../11_mfiles_truecrossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(mI_prev,5) prctile(mI_prev,95)];
    colorplot(pd_prev(idx,1:2),mI_prev(idx),'hot',Property,Value,cax);    
    caxis(cax);
    colorbar;
    title(sprintf('Previous mean trend for %d %d %d',toplot(i,:)));
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('PrevMeanTrend_%d_%d_%d.pdf',toplot(i,:)));
    print(gcf,'-painters','-dpng','-r600',sprintf('PrevMeanTrend_%d_%d_%d.png',toplot(i,:)));
    
    % plot the new mean trend first
    idx = temp(:,1) == toplot(i,1) & temp(:,2) == toplot(i,2) & temp(:,3) == toplot(i,3);
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../11_mfiles_truecrossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../11_mfiles_truecrossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(mI_prev,5) prctile(mI_prev,95)];
    colorplot(ch_nonXval(idx,1:2),mI_nonXval(idx),'hot',Property,Value,cax);    
    caxis(cax);
    colorbar;
    title(sprintf('New mean trend for fold %d %d %d %d',FOLDIDX,toplot(i,:)));
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('NewMeanTrend_fold%d_%d_%d_%d.pdf',FOLDIDX,toplot(i,:)));
    print(gcf,'-painters','-dpng','-r600',sprintf('NewMeanTrend_fold%d_%d_%d_%d.png',FOLDIDX,toplot(i,:)));
        
end

end