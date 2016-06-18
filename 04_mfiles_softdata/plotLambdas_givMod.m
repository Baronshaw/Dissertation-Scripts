function [] = plotLambdas_givMod(dayWHI,daysWHIdisp,n,deltaT,numbins,minpnts,modplots,...
    negval)
% this function will create plots across the US for lambda1 and lambda2 for
% different fixed modeled values

% parameters subject to change
if nargin < 1, dayWHI = [ 2001 06 15 ]; end
if nargin < 2, daysWHIdisp = 20010615; end    
if nargin < 3, n = 3; end % number of closest stations
if nargin < 4, deltaT = 365; end % number of days in the interval
if nargin < 5, numbins = 10; end % number of bins for each plot
if nargin < 6, minpnts = 150; end % min number of point to each plot
if nargin < 7, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 8, negval = 0; end % 0 = there are no negative predicted values

% loading data
load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat',daysWHIdisp,deltaT,n,minpnts,numbins,negval));
load('../09_mfiles_projections/USAstates5.mat');
load('../09_mfiles_projections/USAcontiguous.mat');
maskcontour=[x y];
temp = jet;
cmap=temp(end:-1:1,:);  
lax = [-127 -65.6 24.87 49.7];
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',5,[0 0 0]};
colormap(cmap)

% country outline
load('../09_mfiles_projections/USAcontiguous.mat');
cd ../09_mfiles_projections
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../04_mfiles_softdata

for i = 1:length(modplots)

    %%% lambda1
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(CTMlocs,meanGivMod(:,i+1),lax,[plotax(:,1) plotax(:,2)]);
    colorbar;
    cax = [prctile(Zg(:),5) prctile(Zg(:),95)];
    caxis(cax);

    xlabel('km');
    ylabel('km');

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../04_softdata
        plot(states(:,1),states(:,2),'k-');
    end
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);

    title(sprintf('\\lambda_{1} given %d on %d, \\DeltaT=%d days\n n=%d, bins=%d, minpnts=%d',...
        modplots(i),daysWHIdisp,deltaT,n,numbins,minpnts));

    % save figures
    h = gcf;
    print(h,'-painters','-dpdf','-r600',sprintf('../plots/map_lambda1giv%d_%d_%d_%d_%d_%d.pdf', ...
        modplots(i),daysWHIdisp,deltaT,n,minpnts,numbins));

    %%% lambda2
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(CTMlocs,varGivMod(:,i+1),lax,[plotax(:,1) plotax(:,2)]);
    colorbar;
    cax = [prctile(Zg(:),5) prctile(Zg(:),95)];
    caxis(cax);

    xlabel('km');
    ylabel('km');

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../04_mfiles_softdata
        plot(states(:,1),states(:,2),'k-');
    end
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);

    title(sprintf('\\lambda_{2} given %d on %d, \\DeltaT=%d days\n n=%d, bins=%d, minpnts=%d',...
        modplots(i),daysWHIdisp,deltaT,n,numbins,minpnts));

    % save figures
    h = gcf;
    print(h,'-painters','-dpdf','-r600',sprintf('../plots/map_lambda2giv%d_%d_%d_%d_%d_%d.pdf', ...
    modplots(i),daysWHIdisp,deltaT,n,minpnts,numbins));

end

end