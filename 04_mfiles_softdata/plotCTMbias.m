function [] = plotCTMbias(dayWHI,daysWHIdisp,n,deltaT,numbins,minpnts,modplots,...
    negval)
% this function will plot lambda1-CTM

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
% temp = jet;
% cmap=temp(end:-1:1,:);  
lax = [-127 -65.6 24.87 49.7];
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',5,[0 0 0]};
cmap = jet;
colormap(cmap)

% country outline
load('../09_mfiles_projections/USAcontiguous.mat');
cd ../09_mfiles_projections
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../04_mfiles_softdata

% loading CTM data
load(sprintf('../matfiles/prepCTM_%d.mat',dayWHI(1)));
zh = dailyCTMv;
ch = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];
idx = ch(:,3) == datenum(dayWHI);

if sum(idx) ~= 0 % if soft data exist
    
    % matching the CTM with the lambdas
    if size(CTMlocs,1) == size(meanGivMod(:,1),1)
        [lia lib] = ismember(ch(idx,1:2),CTMlocs,'rows');
    else
        [lia lib] = ismember(ch(idx,1:2),cMSCTMv,'rows');
    end
    lib(lib==0) = [];

    %%% lambda1 - CTM
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(ch(idx,1:2),meanGivMod(lib,1)-zh(idx),lax,[plotax(:,1) plotax(:,2)]);
    colorbar;
    cax = [-10 15];
    caxis(cax);

    xlabel('km');
    ylabel('km');

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for i = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
        cd ../04_mfiles_softdata
        plot(states(:,1),states(:,2),'k-');
    end
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);

    title(sprintf('\\lambda_{1}-CTM on %d, \\DeltaT=%d days\n n=%d, bins=%d, minpnts=%d',...
        daysWHIdisp,deltaT,n,numbins,minpnts));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/map_lambda1bias_%d.pdf',daysWHIdisp));

end

end