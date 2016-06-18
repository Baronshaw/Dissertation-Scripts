function [] = plotLambdasavg(dayWHI,daysWHIdisp,avgval,n,deltaT,numbins,minpnts,modplots,...
    negval)
% the function will plot the average lambda1s and lambda2s over a month or
% a year


% parameters subject to change
if nargin < 1, dayWHI = [ 2001 06 15 ]; end
if nargin < 2, daysWHIdisp = 20010615; end 
if nargin < 3, avgval = 'mon'; end
if nargin < 4, n = 3; end % number of closest stations
if nargin < 5, deltaT = 365; end % number of days in the interval
if nargin < 6, numbins = 10; end % number of bins for each plot
if nargin < 7, minpnts = 150; end % min number of point to each plot
if nargin < 8, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 9, negval = 0; end % 0 = there are no negative predicted values

% loading data
load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',dayWHI(1)));
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

% getting days of interest
if strcmp(avgval,'mon') == 1
    idx = css(:,3) >= datenum([dayWHI(1:2) 1]) & css(:,3) <= datenum([dayWHI(1) dayWHI(2)+1 1]);
elseif strcmp(avgval,'yr') == 1
    idx = css(:,3) >= datenum([dayWHI(1) 1 1]) & css(:,3) <= datenum([dayWHI(1) 12 31]);
end

if sum(idx) ~= 0

    % avgeraging days
    csssub = css(idx,:);
    lambda1sub = lambda1(idx);
    lambda2sub = lambda2(idx);
    unilocs = unique(csssub(:,1:2),'rows');
    lambda1avg = NaN*ones(size(unilocs,1),1);
    lambda2avg = NaN*ones(size(unilocs,1),1);
    for i = 1:size(unilocs,1)
        if mod(i,100) == 0, disp(i); end
        idx = unilocs(i,1) == csssub(:,1) & unilocs(i,2) == csssub(:,2);
       lambda1avg(i) = nanmean(lambda1sub(idx)); % not exactly correct, since not IID
       lambda2avg(i) = (1/sum(idx)).*nanmean(lambda2sub(idx)); % not exactly correct, since not IID
    end

    %%% lambda1
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(unilocs,lambda1avg,lax,[plotax(:,1) plotax(:,2)]);
    colorbar;
    cax = [2 30];
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

    if strcmp(avgval,'mon')
        daydispsave = dayWHI(1)*100 + dayWHI(2);
        str1 = 'XX';
        str2 = 'mon';
    elseif strcmp(avgval,'yr')
        daydispsave = dayWHI(1);
        str1 = 'XXXX';
        str2 = 'yr';
    end

    title(sprintf('\\lambda_{1} on %d%s %s, \\DeltaT=%d days\n n=%d, bins=%d, minpnts=%d',...
        daydispsave,str1,str2,deltaT,n,numbins,minpnts));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/map_lambda1_%d%s_%s.pdf',daydispsave,str1,str2));

    %%% lambda2
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(unilocs,lambda2avg,lax,[plotax(:,1) plotax(:,2)]);
    colorbar;
    cax = [prctile(Zg(:),5) prctile(Zg(:),95)];
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

    title(sprintf('\\lambda_{2} on %d%s %s, \\DeltaT=%d days\n n=%d, bins=%d, minpnts=%d',...
        daydispsave,str1,str2,deltaT,n,numbins,minpnts));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/map_lambda2_%d%s_%s.pdf',daydispsave,str1,str2));

end

end