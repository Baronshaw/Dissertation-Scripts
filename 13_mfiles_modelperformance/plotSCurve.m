function [] = plotSCurve(dayWHI,daysWHIdisp,CTMsearch,n,deltaT,numbins,minpnts,modplots,negval)
% this function will create plots for given mean and variances

% parameters subject to change
if nargin < 1, dayWHI = [ 2001 1 1 ]; end
if nargin < 2, daysWHIdisp = 20010101; end  
if nargin < 3, CTMsearch = [ -990000 486000 ; -1422000 -378000 ; -2214000 990000 ; -1530000 -702000 ]; end
if nargin < 4, n = 3; end % number of closest stations
if nargin < 5, deltaT = 365; end % number of days in the interval
if nargin < 6, numbins = 10; end % number of bins for each plot
if nargin < 7, minpnts = 150; end % min number of point to each plot
if nargin < 8, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 9, negval = 0; end % 0 = there are no negative predicted values

% loading data
load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat',daysWHIdisp,deltaT,n,minpnts,numbins,negval));

[r c] = size(CTMsearch);

% make maps of these locations
figure; hold on;
load('../09_mfiles_projections/USAstates5.mat');
for k = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
    cd ../13_mfiles_modelperformance
    plot(states(:,1),states(:,2),'k-');
end
plot(CTMsearch(:,1),CTMsearch(:,2),'b.');
title('SCurves of interest');    
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600',sprintf('figures/SCurvesOfInterest%d.png',daysWHIdisp));

for i = 1:r
    
    finddist = sqrt( (CTMsearch(i,1)-CTMlocs(:,1)).^2 + (CTMsearch(i,2)-CTMlocs(:,2)).^2 );
    idx = finddist == min(finddist);
    
    figure; hold on;

    plot(valMod(idx,:),valObs(idx,:),'b.');
    plot([get(gca,'XLim')],[get(gca,'XLim')],'k--');
    xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
    ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
    title(sprintf('PM_{2.5} (\\mug/m^3) on %d \nx=%0.0f km, y=%0.0f km',...        
        daysWHIdisp,CTMlocs(idx,1)./1000,CTMlocs(idx,2)./1000));

    % add bins
    for j = 1:numbins+1
        plot([perctile_data(idx,j) perctile_data(idx,j)],get(gca,'YLim'),'r-');
    end

    % calculate means in each bin
    plot(mean_Mod(idx,:),mean_Obs(idx,:),'ro','MarkerFaceColor','r');

    % linear interpolation in each bin
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    interpx = linspace(xlim(1),xlim(2));
    interpy = interp1(mean_Mod(idx,:),mean_Obs(idx,:),interpx,'linear','extrap'); 
    plot(interpx,interpy,'r-');

    % add modeled value
    plot([dailyCTMg(idx) dailyCTMg(idx)],[ylim(1) meanGivMod(idx,1)], ...
        'c--','LineWidth',2);
    plot([xlim(1) dailyCTMg(idx)],[meanGivMod(idx,1) meanGivMod(idx,1)], ...
        'c--','LineWidth',2);

%     % add different modeled values
%     for j = 1:length(modplots)
%         plot([modplots(j) modplots(j)],[ylim(1) meanGivMod(idx,j+1)], ...
%         'k--','LineWidth',2);
%         plot([xlim(1) modplots(j)],[meanGivMod(idx,j+1) meanGivMod(idx,j+1)], ...
%             'k--','LineWidth',2);
%     end

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/Scurve_%d_x%dkm_y%dkm.png', ...
        daysWHIdisp,floor(CTMsearch(i,1)./1000),floor(CTMsearch(i,2)./1000)));

end
close all;

end