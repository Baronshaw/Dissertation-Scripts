function [] = plotScurve_deNazelle(dayWHI,daysWHIdisp,n,deltaT,numbins,minpnts,modplots,...
    negval,plotnum)
% this function will create plots for given mean and variances

% parameters subject to change
if nargin < 1, dayWHI = [ 2001 07 1 ]; end
if nargin < 2, daysWHIdisp = 20010701; end    
if nargin < 3, n = 3; end % number of closest stations
if nargin < 4, deltaT = 365; end % number of days in the interval
if nargin < 5, numbins = 10; end % number of bins for each plot
if nargin < 6, minpnts = 150; end % min number of point to each plot
if nargin < 7, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 8, negval = 0; end % 0 = there are no negative predicted values
if nargin < 9, plotnum = 5; end % number of plots displayed

% loading data
load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d_deNazelle.mat',daysWHIdisp,deltaT,n,minpnts,numbins,negval));
rand('seed',0);
idxMod = randsample(length(meanGivMod),plotnum);

for i = 1:plotnum
    
    figure; hold on;

    plot(valMod(idxMod(i),:),valObs(idxMod(i),:),'b.');
    plot([get(gca,'XLim')],[get(gca,'YLim')],'k--');
    xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
    ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
    title(sprintf('PM_{2.5} (\\mug/m^3) on %d, \\DeltaT=%d days, n=%d\nx=%0.0f km, y=%0.0f km',...        
        daysWHIdisp,deltaT,n,CTMlocs(idxMod(i),1)./1000,CTMlocs(idxMod(i),2)./1000));

    % add bins
    lenbins = length(perctile_data{idxMod(i)});
    for j = 1:lenbins
        plot([perctile_data{idxMod(i)}(j) perctile_data{idxMod(i)}(j)],get(gca,'YLim'),'r-');
    end

    % calculate means in each bin
    plot(mean_Mod{idxMod(i)},mean_Obs{idxMod(i)},'ro','MarkerFaceColor','r');

    % linear interpolation in each bin
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    interpx = linspace(xlim(1),xlim(2));
    if length(mean_Mod{idxMod(i)}) > 1
        interpy = interp1(mean_Mod{idxMod(i)},mean_Obs{idxMod(i)},interpx,'linear','extrap'); 
        plot(interpx,interpy,'r-');
    end

    % add modeled value
    plot([dailyCTMg(idxMod(i)) dailyCTMg(idxMod(i))],[ylim(1) meanGivMod(idxMod(i),1)], ...
        'c--','LineWidth',2);
    plot([xlim(1) dailyCTMg(idxMod(i))],[meanGivMod(idxMod(i),1) meanGivMod(idxMod(i),1)], ...
        'c--','LineWidth',2);

    % add different modeled values
    for j = 1:length(modplots)
        plot([modplots(j) modplots(j)],[ylim(1) meanGivMod(idxMod(i),j+1)], ...
        'k--','LineWidth',2);
        plot([xlim(1) modplots(j)],[meanGivMod(idxMod(i),j+1) meanGivMod(idxMod(i),j+1)], ...
            'k--','LineWidth',2);
    end

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/Scurve_%d_%d_deNazelle.pdf',idxMod(i),daysWHIdisp));

end

end