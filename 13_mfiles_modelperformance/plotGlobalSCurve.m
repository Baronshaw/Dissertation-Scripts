function [] = plotGlobalSCurve
% this function will create a plot of the global SCurve

% load CMAQ paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));

% load data
load('matfiles/globalSCurve.mat');

% perciles 
perctile_data = prctile(Mod,linspace(0,100,numbins+1))';

figure; hold on;

plot(Mod,Obs,'b.');
plot([0 100],[0 100],'k--'); % note: these values are hard-coded
axis([0 100 0 100]); % note: these values are hard-coded
xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
title(sprintf('Global SCurve: PM_{2.5} (\\mug/m^3)'));

% add bins
for j = 1:numbins+1
    plot([perctile_data(j) perctile_data(j)],get(gca,'YLim'),'r-');
end

% calculate means in each bin
plot(mMod,mObs,'ro','MarkerFaceColor','r');

% linear interpolation in each bin
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
interpx = linspace(xlim(1),xlim(2));
interpy = interp1(mMod,mObs,interpx,'linear','extrap'); 
plot(interpx,interpy,'r-');

% add different modeled values
for j = 1:length(modplots)-1
    plot([modplots(j) modplots(j)],[ylim(1) mModGivMod(j)],'k--','LineWidth',2);
    plot([xlim(1) modplots(j)],[mModGivMod(j+1) mModGivMod(j+1)],'k--','LineWidth',2);
end

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpdf','-r600','figures/GlobalScurve.pdf');

end