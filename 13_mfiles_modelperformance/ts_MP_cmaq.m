function [] = ts_MP_cmaq(numplot)
% create time series of bias and error across different locations

if nargin < 1, numplot = 10; end % plot time series of 10 locations

% load files
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
load('matfiles/traditional_performance.mat');
load('matfiles/info.mat');

% over all bias/error
figure; hold on;
plot(unidayz,mModD,'b-',unidayz,mObsD,'r-',unidayz,BiasD,'g-',unidayz,ErrD,'k-');
legend('Mod','Obs','Bias','Err');
ylabel('PM2.5 concentration');
datetick('x',2);
title('Overall Mod and Obs for 2001'); 

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','figures/overall_TS.png');

% by locations
dayz = datenum(yr,mo,da);
uniidx = randsample(length(uniPMM),numplot);
for i = 1:numplot
    
    figure; hold on;
    idx = uniPMM(uniidx(i),1) == coordObs(:,1) & uniPMM(uniidx(i),2) == coordObs(:,2);
    plot(dayz(idx),Mod(idx),'b-',dayz(idx),Obs(idx),'r-', ...
        dayz(idx),Mod(idx)-Obs(idx),'g-',dayz(idx),abs(Mod(idx)-Obs(idx)),'k-');
    legend('Mod','Obs','Bias','Err');
    ylabel('PM2.5 concentration');
    datetick('x',2);
    title(sprintf('Mod and Obs for (%.0f km %0.f km) for 2001',uniPMM(uniidx(i),1),uniPMM(uniidx(i),2))); 
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/overall_TS%d.png',i));
       
end

close all;

end