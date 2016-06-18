function [] = plotForcedIso()
% this function will plot r2 from the cross validation results of the
% forced isolation between hard and soft data as a function of increasing
% radius

constant = 0; 
gauss = 1; 
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

forceddist = 0:100000:1000000;
yearz = 1999:2010;

% load results
for i = 1:length(forceddist)
    
    load(sprintf('Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_results.mat', ...
        '_nosoft',constr,gaussstr,floor(forceddist(i)./1000)));
    r2_yearz_nosoft{i,1} = r2_year';
    
    load(sprintf('Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_results.mat', ...
        '_soft',constr,gaussstr,floor(forceddist(i)./1000)));
    r2_yearz_soft{i,1} = r2_year';

end
r2_yearz_nosoft = cell2mat(r2_yearz_nosoft);
r2_yearz_soft = cell2mat(r2_yearz_soft);

temp = repmat(floor(forceddist/1000),length(n_sites_year),1);
allforceddist = temp(:);
allyearz = repmat(yearz',length(forceddist),1);

% plot for each year soft data exist (i.e. 01, 02, 05, 06, 07)
soft_years = [2001:2002 2005:2007];
for i = 1:length(soft_years)
    
    figure; hold on;
    idx = allyearz == soft_years(i);
    plot(floor(forceddist./1000),r2_yearz_nosoft(idx),'ro-', ...
        floor(forceddist./1000),r2_yearz_soft(idx),'bo-');
    legend('no soft','soft');
    xlim([0 1100]);
    ylim([0 1]);
    title(sprintf('r^{2} of increasing radius for %d',soft_years(i)));
    xlabel('km');
    ylabel(sprintf('r^{2}'));
    
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('plot_forcediso_%d.pdf',soft_years(i)));
%     print(gcf,'-painters','-dpng','-r600',sprintf('plot_forcediso_%d.png',soft_years(i)));
%     
end

% plot for all the years together
c = cell(length(soft_years)*2,1);
stylezh = {'ro-','r.-','-r+','-r*','rx-'};
stylezs = {'bo-','b.-','-b+','-b*','bx-'};
figure; hold on;
for i = 1:length(soft_years)
    
    idx = allyearz == soft_years(i);
    plot(floor(forceddist./1000),r2_yearz_nosoft(idx),stylezh{i}, ...
        floor(forceddist./1000),r2_yearz_soft(idx),stylezs{i});
    c{i*2-1} = sprintf('no soft %d',soft_years(i));
    c{i*2} = sprintf('soft %d',soft_years(i));     
    
end

legend(c);
xlim([0 1100]);
ylim([0 1]);
title('r^{2} of increasing radius for all years');
xlabel('km');
ylabel(sprintf('r^{2}'));

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpdf','-r600','plot_forcediso_allyears.pdf');
print(gcf,'-painters','-dpng','-r600','plot_forcediso_allyears.png');

close all

end