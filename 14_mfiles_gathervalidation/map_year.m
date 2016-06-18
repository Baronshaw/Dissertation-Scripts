function [] = map_year(forceddist)
% this function will plot for each statistic, values of increasing years
% for krig and RAMP

if nargin < 1, forceddist = 0; end

% measure names/values
strnm = {'num of paired modeled and obs';'mean obs';'mean mod';'mean error'; ...
    'mean normalized error';'normalized mean error';'fractional error';'mean absolute error'; ...
    'mean normalized absolute error';'normalized mean absolute error';'fractional absolute error'; ...
    'correlation';'correlation squared';'standard error';'mean squared error'; ...
    'root mean squared error';'normalized root mean squared error'; ...
    'mean error DIV standard error';'mean error squared DIV mean squared error'; ...
    'variance of errors DIV mean squared error';'beta1';'variance of  obs';'varianced of mod'};

for i = 1:length(strnm)
    
    load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','krig',floor(forceddist./1000)));
    toplot_krig = years{i};
    load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','RAMP',floor(forceddist./1000)));
    toplot_RAMP = years{i};
    
    figure; hold on;
    plot(uniyr,toplot_krig,'bo-',uniyr,toplot_RAMP,'ro-');
    title(sprintf('%s across years for %d km',strnm{i},floor(forceddist./1000)));
    xlabel('years');
    ylabel(sprintf('%s',strnm{i}));
    legend('krig','RAMP','Location','best');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/years/%s_%dkm.png', ...
        strnm{i},floor(forceddist./1000)));
        
end
close all;

end