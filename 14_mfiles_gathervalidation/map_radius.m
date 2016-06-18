function [] = map_radius()
% this function will plot for each statistic, values of increasing radius
% for krig and RAMP

% show all these variables for 2001
yrznum = 2001;

% measure names/values
radbuff = 0:100000:900000;
strnm = {'num of paired modeled and obs';'mean obs';'mean mod';'mean error'; ...
    'mean normalized error';'normalized mean error';'fractional error';'mean absolute error'; ...
    'mean normalized absolute error';'normalized mean absolute error';'fractional absolute error'; ...
    'correlation';'correlation squared';'standard error';'mean squared error'; ...
    'root mean squared error';'normalized root mean squared error'; ...
    'mean error DIV standard error';'mean error squared DIV mean squared error'; ...
    'variance of errors DIV mean squared error';'beta1';'variance of  obs';'varianced of mod'};

for i = 1:length(strnm)
    
    for j = 1:length(radbuff)
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','krig',floor(radbuff(j)./1000)));
        yrz = yrznum == uniyr;
        toplot_krig(j,1) = years{i}(yrz);
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','RAMP',floor(radbuff(j)./1000)));
        yrz = yrznum == uniyr;
        toplot_RAMP(j,1) = years{i}(yrz);
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','staticDS',floor(radbuff(j)./1000)));
        yrz = yrznum == uniyr;
        toplot_staticDS(j,1) = years{i}(yrz);
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','stDS_add_ind_muli_ind',floor(radbuff(j)./1000)));
        yrz = yrznum == uniyr;
        toplot_stDS_indind(j,1) = years{i}(yrz);
        load(sprintf('matfiles/allXval_LOO_%s_dist%d.mat','stDS_add_dyn_muli_ind',floor(radbuff(j)./1000)));
        yrz = yrznum == uniyr;
        toplot_stDS_dynind(j,1) = years{i}(yrz);
    end
    
    figure; hold on;
    plot(floor(radbuff./1000),toplot_krig,'bo-',floor(radbuff./1000),toplot_RAMP,'ro-', ...
        floor(radbuff./1000),toplot_staticDS,'co-',floor(radbuff./1000),toplot_stDS_indind,'ko-',...
        floor(radbuff./1000),toplot_stDS_dynind,'yo-');
    title(sprintf('%s of increasing radius for %d',strnm{i},yrznum));
    xlabel('km');
    ylabel(sprintf('%s',strnm{i}));
    legend('krig','RAMP','static DS','stDS add ind muli ind','stDS add dyn muli ind','Location','best');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/radius/%s_%d.png',strnm{i},yrznum));
        
end
close all;

end