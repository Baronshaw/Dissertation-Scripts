function [] = plotDominance()
% this function will plot the variance versus range of given models across
% the different mean trends (aka offsets)

% models and the mean trends
load('../matfiles/covmod_modinfo_joint.mat');
meanNow = {'short','intermediate','long','very_long'};

for i = 1:length(modstr) % there's a figure for every model

    figure; hold on;
    markerOrder = {'o' ; 's' ; 'x' ; '*' };
    plotcolor = cool(6);
    
    % space
    load('../matfiles/covmod_modinfo_joint.mat');
    for j = 1:length(meanNow)
        load(sprintf('../matfiles/covmod_r_%s_%s_joint.mat',meanNow{j},modstr{j,i}));
        
        % getting weighted range
        if sum(strcmp(coeffnames(f),'ar2')) > 0
            ars(j) = f.alp*f.ar1 + (1-f.alp)*f.ar2; % range
        else
            ars(j) = f.ar1;
        end

        Crtests(j) = Crtest(1); % variance        
        subplot(2,1,1); hold on;
        plot(Crtests(j),ars(j),markerOrder{j},'Color',plotcolor(j+1,:),...
            'LineWidth',2);   
    end
    title(sprintf('variance vs range in space/time for %s',modstr{1,i}));
    legend(meanNow,'Location','best');
    dummy = get(gca,'YLim'); set(gca,'YLim',[0 dummy(2)]);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    ylabel(sprintf('spatial range (km)'));
    xlabel('variance r (km^{2})');
    
    % time
    for j = 1:length(meanNow)
        load(sprintf('../matfiles/covmod_t_%s_%s_joint.mat',meanNow{j},modstr{1,i}));
        
        % getting weighted range
        if sum(strcmp(coeffnames(f),'at2')) > 0
            ats(j) = f.alp*f.at1 + (1-f.alp)*f.at2; % range
        else
            ats(j) = f.at1;
        end

        Cttests(j) = Crtest(1); % variance        
        subplot(2,1,2); hold on;
        plot(Cttests(j),ats(j),markerOrder{j},'Color',plotcolor(j+1,:),...
            'LineWidth',2);   
    end
    legend(meanNow,'Location','best');
    ylabel(sprintf('temporal range (days)'));
    xlabel('variance t (days^{2})');
    
    % save figure
    h = gcf;
    print(h,'-painters','-dpdf','-r600',sprintf('../plots_covModel/varvrange_%s_joint.pdf',modstr{1,i}));

end

close all;

end