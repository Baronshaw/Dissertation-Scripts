function [] = SCurve_will()
% this function will investigate SCurves at space/time locations in which a
% nearby location has a very different value for a given model performance
% statistics. These plots will demonstrate if a model performance
% statistics is an artifact of the SCurve (i.e. pulling the nearest three
% stations together).

%%% plotting ME2/MSE with select locations noted

% load data
load('matfiles/traditional_performance_grid.mat')
yr = floor(yrNday./10^4); uniyr = unique(yr);
mo = floor((yrNday - yr*10^4)./10^2);
da = yrNday - yr*10^4 - mo.*10^2;
yrznum = datevec(datenum(uniyr,7,1)); [u v] = size(yrznum);

% measure names/values
strnm = { 'mean bias DIV standard bias' };
valplot = { m2DmsBias }; 

% % loop through each day on the month
% for i = 1:length(valplot)
%     for j = 1:u
%         disp([i j]);
%         idx = yrznum(j,1) == yr & yrznum(j,2) == mo & yrznum(j,3) == da;
%         figure; hold on;
% 
%         load('../09_mfiles_projections/USAcontiguous.mat');
%         cd ../09_mfiles_projections
%         plotax = ell2lambertcc([x,y],'whiproj2001');
%         cd ../13_mfiles_modelperformance
% 
%         % plotting 
%         lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%         [xg yg Zg] = plotField(CTMlocs,valplot{i}(:,idx),lax,[plotax(:,1) plotax(:,2)]);
%         caxis([prctile(valplot{i}(:),5) prctile(valplot{i}(:),90)]);   
%         colorbar;
%         axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%         % setting axis
%         xlabel('km');
%         ylabel('km');
% 
%         % overlaying the states
%         load('../09_mfiles_projections/USAstates5.mat');
%         for k = 1:length(X)
%             cd ../09_mfiles_projections
%             states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
%             cd ../13_mfiles_modelperformance
%             plot(states(:,1),states(:,2),'k-');
%         end
% 
%         % title 
%         title(sprintf('%s of PM_{2.5} (\\mug/m^3) on %s',strnm{i},datestr(yrznum(j,:))));    
% 
%         % save figure 
%         set(gcf,'Position',[0 0 800 600]);
%         set(gcf,'PaperUnits','inches');    
%         set(gcf,'PaperPosition',[0 0 800 600]./100);
%         set(gcf,'PaperPositionMode','manual');
%         
%     end
% end
% 
% finding the exact locations to do the comparisons
SCurveLoc = { [-475000,470000],[-650000,480000] ; [1100000,-780000],[975000,-700000] ; ...
    [-170000,-1350000],[-170000,-1200000] ; [-550000,-1040000],[-800000,-940000] ; ...
    [-1425000,-750000],[-1425000,-600000] ; [-1950000,-400000],[-1800000,-400000] };
% for i = 1:length(SCurveLoc)
%     text(SCurveLoc{i,1}(1),SCurveLoc{i,1}(2),num2str(i),'Color','yellow','FontSize',15);
%     text(SCurveLoc{i,2}(1),SCurveLoc{i,2}(2),num2str(i),'Color','yellow','FontSize',15);
% end
% set(gca,'XTickLabel',get(gca,'XTick')/1000);
% set(gca,'YTickLabel',get(gca,'YTick')/1000);
% print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_%s_ExploreSCurve.png',strnm{1},datestr(yrznum(1,:))));
% 
%%% locations overlaid with monitoring locations
figure; hold on;
load('../09_mfiles_projections/USAcontiguous.mat');
cd ../09_mfiles_projections
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../13_mfiles_modelperformance

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for k = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
    cd ../13_mfiles_modelperformance
    plot(states(:,1),states(:,2),'k-');
end

load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
unimon = unique(coordObs,'rows');
plot(unimon(:,1),unimon(:,2),'k.');

load('temp.mat');
for i = 1:length(SCurveLoc)
    text(SCurveLoc{i,1}(1),SCurveLoc{i,1}(2),num2str(i),'Color','blue','FontSize',15);
    plot(closeMS{i}{1}(:,1),closeMS{i}{1}(:,2),'c.');
    text(SCurveLoc{i,2}(1),SCurveLoc{i,2}(2),num2str(i),'Color','blue','FontSize',15);
    plot(closeMS{i}{2}(:,1),closeMS{i}{2}(:,2),'r.');
end

% save figure 
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
print(gcf,'-painters','-dpng','-r600','figures/Monitoring_ExploreSCurve.png');

%%% plotting these SCurves
for i = 1:length(SCurveLoc)
    for j = 1:12
        daysWHIdisp = 2001*10^4 + j*10^2 + 1;
        plotSCurve([2001 j 01],daysWHIdisp,SCurveLoc{i,1});
        plotSCurve([2001 j 01],daysWHIdisp,SCurveLoc{i,2});
    end
end

% %%% time series for each of these locations
% valplot = { m2DmsBias };
% for i = 1:length(SCurveLoc)
%     figure; hold on;
%     finddist = sqrt( (SCurveLoc{i,1}(1)-CTMlocs(:,1)).^2 + (SCurveLoc{i,1}(2)-CTMlocs(:,2)).^2 );
%     idx = finddist == min(finddist);
%     finddist = sqrt( (SCurveLoc{i,2}(1)-CTMlocs(:,1)).^2 + (SCurveLoc{i,2}(2)-CTMlocs(:,2)).^2 );
%     idx2 = finddist == min(finddist);
%     plot(datenum([yr mo da]),valplot{1}(idx,:),'b-',datenum([yr mo da]),valplot{1}(idx2,:),'r-');
%     set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd mmm'));
%     title(sprintf('TS of %s for location\n(%d km %d km) and (%d km %d km)', ...
%         strnm{1},floor(SCurveLoc{i,1}(1)/1000),floor(SCurveLoc{i,1}(1)/1000), ...
%         floor(SCurveLoc{i,2}(1)/1000),floor(SCurveLoc{i,2}(1)/1000) ));
%     print(gcf,'-painters','-dpng','-r600',sprintf('figures/TS_x%dkm_y%dkm_x%dkm_y%dkm.png', ...
%         floor(SCurveLoc{i,1}(1)/1000),floor(SCurveLoc{i,1}(1)/1000), ...
%         floor(SCurveLoc{i,2}(1)/1000),floor(SCurveLoc{i,2}(1)/1000) ));
% end
% close all;
% 
% %%% for each location, histogram of distances to closest monitors
% for i = 1:length(SCurveLoc)
%     figure; hold on;
%     finddist = sqrt( (SCurveLoc{i,1}(1)-unimon(:,1)).^2 + (SCurveLoc{i,1}(2)-unimon(:,2)).^2 );
%     finddist2 = sqrt( (SCurveLoc{i,2}(1)-unimon(:,1)).^2 + (SCurveLoc{i,2}(2)-unimon(:,2)).^2 );
%     
%     subplot(2,1,1);
%     hist(finddist./1000,100);
%     title(sprintf('Distances to nearest montiors for location\n(%d km %d km)', ...
%         floor(SCurveLoc{i,1}(1)/1000),floor(SCurveLoc{i,1}(1)/1000)));
%     subplot(2,1,2);
%     hist(finddist2./1000,100);    
%     title(sprintf('(%d km %d km)',floor(SCurveLoc{i,2}(1)/1000),floor(SCurveLoc{i,2}(1)/1000) ));
%     
%     print(gcf,'-painters','-dpng','-r600',sprintf('figures/HistDist2Mon_x%dkm_y%dkm_x%dkm_y%dkm.png', ...
%         floor(SCurveLoc{i,1}(1)/1000),floor(SCurveLoc{i,1}(1)/1000), ...
%         floor(SCurveLoc{i,2}(1)/1000),floor(SCurveLoc{i,2}(1)/1000) ));
%     
% end
% close all;

end