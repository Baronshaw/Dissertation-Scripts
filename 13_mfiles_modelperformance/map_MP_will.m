function [] = map_MP_will()
% maps to show about model performance

% load data
load('matfiles/traditional_performance_grid.mat')
yr = floor(yrNday./10^4); uniyr = unique(yr);
mo = floor((yrNday - yr*10^4)./10^2);
da = yrNday - yr*10^4 - mo.*10^2;
yrznum = datevec(datenum(uniyr,1:12,1));

% % average over season
% [r c] = size(num);
% m2DmsBias_season = NaN*ones(r,4);
% s2DmsBias_season = NaN*ones(r,4);
% for i = 1:4
%     
%     if i == 1
%         idx = mo == 1 | mo == 2 | mo == 12;
%     elseif i == 2
%         idx = mo == 3 | mo == 4 | mo == 5;
%     elseif i == 3
%         idx = mo == 6 | mo == 7 | mo == 8;
%     else
%         idx = mo == 9 | mo == 10 | mo == 11;
%     end
% 
%     m2DmsBias_season(:,i) = nanmean(m2DmsBias(:,idx),2);
%     s2DmsBias_season(:,i) = nanmean(s2DmsBias(:,idx),2);
%     
% end
% 
% % average over the year
% m2DmsBias_year = nanmean(m2DmsBias,2);
% s2DmsBias_year = nanmean(s2DmsBias,2);
% 
% % save
% save('matfiles/maps_will.mat', ...
%     'm2DmsBias','s2DmsBias','m2DmsBias_season','s2DmsBias_season', ...
%     'm2DmsBias_year','s2DmsBias_year','CTMlocs','yrNday'); 
% load('matfiles/maps_will.mat');
% 
% % measure names/values
% strnm = { 'mean bias squared DIV mean squared bias - Winter' ; 'variance of bias DIV mean squared bias - Winter' ; ...
%     'mean bias squared DIV mean squared bias - Spring' ; 'variance of bias DIV mean squared bias - Spring' ; ...
%     'mean bias squared DIV mean squared bias - Summer' ; 'variance of bias DIV mean squared bias - Summer' ; ...
%     'mean bias squared DIV mean squared bias - Fall' ; 'variance of bias DIV mean squared bias - Fall' ; ...
%     'mean bias squared DIV mean squared bias - Annual' ; 'variance of bias DIV mean squared bias - Annual' };
% a = m2DmsBias_season(:,1); b = s2DmsBias_season(:,1);
% c = m2DmsBias_season(:,2); d = s2DmsBias_season(:,2); 
% e = m2DmsBias_season(:,3); f = s2DmsBias_season(:,3);
% g = m2DmsBias_season(:,4); h = s2DmsBias_season(:,4);
% valplot = { a ; b ; c ; d ; e ; f ; g ; h ; m2DmsBias_year ; s2DmsBias_year }; 
% 
% % plotting
% for i = 1:length(valplot)
% 
%     figure; hold on;
% 
%     load('../09_mfiles_projections/USAcontiguous.mat');
%     cd ../09_mfiles_projections
%     plotax = ell2lambertcc([x,y],'whiproj2001');
%     cd ../13_mfiles_modelperformance
% 
%     % plotting 
%     lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%     [xg yg Zg] = plotField(CTMlocs,valplot{i},lax,[plotax(:,1) plotax(:,2)]);
%     caxis([0 1]);   
%     colorbar;
%     axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%     % setting axis
%     set(gca,'XTickLabel',get(gca,'XTick')/1000);
%     set(gca,'YTickLabel',get(gca,'YTick')/1000);
%     xlabel('km');
%     ylabel('km');
% 
%     % overlaying the states
%     load('../09_mfiles_projections/USAstates5.mat');
%     for k = 1:length(X)
%         cd ../09_mfiles_projections
%         states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
%         cd ../13_mfiles_modelperformance
%         plot(states(:,1),states(:,2),'k-');
%     end
% 
%     % title 
%     title(sprintf('%s',strnm{i}));    
% 
%     % save figure 
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_grid.png',strnm{i}));
% 
% end
% 
% close all;
% 
%%% look at increasing modeled values
tic
load('matfiles/traditional_performance_extendbin.mat');
toc

% measure names/values
valplot = { m2DmsBiasGivMod ; s2DmsBiasGivMod };
strnm = { 'mean bias squared DIV mean squared bias' ; 'variance of bias DIV mean squared bias' };
yrznum = [2001 1 15 ; 2001 6 15]; [u v] = size(yrznum);
yrznum = [2001 9 13]; [u v] = size(yrznum);

% loop through a few days 
for l = 1:length(valplot) % measure
    for i = 1:u % day
        disp(i);
        idx = yrNday == yrznum(i,1)*10000+yrznum(i,2)*100+yrznum(i,3);
        for j = 1:length(modplots) % modeled values

            figure; hold on;

            load('../09_mfiles_projections/USAcontiguous.mat');
            cd ../09_mfiles_projections
            plotax = ell2lambertcc([x,y],'whiproj2001');
            cd ../13_mfiles_modelperformance

            % plotting 
            lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
            [xg yg Zg] = plotField(CTMlocs,valplot{l}{idx}(:,j),lax,[plotax(:,1) plotax(:,2)]);
            caxis([0 1]);   
            colorbar;
            axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

            % setting axis
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            xlabel('km');
            ylabel('km');

            % overlaying the states
            load('../09_mfiles_projections/USAstates5.mat');
            for k = 1:length(X)
                cd ../09_mfiles_projections
                states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
                cd ../13_mfiles_modelperformance
                plot(states(:,1),states(:,2),'k-');
            end

            % title 
            title(sprintf('%s on %s, giv Mod = %d',strnm{l},datestr(datenum(yrznum(i,:))),modplots(j)));

            % save figure 
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/will_%s_%s_Mod%d_grid.png',strnm{l},datestr(datenum(yrznum(i,:))),modplots(j)));

        end
        close all;
    end
end

% time series
% pick 5 random locations
randidx = randsample(1:length(CTMlocs),5);
vals = [ -2000000 -350000 ; 
    1500000 -300000 ; 
    -1200000 200000 ; 
    1000000 -650000 ; 
    -100000 -1300000 ];
for i = 1:length(vals)
    distz = sqrt( (CTMlocs(:,1)-vals(i,1)).^2 + (CTMlocs(:,2)-vals(i,2)).^2 );
    fidx(i,1) = find(distz==min(distz));
end

figure; hold on;
% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for k = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
    cd ../13_mfiles_modelperformance
    plot(states(:,1),states(:,2),'k-');
end
plot(CTMlocs(fidx,1),CTMlocs(fidx,2),'bo');
% title 
title('Random Times Series');
% save figure 
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
% setting axis
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
xlabel('km');
ylabel('km');
print(gcf,'-painters','-dpng','-r600','figures/will_randlocs_grid.png');

% time series
% measure names/values
valplot = { m2DmsBias ; s2DmsBias };
for j = 1:length(fidx)
    figure; hold on;
    plot(yrNday,valplot{1}(fidx,:),'b-',yrNday,valplot{2}(fidx,:),'r-');
    legend('%% ME','%% SE');
    title(sprintf('TS for location (%f km %f km)',floor(CTMlocs(fidx,1)/1000),floor(CTMlocs(fidx,2))/1000));
end

end