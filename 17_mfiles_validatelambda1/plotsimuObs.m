function [] = plotsimuObs(yr2simu)
% this function will plot the simulated observed values

if nargin < 1, yr2simu = 2001; end

% load data
load(sprintf('simulateObs_%d.mat',yr2simu)); 
% use Mod, mObs_calc, Obssimu_calc

yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
yrmoda = datenum(yr,mo,da);

% days to plot
unidays = datenum(yr2simu,1:12,1);

% looping through each day for the mod values at obs locations
for i = 1:length(unidays) 
    
    figure; hold on;
    
    idx = yrmoda == unidays(i);
    yrmodaSub = yrmoda(idx);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};   
    cax = [0 15];
    colorplot(coordObs(idx,:),Mod(idx),'redpink',Property,Value,cax); 
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Mod on %s at Obs location',datestr(unidays(i))));

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
        {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    for k = 1:length(allstates)
        cd ../09_mfiles_projections
        states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
        cd ../17_mfiles_validatelambda1
        plot(states(:,1),states(:,2),'k-');
    end 

    % save figure
    set(gcf,'Position',[0 0 800 500]);       
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 500]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/recalcMod_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

% looping through each day for the recalculated obs for the CMAQ
% modeled value
for i = 1:length(unidays) 
    
    figure; hold on;
    
    idx = yrmoda == unidays(i);
    yrmodaSub = yrmoda(idx);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};   
    cax = [0 15];
    colorplot(coordObs(idx,:),mObs_calc(idx),'redpink',Property,Value,cax); 
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Recalculated Obs on %s given CMAQ modeled value',datestr(unidays(i))));

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
        {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    for k = 1:length(allstates)
        cd ../09_mfiles_projections
        states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
        cd ../17_mfiles_validatelambda1
        plot(states(:,1),states(:,2),'k-');
    end 

    % save figure
    set(gcf,'Position',[0 0 800 500]);       
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 500]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/recalcObs_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

% looping through each day for the simulated obs for the CMAQ
% modeled value
for i = 1:length(unidays) 
    
    figure; hold on;
    
    idx = yrmoda == unidays(i);
    yrmodaSub = yrmoda(idx);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};   
    cax = [0 15];
    colorplot(coordObs(idx,:),Obssimu_calc(idx),'redpink',Property,Value,cax); 
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Simulated Obs on %s given CMAQ modeled value',datestr(unidays(i))));

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
        {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    for k = 1:length(allstates)
        cd ../09_mfiles_projections
        states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
        cd ../17_mfiles_validatelambda1
        plot(states(:,1),states(:,2),'k-');
    end 

    % save figure
    set(gcf,'Position',[0 0 800 500]);       
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 500]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/SimuObs_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

end