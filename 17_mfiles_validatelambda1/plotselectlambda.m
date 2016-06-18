function [] = plotselectlambda(yr2select)
% this function will plot lambda1 across select days for fixed modeled
% values

if nargin < 1, yr2select = 2001; end

% load data
load(sprintf('PM2p5_meanGivMod_yr%d.mat',yr2select));

% days to plot
unidays = datenum(yr2select,1:12,1);

% looping through each day for lambda1
for i = 1:length(unidays)
    [a b] = size(modplots_all);
    for j = 1:b
        disp(i);      
        idx = css(:,3) == unidays(i);
        cssSub = css(idx,1:2);
        meanGivMod_allSub = meanGivMod_all(idx,j);

        % country outline
        cd ../09_mfiles_projections
        load('USAcontiguous.mat');
        plotax = ell2lambertcc([x,y],'whiproj2001');
        cd ../17_mfiles_validatelambda1

        lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
        [xg yg Zg] = plotField(cssSub,meanGivMod_allSub,lax,[plotax(:,1) plotax(:,2)],redpink);
        cax = [0 30];
        caxis(cax);
        colorbar;
        axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

        % setting axis
        set(gca,'XTickLabel','')
        set(gca,'YTickLabel','')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        title(sprintf('Selected \\lambda_1 on %s gived modeled %d',datestr(unidays(i)),modplots_all(1,j)));

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
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/selectlambda1_%s_mod_%0.2d.png',datestr(unidays(i)),modplots_all(1,j)));

    end 
    close all
end

% looping through each day for lambda2
for i = 1:length(unidays)
    [a b] = size(modplots_all);
    for j = 1:b
        disp(i);      
        idx = css(:,3) == unidays(i);
        cssSub = css(idx,1:2);
        varGivMod_allSub = varGivMod_all(idx,j);

        % country outline
        cd ../09_mfiles_projections
        load('USAcontiguous.mat');
        plotax = ell2lambertcc([x,y],'whiproj2001');
        cd ../17_mfiles_validatelambda1

        lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
        [xg yg Zg] = plotField(cssSub,varGivMod_allSub,lax,[plotax(:,1) plotax(:,2)],redpink);
        cax = [0 50];
        caxis(cax);
        colorbar;
        axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

        % setting axis
        set(gca,'XTickLabel','')
        set(gca,'YTickLabel','')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        title(sprintf('Selected \\lambda_2 on %s gived modeled %d',datestr(unidays(i)),modplots_all(1,j)));

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
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/selectlambda2_%s_mod_%0.2d.png',datestr(unidays(i)),modplots_all(1,j)));

    end 
    close all
end

% looping through each day and plotting lambda1 given CMAQ concentration 
% for select days 
for i = 1:length(unidays)
    disp(i);      
    idx = css(:,3) == unidays(i);
    cssSub = css(idx,1:2);
    meanGivCMAQ_allSub = meanGivCMAQ_all(idx);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cssSub,meanGivCMAQ_allSub,lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [0 15];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Selected \\lambda_1 on %s gived CMAQ modeled',datestr(unidays(i))));

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
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/selectlambda1_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

% looping through each day and plotting lambda2 given CMAQ concentration 
% for select days 
for i = 1:length(unidays)
    disp(i);      
    idx = css(:,3) == unidays(i);
    cssSub = css(idx,1:2);
    varGivCMAQ_allSub = varGivCMAQ_all(idx);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cssSub,varGivCMAQ_allSub,lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [0 50];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Selected \\lambda_2 on %s gived CMAQ modeled',datestr(unidays(i))));

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
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/selectlambda2_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

end