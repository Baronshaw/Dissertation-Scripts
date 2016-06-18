function [] = dispRF_2001(dec2simu)
% this function will plot the randomly generated fields for the first days
% of the month for 2001

if nargin < 1, dec2simu = 1; end % percentile to plot

% load data
load(sprintf('RFlambda1_2001_decile_%0.2d.mat',dec2simu));
FcGen = FcGen + 16;

% days to plot
unidays = datenum(2001,1:12,1);

% looping through each day
for i = 1:length(unidays)
    disp(i);      
    idx = cGen(:,3) == unidays(i);
    cGenSub = cGen(idx,1:2);
    FcGenSub = FcGen(idx);
    
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1
    
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cGenSub,FcGenSub,lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [5 25];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Randomly Generated Field of \\lambda_1 on %s for Decile %0.2d',datestr(unidays(i)),dec2simu));

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
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/genlambda1_%s_RF_%0.2d.png',datestr(unidays(i)),dec2simu));
      
end
close all

end