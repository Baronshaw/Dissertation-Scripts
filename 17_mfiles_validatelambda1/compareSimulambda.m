function [] = compareSimulambda()
% this function will compare the recalculated/simulated lambda1 to the 
% estimated lambda1

% load recalculated lambda1 values
load('recalculatelambdaGivMod.mat');
% css, modplots, lambda1_recalcMod, lambda2_recalcModGrid, modplots, 
% lambda2_recalcsimuMod, lambda1_recalcsimuModGrid
load('recalculate_dailyvCTMorder.mat');
% dailyvCTMorder

% days to plot
unidays = datenum(2001,1:12,1);

% % looping through each day for the recalculated lambda1
% for i = 1:length(modplots)
%     for j = 1:length(unidays)
%         disp(j);      
%         idx = css(:,3) == unidays(j);
%         cssSub = css(idx,1:2);
% 
%         % country outline
%         cd ../09_mfiles_projections
%         load('USAcontiguous.mat');
%         plotax = ell2lambertcc([x,y],'whiproj2001');
%         cd ../17_mfiles_validatelambda1
% 
%         lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%         [xg yg Zg] = plotField(cssSub,lambda1_recalcMod(idx,i),lax,[plotax(:,1) plotax(:,2)],redpink);
%         cax = [0 30];
%         caxis(cax);
%         colorbar;
%         axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%         % setting axis
%         set(gca,'XTickLabel','')
%         set(gca,'YTickLabel','')
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         title(sprintf('Recalculated Field of \\lambda_1 on %s given modeled %d',datestr(unidays(j)),modplots(i)));
% 
%         % overlaying the states
%         load('../09_mfiles_projections/USAstates5.mat');
%         allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
%             {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
%         for k = 1:length(allstates)
%             cd ../09_mfiles_projections
%             states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
%             cd ../17_mfiles_validatelambda1
%             plot(states(:,1),states(:,2),'k-');
%         end 
% 
%         % save figure
%         set(gcf,'Position',[0 0 800 500]);       
%         set(gca,'YTickLabel',get(gca,'YTick')/1000);
%         set(gca,'XTickLabel',get(gca,'XTick')/1000);
%         set(gcf,'PaperUnits','inches');    
%         set(gcf,'PaperPosition',[0 0 800 500]./100);
%         set(gcf,'PaperPositionMode','manual');
%         print(gcf,'-painters','-dpng','-r600',sprintf('figures/recalclambda1_%s_mod_%0.2d.png',datestr(unidays(j)),modplots(i)));
% 
%     end
%     close all
% end
% 
% % looping through each day for the recalculated lambda2
% for i = 1:length(modplots)
%     for j = 1:length(unidays)
%         disp(j);      
%         idx = css(:,3) == unidays(j);
%         cssSub = css(idx,1:2);
% 
%         % country outline
%         cd ../09_mfiles_projections
%         load('USAcontiguous.mat');
%         plotax = ell2lambertcc([x,y],'whiproj2001');
%         cd ../17_mfiles_validatelambda1
% 
%         lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%         [xg yg Zg] = plotField(cssSub,lambda2_recalcMod(idx,i),lax,[plotax(:,1) plotax(:,2)],redpink);
%         cax = [0 50];
%         caxis(cax);
%         colorbar;
%         axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%         % setting axis
%         set(gca,'XTickLabel','')
%         set(gca,'YTickLabel','')
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         title(sprintf('Recalculated Field of \\lambda_2 on %s given modeled %d',datestr(unidays(j)),modplots(i)));
% 
%         % overlaying the states
%         load('../09_mfiles_projections/USAstates5.mat');
%         allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
%             {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
%         for k = 1:length(allstates)
%             cd ../09_mfiles_projections
%             states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
%             cd ../17_mfiles_validatelambda1
%             plot(states(:,1),states(:,2),'k-');
%         end 
% 
%         % save figure
%         set(gcf,'Position',[0 0 800 500]);       
%         set(gca,'YTickLabel',get(gca,'YTick')/1000);
%         set(gca,'XTickLabel',get(gca,'XTick')/1000);
%         set(gcf,'PaperUnits','inches');    
%         set(gcf,'PaperPosition',[0 0 800 500]./100);
%         set(gcf,'PaperPositionMode','manual');
%         print(gcf,'-painters','-dpng','-r600',sprintf('figures/recalclambda2_%s_mod_%0.2d.png',datestr(unidays(j)),modplots(i)));
% 
%     end
%     close all
% end
% 
% % looping through each day for the recalculated lambda1 for the CMAQ
% % modeled value
% for i = 1:length(unidays)    
%     idx = css(:,3) == unidays(i);
%     cssSub = css(idx,1:2);
% 
%     % country outline
%     cd ../09_mfiles_projections
%     load('USAcontiguous.mat');
%     plotax = ell2lambertcc([x,y],'whiproj2001');
%     cd ../17_mfiles_validatelambda1
% 
%     lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%     [xg yg Zg] = plotField(cssSub,lambda1_recalcModGrid(idx),lax,[plotax(:,1) plotax(:,2)],redpink);
%     cax = [0 15];
%     caxis(cax);
%     colorbar;
%     axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%     % setting axis
%     set(gca,'XTickLabel','')
%     set(gca,'YTickLabel','')
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title(sprintf('Recalculated Field of \\lambda_1 on %s given CMAQ modeled value',datestr(unidays(i))));
% 
%     % overlaying the states
%     load('../09_mfiles_projections/USAstates5.mat');
%     allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
%         {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
%     for k = 1:length(allstates)
%         cd ../09_mfiles_projections
%         states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
%         cd ../17_mfiles_validatelambda1
%         plot(states(:,1),states(:,2),'k-');
%     end 
% 
%     % save figure
%     set(gcf,'Position',[0 0 800 500]);       
%     set(gca,'YTickLabel',get(gca,'YTick')/1000);
%     set(gca,'XTickLabel',get(gca,'XTick')/1000);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 500]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpng','-r600',sprintf('figures/recalclambda1_%s_CMAQ_mod.png',datestr(unidays(i))));
% 
% end
% close all
% 
% % looping through each day for the recalculated lambda2 for the CMAQ
% % modeled value
% for i = 1:length(unidays)    
%     idx = css(:,3) == unidays(i);
%     cssSub = css(idx,1:2);
% 
%     % country outline
%     cd ../09_mfiles_projections
%     load('USAcontiguous.mat');
%     plotax = ell2lambertcc([x,y],'whiproj2001');
%     cd ../17_mfiles_validatelambda1
% 
%     lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%     [xg yg Zg] = plotField(cssSub,lambda2_recalcModGrid(idx),lax,[plotax(:,1) plotax(:,2)],redpink);
%     cax = [0 15];
%     caxis(cax);
%     colorbar;
%     axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%     % setting axis
%     set(gca,'XTickLabel','')
%     set(gca,'YTickLabel','')
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title(sprintf('Recalculated Field of \\lambda_2 on %s given CMAQ modeled value',datestr(unidays(i))));
% 
%     % overlaying the states
%     load('../09_mfiles_projections/USAstates5.mat');
%     allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
%         {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
%     for k = 1:length(allstates)
%         cd ../09_mfiles_projections
%         states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
%         cd ../17_mfiles_validatelambda1
%         plot(states(:,1),states(:,2),'k-');
%     end 
% 
%     % save figure
%     set(gcf,'Position',[0 0 800 500]);       
%     set(gca,'YTickLabel',get(gca,'YTick')/1000);
%     set(gca,'XTickLabel',get(gca,'XTick')/1000);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 500]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpng','-r600',sprintf('figures/recalclambda2_%s_CMAQ_mod.png',datestr(unidays(i))));
% 
% end
% close all
% 
% % looping through each day for the simulated lambda1
% for i = 1:length(modplots)
%     for j = 1:length(unidays)
%         disp([i,j]);      
%         idx = css(:,3) == unidays(j);
%         cssSub = css(idx,1:2);
% 
%         % country outline
%         cd ../09_mfiles_projections
%         load('USAcontiguous.mat');
%         plotax = ell2lambertcc([x,y],'whiproj2001');
%         cd ../17_mfiles_validatelambda1
% 
%         lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%         [xg yg Zg] = plotField(cssSub,real(lambda1_recalcsimuMod(idx,i)),lax,[plotax(:,1) plotax(:,2)],redpink);
%         cax = [0 30];
%         caxis(cax);
%         colorbar;
%         axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%         % setting axis
%         set(gca,'XTickLabel','')
%         set(gca,'YTickLabel','')
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         title(sprintf('Simulated Field of \\lambda_1 on %s given modeled %d',datestr(unidays(j)),modplots(i)));
% 
%         % overlaying the states
%         load('../09_mfiles_projections/USAstates5.mat');
%         allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
%             {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
%         for k = 1:length(allstates)
%             cd ../09_mfiles_projections
%             states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
%             cd ../17_mfiles_validatelambda1
%             plot(states(:,1),states(:,2),'k-');
%         end 
% 
%         % save figure
%         set(gcf,'Position',[0 0 800 500]);       
%         set(gca,'YTickLabel',get(gca,'YTick')/1000);
%         set(gca,'XTickLabel',get(gca,'XTick')/1000);
%         set(gcf,'PaperUnits','inches');    
%         set(gcf,'PaperPosition',[0 0 800 500]./100);
%         set(gcf,'PaperPositionMode','manual');
%         print(gcf,'-painters','-dpng','-r600',sprintf('figures/simulambda1_%s_mod_%0.2d.png',datestr(unidays(j)),modplots(i)));
% 
%     end
%     close all
% end
% 
% % looping through each day for the simulated lambda2
% for i = 1:length(modplots)
%     for j = 1:length(unidays)
%         disp([i,j]);      
%         idx = css(:,3) == unidays(j);
%         cssSub = css(idx,1:2);
% 
%         % country outline
%         cd ../09_mfiles_projections
%         load('USAcontiguous.mat');
%         plotax = ell2lambertcc([x,y],'whiproj2001');
%         cd ../17_mfiles_validatelambda1
% 
%         lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
%         [xg yg Zg] = plotField(cssSub,real(lambda2_recalcsimuMod(idx,i)),lax,[plotax(:,1) plotax(:,2)],redpink);
%         cax = [0 30];
%         caxis(cax);
%         colorbar;
%         axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);
% 
%         % setting axis
%         set(gca,'XTickLabel','')
%         set(gca,'YTickLabel','')
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         title(sprintf('Simulated Field of \\lambda_2 on %s given modeled %d',datestr(unidays(j)),modplots(i)));
% 
%         % overlaying the states
%         load('../09_mfiles_projections/USAstates5.mat');
%         allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
%             {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
%         for k = 1:length(allstates)
%             cd ../09_mfiles_projections
%             states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
%             cd ../17_mfiles_validatelambda1
%             plot(states(:,1),states(:,2),'k-');
%         end 
% 
%         % save figure
%         set(gcf,'Position',[0 0 800 500]);       
%         set(gca,'YTickLabel',get(gca,'YTick')/1000);
%         set(gca,'XTickLabel',get(gca,'XTick')/1000);
%         set(gcf,'PaperUnits','inches');    
%         set(gcf,'PaperPosition',[0 0 800 500]./100);
%         set(gcf,'PaperPositionMode','manual');
%         print(gcf,'-painters','-dpng','-r600',sprintf('figures/simulambda2_%s_mod_%0.2d.png',datestr(unidays(j)),modplots(i)));
% 
%     end
%     close all
% end

% looping through each day for the simulated lambda1 for the CMAQ
% modeled value
for i = 1:length(modplots)    
    idx = css(:,3) == unidays(i);
    cssSub = css(idx,1:2);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cssSub,real(lambda1_recalcsimuModGrid(idx)),lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [0 15];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Simulated Field of \\lambda_1 on %s given CMAQ modeled value',datestr(unidays(i))));

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
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/simulambda1_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

% looping through each day for the simulated CMAQ-lambda1 for the CMAQ
% modeled value
for i = 1:length(modplots)    
    idx = css(:,3) == unidays(i);
    cssSub = css(idx,1:2);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cssSub,real(dailyCTMvorder(idx)-lambda1_recalcsimuModGrid(idx)),lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [0 15];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Simulated Field of CMAQ-\\lambda_1 on %s given CMAQ modeled value',datestr(unidays(i))));

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
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/CMAQ-simulambda1_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

% looping through each day for the simulated lambda2 for the CMAQ
% modeled value
for i = 1:length(modplots)    
    idx = css(:,3) == unidays(i);
    cssSub = css(idx,1:2);

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../17_mfiles_validatelambda1

    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(cssSub,real(lambda2_recalcsimuModGrid(idx)),lax,[plotax(:,1) plotax(:,2)],redpink);
    cax = [0 50];
    caxis(cax);
    colorbar;
    axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

    % setting axis
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(sprintf('Simulated Field of \\lambda_2 on %s given CMAQ modeled value',datestr(unidays(i))));

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
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/simulambda2_%s_CMAQ_mod.png',datestr(unidays(i))));

end
close all

end