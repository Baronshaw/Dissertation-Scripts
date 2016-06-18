function [] = plotCTM(daydisp,cax)
% this function will plot figures of the CTM data for a given day across
% the US

if nargin < 1, daydisp = [2001 7 1]; end
if nargin < 2, cax = [2 30]; end

% loading CTM data
load(sprintf('../matfiles/prepCTM_%d.mat',daydisp(1)));
zh = dailyCTMv;
ch = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];

% getting day of interest
idx = ch(:,3) == datenum(daydisp);

% country outline
cd ../09_mfiles_projections
load('USAcontiguous.mat');
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../04_mfiles_softdata

% plot results
if sum(idx) ~= 0 % if soft data exists
    
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(ch(idx,:),zh(idx),lax,[plotax(:,1) plotax(:,2)]);
    caxis(cax);
    colorbar;

    % setting axis
    xlabel('km');
    ylabel('km');

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for i = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
        cd ../04_mfiles_softdata
        plot(states(:,1),states(:,2),'k-');
    end

    % title
    title(sprintf('CTM PM_{2.5} (\\mug/m^3) on %s',datestr(datenum(daydisp))));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    daydispsave = daydisp(1)*10000 + daydisp(2)*100 + daydisp(3);
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/map_CTM_%d.pdf',daydispsave));

end

end