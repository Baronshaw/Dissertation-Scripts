function [] = plotCTMavg(daydisp,avgval,cax)
% this function will average CTM values over a month or year

if nargin < 1, daydisp = [2001 7 1]; end
if nargin < 2, avgval = 'mon'; end
if nargin < 3, cax = [2 30]; end

% loading CTM data
load(sprintf('../matfiles/prepCTM_%d.mat',daydisp(1)));
zh = dailyCTMv;
ch = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];

% getting days of interest
if strcmp(avgval,'mon') == 1
    idx = ch(:,3) >= datenum([daydisp(1:2) 1]) & ch(:,3) <= datenum([daydisp(1) daydisp(2)+1 1]);
elseif strcmp(avgval,'yr') == 1
    idx = ch(:,3) >= datenum([daydisp(1) 1 1]) & ch(:,3) <= datenum([daydisp(1) 12 31]);
end

% avgeraging days of interest
if sum(idx) ~= 0
    
    chsub = ch(idx,:);
    zhsub = zh(idx);
    unilocs = unique(chsub(:,1:2),'rows');
    zhavg = NaN*ones(size(unilocs,1),1);
    for i = 1:size(unilocs,1)
        if mod(i,100) == 0, disp(i); end
        idx = unilocs(i,1) == chsub(:,1) & unilocs(i,2) == chsub(:,2);
        zhavg(i) = nanmean(zhsub(idx));
    end

    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../04_mfiles_softdata

    % plot results
    lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
    [xg yg Zg] = plotField(unilocs,zhavg,lax,[plotax(:,1) plotax(:,2)]);
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

    if strcmp(avgval,'mon')
        daydispsave = daydisp(1)*100 + daydisp(2);
        str1 = 'XX';
        str2 = 'mon';
    elseif strcmp(avgval,'yr')
        daydispsave = daydisp(1);
        str1 = 'XXXX';
        str2 = 'yr';
    end

    % title
    title(sprintf('CTM PM_{2.5} (\\mug/m^3) on %d%s %s',daydispsave,str1,str2));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/map_CTM_%d%s_%s_timezone.pdf',daydispsave,str1,str2));
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/map_CTM_%d%s_%s_timezone.png',daydispsave,str1,str2));

end

end