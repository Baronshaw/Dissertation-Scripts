function [] = plotMeanTrend(years)
% this function will go through all the calculated mean trends and plot
% time series and maps

if nargin < 1, years = 2001; end

% loading 2001 data
load(sprintf('../matfiles/prepObs_%d.mat',years));
yrs = floor(yrmodaObs./10000);
mos = floor( (yrmodaObs-yrs*10000) ./ 100 );
das = yrmodaObs - yrs*10000 - mos*100;
pd = [coordObs datenum(yrs,mos,das)];
zd = Obs;

% find the maximum of closest neighbors
[Zd,MSd,MEd,nanratio]=valstv2stg(pd,zd); MSd = [MSd(:,2) MSd(:,1)];
DMS = sqrt(bsxfun(@plus,dot(MSd,MSd,2),dot(MSd,MSd,2)')-2*(MSd*MSd'));
DME = abs(bsxfun(@minus,MEd,MEd'))';
temp = sort(DMS,2);
maxMinNei = max(temp(:,2));

% all the smoothing parameters to try
ars = [ 20 40 20000 40000 50000:50000:2500000 ];
ats = [ 4:4:20 25:25:300 ];

spatRad = max([2*ars' maxMinNei*ones(length(ars),1)],[],2);
atconstant = 30;
arconstant = 300000;
smoothingParamsSpace = [spatRad ars' 2*atconstant*ones(length(ars),1) atconstant*ones(length(ars),1)];
smoothingParamsTime = [spatRad(1)*ones(length(ats),1) arconstant*ones(length(ats),1) 2*ats' ats'];

% plotting maps of the data for that day

% load data
load(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,smoothingParamsSpace(1,:)));

% getting unique days
uniday = unique(pI(:,3));

% country outline
cd ../09_mfiles_projections
load('USAcontiguous.mat');
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../02_mfiles_offset
lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];

for j = 1:size(uniday,1)

    % getting correct day 
    figure; hold on;
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',6,[0 0 0]};
    cax = [0 20];
    idx = pd(:,3) == uniday(j) & pd(:,1) > lax(1) & pd(:,1) < lax(2) & ...
        pd(:,2) > lax(3) & pd(:,2) < lax(4);
    colorplot([pd(idx,1) pd(idx,2)],zd(idx),jet,Property,Value,cax); 
    colorbar;

    % setting axis
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    xlabel('km');
    ylabel('km');

    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for k = 1:length(X)
        cd ../09_mfiles_projections
        coordObs = ell2lambertcc([X{k},Y{k}],'whiproj2001');
        cd ../02_mfiles_offset
        plot(coordObs(:,1),coordObs(:,2),'k-');
    end

    daystr = datestr(datevec(uniday(j)),'mm/dd/yy');
    title(sprintf('PM_{2.5} (\\mug/m^3) data on %s',daystr));

    % save figures
    daystr = datestr(datevec(uniday(j)),'yyyymmdd');
    h = gcf;
    print(h,'-painters','-dpdf','-r600',sprintf('../plots_meanTrend/meanTrend_map_%s_data.pdf',daystr));

end

% first look at maps
for i = 1:size(smoothingParamsSpace,1)  

    % load data
    load(sprintf('matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,smoothingParamsSpace(i,:)));

    % getting unique days
    uniday = unique(pI(:,3));
    
    for j = 1:size(uniday,1)

        % getting correct day 
        figure; hold on;
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',6,[0 0 0]};
        cax = [0 20];
        idx = pI(:,3) == uniday(j) & pI(:,1) > lax(1) & pI(:,1) < lax(2) & ...
        pI(:,2) > lax(3) & pI(:,2) < lax(4);
        colorplot([pI(idx,1) pI(idx,2)],mI(idx),jet,Property,Value,cax); 
        colorbar;

        % setting axis
        set(gca,'XTickLabel',get(gca,'XTick')/1000);
        set(gca,'YTickLabel',get(gca,'YTick')/1000);
        xlabel('km');
        ylabel('km');

        % overlaying the states
        load('../09_mfiles_projections/USAstates5.mat');
        for k = 1:length(X)
            cd ../09_mfiles_projections
            coordObs = ell2lambertcc([X{k},Y{k}],'whiproj2001');
            cd ../02_mfiles_offset
            plot(coordObs(:,1),coordObs(:,2),'k-');
        end

        daystr = datestr(datevec(uniday(j)),'mm/dd/yy');
        title(sprintf('PM_{2.5} (\\mug/m^3) on %s, R=%0.2f (km), T=%0d (days)',...
            daystr,smoothingParamsSpace(i,2)./1000,smoothingParamsSpace(i,4)));

        % save figures
        daystr = datestr(datevec(uniday(j)),'yyyymmdd');
        h = gcf;
        print(h,'-painters','-dpdf','-r600',sprintf('../plots_meanTrend/meanTrend_map_%s_%d_%d_%d_%d.pdf', ...
            daystr,smoothingParamsSpace(i,:)));
    
    end
    close all
    
end

% second look at time series
for i = 1:size(smoothingParamsTime,1)
    
    % loading file
    load(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,smoothingParamsTime(i,:)));
    
    % getting unique locations
    unilocs = unique(pI(:,1:2),'rows');
    
    for j = 1:size(unilocs,1)
        
       figure; hold on;
        
        % the data
        idx = pd(:,1) == unilocs(j,1) & pd(:,2) == unilocs(j,2);
        plot(pd(idx,3),zd(idx),'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
            'MarkerSize',4);
        
        % plot temporal location
        idx = pI(:,1) == unilocs(j,1) & pI(:,2) == unilocs(j,2);
        temp1 = pI(idx,3);
        temp2 = mI(idx);
        [sorted sortidx] = sort(temp1);
        plot(temp1(sortidx),temp2(sortidx),'b-','LineWidth',2,'MarkerFaceColor','b');
        
        % legend    
        legendCell = cellstr(num2str([smoothingParamsTime(i,2)'./1000 smoothingParamsTime(i,4)'], 'R=%d, T=%d'));
        legendCell = ['data'; legendCell];
        legend(legendCell,'Location','Best');

        % changing axis
        xvals1 = get(gca,'XTick');
        xvals2 = datestr(datevec(xvals1),'mm/dd/yy');
        set(gca,'XTickLabel',xvals2)
        xlabel('dates');
        ylabel(sprintf('PM2.5 \\mug/m^{3}'));
        
        title(sprintf('Time Series, R=%d (km), T=%0d (days)\nlocation (%0.2f km,%0.2f km)',...
            smoothingParamsTime(i,2)/1000,smoothingParamsTime(i,4), ...
            unilocs(j,1)./1000,unilocs(j,2)));   
    
        % save figure
        h = gcf;
        print(h,'-painters','-dpdf','-r600',sprintf('../plots_meanTrend/meanTrend_timeseries_yr%d_loc%0.2f_%0.2f_%d_%d_%d_%d.pdf',...
            years,unilocs(j,1),unilocs(j,2),smoothingParamsTime(i,:)));
        
    end
    close all
      
end

% % final smoothing parameters
% ar2 = [20000 50000 300000 1000000];
% at2 = [10 20 50 200];

end