function [] = getplotMeanTrendfinal(years)
% this function will calculate the mean trend and plot maps and time series 
% of the final mean trend parameters. Note: the time series plots are 
% updated to overlay all the temporal smoothing parameters

if nargin  < 1, years = 2001; end

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

% final smoothing parameters
ars = [20000 50000 300000 1000000];
ats = [10 20 50 200];

% for space, I'll look at one day each month. First see which day in each
% month measured the most
temp = histc(pd(:,3),[datenum(years,1,1):datenum(years,12,31)]');
daysinyear = datevec([datenum(years,1,1):datenum(years,12,31)]');
for i = 1:12
    idx = daysinyear(:,2) == i;
    dummy = find(temp(idx)==max(temp(idx)));
    daysinmonth = daysinyear(idx,:);
    maxshow(i,:) = daysinmonth(dummy(1),1:3);
end
[lia lib] = ismember(pd(:,3),datenum(maxshow)); 
pIr = pd(lia,:);
% for time, I'll look at 10 locations with the most amount of data
uniloc = unique(pd(:,1:2),'rows');
moncnt = NaN*ones(size(uniloc,1),1);
for i = 1:size(uniloc,1)
    idx = uniloc(i,1) == pd(:,1) & uniloc(i,2) == pd(:,2);
    moncnt(i) = sum(idx);
end
[sortm sortidx] = sort(moncnt,'descend');
[lia lib] = ismember(pd(:,1:2),uniloc(sortidx(1:10),1:2),'rows'); 
pIt = pd(lia,:);
pI = [pIr;pIt];

% first calculating the final smoothing parameters
for i = 1:length(ars)
    cd ../10_mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv(pd,zd,[maxMinNei ars(i) 2*ats(i) ats(i)],pI);
    toc
    cd ..    
    % saving results
    save(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,[maxMinNei ars(i) 2*ats(i) ats(i)]), ...
        'pd','zd','pI','mI');
end

% plotting maps of the data for that day

% load data
load(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,[maxMinNei ars(1) 2*ats(1) ats(1)]));

% getting unique days
uniday = unique(pIr(:,3));

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
for i = 1:length(ars)  

    % load data
    load(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',years,[maxMinNei ars(i) 2*ats(i) ats(i)]));

    % getting unique days
    uniday = unique(pIr(:,3));
    
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
            daystr,ars(i)./1000,ats(i)));

        % save figures
        daystr = datestr(datevec(uniday(j)),'yyyymmdd');
        h = gcf;
        print(h,'-painters','-dpdf','-r600',sprintf('../plots_meanTrend/meanTrend_map_%s_%d_%d_%d_%d.pdf', ...
            daystr,[maxMinNei ars(i) 2*ats(i) ats(i)]));
    
    end
    close all
    
end

% finally, visualizing in time
styleOrder = {'-o' ; '-.' ; '--' ; ':'};
plotcolor = cool(length(ats));
% getting unique locations
unilocs = unique(pIt(:,1:2),'rows');

for i = 1:size(unilocs,1)
    
    figure; hold on; 

    % the data
    idx = pd(:,1) == unilocs(i,1) & pd(:,2) == unilocs(i,2);
    plot(pd(idx,3),zd(idx),'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
        'MarkerSize',4);
        
    for j = 1:length(ats) % going through all T's        
        % load data
        load(sprintf('../matfiles/meanTrend_test_%d_%d_%d_%d_%d.mat',2001,[maxMinNei ars(j) 2*ats(j) ats(j)]));

        % plot temporal location
        idx = pI(:,1) == unilocs(j,1) & pI(:,2) == unilocs(j,2);
        temp1 = pI(idx,3);
        temp2 = mI(idx);
        [sorted sortidx] = sort(temp1);
        plot(temp1(sortidx),temp2(sortidx),styleOrder{j},'Color',plotcolor(j,:),'LineWidth',2,'MarkerFaceColor',plotcolor(j,:));
          
    end
    
    % legend    
    legendCell = cellstr(num2str([ars'./1000 ats'], 'R=%d, T=%d'));
    legendCell = ['data'; legendCell];
    legend(legendCell,'Location','Best');

    % changing axis
    xvals1 = get(gca,'XTick');
    xvals2 = datestr(datevec(xvals1),'mm/dd/yy');
    set(gca,'XTickLabel',xvals2)
    xlabel('dates');
    ylabel(sprintf('PM2.5 \\mug/m^{3}'));

    title(sprintf('Time Series, Changing R (km), Changing T (days)\nlocation (%0.2f km,%0.2f km)',...
        unilocs(j,1)./1000,unilocs(j,2)));   

    % save figure
    h = gcf;
    print(h,'-painters','-dpdf','-r600',sprintf('../plots_meanTrend/meanTrend_timeseries_yr%d_loc%0.2f_%0.2f_%d_%d_%d_%d.pdf',...
        years,unilocs(j,1),unilocs(j,2),[maxMinNei ars(j) 2*ats(j) ats(j)]));

end
close all

end