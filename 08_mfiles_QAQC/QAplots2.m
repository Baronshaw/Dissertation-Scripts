function [] = QAplots2(monthsdate,testrun)
% this function will perform some QA analysis on estimation locations,
% namely a maps of the BME mean with observations values, variance, and
% maps of min and max

if nargin < 1, monthsdate = 13; end
if nargin < 2, testrun = 2; end

if testrun == 1 
    teststr = 'test_'; 
    pathstr = 'matfiles_QAQC';
elseif testrun == 0
    teststr = ''; 
    pathstr = 'matfiles_est';
elseif testrun == 2
    teststr = 'BMEest'; 
    pathstr = 'matfiles_est';
end

% load data
inputdates1 = datevec(datenum(1998,monthsdate,1));
inputdates2 = datevec(datenum(1998,monthsdate+1,1)-1);
startdate = inputdates1(:,1:3);
enddate = inputdates2(:,1:3);
daydisp1 = startdate(1)*10000 + startdate(2)*100 + startdate(3);
daydisp2 = enddate(1)*10000 + enddate(2)*100 + enddate(3);
load(sprintf('../%s/%sWHIMS_%d_%d.mat',pathstr,teststr,daydisp1,daydisp2));

% load id
inputdates1 = datevec(datenum(1998,monthsdate,1));
inputdates2 = datevec(datenum(1998,monthsdate+1,1)-1);
startdate = inputdates1(:,1:3);
enddate = inputdates2(:,1:3);
midpnt = (datenum(enddate(1),enddate(2),enddate(3))-datenum(startdate(1),startdate(2),startdate(3)))/2 + ...
    datenum(startdate(1),startdate(2),startdate(3));
lower = datevec(midpnt-182); lower = lower(1);
upper = datevec(midpnt+182); upper = upper(1);
dummy = repmat(1999:2011,2,1);
dummy = dummy(:);
dummy1 = [ 1999 ; dummy ]; dummy2 = [ dummy ; 2011 ];
uniyrsOff = [ dummy1 dummy2 ];
[ dummy loadnum ] = ismember([lower upper],[uniyrsOff],'rows');
if loadnum < 1, loadnum = 1; end

% getting hard data 
load('../matfiles/distprepsoft_newall.mat');
zh = Obsg{loadnum}(:);
cht = repmat(tMEO{loadnum},length(cMSObs{loadnum}),1); cht = cht(:);
chty = floor(cht./10000); chtm = floor( (cht-(chty*10000))./100 );
chtd = cht - chty*10000 - chtm*100;
cht = datenum(chty,chtm,chtd);
ch = [ repmat(cMSObs{loadnum},length(tMEO{loadnum}),1) cht ];
idxs = isnan(zh);
zh(idxs) = []; ch(idxs,:) = [];

% getting best day
testdays = datenum(startdate(1),startdate(2),12:18);
bestday = cell2mat( arrayfun(@(x) sum(ch(:,3)==x),testdays,'UniformOutput',false) );
ckt = testdays(bestday==max(bestday)); ckt = ckt(1); % if repeates

load('../09_mfiles_projections/USAcontiguous.mat');
cd ../09_mfiles_projections
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../08_mfiles_QAQC

%%% figure 1, plotting BME mean
lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
zk(idxoor) = 0;
zk(idxoor) = zkoor;
zdisplay = zk + mItest;

if testrun == 2
    [xg yg Zg] = plotField(ck,zdisplay,lax,[plotax(:,1) plotax(:,2)]);
elseif testrun == 1 | testrun == 0
    figure; hold on; box on;
    idxm = ck(:,3) == ckt;
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    cax = [2 30];
    colorplot([ck(idxm,1) ck(idxm,2)],zdisplay(idxm),'jet',Property,{'s',6,[0 0 0]
    plotFieldmaskonly(lax,[plotax(:,1) plota
end
caxis([2 30]);
colorbar;
axis(lax);

% setting axis
xlabel('km');
ylabel('km');

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for i = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
    cd ../08_mfiles_QAQC
    plot(states(:,1),states(:,2),'k-');
end

% title
title(sprintf('PM_{2.5} (\\mug/m^3) on %s',datestr(ckt)));

% overlaying the hard data
idxh = ch(:,3) == ckt;
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'o',7,[0 0 0]};
cax = [2 30];
colorplot([ch(idxh,1) ch(idxh,2)],zh(idxh),'jet',Property,Value,cax);

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
temp = datevec(ckt);
daydispsave = temp(1).*10000 + temp(2).*100 + temp(3);
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_QAQC/%smap_%d_mean.pdf',teststr,daydispsave));
print(gcf,'-painters','-dpng','-r600',sprintf('../plots_QAQC/%smap_%d_mean.png',teststr,daydispsave));

%%% figure 2, plotting variance
lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
vk(idxoor) = vkoor;

if testrun == 2
    [xg yg Zg] = plotField(ck,vk,lax,[plotax(:,1) plotax(:,2)]);
elseif testrun == 1 | testrun == 0
    figure; hold on; box on;
    idxm = ck(:,3) == ckt;
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    cax = [2 30];
    colorplot([ck(idxm,1) ck(idxm,2)],vk(idxm),'jet',Property,{'s',6,[0 0 0]},cax);
    plotFieldmaskonly(lax,[plotax(:,1) plotax(:,2)]);
end
caxis([0 50]);
colorbar;
axis(lax);

% setting axis
xlabel('km');
ylabel('km');

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for i = 1:length(X)
    cd ../GeodeticTransformations
    states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
    cd ../08_mfiles_QAQC
    plot(states(:,1),states(:,2),'k-');
end

% title
title(sprintf('variance PM_{2.5} (\\mug/m^3) on %s',datestr(ckt)));

% load id
inputdates1 = datevec(datenum(1998,monthsdate,1));
inputdates2 = datevec(datenum(1998,monthsdate+1,1)-1);
startdate = inputdates1(:,1:3);
enddate = inputdates2(:,1:3);
midpnt = (datenum(enddate(1),enddate(2),enddate(3))-datenum(startdate(1),startdate(2),startdate(3)))/2 + ...
    datenum(startdate(1),startdate(2),startdate(3));
lower = datevec(midpnt-182); lower = lower(1);
upper = datevec(midpnt+182); upper = upper(1);
dummy = repmat(1999:2011,2,1);
dummy = dummy(:);
dummy1 = [ 1999 ; dummy ]; dummy2 = [ dummy ; 2011 ];
uniyrsOff = [ dummy1 dummy2 ];
[ dummy loadnum ] = ismember([lower upper],[uniyrsOff],'rows');
if loadnum < 1, loadnum = 1; end

% overlaying the monitoring locations
load('../mfiles/distprepsoft_newall.mat');
if testrun == 2   
    plot(cMSObs{loadnum}(:,1),cMSObs{loadnum}(:,2),'kx','LineWidth',1.75);
elseif testrun == 0 | testrun == 1
    plot(ch(idxh,1),ch(idxh,2),'kx','LineWidth',1.75);
end

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
temp = datevec(ckt);
daydispsave = temp(1).*10000 + temp(2).*100 + temp(3);
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_QAQC/%smap_%d_var.pdf',teststr,daydispsave));
print(gcf,'-painters','-dpng','-r600',sprintf('./plotss_QAQC/%smap_%d_var.png',teststr,daydispsave));

end