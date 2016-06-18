function [] = displayBMEmaps(daydisp,soft,constant,gauss)
% this will take the results of BMEmaps.m and display the results in
% the form of maps

if nargin < 1, daydisp = [2001 1 15]; end
if nargin < 2, soft = 1; end
if nargin < 3, constant = 0; end
if nargin < 4, gauss = 1; end

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

if soft == 0,
    plotstr = '(observed only)';
else
    plotstr = '(observed and modeled)';
end

% load BME
cd ../BMELIB2.0b
startup();
cd ../07_mfiles_map

% loading data
load(sprintf('../matfiles/BMEmaps_%d_%d_%d%s%s%s.mat',daydisp(1),daydisp(2),daydisp(3), ...
        softstr,constr,gaussstr));
    
daydisp = unique(datevec(ck(:,3)),'rows'); daydisp = daydisp(:,1:3);

load('../09_mfiles_projections/USAcontiguous.mat');
cd ../09_mfiles_projections
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../07_mfiles_map

%%% plotting mean
lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
[xg yg Zg] = plotField(ck,zk_madd,lax,[plotax(:,1) plotax(:,2)]);
caxis([2 30]);
colorbar;
axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

% setting axis
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
xlabel('km');
ylabel('km');

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for i = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
    cd ../07_mfiles_map
    plot(states(:,1),states(:,2),'k-');
end

% title
title(sprintf('PM_{2.5} (\\mug/m^3) on %s %s',datestr(datenum(daydisp)),plotstr));

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
daydispsave = daydisp(1)*10000 + daydisp(2)*100 + daydisp(3);
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/BME_maps_%d%s%s%s.pdf',daydispsave,softstr,constr,gaussstr));

%%% plotting the hard data

% overlay the hard data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;
idx = ch(:,3)==datenum(daydisp);
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'o',6,[0 0 0]};
%cax = [min(Zg(:)) max(Zg(:))];
cax = [2 30];
colorplot([ch(idx,1) ch(idx,2)],zh(idx),'jet',Property,Value,cax);
caxis(cax);
colorbar;
axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

% setting axis
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
xlabel('km');
ylabel('km');

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for i = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
    cd ../07_mfiles_map
    plot(states(:,1),states(:,2),'k-');
end

% title
title(sprintf('PM_{2.5} (\\mug/m^3) on %s %s',datestr(datenum(daydisp)),plotstr));

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
daydispsave = daydisp(1)*10000 + daydisp(2)*100 + daydisp(3);
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/Obs_maps_%d%s%s%s.pdf',daydispsave,softstr,constr,gaussstr));
print(gcf,'-painters','-dpng','-r600',sprintf('../plots/Obs_maps_%d%s%s%s.png',daydispsave,softstr,constr,gaussstr));

%%% plotting the variance

figure; hold on;

% plot results
lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
cax = [prctile(sqrt(vk),5) prctile(sqrt(vk),95)];
[xg yg Zg] = plotField(ck,sqrt(vk),lax,[plotax(:,1) plotax(:,2)]);
caxis(cax);
colorbar;
axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

% location of monitors
plot(ch(idx,1),ch(idx,2),'kx');

% setting axis
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
xlabel('km');
ylabel('km');

% overlaying the states
load('../09_mfiles_projections/USAstates5.mat');
for i = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
    cd ../07_mfiles_map
    plot(states(:,1),states(:,2),'k-');
end

% title
title(sprintf('Standard Deviation PM_{2.5} on %s %s',datestr(datenum(daydisp)),plotstr));

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
daydispsave = daydisp(1)*10000 + daydisp(2)*100 + daydisp(3);
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/BME_maps_%d%s%s%s_variance.pdf',daydispsave,softstr,constr,gaussstr));
print(gcf,'-painters','-dpng','-r600',sprintf('../plots/BME_maps_%d%s%s%s_variance.png',daydispsave,softstr,constr,gaussstr));


end