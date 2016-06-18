function [] = BMEtimeseries(nummon,years,soft,constant,gauss)
% this function will create and display time series for randomly selected
% monitoring stations

if nargin < 1, nummon = 20; end % number of monitoring stations
if nargin < 2, years = 2001; end % years for time series
if nargin < 2, soft = 1; end % soft data or not
if nargin < 3, constant = 0; end % constant offset or not
if nargin < 4, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% gathering all the data
load(sprintf('../matfiles/Xvalforcediso_LOOCV__soft_long_gauss_foriso%dkm.mat',0));
zkalls = zk_madd; zhalls = zh_Xval; ckalls = ck; vkalls = vk; 

load(sprintf('../matfiles/Xvalforcediso_LOOCV__nosoft_long_gauss_foriso%dkm.mat',0));
zkallh = zk_madd; zhallh = zh_Xval; ckallh = ck; vkallh = vk;  

zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
zhallh = cell2mat(zhallh); zhalls = cell2mat(zhalls);
idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx); zhallh = zhallh(idx); 
zkalls = zkalls(idx); zhalls = zhalls(idx); 

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckalls = cell(length(cMS),1);
for i = 1:length(ckalls)
    idx2 = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckalls{i,1} = ch(idx2,:);
end
ckalls = cell2mat(ckalls);
ckalls = ckalls(idx,:);

temp = datevec(ckalls(:,3));
yrs = temp(:,1);

% calculate r2 by years and location
R2soft = NaN*ones(length(cMS),1);
R2hard = NaN*ones(length(cMS),1);
for i = 1:length(R2soft)
    idx = yrs == years & cMS(i,1) == ckalls(:,1) & cMS(i,2) == ckalls(:,2);
    if sum(idx) > 0
        R2soft(i) = (corr(zkalls(idx),zhalls(idx))).^2;
        R2hard(i) = (corr(zkallh(idx),zhallh(idx))).^2;
    end
end

% plot time series for a given year by how good the R2 is 
idx = isnan(R2soft);
cMS = cMS(~idx,:);
R2diff = R2soft(~idx) - R2hard(~idx);
[dummy sortidx] = sort(R2diff,'descend');
tempyrs = datevec(ch(:,3)); tempyrs = tempyrs(:,1);
for i = 1:nummon
    figure; hold on;
    idx = tempyrs == years & cMS(sortidx(i),1) == ch(:,1) & cMS(sortidx(i),2) == ch(:,2);
    plot(ch(idx,3),zh(idx),'ks'); % observed values
    idx = yrs == years & cMS(sortidx(i),1) == ckalls(:,1) & cMS(sortidx(i),2) == ckalls(:,2);
    plot(ckalls(idx,3),zkalls(idx),'b--',ckalls(idx,3),zkallh(idx),'r--');
    
    % xlabel
    temp = num2cell( datestr(datevec(get(gca,'XTick'))), 2 );
    set(gca,'XTickLabel',temp);
    % ylabel
    ylabel('PM_{2.5} concentration');
    % title
    title(sprintf('PM_{2.5} concentration for a monitoring station (%d km, %d km)',...
        floor(cMS(sortidx(i),1)./1000),floor(cMS(sortidx(i),2)./1000)));
    legend('Observed','Soft','Hard');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/BMETS_mon%dyr%d_.pdf',i,years));
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/BMETS_mon%dyr%d_.png',i,years));
   
end

end