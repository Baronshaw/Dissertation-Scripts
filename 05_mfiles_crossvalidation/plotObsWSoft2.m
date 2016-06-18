function [] = plotObsWSoft2(best)
% this function will create time series that plots showing closest hard 
% data, closests soft data, CTM, lambda1 with lambda2 error bars, mean 
% trend with variance, observed data, soft estimate, hard estimate from LOO 
% estimation. This will not create S-curves for the corresponding soft 
% data. After plots are made, I will zoom in to the best/worst days to 
% better see what going on

if nargin < 1, best = 1; end % getting the best performances

if best == 1, beststr = 'best'; else beststr = 'worst'; end

forceddist = 0; 
soft_years = 2001; 
numplots = 10;

% gathering all the data
load(sprintf('../matfiles/Xvalforcediso_LOOCV__nosoft_long_gauss_foriso%dkm.mat',floor(forceddist./1000)));
zkallh = zk_madd; zhallh = zh_Xval; ckallh = ck; vkallh = vk;  

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckallh = cell(length(cMS),1);
meantrendh = cell(length(cMS),1);
for i = 1:length(ckallh)
    disp(i);
    idx2 = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckallh{i,1} = ch(idx2,:);
    meantrendh{i,1} = mI(idx2,:);
end
ckallh = cell2mat(ckallh);

% load variance information
load('../matfiles/covmod_constant_joint exponential exponential_joint.mat');
meantrendErr = f.Cr1;

load(sprintf('../matfiles/Xvalforcediso_LOOCV__soft_long_gauss_foriso%dkm_timezone_n6T90.mat',floor(forceddist./1000)));
zkalls = zk_madd; zhalls = zh_Xval; ckalls = ck; vkalls = vk; 

% load soft mean trend
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 300000 100 50],soft_years));
mIsoft = mI;
pIsoft = pI;
pdsoft = pd;
zdsoft = zd;

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckalls = cell(length(cMS),1);
for i = 1:length(ckalls)
    disp(i);
    idx2 = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckalls{i,1} = ch(idx2 & temp2001==2001,:);
end
ckallcell = ckalls;
ckalls = cell2mat(ckalls);

zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
zhallh = cell2mat(zhallh); zhalls = cell2mat(zhalls);

[lia lib] = ismember(ckallh,ckalls,'rows');
ckallh = ckallh(lia,:);
zkallh = zkallh(lia); zhallh = zhallh(lia);

idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx); zhallh = zhallh(idx); 
zkalls = zkalls(idx); zhalls = zhalls(idx); 
ckallh = ckallh(idx,:); ckalls = ckalls(idx,:);

% load CTM values, dailyCTMv and distCTMv
load(sprintf('../matfiles/prepCTM_%d.mat',soft_years));
distCTMv = [distCTMv datenum(yrmodaCTMv)];

% looking at performance by station/year
temp = datevec(ckalls(:,3));
yrs = temp(:,1);
diffR2 = NaN*ones(length(cMS),1);
corrH = NaN*ones(length(cMS),1);
corrS = NaN*ones(length(cMS),1);
for i = 1:length(cMS)
    
    idx = cMS(i,1) == ckalls(:,1) & cMS(i,2) == ckalls(:,2) & yrs == soft_years;
    if sum(idx) > 0
        corrH(i) = (corr(zhallh(idx),zkallh(idx))).^2;
        corrS(i) = (corr(zhalls(idx),zkalls(idx))).^2;
        diffR2(i) = corrS(i) - corrH(i);
    end
    
end

% plotting time series
if best == 1
    [aorig borig] = sort(diffR2,'descend');
    a = aorig(~isnan(aorig)); b = borig(~isnan(aorig)); 
else
    [a b] = sort(diffR2);
end
for i = 1:numplots
    
    figure; hold on;
    
    idx = cMS(b(i),1) == ckalls(:,1) & cMS(b(i),2) == ckalls(:,2) & yrs == soft_years;
    plot(ckalls(idx,3),zhallh(idx),'ks'); % observed data
    plot(ckalls(idx,3),zkallh(idx),'bo-'); % hard estimate
    plot(ckalls(idx,3),zkalls(idx),'ro-'); % soft estimate
    
    % gather all the lambda1's
    tempval = cell2mat(alllambda1{b(i)});
    templ2 = cell2mat(alllambda2{b(i)});
    tempcs = cell2mat(allcs{b(i)});
    s_yrs = datevec(tempcs(:,3)); s_yrs = s_yrs(:,1);
    tempcs = tempcs(s_yrs==2001,:);
    tempval = tempval(s_yrs==2001);
    templ2 = templ2(s_yrs==2001);
    [lia lib] = ismember(round(tempcs),round(pIsoft),'rows');
    lib(lib==0)=[];
    lam1 = tempval + mIsoft(lib);
    lam2 = templ2;
    
    % gather all the CTM values
    [lia lib] = ismember(round(tempcs),round(distCTMv),'rows');
    lib(lib==0) = [];
    ctm = dailyCTMv(lib);
    
    % gather all the hard data
    temphard = cell2mat(allhard{b(i)});
    tempch = cell2mat(allch{b(i)});
    h_yrs = datevec(tempch(:,3)); h_yrs = h_yrs(:,1);
    tempch = tempch(h_yrs==2001,:);
    temphard = temphard(h_yrs==2001);
    [lia lib] = ismember(round(tempch),round(pI),'rows');
    lib(lib==0)=[];
    harddata = temphard + mI(lib);
    
    % add the hard data
    plot(tempch(:,3),harddata,'bx');
    
    % add the mean trend with variance
    [lia lib] = ismember(round(ckalls(idx,:)),round(pI),'rows');
    lib(lib==0)=[];
    errorbar(ckalls(idx,3),mI(lib),sqrt(meantrendErr)*ones(length(lib),1),'gs'); 
    
    % add the lam1's with lam2's
    errorbar(tempcs(:,3),lam1,sqrt(lam2),'rx'); 
    
    % add the ctm values
    plot(tempcs(:,3),ctm,'rs');

    % find the worst performances for monitoring station
    ckallsub = ckalls(idx,:);
    stationErr = abs(zhallh(idx)-zkallh(idx)) - abs(zhalls(idx)-zkalls(idx));
    if best == 1
        [a2orig b2orig] = sort(stationErr,'descend');
        a2 = a2orig(~isnan(a2orig)); b2 = b2orig(~isnan(b2orig));
    else
        [a2 b2] = sort(stationErr);
    end
    
    % the ones I'll be viewing
    idxB = [];
    for j = 1:10
        idxB = [ idxB ; find( tempcs(:,3) == ckallsub(b2(j),3) ) ];
    end
    errorbar(tempcs(idxB,3),lam1(idxB),sqrt(lam2(idxB)),'cx');
    
    legend('obs','est h','est s','hard','mean','lam1','ctm');
    title(sprintf('obs, est, and soft with diffR2=%f-%f=%f',corrS(b(i)),corrH(b(i)),diffR2(b(i))));
    set(gca,'XTickLabel',datestr(datevec(get(gca,'XTick'))));
    xlabel('time');
    ylabel('PM2.5 conc');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/fulldiag_%s_station_%d_n6.pdf',beststr,i));
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/fulldiag_%s_station_%d_n6.png',beststr,i));
    
    % for each of the 10 best/worst days, zoom into specific day
    for j = 1:numplots                
       
        dayWHI = datevec(ckallsub(b2(j),3)); dayWHI = dayWHI(:,1:3);       
        daysWHIdisp = dayWHI(1)*10^4 + dayWHI(2)*10^2 + dayWHI(3);
        
        idxcell = ckallcell{b(i)}(:,3) == ckallsub(b2(j),3);
        idxday = cMS(b(i),1) == ckalls(:,1) & cMS(b(i),2) == ckalls(:,2) & ckallsub(b2(j),3) == ckalls(:,3);

        figure; hold on;
        plot(ckalls(idxday,3),zhallh(idxday),'ks','MarkerSize',10); % observed data
        plot(ckalls(idxday,3),zkallh(idxday),'bo','MarkerSize',10); % hard estimate
        plot(ckalls(idxday,3),zkalls(idxday),'ro','MarkerSize',10); % soft estimate
        
        chsub = cell2mat(allch{b(i)}(idxcell));
        zhsub = cell2mat(allhard{b(i)}(idxcell));
        lenhard = length(zhsub);
        [lia lib] = ismember(chsub,pI,'rows');
        plot(repmat(ckallsub(b2(j),3),lenhard,1),zhsub+mI(lib),'bx','MarkerSize',6); % hard
        
        errorbar(ckalls(idxday,3),mI(idxday),sqrt(meantrendErr),'gx','MarkerSize',6); % mean

        cssub = cell2mat(allcs{b(i)}(idxcell));
        zssub = cell2mat(alllambda1{b(i)}(idxcell));
        vssub = cell2mat(alllambda2{b(i)}(idxcell));
        lensoft = length(zssub);
        [lia lib] = ismember(cssub,round(pIsoft),'rows');
        errorbar(repmat(ckallsub(b2(j),3),lensoft,1),zssub+mIsoft(lib),sqrt(vssub),'rx','MarkerSize',6);% lam1       
        
        [lia lib] = ismember(round(cssub),round(distCTMv),'rows');
        plot(repmat(ckallsub(b2(j),3),lensoft,1),dailyCTMv(lib),'rs','MarkerSize',6); % ctm
        
        legend('obs','est h','est s','hard','mean','lam1','ctm');
        title(sprintf('obs, est, and soft with diffR2=%f on %d',diffR2(b(i)),daysWHIdisp));
        lenx = get(gca,'XLim'); len = lenx(2)-lenx(1);
        set(gca,'XLim',[lenx(1)-len lenx(2)+len]);
        set(gca,'XTickLabel','');
        xlabel('day of interest');
        ylabel('PM2.5 conc');
        
        % save zoomed in day
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/fulldiag_%s_station_%d_%d_n6.pdf', ...
            beststr,i,daysWHIdisp));
        print(gcf,'-painters','-dpng','-r600',sprintf('../plots/fulldiag_%s_station_%d_%d_n6.png', ...
            beststr,i,daysWHIdisp));
       
    end
    
    close all;
       
end

end