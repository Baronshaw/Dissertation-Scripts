function [] = evalMistakes()
% this function will load results and compare performance for krig and
% krigcorrected and by nsmax

soft = 1; 
constant = 0; 
gauss = 1;
forceddist = 0; 
eststr = '';

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% load hard data
load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm.mat', ...
        '_nosoft',constr,gaussstr,floor(forceddist./1000)));
zkall = zk;
zhall = zh_Xval;
vkall = vk;     
zkallh = cell2mat(zkall);
zhallh = cell2mat(zhall); 
vkallh = cell2mat(vkall); 

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckallC = cell(size(cMS,1),1);
ckallH = cell(size(cMS,1),1);
zhallC = cell(size(cMS,1),1);
zkHC = cell(size(cMS,1),1);
vkHC = cell(size(cMS,1),1);
for i = 1:500 % size(cMS,1)
    disp(i);
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckallC{i,1} = ch(idx & temp2001==2001,:);
    ckallH{i,1} = ch(idx,:);
    zhallC{i,1} = zh(idx & temp2001==2001) - mI(idx & temp2001==2001);
end

ckallh = cell2mat(ckallH);
    
for i = 1:500 % size(cMS,1)
    disp(i);
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    [lia lib] = ismember(ch(idx & temp2001==2001,:),ckallh,'rows');
    zkHC{i,1} = zkallh(lib);
    vkHC{i,1} = vkallh(lib);
end

for i = 3:3 % loop for each nsmax
    
    % krigingME2
    load(sprintf('../matfiles/testingkrigME%d_nsmax%d_Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s.mat', ...
            1,i,softstr,constr,gaussstr,floor(forceddist./1000),eststr))
    zkME = zk; vkME = vk; cellLamWME = cellLamW; lam2ME = alllambda2; 
        
    % krigingME2_correct
    load(sprintf('../matfiles/testingkrigME%d_nsmax%d_Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s.mat', ...
            0,i,softstr,constr,gaussstr,floor(forceddist./1000),eststr))
    zkMEW = zk; vkMEW = vk; cellLamWMEW = cellLamW; lam2MEW = alllambda2;
    
    % load small lambda2 krig
    load(sprintf('../matfiles/testingkrigME%d_nsmax%d_Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s_lam20.mat', ...
        1,i,softstr,constr,gaussstr,floor(forceddist./1000),eststr))
    zkMEl2 = zk; vkMEl2 = vk; cellLamWMEl2 = cellLamW; 
    
    % load small lambda2 krigcorrected
    load(sprintf('../matfiles/testingkrigME%d_nsmax%d_Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s_lam20.mat', ...
            0,i,softstr,constr,gaussstr,floor(forceddist./1000),eststr))
    zkMEWl2 = zk; vkMEWl2 = vk; cellLamWMEWl2 = cellLamW; 
    
%     % BME
%     load(sprintf('../matfiles/testingkrigME%d_nsmax%d_Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s.mat', ...
%             1,i,softstr,constr,'_nongauss',floor(forceddist./1000),eststr))
    
    idx = cellfun(@isempty,zk);
    zkME = zkME(~idx); vkME = vkME(~idx); cellLamWME = cellLamWME(~idx); lam2ME = lam2ME(~idx);
    zkMEW = zkMEW(~idx); vkMEW = vkMEW(~idx); cellLamWMEW = cellLamWMEW(~idx); lam2MEW = lam2MEW(~idx);
    zkMEl2 = zkMEl2(~idx); vkMEl2 = vkMEl2(~idx); cellLamWMEl2 = cellLamWMEl2(~idx); 
    zkMEWl2 = zkMEWl2(~idx); vkMEWl2 = vkMEWl2(~idx); cellLamWMEWl2 = cellLamWMEWl2(~idx); 
%     moments = moments(~idx); momentsA = cell2mat(moments);
    zhall = zhallC(~idx); ckall = ckallC(~idx);
    zkH = zkHC(~idx); vkH = vkHC(~idx);
    
    r2H(i) = (corr(cell2mat(zkH),cell2mat(zhall),'type','Pearson')).^2;
    r2ME(i) = (corr(cell2mat(zkME),cell2mat(zhall),'type','Pearson')).^2;
    r2MEW(i) = (corr(cell2mat(zkMEW),cell2mat(zhall),'type','Pearson')).^2;
    r2MEl2(i) = (corr(cell2mat(zkMEl2),cell2mat(zhall),'type','Pearson')).^2;
    r2MEWl2(i) = (corr(cell2mat(zkMEWl2),cell2mat(zhall),'type','Pearson')).^2;
%     r2BME(i) = (corr(momentsA(:,1),cell2mat(zhall),'type','Pearson')).^2;
    
    % for every estimation location, look at % kriging weight of soft
    PerME = cell(length(zhall),1);
    PerMEW = cell(length(zhall),1);
    PerMEl2 = cell(length(zhall),1);
    PerMEWl2 = cell(length(zhall),1);
    meanlam2ME = cell(length(zhall),1);
    meanlam2MEW = cell(length(zhall),1);
    for j = 1:length(PerME)
        temp = cellLamWME{j};
        len = length(temp{1});
        PerME{j} = cell2mat( cellfun( @(x) 100*sum(x(8:len))/sum(x), temp, 'UniformOutput',false) );
        tempW = cellLamWMEW{j};
        lenW = length(tempW{1});
        PerMEW{j} = cell2mat( cellfun( @(x) 100*sum(x(8:lenW))/sum(x), tempW, 'UniformOutput',false) );
        templ2 = cellLamWMEl2{j};
        lenl2 = length(templ2{1});
        PerMEl2{j} = cell2mat( cellfun( @(x) 100*sum(x(8:lenl2))/sum(x), templ2, 'UniformOutput',false) );
        tempWl2 = cellLamWMEWl2{j};
        lenWl2 = length(tempWl2{1});
        PerMEWl2{j} = cell2mat( cellfun( @(x) 100*sum(x(8:lenWl2))/sum(x), tempWl2, 'UniformOutput',false) );
        temp2ME = lam2ME{j};
        meanlam2ME{j} = cell2mat( cellfun( @(x) mean(x), temp2ME, 'UniformOutput',false) );
        temp2MEW = lam2MEW{j};
        meanlam2MEW{j} = cell2mat( cellfun( @(x) mean(x), temp2MEW, 'UniformOutput',false) );
    end
    
    % scaled histogram of kriging weights
    figure; hold on;
    title(sprintf('%% soft weights KriginME, nsmax = %d',i));
    bar(hist(cell2mat(PerME),100) ./ sum(hist(cell2mat(PerME),100)));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/hist_krigME_nsmax%d.png',i));
    figure; hold on;
    title(sprintf('%% soft weights KriginME corrected, nsmax = %d',i));
    bar(hist(cell2mat(PerMEW),100) ./ sum(hist(cell2mat(PerMEW),100)));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/hist_krigMEW_nsmax%d.png',i));
    figure; hold on;
    title(sprintf('%% soft weights KriginME l2, nsmax = %d',i));
    bar(hist(cell2mat(PerMEl2),100) ./ sum(hist(cell2mat(PerMEl2),100)));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/hist_krigME_nsmax%dl2.png',i));
    figure; hold on;
    title(sprintf('%% soft weights KriginME corrected l2, nsmax = %d',i));
    bar(hist(cell2mat(PerMEWl2),100) ./ sum(hist(cell2mat(PerMEWl2),100)));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/hist_krigMEW_nsmax%dl2.png',i));
    
    % compare weights with bias    
    temp1a = abs(cell2mat(zkME) - cell2mat(zhall));
    temp1b = abs(cell2mat(zkMEW) - cell2mat(zhall));
    temp1c = abs(cell2mat(zkMEl2) - cell2mat(zhall));
    temp1d = abs(cell2mat(zkMEWl2) - cell2mat(zhall));
    temp2a = cell2mat(PerME);
    temp2b = cell2mat(PerMEW);
    temp2c = cell2mat(PerMEl2);
    temp2d = cell2mat(PerMEWl2);
    temp3a = abs(cell2mat(zkME) - cell2mat(zhall)) - abs(cell2mat(zkH) - cell2mat(zhall));
    temp3b = abs(cell2mat(zkMEW) - cell2mat(zhall)) - abs(cell2mat(zkH) - cell2mat(zhall));
    temp3c = abs(cell2mat(zkMEl2) - cell2mat(zhall)) - abs(cell2mat(zkH) - cell2mat(zhall));
    temp3d = abs(cell2mat(zkMEWl2) - cell2mat(zhall)) - abs(cell2mat(zkH) - cell2mat(zhall));
    figure; hold on;
    plot(temp2a,temp1a,'b.',temp2b,temp1b,'r.');
    title(sprintf('%% soft weights vs bias, nsmax = %d',i));
    legend('krig','krigcorrected');
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvbias_krigME_nsmax%d.png',i));
    figure; hold on;
    plot(temp2c,temp1c,'b.',temp2d,temp1d,'r.');
    title(sprintf('%% soft weights vs bias l2, nsmax = %d',i));
    legend('krig','krigcorrected');
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvbias_krigME_nsmax%dl2.png',i));
    figure; hold on;
    plot(temp2a,temp3a,'b.',temp2b,temp3b,'r.');
    title(sprintf('%% soft weights vs bias with hard, nsmax = %d',i));
    legend('krig','krigcorrected');
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvbiashard_krigME_nsmax%d.png',i));
    
    % weights and bias by station
    ckallM = cell2mat(ckall);
    temp = cell2mat(meanlam2ME);
    uniID = unique(ckallM(:,1:2),'rows');
    for j = 1:size(uniID,1)
        idx = ckallM(:,1) == uniID(j,1) & ckallM(:,2) == uniID(j,2);
        unitemp2a(j) = mean(temp2a(idx));
        unitemp3a(j) = mean(temp3a(idx));
        unimeanlam2ME(j) = mean(temp(idx));
    end
    figure; hold on; plot(unitemp2a,unitemp3a,'r.');
    figure; hold on; plot(unimeanlam2ME,unitemp3a,'r.');
    
    % colorplot
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for k = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(unitemp2a,5) prctile(unitemp2a,95)];
    colorplot(uniID,unitemp2a',hot(10),Property,Value,cax);       
    caxis(cax);
    colorbar;
    
    % colorplot
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for k = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(unitemp3a,5) prctile(unitemp3a,95)];
    colorplot(uniID,unitemp3a',hot(10),Property,Value,cax);       
    caxis(cax);
    colorbar;
    
    % colorplot
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for k = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(unimeanlam2ME,5) prctile(unimeanlam2ME,95)];
    colorplot(uniID,unimeanlam2ME',hot(10),Property,Value,cax);       
    caxis(cax);
    colorbar;
    
    % % bias by station
    n = 0;
    for j = 0:10:90
        n = n + 1;
        idxMEh = unitemp3a >= prctile(unitemp3a,j) & unitemp3a < prctile(unitemp3a,j+10);
        uniweightbiasMEh(n) = median(unitemp2a(idxMEh));
        uniweightbiasMElam2h(n) = median(unimeanlam2ME(idxMEh));
    end
    figure; hold on; plot(1:10,uniweightbiasMEh,'r.');
    figure; hold on; plot(1:10,uniweightbiasMElam2h,'r.');
    
    % compare lambda2 with bias
    figure; hold on;
    plot(cell2mat(meanlam2ME),temp1a,'b.',cell2mat(meanlam2MEW),temp1b,'r.');
    title(sprintf('lambda2 v bias, nsmax = %d',i));
    legend('krig','krigcorrected');
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/Lam2vbias_krigME_nsmax%d.png',i));
    figure; hold on;
    plot(cell2mat(meanlam2ME),temp3a,'b.',cell2mat(meanlam2MEW),temp3b,'r.');
    title(sprintf('lambda2 vs bias with hard, nsmax = %d',i));
    legend('krig','krigcorrected');
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/Lam2vbiashard_krigME_nsmax%d.png',i));
    
    % which weight contributes the most bias?
    temp2alam2 = cell2mat(meanlam2ME);
    temp2blam2 = cell2mat(meanlam2MEW);
    weightbiasME = NaN*ones(100,1);
    weightbiasMEW = NaN*ones(100,1);
    weightbiasMEl2 = NaN*ones(100,1);
    weightbiasMEWl2 = NaN*ones(100,1);
    weightbiasMEh = NaN*ones(100,1);
    weightbiasMEWh = NaN*ones(100,1);
    weightbiasMEl2h = NaN*ones(100,1);
    weightbiasMEWl2h = NaN*ones(100,1);
    weightbiasMElam2 = NaN*ones(100,1);
    weightbiasMEWlam2 = NaN*ones(100,1);
    weightbiasMElam2h = NaN*ones(100,1);
    weightbiasMEWlam2h = NaN*ones(100,1);
    for j = 0:100-1
        idxME = temp1a >= prctile(temp1a,j) & temp1a < prctile(temp1a,j+1);
        idxMEh = temp3a >= prctile(temp3a,j) & temp3a < prctile(temp3a,j+1);
        weightbiasME(j+1) = median(temp2a(idxME));
        weightbiasMEh(j+1) = median(temp2a(idxMEh));
        idxMEW = temp1b >= prctile(temp1b,j) & temp1b < prctile(temp1b,j+1);
        idxMEWh = temp3b >= prctile(temp3b,j) & temp3b < prctile(temp3b,j+1);
        weightbiasMEW(j+1) = median(temp2b(idxMEW));
        weightbiasMEWh(j+1) = median(temp2b(idxMEWh));
        idxMEl2 = temp1c >= prctile(temp1c,j) & temp1c < prctile(temp1c,j+1);
        idxMEl2h = temp3c >= prctile(temp3c,j) & temp3c < prctile(temp3c,j+1);
        weightbiasMEl2(j+1) = median(temp2c(idxMEl2));
        weightbiasMEl2h(j+1) = median(temp2c(idxMEl2h));
        idxMEWl2 = temp1d >= prctile(temp1d,j) & temp1d < prctile(temp1d,j+1);
        idxMEWl2h = temp3d >= prctile(temp3d,j) & temp3d < prctile(temp3d,j+1);
        weightbiasMEWl2(j+1) = median(temp2d(idxMEWl2));
        weightbiasMEWl2h(j+1) = median(temp2d(idxMEWl2h));
        weightbiasMElam2(j+1) = median(temp2alam2(idxME));
        weightbiasMEWlam2(j+1) = median(temp2blam2(idxMEW));
        weightbiasMElam2h(j+1) = median(temp2alam2(idxMEh));
        weightbiasMEWlam2h(j+1) = median(temp2blam2(idxMEWh));
    end
    figure; hold on;
    plot(0:99,weightbiasME,'b.',0:99,weightbiasMEW,'r.');
    legend('krig','krigcorrected');
    xlabel('percentile of bias');
    ylabel('median weights');
    title(sprintf('kriging weights by bias percentile, nsmax = %d',i));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvpercent_krigME_nsmax%d.png',i));
    figure; hold on;
    plot(0:99,weightbiasMEh,'b.',0:99,weightbiasMEWh,'r.');
    legend('krig','krigcorrected');
    xlabel('percentile of bias');
    ylabel('median weights');
    title(sprintf('kriging weights by hard bias percentile, nsmax = %d',i));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvpercenthard_krigME_nsmax%d.png',i));
    figure; hold on;
    plot(0:99,weightbiasMEl2,'b.',0:99,weightbiasMEWl2,'r.');
    legend('krig','krigcorrected');
    xlabel('percentile of bias');
    ylabel('median weights');
    title(sprintf('kriging weights by bias percentile l2, nsmax = %d',i));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvpercent_krigME_nsmax%dl2.png',i));
    figure; hold on;
    plot(0:99,weightbiasMEl2h,'b.',0:99,weightbiasMEWl2h,'r.');
    legend('krig','krigcorrected');
    xlabel('percentile of bias');
    ylabel('median weights');
    title(sprintf('kriging weights by hard bias percentile l2, nsmax = %d',i));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/weightvpercenthard_krigME_nsmax%dl2.png',i));
    
    % which median lambda2 contributes most to bias?
    figure; hold on;
    plot(0:99,weightbiasMElam2,'b.',0:99,weightbiasMEWlam2,'r.');
    legend('krig','krigcorrected');
    xlabel('percentile of bias');
    ylabel('mean lambda2');
    title(sprintf('ambda2 by bias, nsmax = %d',i));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/lam2vpercent_krigME_nsmax%d.png',i));
    figure; hold on;
    plot(0:99,weightbiasMElam2h,'b.',0:99,weightbiasMEWlam2h,'r.');
    legend('krig','krigcorrected');
    xlabel('percentile of bias hard');
    ylabel('mean lambda2');
    title(sprintf('lambda2 by bias hard, nsmax = %d',i));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/lam2vpercenthard_krigME_nsmax%d.png',i));

    close all;
        
end

end