function [] = plotXvalStats(soft,constant,gauss)
% this function will plot the differences in different crossvalidation
% statistics across the United States for 2001

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

forceddist = 0;

% looking at lambda2 and % soft kriging weights
load(sprintf('../matfiles/testingkrigME%d_nsmax%d_Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone%s.mat', ...
    0,3,softstr,constr,gaussstr,floor(forceddist./1000),''));
len = length(cellLamW);
matLamW = cell(len,1);
matLam2 = cell(len,1);
for i = 1:len
    temp1 = cellLamW{i};
    if ~isempty(temp1)
        matLamW{i} = cell2mat( cellfun( @(x) 100*sum(x(8:10))/sum(x), temp1, 'UniformOutput',false) );
    end
    temp2 = alllambda2{i};
    if ~isempty(temp2)   
        matLam2{i} = cell2mat( cellfun( @(x) mean(x), temp2, 'UniformOutput',false) );
    end
end

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckall = cell(size(cMS,1),1);
ckallH = cell(size(cMS,1),1);
for i = 1:size(cMS,1)   
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckall{i,1} = ch(idx & temp2001==2001,:);
    ckallH{i,1} = ch(idx,:);
end

% calculate distance to closest neighbor
for i = 1:size(cMS,1)
    distMS = sqrt( (cMS(i,1)-cMS(:,1)).^2 + (cMS(i,2)-cMS(:,2)).^2 );
    [a b] = sort(distMS);
    closeneib(i) = a(2)./1000;
end
idx = cellfun( @isempty, ckall );
closeneib(idx) = NaN;

% load hard data
load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm.mat', ...
    '_nosoft',constr,gaussstr,floor(forceddist./1000)));
zkall = zk_madd;
zhall = zh_Xval;
vkall = vk;     
zkallh = cell2mat(zkall);
zhallh = cell2mat(zhall); 
ckallh = cell2mat(ckallH); 
vkallh = cell2mat(vkall); 

% load soft data
load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone.mat', ...
        softstr,constr,gaussstr,floor(forceddist./1000)));
zkall = zk_madd;
zhall = zh_Xval;
vkall = vk;     
zkalls = cell2mat(zkall);
zhalls = cell2mat(zhall); 
ckalls = cell2mat(ckall); 
vkalls = cell2mat(vkall); 

temp = cell2mat( cellfun( @(x,y) ~isempty(x) & isempty(y), zkall, matLamW, 'UniformOutput', false ) );
temp = find(temp == 1);
for i = 1:length(temp)
    a = length(zkall{temp(i)});
    matLamW{temp(i)} = NaN*ones(a,1);
    matLam2{temp(i)} = NaN*ones(a,1);
end
matLamW = cell2mat(matLamW);
matLam2 = cell2mat(matLam2);

[lia lib] = ismember(cell2mat(ckall),cell2mat(ckallH),'rows');
lib(lib==0) = [];
zkallh = zkallh(lib); zhallh = zhallh(lib); ckallh = ckallh(lib,:); vkallh = vkallh(lib);

idx = ~isnan(zkallh) & ~isnan(zkalls);
zkallh = zkallh(idx); zhallh = zhallh(idx); ckallh = ckallh(idx,:); vkallh = vkallh(idx);
zkalls = zkalls(idx); zhalls = zhalls(idx); ckalls = ckalls(idx,:); vkalls = vkalls(idx);
matLamW = matLamW(idx); matLam2 = matLam2(idx);

lenvkh = (zkallh + 2*sqrt(vkallh)) - (zkallh - 2*sqrt(vkallh));
lenvks = (zkalls + 2*sqrt(vkalls)) - (zkalls - 2*sqrt(vkalls));
disp(sum(lenvks<=lenvkh)/length(lenvks));
disp([mean(lenvkh) mean(lenvks)]);

coverageh = zkallh - 2*sqrt(vkallh) <= zkallh & zhallh <= zkallh + 2*sqrt(vkallh);
coverages = zkalls - 2*sqrt(vkalls) <= zkalls & zhalls <= zkalls + 2*sqrt(vkalls);
disp(sum(coverageh)/length(coverageh));
disp(sum(coverages)/length(coverages));

% comparing to downscaler
len = length(unique(ckallh(:,3)));
unidays = unique(ckallh(:,3));
for i = 1:len
    idx = ckallh(:,1)>0 & ckallh(:,3) == unidays(i);
    zhwest(i) = mean(zhallh(idx));
end
unitemp = datevec(unidays);
[a b] = sort(zhwest);
daysDS = unitemp(b(end-2:end),:);
for i = 1:size(daysDS,1)
    idx = ckallh(:,1)>0 & ckallh(:,3)==datenum(daysDS(i,:));
    errorsh = zkallh(idx) - zhallh(idx);
    DSMSEh(i) = mean(errorsh.^2);
    DSMAEh(i) = mean(abs(errorsh));
    DScovh(i) = 100*sum(coverageh(idx))/length(coverageh(idx));
    DSlenvkh(i) = mean(lenvkh(idx));
    DSavgvkh(i) = mean(vkallh(idx));
    errorss = zkalls(idx) - zhalls(idx);
    DSMSEs(i) = mean(errorss.^2);
    DSMAEs(i) = mean(abs(errorss));
    DScovs(i) = 100*sum(coverages(idx))/length(coverages(idx));
    DSlenvks(i) = mean(lenvks(idx));
    DSavgvks(i) = mean(vkalls(idx));
end
statsDS = [ DSMSEh(1) DSMAEh(1) DScovh(1) DSlenvkh(1) DSavgvkh(1) ; ...
    DSMSEs(1) DSMAEs(1) DScovs(1) DSlenvks(1) DSavgvks(1) ; ...
    DSMSEh(2) DSMAEh(2) DScovh(2) DSlenvkh(2) DSavgvkh(2) ; ...
    DSMSEs(2) DSMAEs(2) DScovs(2) DSlenvks(2) DSavgvks(2) ; ...
    DSMSEh(3) DSMAEh(3) DScovh(3) DSlenvkh(3) DSavgvkh(3) ; ...
    DSMSEs(3) DSMAEs(3) DScovs(3) DSlenvks(3) DSavgvks(3) ];

% calculate statistics for hard data for each station
for i = 1:size(cMS,1)
 
    idx = ckalls(:,1) == cMS(i,1) & ckalls(:,2) == cMS(i,2);
    if sum(idx) == 0
        n_siteh(i) = NaN; n_h(i) = NaN; RMSEh(i) = NaN; MAEh(i) = NaN;
        MEh(i) = NaN; r2h(i) = NaN; MSh(i) = NaN; RMSSh(i) = NaN;
        MRh(i) = NaN; std_esth(i) = NaN; std_obsh(i) = NaN; QAr2h(i) = NaN;
        n_sites(i) = NaN; n_s(i) = NaN; RMSEs(i) = NaN; MAEs(i) = NaN;
        MEs(i) = NaN; r2s(i) = NaN; MSs(i) = NaN; RMSSs(i) = NaN;
        MRs(i) = NaN; std_ests(i) = NaN; std_obss(i) = NaN; QAr2s(i) = NaN;
        matLamWs(i) = NaN; matLam2s(i) = NaN;
    else 
        n_siteh(i) = size(unique(ckallh(idx,1:2),'rows'),1);
        n_h(i) = size(ckallh(idx,:),1);
        errorsh = zkallh(idx) - zhallh(idx);
        RMSEh(i) = sqrt( mean(errorsh.^2) );
        MAEh(i) = mean(abs(errorsh));
        MEh(i) = mean(errorsh);
        r2h(i) = (corr(zkallh(idx),zhallh(idx),'type','Pearson')).^2;
        MSh(i) = mean(errorsh./sqrt(vkallh(idx)));
        RMSSh(i) = std(errorsh./sqrt(vkallh(idx)));
        MRh(i) = mean(sqrt(vkallh(idx)));
        std_esth(i) = std(zkallh(idx));
        std_obsh(i) = std(zhallh(idx));
        QAr2h(i) = ( ( std_esth(i).^2 + std_obsh(i).^2 - (RMSEh(i).^2-MEh(i).^2) )./(2.*std_esth(i).*std_obsh(i)) ).^2;
        n_sites(i) = size(unique(ckalls(idx,1:2),'rows'),1);
        n_s(i) = size(ckalls(idx,:),1);
        errorss = zkalls(idx) - zhalls(idx);
        RMSEs(i) = sqrt( mean(errorss.^2) );
        MAEs(i) = mean(abs(errorss));
        MEs(i) = mean(errorss);
        r2s(i) = (corr(zkalls(idx),zhalls(idx),'type','Pearson')).^2;
        MSs(i) = mean(errorss./sqrt(vkalls(idx)));
        RMSSs(i) = std(errorss./sqrt(vkalls(idx)));
        MRs(i) = mean(sqrt(vkalls(idx)));
        std_ests(i) = std(zkalls(idx));
        std_obss(i) = std(zhalls(idx));
        QAr2s(i) = ( ( std_ests(i).^2 + std_obss(i).^2 - (RMSEs(i).^2-MEs(i).^2) )./(2.*std_ests(i).*std_obss(i)) ).^2;
        matLamWs(i) = mean(matLamW(idx));
        matLam2s(i) = mean(matLam2(idx));
    end
end

% save results
save(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s_foriso%dkm_timezone_Xvalresults.mat', ...
    constr,gaussstr,floor(forceddist./1000)),'cMS','n_siteh','n_h',...
    'RMSEh','MAEh','MEh','r2h','MSh','RMSSh','MRh','std_esth','std_obsh',...
    'QAr2h','n_sites','n_s','RMSEs','MAEs','MEs','r2s','MSs','RMSSs','MRs',...
    'std_ests','std_obss','QAr2s','matLamWs','matLam2s');

% count to see number of times each station was better

% have maps showing bias, r2, rmse, mr, rmss, ms
n = 1;
for i = 1:9
    
    if i == 1 
        toplot1 = MEh; toplot2 = MEs; toplot3 = MEs - MEh; toplotstr = 'Mean Error';
    elseif i == 2
        toplot1 = r2h; toplot2 = r2s; toplot3 = r2s - r2h; toplotstr = 'R2';
    elseif i == 3
        toplot1 = RMSEh; toplot2 = RMSEs; toplot3 = RMSEs - RMSEh; toplotstr = 'RMSE';
    elseif i == 4
        toplot1 = MRh; toplot2 = MRs; toplot3 = MRs - MRh; toplotstr = 'MR';
    elseif i == 5
        toplot1 = RMSSh; toplot2 = RMSSs; toplot3 = RMSSs - RMSSh; toplotstr = 'RMSS';
    elseif i == 6
        toplot1 = MSh; toplot2 = MSs; toplot3 = MSs - MSh; toplotstr = 'MS';
    elseif i == 7
        toplot1 = matLamWs; toplot2 = matLamWs; toplot3 = matLamWs; toplotstr = 'krig weights'; 
    elseif i == 8;
        toplot1 = matLam2s; toplot2 = matLam2s; toplot3 = matLam2s; toplotstr = 'lambda2';
    else
        toplot1 = closeneib; toplot2 = closeneib; toplot3 = closeneib; toplotstr = 'closest neib';
    end
    
    for j = 1:3
        
        if j == 1, toplot = toplot1; toplotstrsub = 'hard';
        elseif j == 2, toplot = toplot2; toplotstrsub = 'soft';
        elseif j == 3, toplot = toplot3; toplotstrsub = 'difference'; end
        if i == 7 | i == 8 | i == 9, toplotstrsub = ''; end
        
        figure; hold on;
        hist(toplot,100);
        title(sprintf('histogram of %s %s by station',toplotstrsub,toplotstr));
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpng','-r600',sprintf('hist_%s_%s.png',toplotstrsub,toplotstr));

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
        if j == 1 | j == 2
            fixed = prctile(toplot1,0:10:100);
            cax = [prctile(toplot1,0) prctile(toplot1,100)];
        else
            fixed = prctile(toplot,0:10:100);
            cax = [prctile(toplot,0) prctile(toplot,100)];
        end
        colorplot_fixed(cMS,toplot',hot(10),Property,Value,cax,fixed);       
        caxis(cax);
        hcb = colorbar;
        set(hcb,'ytick',linspace(cax(1),cax(2),length(fixed)));
        set(hcb,'yticklabel',arrayfun(@num2str,fixed,'uni',false));
        title(sprintf('map of %s %s in 2001',toplotstrsub,toplotstr));
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpng','-r600',sprintf('map_%s_%s.png',toplotstrsub,toplotstr)); 
        
        prcStats(:,n) = prctile(toplot,0:5:100);
        n = n + 1;
    
    end

end

blah = 5;

end