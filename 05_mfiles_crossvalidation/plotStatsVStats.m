function [] = plotStatsVStats()
% this function will plot all the statistics by station one versus the
% other in hopes of seeing some sort of performance pattern

soft = 1; constant = 0; gauss = 1;

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

% save results by stations
save('../matfiles/Xvalresultsbystation.mat','cMS', ...
    'n_siteh','n_h','RMSEh','MAEh','MEh','r2h','MSh','RMSSh', ...
    'MRh','std_esth','std_obsh','QAr2h','n_sites','n_s','RMSEs', ...
    'MAEs','MEs','r2s','MSs','RMSSs','MRs','std_ests','std_obss', ...
    'QAr2s','matLamWs','matLam2s');

% plot statistics versus statistics
stattry = [ RMSEs'-RMSEh' MAEs'-MAEh' MEs'-MEh' r2s'-r2h' MSs'-MSh' ...
    RMSSs'-RMSSh' MRs'-MRh' matLam2s' closeneib' matLamWs' ];
stattrystr = {'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'lambda2' ; 'closest' ; 'krig weight'};
idx1 = cMS(:,1) <= -1500000; idx2 = cMS(:,1) > -1500000 & cMS(:,1) <= -500000;
idx3 = cMS(:,1) > -500000 & cMS(:,1) <= 500000; idx4 = cMS(:,1) > 500000;

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
plot(cMS(idx1,1),cMS(idx1,2),'r.'); plot(cMS(idx2,1),cMS(idx2,2),'b.');
plot(cMS(idx3,1),cMS(idx3,2),'g.'); plot(cMS(idx4,1),cMS(idx4,2),'c.');
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','country_divide.png'); 

for i = 1:length(stattrystr)-1
    for j = i+1:length(stattrystr)
        
        figure; hold on;
        % plot(stattry(:,i),stattry(:,j),'r.');
        plot(stattry(idx1,i),stattry(idx1,j),'r.');
        plot(stattry(idx2,i),stattry(idx2,j),'b.');
        plot(stattry(idx3,i),stattry(idx3,j),'g.');
        plot(stattry(idx4,i),stattry(idx4,j),'c.');
        title(sprintf('%s versus %s',stattrystr{i},stattrystr{j}));
        xlabel(sprintf('%s',stattrystr{i}));
        ylabel(sprintf('%s',stattrystr{j}));
        legend('first','second','third','fourth');
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpng','-r600',sprintf('%s_vs_%s.png',stattrystr{i},stattrystr{j})); 
        
    end
end

end