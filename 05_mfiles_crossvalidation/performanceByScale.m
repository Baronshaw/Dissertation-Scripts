function [] = performanceByScale()
% this function will see how soft data estimates perform by the scale of
% the monitors used to construct the soft data

% bring in scale information for each monitor
load('../datafiles/Observed_PM2p5/MasterDaily_PM2p5_2001.mat');
toscale = unique([longitude latitude],'rows');
[lia lib] = ismember(toscale,[longitude latitude],'rows');
% projection location data
cd ../09_mfiles_projections
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load Ellipsoids.mat
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
coordObs = ell2lambertcc(toscale,'whiproj2001'); 
cd ../05_mfiles_crossvalidation
toscale = [scale(lib) coordObs];

% scale info
scaleinfo = [NaN NaN;NaN NaN;100 500;500 4000;4000 50000;50000 10000000];
scaleinfonum = 1:6;

% look at all distance from each monitor to the three closest stations
% what is the scale induced? were there other monitors within range that
% were excluded?
X = toscale(:,2:3)'; Y = X;
D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
[sorted sortidx] = sort(D,2);
inducedScale = NaN*ones(length(toscale),1);
monitorsInRange = NaN*ones(length(toscale),1);
for i = 1:length(toscale) 
    if toscale(i,1) ~= 1
        a = scaleinfo(toscale(i,1),:);
        b = sorted(i,3); % furthest monitor
        idx = b >= scaleinfo(:,1) & b <= scaleinfo(:,2);
        inducedScale(i) = scaleinfonum(idx);
        monitorsInRange(i) = sum(sorted(i,:)<=a(2));        
    end    
end

% look at the performance of the hard and soft estimates by scale/induced
% scale
softstr = '_soft';  constr = '_long';  gaussstr = '_gauss'; forceddist = 0;

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

[lia lib] = ismember(cell2mat(ckall),cell2mat(ckallH),'rows');
lib(lib==0) = [];
zkallh = zkallh(lib); zhallh = zhallh(lib); ckallh = ckallh(lib,:); vkallh = vkallh(lib);

idx = ~isnan(zkallh) & ~isnan(zkalls);
zkallh = zkallh(idx); zhallh = zhallh(idx); ckallh = ckallh(idx,:); vkallh = vkallh(idx);
zkalls = zkalls(idx); zhalls = zhalls(idx); ckalls = ckalls(idx,:); vkalls = vkalls(idx);

% calculate statistics for each station
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
    end
end

% fill in
toscaleall = NaN*ones(length(ckallh),1);
inducedScaleall = NaN*ones(length(ckallh),1);
for i = 1:length(toscale)
    idx = toscale(i,2)==ckalls(:,1) & toscale(i,3)==ckalls(:,2);
    toscaleall(idx) = toscale(i,1);
    inducedScaleall(idx) = inducedScale(i);
end

% stats overall
for i = 1:6
 
    idx1 = toscaleall == i;
    idx2 = inducedScaleall == i;
    if sum(idx1) == 0
        n_siteh1(i) = NaN; n_h1(i) = NaN; RMSEh1(i) = NaN; MAEh1(i) = NaN;
        MEh1(i) = NaN; r2h1(i) = NaN; MSh1(i) = NaN; RMSSh1(i) = NaN;
        MRh1(i) = NaN; std_esth1(i) = NaN; std_obsh1(i) = NaN; QAr2h1(i) = NaN;
        n_sites1(i) = NaN; n_s1(i) = NaN; RMSEs1(i) = NaN; MAEs1(i) = NaN;
        MEs1(i) = NaN; r2s1(i) = NaN; MSs1(i) = NaN; RMSSs1(i) = NaN;
        MRs1(i) = NaN; std_ests1(i) = NaN; std_obss1(i) = NaN; QAr2s1(i) = NaN;
        matLamWs1(i) = NaN; matLam2s1(i) = NaN;
    else 
        n_siteh1(i) = size(unique(ckallh(idx1,1:2),'rows'),1);
        n_h1(i) = size(ckallh(idx1,:),1);
        errorsh1 = zkallh(idx1) - zhallh(idx1);
        RMSEh1(i) = sqrt( mean(errorsh1.^2) );
        MAEh1(i) = mean(abs(errorsh1));
        MEh1(i) = mean(errorsh1);
        r2h1(i) = (corr(zkallh(idx1),zhallh(idx1),'type','Pearson')).^2;
        MSh1(i) = mean(errorsh1./sqrt(vkallh(idx1)));
        RMSSh1(i) = std(errorsh1./sqrt(vkallh(idx1)));
        MRh1(i) = mean(sqrt(vkallh(idx1)));
        std_esth1(i) = std(zkallh(idx1));
        std_obsh1(i) = std(zhallh(idx1));
        QAr2h1(i) = ( ( std_esth1(i).^2 + std_obsh1(i).^2 - (RMSEh1(i).^2-MEh1(i).^2) )./(2.*std_esth1(i).*std_obsh1(i)) ).^2;
        n_sites1(i) = size(unique(ckalls(idx1,1:2),'rows'),1);
        n_s1(i) = size(ckalls(idx1,:),1);
        errorss1 = zkalls(idx1) - zhalls(idx1);
        RMSEs1(i) = sqrt( mean(errorss1.^2) );
        MAEs1(i) = mean(abs(errorss1));
        MEs1(i) = mean(errorss1);
        r2s1(i) = (corr(zkalls(idx1),zhalls(idx1),'type','Pearson')).^2;
        MSs1(i) = mean(errorss1./sqrt(vkalls(idx1)));
        RMSSs1(i) = std(errorss1./sqrt(vkalls(idx1)));
        MRs1(i) = mean(sqrt(vkalls(idx1)));
        std_ests1(i) = std(zkalls(idx1));
        std_obss1(i) = std(zhalls(idx1));
        QAr2s1(i) = ( ( std_ests1(i).^2 + std_obss1(i).^2 - (RMSEs1(i).^2-MEs1(i).^2) )./(2.*std_ests1(i).*std_obss1(i)) ).^2;
    end
    if sum(idx2) == 0
        n_siteh2(i) = NaN; n_h2(i) = NaN; RMSEh2(i) = NaN; MAEh2(i) = NaN;
        MEh2(i) = NaN; r2h2(i) = NaN; MSh2(i) = NaN; RMSSh2(i) = NaN;
        MRh2(i) = NaN; std_esth2(i) = NaN; std_obsh2(i) = NaN; QAr2h2(i) = NaN;
        n_sites2(i) = NaN; n_s2(i) = NaN; RMSEs2(i) = NaN; MAEs2(i) = NaN;
        MEs2(i) = NaN; r2s2(i) = NaN; MSs2(i) = NaN; RMSSs2(i) = NaN;
        MRs2(i) = NaN; std_ests2(i) = NaN; std_obss2(i) = NaN; QAr2s2(i) = NaN;
        matLamWs2(i) = NaN; matLam2s2(i) = NaN;
    else 
        n_siteh2(i) = size(unique(ckallh(idx2,1:2),'rows'),1);
        n_h2(i) = size(ckallh(idx2,:),1);
        errorsh2 = zkallh(idx2) - zhallh(idx2);
        RMSEh2(i) = sqrt( mean(errorsh2.^2) );
        MAEh2(i) = mean(abs(errorsh2));
        MEh2(i) = mean(errorsh2);
        r2h2(i) = (corr(zkallh(idx2),zhallh(idx2),'type','Pearson')).^2;
        MSh2(i) = mean(errorsh2./sqrt(vkallh(idx2)));
        RMSSh2(i) = std(errorsh2./sqrt(vkallh(idx2)));
        MRh2(i) = mean(sqrt(vkallh(idx2)));
        std_esth2(i) = std(zkallh(idx2));
        std_obsh2(i) = std(zhallh(idx2));
        QAr2h2(i) = ( ( std_esth2(i).^2 + std_obsh2(i).^2 - (RMSEh2(i).^2-MEh2(i).^2) )./(2.*std_esth2(i).*std_obsh2(i)) ).^2;
        n_sites2(i) = size(unique(ckalls(idx2,1:2),'rows'),1);
        n_s2(i) = size(ckalls(idx2,:),1);
        errorss2 = zkalls(idx2) - zhalls(idx2);
        RMSEs2(i) = sqrt( mean(errorss2.^2) );
        MAEs2(i) = mean(abs(errorss2));
        MEs2(i) = mean(errorss2);
        r2s2(i) = (corr(zkalls(idx2),zhalls(idx2),'type','Pearson')).^2;
        MSs2(i) = mean(errorss2./sqrt(vkalls(idx2)));
        RMSSs2(i) = std(errorss2./sqrt(vkalls(idx2)));
        MRs2(i) = mean(sqrt(vkalls(idx2)));
        std_ests2(i) = std(zkalls(idx2));
        std_obss2(i) = std(zhalls(idx2));
        QAr2s2(i) = ( ( std_ests2(i).^2 + std_obss2(i).^2 - (RMSEs2(i).^2-MEs2(i).^2) )./(2.*std_ests2(i).*std_obss2(i)) ).^2;
    end
end

% stats overall by changeing scale
for i = 1:6
 for j = 1:6
    idx = toscaleall == i & inducedScaleall == j;
    if sum(idx) == 0
        n_sitehij(i,j) = NaN; n_hij(i,j) = NaN; RMSEhij(i,j) = NaN; MAEhij(i,j) = NaN;
        MEhij(i,j) = NaN; r2hij(i,j) = NaN; MShij(i,j) = NaN; RMSShij(i,j) = NaN;
        MRhij(i,j) = NaN; std_esthij(i,j) = NaN; std_obshij(i,j) = NaN; QAr2hij(i,j) = NaN;
        n_sitesij(i,j) = NaN; n_sij(i,j) = NaN; RMSEsij(i,j) = NaN; MAEsij(i,j) = NaN;
        MEsij(i,j) = NaN; r2sij(i,j) = NaN; MSsij(i,j) = NaN; RMSSsij(i,j) = NaN;
        MRsij(i,j) = NaN; std_estsij(i,j) = NaN; std_obssij(i,j) = NaN; QAr2sij(i,j) = NaN;
        matLamWsij(i,j) = NaN; matLam2sij(i,j) = NaN;
    else 
        n_sitehij(i,j) = size(unique(ckallh(idx,1:2),'rows'),1);
        n_hij(i,j) = size(ckallh(idx,:),1);
        errorshij = zkallh(idx) - zhallh(idx);
        RMSEhij(i,j) = sqrt( mean(errorshij.^2) );
        MAEhij(i,j) = mean(abs(errorshij));
        MEhij(i,j) = mean(errorshij);
        r2hij(i,j) = (corr(zkallh(idx),zhallh(idx),'type','Pearson')).^2;
        MShij(i,j) = mean(errorshij./sqrt(vkallh(idx)));
        RMSShij(i,j) = std(errorshij./sqrt(vkallh(idx)));
        MRhij(i,j) = mean(sqrt(vkallh(idx)));
        std_esthij(i,j) = std(zkallh(idx));
        std_obshij(i,j) = std(zhallh(idx));
        QAr2hij(i,j) = ( ( std_esthij(i,j).^2 + std_obshij(i,j).^2 - (RMSEhij(i,j).^2-MEhij(i,j).^2) )./(2.*std_esthij(i,j).*std_obshij(i,j)) ).^2;
        n_sitesij(i,j) = size(unique(ckalls(idx,1:2),'rows'),1);
        n_sij(i,j) = size(ckalls(idx,:),1);
        errorssij = zkalls(idx) - zhalls(idx);
        RMSEsij(i,j) = sqrt( mean(errorssij.^2) );
        MAEsij(i,j) = mean(abs(errorssij));
        MEsij(i,j) = mean(errorssij);
        r2sij(i,j) = (corr(zkalls(idx),zhalls(idx),'type','Pearson')).^2;
        MSsij(i,j) = mean(errorssij./sqrt(vkalls(idx)));
        RMSSsij(i,j) = std(errorssij./sqrt(vkalls(idx)));
        MRsij(i,j) = mean(sqrt(vkalls(idx)));
        std_estsij(i,j) = std(zkalls(idx));
        std_obssij(i,j) = std(zhalls(idx));
        QAr2sij(i,j) = ( ( std_estsij(i,j).^2 + std_obssij(i,j).^2 - (RMSEsij(i,j).^2-MEsij(i,j).^2) )./(2.*std_estsij(i,j).*std_obssij(i,j)) ).^2;
    end
 end
end

% put results in a text file
allstats1 = [n_h1' n_s1' n_siteh1' n_sites1' RMSEh1' RMSEs1' MAEh1' MAEs1' MEh1' MEs1' r2h1' r2s1' ...
    MSh1' MSs1' RMSSh1' RMSSs1' MRh1' MRs1' std_esth1' std_ests1' std_obsh1' std_obss1' QAr2h1' QAr2s1'];
allstats2 = [n_h2' n_s2' n_siteh2' n_sites2' RMSEh2' RMSEs2' MAEh2' MAEs2' MEh2' MEs2' r2h2' r2s2' ...
    MSh2' MSs2' RMSSh2' RMSSs2' MRh2' MRs2' std_esth2' std_ests2' std_obsh2' std_obss2' QAr2h2' QAr2s2'];

% plot results by statistic
stattrystr = {'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR'};
plotstr1 = {'b.-' ; 'r.-' ; 'g.-' ; 'm.-' ; 'y.-' ; 'k.-'};
plotstr2 = {'bo-' ; 'ro-' ; 'go-' ; 'mo-' ; 'yo-' ; 'ko-'};
for i = 1:length(stattrystr)
    figure; hold on;
    if i == 1, toplot1 = RMSEhij; toplot2 = RMSEsij; 
    elseif i == 2, toplot1 = MAEhij; toplot2 = MAEsij; 
    elseif i == 3, toplot1 = MEhij; toplot2 = MEsij;
    elseif i == 4, toplot1 = r2hij; toplot2 = r2sij;
    elseif i == 5, toplot1 = MShij; toplot2 = MSsij;
    elseif i == 6, toplot1 = RMSShij; toplot2 = RMSSsij;
    elseif i == 7, toplot1 = MRhij; toplot2 = MRsij; end
    for j = 1:6
        plot(1:6,toplot1(j,:),plotstr1{j});
        plot(1:6,toplot2(j,:),plotstr2{j});
    end
    legend('1h','1s','2h','2s','3h','3s','4h','4s','5h','5s','6h','6s','Location','best');
    title(sprintf('Scale vs Induced Scale, %s',stattrystr{i}));
    xlabel('Induced Scale');
    ylabel(sprintf('%s',stattrystr{i}));
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('OverallScale_%s.png',stattrystr{i})); 
end

blah = 5;

end