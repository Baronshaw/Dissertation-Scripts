function [] = testingCovariance()
% this function will test the stationarity assumption of the covariance
% model. This will look at covariance models regionally (by EPA region). 
% I will also test to see if there is connection between variance of a 
% station and RMSS

% load observed data (with mean trend removed)
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
zh_rmv = zd - mI;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckall = cell(size(cMS,1),1);
ckallH = cell(size(cMS,1),1);
for i = 1:size(cMS,1)   
    idx = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckall{i,1} = ch(idx & temp2001==2001,:);
    ckallH{i,1} = ch(idx,:);
end

% assign EPA regions to the data
temp = datevec(ch(:,3));
uniyr = unique(temp(:,1));
cd Observed_PM2p5
for i = 1:length(uniyr)
    load(sprintf('MasterDaily_PM2p5_%d.mat',uniyr(i)));
    alllon{i,1} = longitude;
    alllat{i,1} = latitude;
    allID{i,1} = location;
end
alllon = cell2mat(alllon);
alllat = cell2mat(alllat);
allID = cell2mat(allID);
[uniloc uniidx] = unique([alllon alllat],'rows');
unilon = uniloc(:,1);
unilat = uniloc(:,2);
uniID = allID(uniidx);
cd ..

% convert them to coordinate system
load projexample.mat;
cd 09_mfiles_projections
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load Ellipsoids.mat
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
projuniID = ell2lambertcc([unilon unilat],'whiproj2001'); 
cd ..
% match estimations locations with station id
sIDall = NaN*ones(size(ch,1),1);
for i = 1:length(projuniID)
    idx = projuniID(i,1) == ch(:,1) & projuniID(i,2) == ch(:,2);
    sIDall(idx) = uniID(i);
end
stateall = floor(sIDall./10^7);
% match id's with EPA region
% from https://aqs.epa.gov/aqsweb/codes/data/StateCountyCodes.csv
M = csvread('StateCountyCodes_mod.csv');
temp = unique(M,'rows');
stateIDs = temp(:,1);
EPAregionsIDs = temp(:,2);
EPAall = NaN*ones(size(cell2mat(ckallH),1),1);
for i = 1:length(stateIDs)
    idx = stateIDs(i) == stateall;
    EPAall(idx) = EPAregionsIDs(i);
end
uniEPA = unique(EPAall);

% getting EPA region for each unique ID
stateuni = floor(uniID./10^7);
EPAuni = NaN*ones(size(ckallH,1),1);
for i = 1:length(stateIDs)
    idx = stateIDs(i) == stateuni;
    EPAuni(idx) = EPAregionsIDs(i);
end

% call the automatic covariance fitting function, assume that the
% covariance model is fixed and the parameters (range and sill) are
% variable
for i = 1:length(uniEPA)
    idx = uniEPA(i) == EPAall;
    [Z_mrmvd,cMS,tME,nanratio]=valstv2stg(ch(idx,:),zh_rmv(idx));
    DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));
    DME = abs(bsxfun(@minus,tME,tME'))';
    rLag = [0 prctile(unique(DMS(:)),[0.25 0.5 0.75 1 1.5 2:10 12.5 15:5:50])];
    rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];
    tLag = [0:10 15:5:40 50:25:150];
    tTol = [0 repmat(0.5,1,10) repmat(2.5,1,6) repmat(25,1,5)];
    [Crtest nprtest]=stcov(Z_mrmvd,cMS,tME,Z_mrmvd,cMS,tME,rLag,rTol,0,0);
    [Cttest npttest]=stcov(Z_mrmvd,cMS,tME,Z_mrmvd,cMS,tME,0,0,tLag,tTol);
    Cr1 = Crtest(1); 
    x = [rLag' ; zeros(length(tLag),1)];
    y = [zeros(length(rLag),1) ; tLag'];
    z = [Crtest ; Cttest'];
    s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
        'Upper',[1,5000000,5000000,Inf,Inf],'Startpoint',[0.5,3000000,3000000,500,500]);
    cd ../03_mfiles_covariance/covariance_models
    g = fittype( 'jointexpexp(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s,...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'} );    
    [f{i} gof{i} output{i}] = fit([x,y],z,g,'problem',Cr1);
    cd ../../05_mfiles_crossvalidation
    disp(f{i});
    disp(gof{i}.rsquare);
end

% load hard data
softstr = '_soft';  constr = '_long';  gaussstr = '_gauss'; forceddist = 0;
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

% getting the EPA region for the soft data
[lia lib] = ismember(cell2mat(ckall),ch,'rows');
EPAall = EPAall(lib);
uniEPA = unique(EPAall);

% calculate the RMSS for each EPA region
for i = 1:length(uniEPA)
 
    idx = uniEPA(i) == EPAall;
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

% create table showing parameter values by region with RMSS, is there a
% connection between covariance parameters and RMSS?
for i = 1:10
    params(i,:) = [f{i}.Cr1 f{i}.alp f{i}.ar1 f{i}.ar2 f{i}.at1 f{i}.at2];
    weightedranges(i,:) = [f{i}.alp*f{i}.ar1+(1-f{i}.alp)*f{i}.ar2 f{i}.alp*f{i}.at1+(1-f{i}.alp)*f{i}.at2];
end
allstats = [uniEPA params weightedranges RMSSh' RMSSs'];
allstats1 = [n_h' n_s' n_siteh' n_sites' RMSEh' RMSEs' MAEh' MAEs' MEh' MEs' r2h' r2s' ...
    MSh' MSs' RMSSh' RMSSs' MRh' MRs' std_esth' std_ests' std_obsh' std_obss' QAr2h' QAr2s'];

figure; hold on;
plot(params(:,1),RMSSh','bo',params(:,1),RMSSs','ro');
legend('hard','soft');
title('Variance vs RMSS');
xlabel('variance');
ylabel('RMSS');
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','VarvRMSS.png'); 

for i = 1:2
    figure; hold on;
    if i == 1, plotstr = 'ar'; else plotstr = 'at'; end
    plot(weightedranges(:,i),RMSSh','bo',weightedranges(:,i),RMSSs','ro');
    legend('hard','soft');
    title(sprintf('%s vs RMSS',plotstr));
    xlabel(sprintf('%s',plotstr));
    ylabel('RMSS');
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('%svRMSS.png',plotstr)); 
end

% saving results
cMS = unique(ch(:,1:2),'rows');
save('RegionByID.mat','cMS','EPAuni','f');

blah = 5;

end