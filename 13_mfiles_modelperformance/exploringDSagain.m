function [] = exploringDSagain()
% this function will examine the old downscaler results to refresh my
% memory

% load data
load('matfiles/allInfo.mat');

% all the variables
allkrig = {zk_hard_0km,zk_hard_100km,zk_hard_200km,zk_hard_300km, ...
    zk_hard_400km,zk_hard_500km,zk_hard_600km,zk_hard_700km,zk_hard_800km, ...
    zk_hard_900km};
allBME = {zk_soft_0km,zk_soft_100km,zk_soft_200km,zk_soft_300km, ...
    zk_soft_400km,zk_soft_500km,zk_soft_600km,zk_soft_700km,zk_soft_800km, ...
    zk_soft_900km};
allCAMP = {zk_camp_0km,zk_camp_100km,zk_camp_200km,zk_camp_300km, ...
    zk_camp_400km,zk_camp_500km,zk_camp_600km,zk_camp_700km,zk_camp_800km, ...
    zk_camp_900km};
allDS = {zk_STDSI_0km,zk_STDSI_100km,zk_STDSI_200km,zk_STDSI_300km, ...
    zk_STDSI_400km,zk_STDSI_500km,zk_STDSI_600km,zk_STDSI_700km,zk_STDSI_800km, ...
    zk_STDSI_900km};

% getting 2001
[yr mo da] = datevec(ck_base(:,3));

% calculcate statistics
% first column krig, second column BME, third column CAMP, fourth column DS
MSE_4methods = NaN*ones(10,4);
for i = 1:10
    
    idxnan3 = ~isnan(allDS{i});
    MSE_4methods(i,1) = mean((allkrig{i}(idxnan3)-zh(idxnan3)).^2);
    MSE_4methods(i,2) = mean((allBME{i}(idxnan3)-zh(idxnan3)).^2);
    MSE_4methods(i,3) = mean((allCAMP{i}(idxnan3)-zh(idxnan3)).^2);
    MSE_4methods(i,4) = mean((allDS{i}(idxnan3)-zh(idxnan3)).^2);
    
end

% make the increasing radius plots
figure; hold on;
plot(0:100:900,MSE_4methods(:,1),'bo-',0:100:900,MSE_4methods(:,2),'ro-', ...
    0:100:900,MSE_4methods(:,3),'ko-',0:100:900,MSE_4methods(:,4),'go-');
legend('Kriging','BME','CAMP','DS','location','best');
title('MSE for increasing radius for three methods');
xlabel('km');

% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures/DS_4methods.png');

end