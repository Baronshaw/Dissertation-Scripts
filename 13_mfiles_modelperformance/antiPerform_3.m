function [] = antiPerform_3()
% this function will create figures that show statistics by increasing
% radius with and without anticipating results

% increasing radius for R2 and MSE for ME2>90th% & in{22}
forceddistz = 0:100000:900000;


for i = 1:length(forceddistz)
    
    %%% getting info
    
    forceddist = forceddistz(i);
    
    % loading info
    load(sprintf('matfiles/allInfo_%dkm.mat',floor(forceddist./1000)));

    % get region
    if ~exist('matfiles/inregion.mat')
        allregions = shaperead('FWS_LCC/FWS_LCC.shp');
        for k = 1:length(allregions)
            tic
            cd ../09_mfiles_projections
            allregions_p{k,1} = ell2lambertcc([allregions(k).X',allregions(k).Y'],'whiproj2001');
            cd ../13_mfiles_modelperformance
            in{k,1} = inpolygon(ck_base(:,1),ck_base(:,2),allregions_p{k,1}(:,1),allregions_p{k,1}(:,2));
            toc
        end
        save('matfiles/inregion.mat','in','allregions');
    else
        load('matfiles/inregion.mat')
    end

    % get season
    [yr mo da] = datevec(ck_base(:,3));
    IsWinter = NaN*ones(length(ck_base),1); IsWinter(mo==1|mo==2|mo==12) = 1;
    IsSpring = NaN*ones(length(ck_base),1); IsSpring(mo==3|mo==4|mo==5) = 1;
    IsSummer = NaN*ones(length(ck_base),1); IsSummer(mo==6|mo==7|mo==8) = 1;
    IsFall = NaN*ones(length(ck_base),1); IsFall(mo==9|mo==10|mo==11) = 1;

    % get longitude
    cd ../09_mfiles_projections
    testing = [-125:5:-90]'; 
    M = ell2lambertcc([testing 40*ones(length(testing),1)],'whiproj2001');
    cd ../13_mfiles_modelperformance
    Is125 = NaN*ones(length(ck_base),1); Is125(ck_base(:,1)<M(1,1)) = 1;
    Is120 = NaN*ones(length(ck_base),1); Is120(ck_base(:,1)<M(2,1)) = 1;
    Is115 = NaN*ones(length(ck_base),1); Is115(ck_base(:,1)<M(3,1)) = 1;
    Is110 = NaN*ones(length(ck_base),1); Is110(ck_base(:,1)<M(4,1)) = 1;
    Is105 = NaN*ones(length(ck_base),1); Is105(ck_base(:,1)<M(5,1)) = 1;
    Is100 = NaN*ones(length(ck_base),1); Is100(ck_base(:,1)<M(6,1)) = 1;
    Is95 = NaN*ones(length(ck_base),1); Is95(ck_base(:,1)<M(7,1)) = 1;
    Is90 = NaN*ones(length(ck_base),1); Is90(ck_base(:,1)<M(8,1)) = 1; 

    %%% calculate with and without subsetting
    
    % overall
    overall_r2_kriging(i,1) = (corr(zk_hard,zh_base)).^2; % R2 kriging
    overall_r2_ramp(i,1) = (corr(zk_soft,zh_base)).^2; % R2 RAMP
    overall_mse_kriging(i,1) = mean((zk_hard-zh_base).^2); % MSE kriging
    overall_mse_ramp(i,1) = mean((zk_soft-zh_base).^2); % MSE RAMP
    
    % subset 1
    idx = in{14}==1 & ~isnan(Mod_2);
    sub1_r2_kriging(i,1) = (corr(zk_hard(idx),zh_base(idx))).^2; % R2 kriging
    sub1_r2_ramp(i,1) = (corr(zk_soft(idx),zh_base(idx))).^2; % R2 RAMP
    sub1_mse_kriging(i,1) = mean((zk_hard(idx)-zh_base(idx)).^2); % MSE kriging
    sub1_mse_ramp(i,1) = mean((zk_soft(idx)-zh_base(idx)).^2); % MSE RAMP
    
    % subset 2
    idx = (mBias_2).^2 > prctile((mBias_2).^2,90) & in{22}==1 & ~isnan(Mod_2);
    sub2_r2_kriging(i,1) = (corr(zk_hard(idx),zh_base(idx))).^2; % R2 kriging
    sub2_r2_ramp(i,1) = (corr(zk_soft(idx),zh_base(idx))).^2; % R2 RAMP
    sub2_mse_kriging(i,1) = mean((zk_hard(idx)-zh_base(idx)).^2); % MSE kriging
    sub2_mse_ramp(i,1) = mean((zk_soft(idx)-zh_base(idx)).^2); % MSE RAMP
        
end

% plot of R2
figure; hold on;
plot(floor(forceddistz./1000),overall_r2_kriging,'bo-',floor(forceddistz./1000),overall_r2_ramp,'ro-', ...
    floor(forceddistz./1000),sub1_r2_kriging,'bo--',floor(forceddistz./1000),sub1_r2_ramp,'ro--', ...
    floor(forceddistz./1000),sub2_r2_kriging,'bo:',floor(forceddistz./1000),sub2_r2_ramp,'ro:');
title('R2 of increasing radius for 2001');
xlabel('km');
ylabel('R2');
legend('krig','RAMP','sub1 krig','sub1 RAMP','sub2 krig','sub2 RAMP','Location','best');

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','figures/radius_r2_antiperform.png');

% plot of MSE
figure; hold on;
plot(floor(forceddistz./1000),overall_mse_kriging,'bo-',floor(forceddistz./1000),overall_mse_ramp,'ro-', ...
    floor(forceddistz./1000),sub1_mse_kriging,'bo--',floor(forceddistz./1000),sub1_mse_ramp,'ro--', ...
    floor(forceddistz./1000),sub2_mse_kriging,'bo:',floor(forceddistz./1000),sub2_mse_ramp,'ro:');
title('MSE of increasing radius for 2001');
xlabel('km');
ylabel('MSE');
legend('krig','RAMP','sub1 krig','sub1 RAMP','sub2 krig','sub2 RAMP','Location','best');

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','figures/radius_mse_antiperform.png');

end