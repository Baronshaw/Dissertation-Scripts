function [] = displayCovariancePlots(FOLDIDX)
% this function will display the new covariance plots for each fold on top
% of the previous covariance function

if nargin < 1, FOLDIDX = 1; end
    
covmodel = {'exponentialC/exponentialC','exponentialC/exponentialC'};

% load old covariance parameters
load(sprintf('../matfiles/expcov_%s.mat','long')); 
rLag_old = rLag; Crtest_old = Crtest; tLag_old = tLag; Cttest_old = Cttest;     
load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
covparam_old = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};   

% load new covariance parameters
load(sprintf('../matfiles/covmod_true10fold_%d.mat',FOLDIDX));
rLag_new = rLag; Crtest_new = Crtest; tLag_new = tLag; Cttest_new = Cttest;   
covparam_new = {[f.Cr1*f.alp f.ar1 f.at1] [f.Cr1*(1-f.alp) f.ar2 f.at2]};

% assign x and y values
xr = 0:0.01:rLag_old(end);
yr_old = covparam_old{1}(1).*exp(-3*xr./covparam_old{1}(2)) + ...
    covparam_old{2}(1).*exp(-3*xr./covparam_old{2}(2));
yr_new = covparam_new{1}(1).*exp(-3*xr./covparam_new{1}(2)) + ...
    covparam_new{2}(1).*exp(-3*xr./covparam_new{2}(2));
xt = 0:0.01:tLag_old(end);
yt_old = covparam_old{1}(1).*exp(-3*xt./covparam_old{1}(3)) + ...
    covparam_old{2}(1).*exp(-3*xt./covparam_old{2}(3));
yt_new = covparam_new{1}(1).*exp(-3*xt./covparam_new{1}(3)) + ...
    covparam_new{2}(1).*exp(-3*xt./covparam_new{2}(3));

% diplaying results
figure; hold on;
plot(rLag_old,Crtest_old,'bo');
plot(xr,yr_old,'b-');
plot(rLag_new,Crtest_new,'ro');
plot(xr,yr_new,'r-');
legend('Prev Exp','Pre Mod','Fold Exp','Fold Mod');
title(sprintf('Experimental and Modeled Covariance Old and New Space: fold %d',FOLDIDX));
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpdf','-r600',sprintf('CovR_true10fold_fold%d.pdf',FOLDIDX));
print(gcf,'-painters','-dpng','-r600',sprintf('CovR_true10fold_fold%d.png',FOLDIDX));

figure; hold on;
plot(tLag_old,Cttest_old,'bo');
plot(xt,yt_old,'b-');
plot(tLag_new,Cttest_new,'ro');
plot(xt,yt_new,'r-');
legend('Prev Exp','Pre Mod','Fold Exp','Fold Mod');
title(sprintf('Experimental and Modeled Covariance Old and New Time: fold %d',FOLDIDX));
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpdf','-r600',sprintf('CovT_true10fold_fold%d.pdf',FOLDIDX));
print(gcf,'-painters','-dpng','-r600',sprintf('CovT_true10fold_fold%d.png',FOLDIDX));
   
end