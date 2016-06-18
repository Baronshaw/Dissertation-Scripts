function [] = RbetweenMethods()
% this function will calculcate the R2 between the true lambda1/2 and the
% simulated lambda1/2 for each of the methods

% pick day
dayz = datenum(2001,7,1);

% load CMAQ
load('recalculate_dailyvCTMorder.mat')

% load true
load('PM2p5_meanGivMod_yr2001.mat');

% load constant
load('recalculatelambdaGivMod_Constant.mat');
% R between true lambda1 and simulated lambda1
R_lambda1_con = real(corr(meanGivCMAQ_all(css(:,3)==dayz),lambda1_recalcsimuModGrid(css(:,3)==dayz))); 
% R between true lambda2 and simulated lambda2
R_lambda2_con = real(corr(varGivCMAQ_all(css(:,3)==dayz),lambda2_recalcsimuModGrid(css(:,3)==dayz)));
% R between true lambda1-CMAQ and simulated lambda1
R_error_con = real(corr(meanGivCMAQ_all(css(:,3)==dayz)-dailyCTMvorder(css(:,3)==dayz), ...
    lambda1_recalcsimuModGrid(css(:,3)==dayz)-dailyCTMvorder(css(:,3)==dayz))); 

% load camp
load('recalculatelambdaGivMod_CAMP.mat');
% R between true lambda1 and simulated lambda1
R_lambda1_camp = real(corr(meanGivCMAQ_all(css(:,3)==dayz),lambda1_recalcsimuModGrid(css(:,3)==dayz))); 
% R between true lambda2 and simulated lambda2
R_lambda2_camp = real(corr(varGivCMAQ_all(css(:,3)==dayz),lambda2_recalcsimuModGrid(css(:,3)==dayz)));
% R between true lambda1-CMAQ and simulated lambda1
R_error_camp = real(corr(meanGivCMAQ_all(css(:,3)==dayz)-dailyCTMvorder(css(:,3)==dayz), ...
    lambda1_recalcsimuModGrid(css(:,3)==dayz)-dailyCTMvorder(css(:,3)==dayz))); 

% load ramp
load('recalculatelambdaGivMod.mat');
% R between true lambda1 and simulated lambda1
R_lambda1_ramp = corr(meanGivCMAQ_all(css(:,3)==dayz),lambda1_recalcsimuModGrid(css(:,3)==dayz)); 
% R between true lambda2 and simulated lambda2
R_lambda2_ramp = corr(varGivCMAQ_all(css(:,3)==dayz),lambda2_recalcsimuModGrid(css(:,3)==dayz));
% R between true lambda1-CMAQ and simulated lambda1
R_error_ramp = corr(meanGivCMAQ_all(css(:,3)==dayz)-dailyCTMvorder(css(:,3)==dayz), ...
    lambda1_recalcsimuModGrid(css(:,3)==dayz)-dailyCTMvorder(css(:,3)==dayz)); 

% save results
save('RbetweenMethods.mat','R_lambda1_con','R_lambda2_con','R_error_con', ...
    'R_lambda1_camp','R_lambda2_camp','R_error_camp','R_lambda1_ramp','R_lambda2_ramp','R_error_ramp');

end