function [] = stratCMAQRAMP()
% this function will calculate different statistics by different
% stratifications

% load all characteristic info
load('performanceInfo.mat')

% calculate statistics for each stratification using the following
% statistics: ME, SE, ME2, SE2, MSE, R, R2, RMSE, MNE, NME, MNAE, NMAE, ME,
% SE2, RMSS, MR

% overall 

% CMAQ
ME_cmaq = mean(Mod-Obs);
SE_cmaq = std(Mod-Obs);
ME2_cmaq = ME_cmaq.^2;
SE2_cmaq = var(Mod-Obs);
MSE_cmaq = mean((Mod-Obs).^2);
R_cmaq = corr(Mod,Obs);
R2_cmaq = R_cmaq.^2;
RMSE_cmaq = sqrt(MSE_cmaq);
MNE_cmaq = mean((Mod(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NME_cmaq = ME_cmaq/mean(Obs);
MAE_cmaq = mean(abs(Mod-Obs));
MNAE_cmaq = mean(abs(Mod(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NMAE_cmaq = MAE_cmaq/mean(Obs);
ME_cmaq2 = NaN;
SE2_cmaq2 = NaN;
RMSS_cmaq = NaN;
MR_cmaq = NaN;
MSNS_cmaq = NaN;

% RAMP
ME_ramp = mean(lambda1-Obs);
SE_ramp = std(lambda1-Obs);
ME2_ramp = ME_ramp.^2;
SE2_ramp = var(lambda1-Obs);
MSE_ramp = mean((lambda1-Obs).^2);
R_ramp = corr(lambda1,Obs);
R2_ramp = R_ramp.^2;
RMSE_ramp = sqrt(MSE_ramp);
MNE_ramp = mean((lambda1(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NME_ramp = ME_ramp/mean(Obs);
MAE_ramp = mean(abs(lambda1-Obs));
MNAE_ramp = mean(abs(lambda1(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NMAE_ramp = MAE_ramp/mean(Obs);
ME_ramp2 = mean((Obs-lambda1)./sqrt(lambda2));
SE2_ramp2 = var((Obs-lambda1)./sqrt(lambda2));
RMSS_ramp = sqrt(mean( ((Obs-lambda1)./sqrt(lambda2)).^2 ));
MR_ramp = mean(sqrt(lambda2));
MSNS_ramp = mean( ( ((Obs-lambda1)./sqrt(lambda2)).^2 - 1 ).^2 );

% CAMP
ME_camp = mean(lambda1CAMP-Obs);
SE_camp = std(lambda1CAMP-Obs);
ME2_camp = ME_camp.^2;
SE2_camp = var(lambda1CAMP-Obs);
MSE_camp = mean((lambda1CAMP-Obs).^2);
R_camp = corr(lambda1CAMP,Obs);
R2_camp = R_camp.^2;
RMSE_camp = sqrt(MSE_camp);
MNE_camp = mean((lambda1CAMP(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NME_camp = ME_camp/mean(Obs);
MAE_camp = mean(abs(lambda1CAMP-Obs));
MNAE_camp = mean(abs(lambda1CAMP(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NMAE_camp = MAE_camp/mean(Obs);
ME_camp2 = mean((Obs-lambda1CAMP)./sqrt(lambda2CAMP));
SE2_camp2 = var((Obs-lambda1CAMP)./sqrt(lambda2CAMP));
RMSS_camp = sqrt(mean( ((Obs-lambda1CAMP)./sqrt(lambda2CAMP)).^2 ));
MR_camp = mean(sqrt(lambda2CAMP));
MSNS_camp = mean( ( ((Obs-lambda1CAMP)./sqrt(lambda2CAMP)).^2 - 1 ).^2 );

% Constant
ME_con = mean(lambda1Constant-Obs);
SE_con = std(lambda1Constant-Obs);
ME2_con = ME_con.^2;
SE2_con = var(lambda1Constant-Obs);
MSE_con = mean((lambda1Constant-Obs).^2);
R_con = corr(lambda1Constant,Obs);
R2_con = R_con.^2;
RMSE_con = sqrt(MSE_con);
MNE_con = mean((lambda1Constant(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NME_con = ME_con/mean(Obs);
MAE_con = mean(abs(lambda1Constant-Obs));
MNAE_con = mean(abs(lambda1Constant(Obs>0)-Obs(Obs>0))./Obs(Obs>0));
NMAE_con = MAE_con/mean(Obs);
ME_con2 = mean((Obs-lambda1Constant)./sqrt(lambda2Constant));
SE2_con2 = var((Obs-lambda1Constant)./sqrt(lambda2Constant));
RMSS_con = sqrt(mean( ((Obs-lambda1Constant)./sqrt(lambda2Constant)).^2 ));
MR_con = mean(sqrt(lambda2Constant));
MSNS_con = mean( ( ((Obs-lambda1Constant)./sqrt(lambda2Constant)).^2 - 1 ).^2 );

% season
for i = 1:4
    
    % pick season
    season_str = {'Winter','Spring','Summer','Fall'};
    if i == 1, idx = IsWinter == 1; end
    if i == 2, idx = IsSpring == 1; end
    if i == 3, idx = IsSummer == 1; end
    if i == 4, idx = IsFall == 1; end
    
    % CMAQ
    ME_cmaq_season(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_season(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_season(i,1) = ME_cmaq_season(i,1).^2;
    SE2_cmaq_season(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_season(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_season(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_season(i,1) = R_cmaq_season(i,1).^2;
    RMSE_cmaq_season(i,1) = sqrt(MSE_cmaq_season(i,1));
    MNE_cmaq_season(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_season(i,1) = ME_cmaq_season(i,1)/mean(Obs(idx));
    MAE_cmaq_season(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_season(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_season(i,1) = MAE_cmaq_season(i,1)/mean(Obs(idx));
    ME_cmaq2_season(i,1) = NaN;
    SE2_cmaq2_season(i,1) = NaN;
    RMSS_cmaq_season(i,1) = NaN;
    MR_cmaq_season(i,1) = NaN;

    % RAMP
    ME_ramp_season(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_season(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_season(i,1) = ME_ramp_season(i,1).^2;
    SE2_ramp_season(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_season(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_season(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_season(i,1) = R_ramp_season(i,1).^2;
    RMSE_ramp_season(i,1) = sqrt(MSE_ramp_season(i,1));
    MNE_ramp_season(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_season(i,1) = ME_ramp_season(i,1)/mean(Obs(idx));
    MAE_ramp_season(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_season(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_season(i,1) = MAE_ramp_season(i,1)/mean(Obs(idx));
    ME_ramp2_season(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_season(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_season(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_season(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_season(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_season(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_season(i,1) = ME_camp_season(i,1).^2;
    SE2_camp_season(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_season(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_season(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_season(i,1) = R_camp_season(i,1).^2;
    RMSE_camp_season(i,1) = sqrt(MSE_camp_season(i,1));
    MNE_camp_season(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_season(i,1) = ME_camp_season(i,1)/mean(Obs(idx));
    MAE_camp_season(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_season(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_season(i,1) = MAE_camp_season(i,1)/mean(Obs(idx));
    ME_camp2_season(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_season(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_season(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_season(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % seasonal CAMP 
    ME_campS_season(i,1) = mean(lambda1CAMPS(idx)-Obs(idx));
    SE_campS_season(i,1) = std(lambda1CAMPS(idx)-Obs(idx));
    ME2_campS_season(i,1) = ME_campS_season(i,1).^2;
    SE2_campS_season(i,1) = var(lambda1CAMPS(idx)-Obs(idx));
    MSE_campS_season(i,1) = mean((lambda1CAMPS(idx)-Obs(idx)).^2);
    R_campS_season(i,1) = corr(lambda1CAMPS(idx),Obs(idx));
    R2_campS_season(i,1) = R_campS_season(i,1).^2;
    RMSE_campS_season(i,1) = sqrt(MSE_campS_season(i,1));
    MNE_campS_season(i,1) = mean((lambda1CAMPS(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_campS_season(i,1) = ME_campS_season(i,1)/mean(Obs(idx));
    MAE_campS_season(i,1) = mean(abs(lambda1CAMPS(idx)-Obs(idx)));
    MNAE_campS_season(i,1) = mean(abs(lambda1CAMPS(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_campS_season(i,1) = MAE_campS_season(i,1)/mean(Obs(idx));
    ME_campS2_season(i,1) = mean((Obs(idx)-lambda1CAMPS(idx))./sqrt(lambda2CAMPS(idx)));
    SE2_campS2_season(i,1) = var((Obs(idx)-lambda1CAMPS(idx))./sqrt(lambda2CAMPS(idx)));
    RMSS_campS_season(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMPS(idx))./sqrt(lambda2CAMPS(idx))).^2 ));
    MR_campS_season(i,1) = mean(sqrt(lambda2CAMPS(idx)));
    
    % Constant
    ME_con_season(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_season(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_season(i,1) = ME_con_season(i,1).^2;
    SE2_con_season(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_season(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_season(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_season(i,1) = R_con_season(i,1).^2;
    RMSE_con_season(i,1) = sqrt(MSE_con_season(i,1));
    MNE_con_season(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_season(i,1) = ME_con_season(i,1)/mean(Obs(idx));
    MAE_con_season(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_season(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_season(i,1) = MAE_con_season(i,1)/mean(Obs(idx));
    ME_con2_season(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_season(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_season(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_season(i,1) = mean(sqrt(lambda2Constant(idx)));
    
    % seasonal Constant 
    ME_conS_season(i,1) = mean(lambda1ConstantS(idx)-Obs(idx));
    SE_conS_season(i,1) = std(lambda1ConstantS(idx)-Obs(idx));
    ME2_conS_season(i,1) = ME_conS_season(i,1).^2;
    SE2_conS_season(i,1) = var(lambda1ConstantS(idx)-Obs(idx));
    MSE_conS_season(i,1) = mean((lambda1ConstantS(idx)-Obs(idx)).^2);
    R_conS_season(i,1) = corr(lambda1ConstantS(idx),Obs(idx));
    R2_conS_season(i,1) = R_conS_season(i,1).^2;
    RMSE_conS_season(i,1) = sqrt(MSE_conS_season(i,1));
    MNE_conS_season(i,1) = mean((lambda1ConstantS(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_conS_season(i,1) = ME_conS_season(i,1)/mean(Obs(idx));
    MAE_conS_season(i,1) = mean(abs(lambda1ConstantS(idx)-Obs(idx)));
    MNAE_conS_season(i,1) = mean(abs(lambda1ConstantS(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_conS_season(i,1) = MAE_conS_season(i,1)/mean(Obs(idx));
    ME_conS2_season(i,1) = mean((Obs(idx)-lambda1ConstantS(idx))./sqrt(lambda2ConstantS(idx)));
    SE2_conS2_season(i,1) = var((Obs(idx)-lambda1ConstantS(idx))./sqrt(lambda2ConstantS(idx)));
    RMSS_conS_season(i,1) = sqrt(mean( ((Obs(idx)-lambda1ConstantS(idx))./sqrt(lambda2ConstantS(idx))).^2 ));
    MR_conS_season(i,1) = mean(sqrt(lambda2ConstantS(idx)));
    
end

% network
for i = 1:3
    
    % pick network
    network_str = {'IMPROVE','STN','Other'};
    if i == 1, idx = IsIMPROVE == 1; end
    if i == 2, idx = IsSTN == 1; end
    if i == 3, idx = isnan(IsIMPROVE) & isnan(IsSTN); end
    
    % CMAQ
    ME_cmaq_network(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_network(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_network(i,1) = ME_cmaq_network(i,1).^2;
    SE2_cmaq_network(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_network(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_network(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_network(i,1) = R_cmaq_network(i,1).^2;
    RMSE_cmaq_network(i,1) = sqrt(MSE_cmaq_network(i,1));
    MNE_cmaq_network(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_network(i,1) = ME_cmaq_network(i,1)/mean(Obs(idx));
    MAE_cmaq_network(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_network(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_network(i,1) = MAE_cmaq_network(i,1)/mean(Obs(idx));
    ME_cmaq2_network(i,1) = NaN;
    SE2_cmaq2_network(i,1) = NaN;
    RMSS_cmaq_network(i,1) = NaN;
    MR_cmaq_network(i,1) = NaN;

    % RAMP
    ME_ramp_network(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_network(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_network(i,1) = ME_ramp_network(i,1).^2;
    SE2_ramp_network(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_network(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_network(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_network(i,1) = R_ramp_network(i,1).^2;
    RMSE_ramp_network(i,1) = sqrt(MSE_ramp_network(i,1));
    MNE_ramp_network(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_network(i,1) = ME_ramp_network(i,1)/mean(Obs(idx));
    MAE_ramp_network(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_network(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_network(i,1) = MAE_ramp_network(i,1)/mean(Obs(idx));
    ME_ramp2_network(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_network(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_network(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_network(i,1) = sqrt(mean(lambda2(idx)));
    
    % CAMP
    ME_camp_network(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_network(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_network(i,1) = ME_camp_network(i,1).^2;
    SE2_camp_network(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_network(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_network(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_network(i,1) = R_camp_network(i,1).^2;
    RMSE_camp_network(i,1) = sqrt(MSE_camp_network(i,1));
    MNE_camp_network(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_network(i,1) = ME_camp_network(i,1)/mean(Obs(idx));
    MAE_camp_network(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_network(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_network(i,1) = MAE_camp_network(i,1)/mean(Obs(idx));
    ME_camp2_network(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_network(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_network(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_network(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_network(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_network(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_network(i,1) = ME_con_network(i,1).^2;
    SE2_con_network(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_network(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_network(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_network(i,1) = R_con_network(i,1).^2;
    RMSE_con_network(i,1) = sqrt(MSE_con_network(i,1));
    MNE_con_network(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_network(i,1) = ME_con_network(i,1)/mean(Obs(idx));
    MAE_con_network(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_network(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_network(i,1) = MAE_con_network(i,1)/mean(Obs(idx));
    ME_con2_network(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_network(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_network(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_network(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% FRM/TEOM
for i = 1:2
    
    % pick frm/teom
    frm_str = {'FRM','TEOM'};
    if i == 1, idx = IsFRM == 1; end
    if i == 2, idx = IsTEOM == 1; end

    % CMAQ
    ME_cmaq_frm(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_frm(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_frm(i,1) = ME_cmaq_frm(i,1).^2;
    SE2_cmaq_frm(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_frm(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_frm(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_frm(i,1) = R_cmaq_frm(i,1).^2;
    RMSE_cmaq_frm(i,1) = sqrt(MSE_cmaq_frm(i,1));
    MNE_cmaq_frm(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_frm(i,1) = ME_cmaq_frm(i,1)/mean(Obs(idx));
    MAE_cmaq_frm(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_frm(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_frm(i,1) = MAE_cmaq_frm(i,1)/mean(Obs(idx));
    ME_cmaq2_frm(i,1) = NaN;
    SE2_cmaq2_frm(i,1) = NaN;
    RMSS_cmaq_frm(i,1) = NaN;
    MR_cmaq_frm(i,1) = NaN;

    % RAMP
    ME_ramp_frm(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_frm(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_frm(i,1) = ME_ramp_frm(i,1).^2;
    SE2_ramp_frm(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_frm(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_frm(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_frm(i,1) = R_ramp_frm(i,1).^2;
    RMSE_ramp_frm(i,1) = sqrt(MSE_ramp_frm(i,1));
    MNE_ramp_frm(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_frm(i,1) = ME_ramp_frm(i,1)/mean(Obs(idx));
    MAE_ramp_frm(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_frm(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_frm(i,1) = MAE_ramp_frm(i,1)/mean(Obs(idx));
    ME_ramp2_frm(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_frm(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_frm(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_frm(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_frm(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_frm(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_frm(i,1) = ME_camp_frm(i,1).^2;
    SE2_camp_frm(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_frm(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_frm(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_frm(i,1) = R_camp_frm(i,1).^2;
    RMSE_camp_frm(i,1) = sqrt(MSE_camp_frm(i,1));
    MNE_camp_frm(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_frm(i,1) = ME_camp_frm(i,1)/mean(Obs(idx));
    MAE_camp_frm(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_frm(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_frm(i,1) = MAE_camp_frm(i,1)/mean(Obs(idx));
    ME_camp2_frm(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_frm(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_frm(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_frm(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_frm(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_frm(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_frm(i,1) = ME_con_frm(i,1).^2;
    SE2_con_frm(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_frm(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_frm(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_frm(i,1) = R_con_frm(i,1).^2;
    RMSE_con_frm(i,1) = sqrt(MSE_con_frm(i,1));
    MNE_con_frm(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_frm(i,1) = ME_con_frm(i,1)/mean(Obs(idx));
    MAE_con_frm(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_frm(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_frm(i,1) = MAE_con_frm(i,1)/mean(Obs(idx));
    ME_con2_frm(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_frm(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_frm(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_frm(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% east/west
for i = 1:2
    
    % pick west/east
    west_str = {'West','East'};
    if i == 1, idx = IsWest == 1; end
    if i == 2, idx = IsWest == 0; end

    % CMAQ
    ME_cmaq_west(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_west(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_west(i,1) = ME_cmaq_west(i,1).^2;
    SE2_cmaq_west(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_west(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_west(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_west(i,1) = R_cmaq_west(i,1).^2;
    RMSE_cmaq_west(i,1) = sqrt(MSE_cmaq_west(i,1));
    MNE_cmaq_west(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_west(i,1) = ME_cmaq_west(i,1)/mean(Obs(idx));
    MAE_cmaq_west(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_west(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_west(i,1) = MAE_cmaq_west(i,1)/mean(Obs(idx));
    ME_cmaq2_west(i,1) = NaN;
    SE2_cmaq2_west(i,1) = NaN;
    RMSS_cmaq_west(i,1) = NaN;
    MR_cmaq_west(i,1) = NaN;

    % RAMP
    ME_ramp_west(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_west(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_west(i,1) = ME_ramp_west(i,1).^2;
    SE2_ramp_west(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_west(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_west(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_west(i,1) = R_ramp_west(i,1).^2;
    RMSE_ramp_west(i,1) = sqrt(MSE_ramp_west(i,1));
    MNE_ramp_west(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_west(i,1) = ME_ramp_west(i,1)/mean(Obs(idx));
    MAE_ramp_west(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_west(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_west(i,1) = MAE_ramp_west(i,1)/mean(Obs(idx));
    ME_ramp2_west(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_west(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_west(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_west(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_west(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_west(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_west(i,1) = ME_camp_west(i,1).^2;
    SE2_camp_west(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_west(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_west(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_west(i,1) = R_camp_west(i,1).^2;
    RMSE_camp_west(i,1) = sqrt(MSE_camp_west(i,1));
    MNE_camp_west(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_west(i,1) = ME_camp_west(i,1)/mean(Obs(idx));
    MAE_camp_west(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_west(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_west(i,1) = MAE_camp_west(i,1)/mean(Obs(idx));
    ME_camp2_west(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_west(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_west(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_west(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_west(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_west(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_west(i,1) = ME_con_west(i,1).^2;
    SE2_con_west(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_west(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_west(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_west(i,1) = R_con_west(i,1).^2;
    RMSE_con_west(i,1) = sqrt(MSE_con_west(i,1));
    MNE_con_west(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_west(i,1) = ME_con_west(i,1)/mean(Obs(idx));
    MAE_con_west(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_west(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_west(i,1) = MAE_con_west(i,1)/mean(Obs(idx));
    ME_con2_west(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_west(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_west(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_west(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% urban/rural/suburban
for i = 1:3
    
    % pick urban/rural/suburban
    urban_str = {'Urban','Rurual','Suburban'};
    if i == 1, idx = IsUrban == 1; end
    if i == 2, idx = IsRural == 1; end
    if i == 3, idx = IsSuburban == 1; end

    % CMAQ
    ME_cmaq_urban(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_urban(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_urban(i,1) = ME_cmaq_urban(i,1).^2;
    SE2_cmaq_urban(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_urban(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_urban(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_urban(i,1) = R_cmaq_urban(i,1).^2;
    RMSE_cmaq_urban(i,1) = sqrt(MSE_cmaq_urban(i,1));
    MNE_cmaq_urban(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_urban(i,1) = ME_cmaq_urban(i,1)/mean(Obs(idx));
    MAE_cmaq_urban(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_urban(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_urban(i,1) = MAE_cmaq_urban(i,1)/mean(Obs(idx));
    ME_cmaq2_urban(i,1) = NaN;
    SE2_cmaq2_urban(i,1) = NaN;
    RMSS_cmaq_urban(i,1) = NaN;
    MR_cmaq_urban(i,1) = NaN;

    % RAMP
    ME_ramp_urban(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_urban(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_urban(i,1) = ME_ramp_urban(i,1).^2;
    SE2_ramp_urban(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_urban(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_urban(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_urban(i,1) = R_ramp_urban(i,1).^2;
    RMSE_ramp_urban(i,1) = sqrt(MSE_ramp_urban(i,1));
    MNE_ramp_urban(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_urban(i,1) = ME_ramp_urban(i,1)/mean(Obs(idx));
    MAE_ramp_urban(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_urban(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_urban(i,1) = MAE_ramp_urban(i,1)/mean(Obs(idx));
    ME_ramp2_urban(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_urban(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));  
    RMSS_ramp_urban(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_urban(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_urban(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_urban(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_urban(i,1) = ME_camp_urban(i,1).^2;
    SE2_camp_urban(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_urban(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_urban(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_urban(i,1) = R_camp_urban(i,1).^2;
    RMSE_camp_urban(i,1) = sqrt(MSE_camp_urban(i,1));
    MNE_camp_urban(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_urban(i,1) = ME_camp_urban(i,1)/mean(Obs(idx));
    MAE_camp_urban(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_urban(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_urban(i,1) = MAE_camp_urban(i,1)/mean(Obs(idx));
    ME_camp2_urban(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_urban(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_urban(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_urban(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_urban(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_urban(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_urban(i,1) = ME_con_urban(i,1).^2;
    SE2_con_urban(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_urban(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_urban(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_urban(i,1) = R_con_urban(i,1).^2;
    RMSE_con_urban(i,1) = sqrt(MSE_con_urban(i,1));
    MNE_con_urban(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_urban(i,1) = ME_con_urban(i,1)/mean(Obs(idx));
    MAE_con_urban(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_urban(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_urban(i,1) = MAE_con_urban(i,1)/mean(Obs(idx));
    ME_con2_urban(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_urban(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_urban(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_urban(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% region
for i = 1:length(inregions)
    
    % region
    idx = inregions{i} == 1; 

    if sum(idx) > 0 
        
        % CMAQ
        ME_cmaq_region(i,1) = mean(Mod(idx)-Obs(idx));
        SE_cmaq_region(i,1) = std(Mod(idx)-Obs(idx));
        ME2_cmaq_region(i,1) = ME_cmaq_region(i,1).^2;
        SE2_cmaq_region(i,1) = var(Mod(idx)-Obs(idx));
        MSE_cmaq_region(i,1) = mean((Mod(idx)-Obs(idx)).^2);
        R_cmaq_region(i,1) = corr(Mod(idx),Obs(idx));
        R2_cmaq_region(i,1) = R_cmaq_region(i,1).^2;
        RMSE_cmaq_region(i,1) = sqrt(MSE_cmaq_region(i,1));
        MNE_cmaq_region(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NME_cmaq_region(i,1) = ME_cmaq_region(i,1)/mean(Obs(idx));
        MAE_cmaq_region(i,1) = mean(abs(Mod(idx)-Obs(idx)));
        MNAE_cmaq_region(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NMAE_cmaq_region(i,1) = MAE_cmaq_region(i,1)/mean(Obs(idx));
        ME_cmaq2_region(i,1) = NaN;
        SE2_cmaq2_region(i,1) = NaN;
        RMSS_cmaq_region(i,1) = NaN;
        MR_cmaq_region(i,1) = NaN;

        % RAMP
        ME_ramp_region(i,1) = mean(lambda1(idx)-Obs(idx));
        SE_ramp_region(i,1) = std(lambda1(idx)-Obs(idx));
        ME2_ramp_region(i,1) = ME_ramp_region(i,1).^2;
        SE2_ramp_region(i,1) = var(lambda1(idx)-Obs(idx));
        MSE_ramp_region(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
        R_ramp_region(i,1) = corr(lambda1(idx),Obs(idx));
        R2_ramp_region(i,1) = R_ramp_region(i,1).^2;
        RMSE_ramp_region(i,1) = sqrt(MSE_ramp_region(i,1));
        MNE_ramp_region(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NME_ramp_region(i,1) = ME_ramp_region(i,1)/mean(Obs(idx));
        MAE_ramp_region(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
        MNAE_ramp_region(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NMAE_ramp_region(i,1) = MAE_ramp_region(i,1)/mean(Obs(idx));
        ME_ramp2_region(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
        SE2_ramp2_region(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
        RMSS_ramp_region(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
        MR_ramp_region(i,1) = mean(sqrt(lambda2(idx)));
        
        % CAMP
        ME_camp_region(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
        SE_camp_region(i,1) = std(lambda1CAMP(idx)-Obs(idx));
        ME2_camp_region(i,1) = ME_camp_region(i,1).^2;
        SE2_camp_region(i,1) = var(lambda1CAMP(idx)-Obs(idx));
        MSE_camp_region(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
        R_camp_region(i,1) = corr(lambda1CAMP(idx),Obs(idx));
        R2_camp_region(i,1) = R_camp_region(i,1).^2;
        RMSE_camp_region(i,1) = sqrt(MSE_camp_region(i,1));
        MNE_camp_region(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NME_camp_region(i,1) = ME_camp_region(i,1)/mean(Obs(idx));
        MAE_camp_region(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
        MNAE_camp_region(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NMAE_camp_region(i,1) = MAE_camp_region(i,1)/mean(Obs(idx));
        ME_camp2_region(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
        SE2_camp2_region(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
        RMSS_camp_region(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
        MR_camp_region(i,1) = mean(sqrt(lambda2CAMP(idx)));
        
        % Constant
        ME_con_region(i,1) = mean(lambda1Constant(idx)-Obs(idx));
        SE_con_region(i,1) = std(lambda1Constant(idx)-Obs(idx));
        ME2_con_region(i,1) = ME_con_region(i,1).^2;
        SE2_con_region(i,1) = var(lambda1Constant(idx)-Obs(idx));
        MSE_con_region(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
        R_con_region(i,1) = corr(lambda1Constant(idx),Obs(idx));
        R2_con_region(i,1) = R_con_region(i,1).^2;
        RMSE_con_region(i,1) = sqrt(MSE_con_region(i,1));
        MNE_con_region(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NME_con_region(i,1) = ME_con_region(i,1)/mean(Obs(idx));
        MAE_con_region(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
        MNAE_con_region(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
        NMAE_con_region(i,1) = MAE_con_region(i,1)/mean(Obs(idx));
        ME_con2_region(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
        SE2_con2_region(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
        RMSS_con_region(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
        MR_con_region(i,1) = mean(sqrt(lambda2Constant(idx)));

    else
        
        % CMAQ
        ME_cmaq_region(i,1) = NaN;
        SE_cmaq_region(i,1) = NaN;
        ME2_cmaq_region(i,1) = NaN;
        SE2_cmaq_region(i,1) = NaN;
        MSE_cmaq_region(i,1) = NaN;
        R_cmaq_region(i,1) = NaN;
        R2_cmaq_region(i,1) = NaN;
        RMSE_cmaq_region(i,1) = NaN;
        MNE_cmaq_region(i,1) = NaN;
        NME_cmaq_region(i,1) = NaN;
        MAE_cmaq_region(i,1) = NaN;
        MNAE_cmaq_region(i,1) = NaN;
        NMAE_cmaq_region(i,1) = NaN;
        ME_cmaq2_region(i,1) = NaN;
        SE2_cmaq2_region(i,1) = NaN;
        RMSS_cmaq_region(i,1) = NaN;
        MR_cmaq_region(i,1) = NaN;

        % RAMP
        ME_ramp_region(i,1) = NaN;
        SE_ramp_region(i,1) = NaN;
        ME2_ramp_region(i,1) = NaN;
        SE2_ramp_region(i,1) = NaN;
        MSE_ramp_region(i,1) = NaN;
        R_ramp_region(i,1) = NaN;
        R2_ramp_region(i,1) = NaN;
        RMSE_ramp_region(i,1) = NaN;
        MNE_ramp_region(i,1) = NaN;
        NME_ramp_region(i,1) = NaN;
        MAE_ramp_region(i,1) = NaN;
        MNAE_ramp_region(i,1) = NaN;
        NMAE_ramp_region(i,1) = NaN;
        ME_ramp2_region(i,1) = NaN;
        SE2_ramp2_region(i,1) = NaN;
        RMSS_ramp_region(i,1) = NaN;
        MR_ramp_region(i,1) = NaN;
        
        % CAMP
        ME_camp_region(i,1) = NaN;
        SE_camp_region(i,1) = NaN;
        ME2_camp_region(i,1) = NaN;
        SE2_camp_region(i,1) = NaN;
        MSE_camp_region(i,1) = NaN;
        R_camp_region(i,1) = NaN;
        R2_camp_region(i,1) = NaN;
        RMSE_camp_region(i,1) = NaN;
        MNE_camp_region(i,1) = NaN;
        NME_camp_region(i,1) = NaN;
        MAE_camp_region(i,1) = NaN;
        MNAE_camp_region(i,1) = NaN;
        NMAE_camp_region(i,1) = NaN;
        ME_camp2_region(i,1) = NaN;
        SE2_camp2_region(i,1) = NaN;
        RMSS_camp_region(i,1) = NaN;
        MR_camp_region(i,1) = NaN;
        
        % Constant
        ME_con_region(i,1) = NaN;
        SE_con_region(i,1) = NaN;
        ME2_con_region(i,1) = NaN;
        SE2_con_region(i,1) = NaN;
        MSE_con_region(i,1) = NaN;
        R_con_region(i,1) = NaN;
        R2_con_region(i,1) = NaN;
        RMSE_con_region(i,1) = NaN;
        MNE_con_region(i,1) = NaN;
        NME_con_region(i,1) = NaN;
        MAE_con_region(i,1) = NaN;
        MNAE_con_region(i,1) = NaN;
        NMAE_con_region(i,1) = NaN;
        ME_con2_region(i,1) = NaN;
        SE2_con2_region(i,1) = NaN;
        RMSS_con_region(i,1) = NaN;
        MR_con_region(i,1) = NaN;
        
    end
    
end

% regional CAMP
region6_str = {'northeast';'southeast';'uppermidwest';'lowermidwest';'rockymountains';'pacificcoast'};
for i = 1:6
    
    % region
    idx = inregion6 == i; 
        
    % CMAQ
    ME_cmaq_region6(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_region6(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_region6(i,1) = ME_cmaq_region6(i,1).^2;
    SE2_cmaq_region6(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_region6(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_region6(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_region6(i,1) = R_cmaq_region6(i,1).^2;
    RMSE_cmaq_region6(i,1) = sqrt(MSE_cmaq_region6(i,1));
    MNE_cmaq_region6(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_region6(i,1) = ME_cmaq_region6(i,1)/mean(Obs(idx));
    MAE_cmaq_region6(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_region6(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_region6(i,1) = MAE_cmaq_region6(i,1)/mean(Obs(idx));
    ME_cmaq2_region6(i,1) = NaN;
    SE2_cmaq2_region6(i,1) = NaN;
    RMSS_cmaq_region6(i,1) = NaN;
    MR_cmaq_region6(i,1) = NaN;

    % RAMP
    ME_ramp_region6(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_region6(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_region6(i,1) = ME_ramp_region6(i,1).^2;
    SE2_ramp_region6(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_region6(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_region6(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_region6(i,1) = R_ramp_region6(i,1).^2;
    RMSE_ramp_region6(i,1) = sqrt(MSE_ramp_region6(i,1));
    MNE_ramp_region6(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_region6(i,1) = ME_ramp_region6(i,1)/mean(Obs(idx));
    MAE_ramp_region6(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_region6(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_region6(i,1) = MAE_ramp_region6(i,1)/mean(Obs(idx));
    ME_ramp2_region6(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_region6(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_region6(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_region6(i,1) = mean(sqrt(lambda2(idx)));

    % CAMP
    ME_camp_region6(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_region6(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_region6(i,1) = ME_camp_region6(i,1).^2;
    SE2_camp_region6(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_region6(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_region6(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_region6(i,1) = R_camp_region6(i,1).^2;
    RMSE_camp_region6(i,1) = sqrt(MSE_camp_region6(i,1));
    MNE_camp_region6(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_region6(i,1) = ME_camp_region6(i,1)/mean(Obs(idx));
    MAE_camp_region6(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_region6(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_region6(i,1) = MAE_camp_region6(i,1)/mean(Obs(idx));
    ME_camp2_region6(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_region6(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_region6(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_region6(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % regional CAMP
    ME_camp6_region6(i,1) = mean(lambda1CAMP6(idx)-Obs(idx));
    SE_camp6_region6(i,1) = std(lambda1CAMP6(idx)-Obs(idx));
    ME2_camp6_region6(i,1) = ME_camp6_region6(i,1).^2;
    SE2_camp6_region6(i,1) = var(lambda1CAMP6(idx)-Obs(idx));
    MSE_camp6_region6(i,1) = mean((lambda1CAMP6(idx)-Obs(idx)).^2);
    R_camp6_region6(i,1) = corr(lambda1CAMP6(idx),Obs(idx));
    R2_camp6_region6(i,1) = R_camp6_region6(i,1).^2;
    RMSE_camp6_region6(i,1) = sqrt(MSE_camp6_region6(i,1));
    MNE_camp6_region6(i,1) = mean((lambda1CAMP6(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp6_region6(i,1) = ME_camp6_region6(i,1)/mean(Obs(idx));
    MAE_camp6_region6(i,1) = mean(abs(lambda1CAMP6(idx)-Obs(idx)));
    MNAE_camp6_region6(i,1) = mean(abs(lambda1CAMP6(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp6_region6(i,1) = MAE_camp6_region6(i,1)/mean(Obs(idx));
    ME_camp62_region6(i,1) = mean((Obs(idx)-lambda1CAMP6(idx))./sqrt(lambda2CAMP6(idx)));
    SE2_camp62_region6(i,1) = var((Obs(idx)-lambda1CAMP6(idx))./sqrt(lambda2CAMP6(idx)));
    RMSS_camp6_region6(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP6(idx))./sqrt(lambda2CAMP6(idx))).^2 ));
    MR_camp6_region6(i,1) = mean(sqrt(lambda2CAMP6(idx)));
    
    % Constant
    ME_con_region6(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_region6(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_region6(i,1) = ME_con_region6(i,1).^2;
    SE2_con_region6(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_region6(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_region6(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_region6(i,1) = R_con_region6(i,1).^2;
    RMSE_con_region6(i,1) = sqrt(MSE_con_region6(i,1));
    MNE_con_region6(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_region6(i,1) = ME_con_region6(i,1)/mean(Obs(idx));
    MAE_con_region6(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_region6(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_region6(i,1) = MAE_con_region6(i,1)/mean(Obs(idx));
    ME_con2_region6(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_region6(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_region6(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_region6(i,1) = mean(sqrt(lambda2Constant(idx)));
    
    % regional Constant
    ME_con6_region6(i,1) = mean(lambda1Constant6(idx)-Obs(idx));
    SE_con6_region6(i,1) = std(lambda1Constant6(idx)-Obs(idx));
    ME2_con6_region6(i,1) = ME_con6_region6(i,1).^2;
    SE2_con6_region6(i,1) = var(lambda1Constant6(idx)-Obs(idx));
    MSE_con6_region6(i,1) = mean((lambda1Constant6(idx)-Obs(idx)).^2);
    R_con6_region6(i,1) = corr(lambda1Constant6(idx),Obs(idx));
    R2_con6_region6(i,1) = R_con6_region6(i,1).^2;
    RMSE_con6_region6(i,1) = sqrt(MSE_con6_region6(i,1));
    MNE_con6_region6(i,1) = mean((lambda1Constant6(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con6_region6(i,1) = ME_con6_region6(i,1)/mean(Obs(idx));
    MAE_con6_region6(i,1) = mean(abs(lambda1Constant6(idx)-Obs(idx)));
    MNAE_con6_region6(i,1) = mean(abs(lambda1Constant6(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con6_region6(i,1) = MAE_con6_region6(i,1)/mean(Obs(idx));
    ME_con62_region6(i,1) = mean((Obs(idx)-lambda1Constant6(idx))./sqrt(lambda2Constant6(idx)));
    SE2_con62_region6(i,1) = var((Obs(idx)-lambda1Constant6(idx))./sqrt(lambda2Constant6(idx)));
    RMSS_con6_region6(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant6(idx))./sqrt(lambda2Constant6(idx))).^2 ));
    MR_con6_region6(i,1) = mean(sqrt(lambda2Constant6(idx)));

end

% distance to closest monitor
unilength = unique(Dist2NMon);
lendist = prctile(unilength,0:10:100);
for i = 1:10
    
    % distance to closest monitor
    idx = Dist2NMon >= lendist(i) & Dist2NMon < lendist(i+1);

    % CMAQ
    ME_cmaq_dist(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_dist(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_dist(i,1) = ME_cmaq_dist(i,1).^2;
    SE2_cmaq_dist(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_dist(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_dist(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_dist(i,1) = R_cmaq_dist(i,1).^2;
    RMSE_cmaq_dist(i,1) = sqrt(MSE_cmaq_dist(i,1));
    MNE_cmaq_dist(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_dist(i,1) = ME_cmaq_dist(i,1)/mean(Obs(idx));
    MAE_cmaq_dist(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_dist(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_dist(i,1) = MAE_cmaq_dist(i,1)/mean(Obs(idx));
    ME_cmaq2_dist(i,1) = NaN;
    SE2_cmaq2_dist(i,1) = NaN;
    RMSS_cmaq_dist(i,1) = NaN;
    MR_cmaq_dist(i,1) = NaN;

    % RAMP
    ME_ramp_dist(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_dist(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_dist(i,1) = ME_ramp_dist(i,1).^2;
    SE2_ramp_dist(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_dist(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_dist(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_dist(i,1) = R_ramp_dist(i,1).^2;
    RMSE_ramp_dist(i,1) = sqrt(MSE_ramp_dist(i,1));
    MNE_ramp_dist(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_dist(i,1) = ME_ramp_dist(i,1)/mean(Obs(idx));
    MAE_ramp_dist(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_dist(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_dist(i,1) = MAE_ramp_dist(i,1)/mean(Obs(idx));
    ME_ramp2_dist(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_dist(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_dist(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_dist(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_dist(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_dist(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_dist(i,1) = ME_camp_dist(i,1).^2;
    SE2_camp_dist(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_dist(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_dist(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_dist(i,1) = R_camp_dist(i,1).^2;
    RMSE_camp_dist(i,1) = sqrt(MSE_camp_dist(i,1));
    MNE_camp_dist(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_dist(i,1) = ME_camp_dist(i,1)/mean(Obs(idx));
    MAE_camp_dist(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_dist(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_dist(i,1) = MAE_camp_dist(i,1)/mean(Obs(idx));
    ME_camp2_dist(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_dist(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_dist(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_dist(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_dist(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_dist(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_dist(i,1) = ME_con_dist(i,1).^2;
    SE2_con_dist(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_dist(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_dist(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_dist(i,1) = R_con_dist(i,1).^2;
    RMSE_con_dist(i,1) = sqrt(MSE_con_dist(i,1));
    MNE_con_dist(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_dist(i,1) = ME_con_dist(i,1)/mean(Obs(idx));
    MAE_con_dist(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_dist(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_dist(i,1) = MAE_con_dist(i,1)/mean(Obs(idx));
    ME_con2_dist(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_dist(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_dist(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_dist(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% observed value
lenobs = prctile(Obs,0:10:100);
for i = 1:10
    
    % percentile of obs
    idx = Obs >= lenobs(i) & Obs < lenobs(i+1);

    % CMAQ
    ME_cmaq_obs(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_obs(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_obs(i,1) = ME_cmaq_obs(i,1).^2;
    SE2_cmaq_obs(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_obs(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_obs(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_obs(i,1) = R_cmaq_obs(i,1).^2;
    RMSE_cmaq_obs(i,1) = sqrt(MSE_cmaq_obs(i,1));
    MNE_cmaq_obs(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_obs(i,1) = ME_cmaq_obs(i,1)/mean(Obs(idx));
    MAE_cmaq_obs(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_obs(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_obs(i,1) = MAE_cmaq_obs(i,1)/mean(Obs(idx));
    ME_cmaq2_obs(i,1) = NaN;
    SE2_cmaq2_obs(i,1) = NaN;
    RMSS_cmaq_obs(i,1) = NaN;
    MR_cmaq_obs(i,1) = NaN;

    % RAMP
    ME_ramp_obs(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_obs(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_obs(i,1) = ME_ramp_obs(i,1).^2;
    SE2_ramp_obs(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_obs(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_obs(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_obs(i,1) = R_ramp_obs(i,1).^2;
    RMSE_ramp_obs(i,1) = sqrt(MSE_ramp_obs(i,1));
    MNE_ramp_obs(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_obs(i,1) = ME_ramp_obs(i,1)/mean(Obs(idx));
    MAE_ramp_obs(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_obs(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_obs(i,1) = MAE_ramp_obs(i,1)/mean(Obs(idx));
    ME_ramp2_obs(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_obs(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_obs(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_obs(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_obs(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_obs(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_obs(i,1) = ME_camp_obs(i,1).^2;
    SE2_camp_obs(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_obs(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_obs(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_obs(i,1) = R_camp_obs(i,1).^2;
    RMSE_camp_obs(i,1) = sqrt(MSE_camp_obs(i,1));
    MNE_camp_obs(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_obs(i,1) = ME_camp_obs(i,1)/mean(Obs(idx));
    MAE_camp_obs(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_obs(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_obs(i,1) = MAE_camp_obs(i,1)/mean(Obs(idx));
    ME_camp2_obs(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_obs(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_obs(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_obs(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_obs(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_obs(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_obs(i,1) = ME_con_obs(i,1).^2;
    SE2_con_obs(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_obs(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_obs(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_obs(i,1) = R_con_obs(i,1).^2;
    RMSE_con_obs(i,1) = sqrt(MSE_con_obs(i,1));
    MNE_con_obs(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_obs(i,1) = ME_con_obs(i,1)/mean(Obs(idx));
    MAE_con_obs(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_obs(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_obs(i,1) = MAE_con_obs(i,1)/mean(Obs(idx));
    ME_con2_obs(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_obs(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_obs(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_obs(i,1) = mean(sqrt(lambda2Constant(idx)));
        
end

% modeled value
lenmod = prctile(Mod,0:10:100);
for i = 1:10
    
    % distance to closest monitor
    idx = Mod >= lenmod(i) & Mod < lenmod(i+1);

    % CMAQ
    ME_cmaq_mod(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_mod(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_mod(i,1) = ME_cmaq_mod(i,1).^2;
    SE2_cmaq_mod(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_mod(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_mod(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_mod(i,1) = R_cmaq_mod(i,1).^2;
    RMSE_cmaq_mod(i,1) = sqrt(MSE_cmaq_mod(i,1));
    MNE_cmaq_mod(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_mod(i,1) = ME_cmaq_mod(i,1)/mean(Obs(idx));
    MAE_cmaq_mod(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_mod(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_mod(i,1) = MAE_cmaq_mod(i,1)/mean(Obs(idx));
    ME_cmaq2_mod(i,1) = NaN;
    SE2_cmaq2_mod(i,1) = NaN;
    RMSS_cmaq_mod(i,1) = NaN;
    MR_cmaq_mod(i,1) = NaN;

    % RAMP
    ME_ramp_mod(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_mod(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_mod(i,1) = ME_ramp_mod(i,1).^2;
    SE2_ramp_mod(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_mod(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_mod(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_mod(i,1) = R_ramp_mod(i,1).^2;
    RMSE_ramp_mod(i,1) = sqrt(MSE_ramp_mod(i,1));
    MNE_ramp_mod(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_mod(i,1) = ME_ramp_mod(i,1)/mean(Obs(idx));
    MAE_ramp_mod(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_mod(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_mod(i,1) = MAE_ramp_mod(i,1)/mean(Obs(idx));
    ME_ramp2_mod(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_mod(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_mod(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_mod(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_mod(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_mod(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_mod(i,1) = ME_camp_mod(i,1).^2;
    SE2_camp_mod(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_mod(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_mod(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_mod(i,1) = R_camp_mod(i,1).^2;
    RMSE_camp_mod(i,1) = sqrt(MSE_camp_mod(i,1));
    MNE_camp_mod(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_mod(i,1) = ME_camp_mod(i,1)/mean(Obs(idx));
    MAE_camp_mod(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_mod(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_mod(i,1) = MAE_camp_mod(i,1)/mean(Obs(idx));
    ME_camp2_mod(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_mod(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_mod(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_mod(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_mod(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_mod(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_mod(i,1) = ME_con_mod(i,1).^2;
    SE2_con_mod(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_mod(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_mod(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_mod(i,1) = R_con_mod(i,1).^2;
    RMSE_con_mod(i,1) = sqrt(MSE_con_mod(i,1));
    MNE_con_mod(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_mod(i,1) = ME_con_mod(i,1)/mean(Obs(idx));
    MAE_con_mod(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_mod(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_mod(i,1) = MAE_con_mod(i,1)/mean(Obs(idx));
    ME_con2_mod(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_mod(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_mod(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_mod(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% lambda2
lenlambda2 = prctile(lambda2,0:10:100);
for i = 1:10
    
    % distance to closest monitor
    idx = lambda2 >= lenlambda2(i) & lambda2 < lenlambda2(i+1);

    % CMAQ
    ME_cmaq_lambda2(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_lambda2(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_lambda2(i,1) = ME_cmaq_lambda2(i,1).^2;
    SE2_cmaq_lambda2(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_lambda2(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_lambda2(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_lambda2(i,1) = R_cmaq_lambda2(i,1).^2;
    RMSE_cmaq_lambda2(i,1) = sqrt(MSE_cmaq_lambda2(i,1));
    MNE_cmaq_lambda2(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_lambda2(i,1) = ME_cmaq_lambda2(i,1)/mean(Obs(idx));
    MAE_cmaq_lambda2(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_lambda2(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_lambda2(i,1) = MAE_cmaq_lambda2(i,1)/mean(Obs(idx));
    ME_cmaq2_lambda2(i,1) = NaN;
    SE2_cmaq2_lambda2(i,1) = NaN;
    RMSS_cmaq_lambda2(i,1) = NaN;
    MR_cmaq_lambda2(i,1) = NaN;

    % RAMP
    ME_ramp_lambda2(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_lambda2(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_lambda2(i,1) = ME_ramp_lambda2(i,1).^2;
    SE2_ramp_lambda2(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_lambda2(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_lambda2(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_lambda2(i,1) = R_ramp_lambda2(i,1).^2;
    RMSE_ramp_lambda2(i,1) = sqrt(MSE_ramp_lambda2(i,1));
    MNE_ramp_lambda2(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_lambda2(i,1) = ME_ramp_lambda2(i,1)/mean(Obs(idx));
    MAE_ramp_lambda2(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_lambda2(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_lambda2(i,1) = MAE_ramp_lambda2(i,1)/mean(Obs(idx));
    ME_ramp2_lambda2(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_lambda2(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_lambda2(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_lambda2(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_lambda2(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_lambda2(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_lambda2(i,1) = ME_camp_lambda2(i,1).^2;
    SE2_camp_lambda2(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_lambda2(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_lambda2(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_lambda2(i,1) = R_camp_lambda2(i,1).^2;
    RMSE_camp_lambda2(i,1) = sqrt(MSE_camp_lambda2(i,1));
    MNE_camp_lambda2(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_lambda2(i,1) = ME_camp_lambda2(i,1)/mean(Obs(idx));
    MAE_camp_lambda2(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_lambda2(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_lambda2(i,1) = MAE_camp_lambda2(i,1)/mean(Obs(idx));
    ME_camp2_lambda2(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_lambda2(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_lambda2(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_lambda2(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_lambda2(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_lambda2(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_lambda2(i,1) = ME_con_lambda2(i,1).^2;
    SE2_con_lambda2(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_lambda2(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_lambda2(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_lambda2(i,1) = R_con_lambda2(i,1).^2;
    RMSE_con_lambda2(i,1) = sqrt(MSE_con_lambda2(i,1));
    MNE_con_lambda2(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_lambda2(i,1) = ME_con_lambda2(i,1)/mean(Obs(idx));
    MAE_con_lambda2(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_lambda2(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_lambda2(i,1) = MAE_con_lambda2(i,1)/mean(Obs(idx));
    ME_con2_lambda2(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_lambda2(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_lambda2(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_lambda2(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% by day
uniday = unique(yrmoda);
for i = 1:length(uniday)
    
    % unique day
    idx = uniday(i) == yrmoda;

    % CMAQ
    ME_cmaq_day(i,1) = mean(Mod(idx)-Obs(idx));
    SE_cmaq_day(i,1) = std(Mod(idx)-Obs(idx));
    ME2_cmaq_day(i,1) = ME_cmaq_day(i,1).^2;
    SE2_cmaq_day(i,1) = var(Mod(idx)-Obs(idx));
    MSE_cmaq_day(i,1) = mean((Mod(idx)-Obs(idx)).^2);
    R_cmaq_day(i,1) = corr(Mod(idx),Obs(idx));
    R2_cmaq_day(i,1) = R_cmaq_day(i,1).^2;
    RMSE_cmaq_day(i,1) = sqrt(MSE_cmaq_day(i,1));
    MNE_cmaq_day(i,1) = mean((Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_cmaq_day(i,1) = ME_cmaq_day(i,1)/mean(Obs(idx));
    MAE_cmaq_day(i,1) = mean(abs(Mod(idx)-Obs(idx)));
    MNAE_cmaq_day(i,1) = mean(abs(Mod(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_cmaq_day(i,1) = MAE_cmaq_day(i,1)/mean(Obs(idx));
    ME_cmaq2_day(i,1) = NaN;
    SE2_cmaq2_day(i,1) = NaN;
    RMSS_cmaq_day(i,1) = NaN;
    MR_cmaq_day(i,1) = NaN;

    % RAMP
    ME_ramp_day(i,1) = mean(lambda1(idx)-Obs(idx));
    SE_ramp_day(i,1) = std(lambda1(idx)-Obs(idx));
    ME2_ramp_day(i,1) = ME_ramp_day(i,1).^2;
    SE2_ramp_day(i,1) = var(lambda1(idx)-Obs(idx));
    MSE_ramp_day(i,1) = mean((lambda1(idx)-Obs(idx)).^2);
    R_ramp_day(i,1) = corr(lambda1(idx),Obs(idx));
    R2_ramp_day(i,1) = R_ramp_day(i,1).^2;
    RMSE_ramp_day(i,1) = sqrt(MSE_ramp_day(i,1));
    MNE_ramp_day(i,1) = mean((lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_ramp_day(i,1) = ME_ramp_day(i,1)/mean(Obs(idx));
    MAE_ramp_day(i,1) = mean(abs(lambda1(idx)-Obs(idx)));
    MNAE_ramp_day(i,1) = mean(abs(lambda1(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_ramp_day(i,1) = MAE_ramp_day(i,1)/mean(Obs(idx));
    ME_ramp2_day(i,1) = mean((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    SE2_ramp2_day(i,1) = var((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx)));
    RMSS_ramp_day(i,1) = sqrt(mean( ((Obs(idx)-lambda1(idx))./sqrt(lambda2(idx))).^2 ));
    MR_ramp_day(i,1) = mean(sqrt(lambda2(idx)));
    
    % CAMP
    ME_camp_day(i,1) = mean(lambda1CAMP(idx)-Obs(idx));
    SE_camp_day(i,1) = std(lambda1CAMP(idx)-Obs(idx));
    ME2_camp_day(i,1) = ME_camp_day(i,1).^2;
    SE2_camp_day(i,1) = var(lambda1CAMP(idx)-Obs(idx));
    MSE_camp_day(i,1) = mean((lambda1CAMP(idx)-Obs(idx)).^2);
    R_camp_day(i,1) = corr(lambda1CAMP(idx),Obs(idx));
    R2_camp_day(i,1) = R_camp_day(i,1).^2;
    RMSE_camp_day(i,1) = sqrt(MSE_camp_day(i,1));
    MNE_camp_day(i,1) = mean((lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_camp_day(i,1) = ME_camp_day(i,1)/mean(Obs(idx));
    MAE_camp_day(i,1) = mean(abs(lambda1CAMP(idx)-Obs(idx)));
    MNAE_camp_day(i,1) = mean(abs(lambda1CAMP(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_camp_day(i,1) = MAE_camp_day(i,1)/mean(Obs(idx));
    ME_camp2_day(i,1) = mean((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    SE2_camp2_day(i,1) = var((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx)));
    RMSS_camp_day(i,1) = sqrt(mean( ((Obs(idx)-lambda1CAMP(idx))./sqrt(lambda2CAMP(idx))).^2 ));
    MR_camp_day(i,1) = mean(sqrt(lambda2CAMP(idx)));
    
    % Constant
    ME_con_day(i,1) = mean(lambda1Constant(idx)-Obs(idx));
    SE_con_day(i,1) = std(lambda1Constant(idx)-Obs(idx));
    ME2_con_day(i,1) = ME_con_day(i,1).^2;
    SE2_con_day(i,1) = var(lambda1Constant(idx)-Obs(idx));
    MSE_con_day(i,1) = mean((lambda1Constant(idx)-Obs(idx)).^2);
    R_con_day(i,1) = corr(lambda1Constant(idx),Obs(idx));
    R2_con_day(i,1) = R_con_day(i,1).^2;
    RMSE_con_day(i,1) = sqrt(MSE_con_day(i,1));
    MNE_con_day(i,1) = mean((lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NME_con_day(i,1) = ME_con_day(i,1)/mean(Obs(idx));
    MAE_con_day(i,1) = mean(abs(lambda1Constant(idx)-Obs(idx)));
    MNAE_con_day(i,1) = mean(abs(lambda1Constant(idx&Obs>0)-Obs(idx&Obs>0))./Obs(idx&Obs>0));
    NMAE_con_day(i,1) = MAE_con_day(i,1)/mean(Obs(idx));
    ME_con2_day(i,1) = mean((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    SE2_con2_day(i,1) = var((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx)));
    RMSS_con_day(i,1) = sqrt(mean( ((Obs(idx)-lambda1Constant(idx))./sqrt(lambda2Constant(idx))).^2 ));
    MR_con_day(i,1) = mean(sqrt(lambda2Constant(idx)));
    
end

% save all results
save('stratCMAQRAMP.mat','ME_cmaq','SE_cmaq','ME2_cmaq','SE2_cmaq','MSE_cmaq','R_cmaq', ...
    'R2_cmaq','RMSE_cmaq','MNE_cmaq','NME_cmaq','MAE_cmaq','MNAE_cmaq','NMAE_cmaq','ME_cmaq2','SE2_cmaq2','RMSS_cmaq','MR_cmaq','MSNS_cmaq', ...
    'ME_ramp','SE_ramp','ME2_ramp','SE2_ramp','MSE_ramp','R_ramp', ...
    'R2_ramp','RMSE_ramp','MNE_ramp','NME_ramp','MAE_ramp','MNAE_ramp','NMAE_ramp','ME_ramp2','SE2_ramp2','RMSS_ramp','MR_ramp','MSNS_ramp', ... 
    'ME_camp','SE_camp','ME2_camp','SE2_camp','MSE_camp','R_camp', ...
    'R2_camp','RMSE_camp','MNE_camp','NME_camp','MAE_camp','MNAE_camp','NMAE_camp','ME_camp2','SE2_camp2','RMSS_camp','MR_camp','MSNS_camp', ...
    'ME_con','SE_con','ME2_con','SE2_con','MSE_con','R_con', ...
    'R2_con','RMSE_con','MNE_con','NME_con','MAE_con','MNAE_con','NMAE_con','ME_con2','SE2_con2','RMSS_con','MR_con','MSNS_con', ... 
    'season_str','ME_cmaq_season','SE_cmaq_season','ME2_cmaq_season','SE2_cmaq_season','MSE_cmaq_season','R_cmaq_season', ...
    'R2_cmaq_season','RMSE_cmaq_season','MNE_cmaq_season','NME_cmaq_season','MAE_cmaq_season','MNAE_cmaq_season','NMAE_cmaq_season','ME_cmaq2_season','SE2_cmaq2_season','RMSS_cmaq_season','MR_cmaq_season', ...
    'ME_ramp_season','SE_ramp_season','ME2_ramp_season','SE2_ramp_season','MSE_ramp_season','R_ramp_season', ...
    'R2_ramp_season','RMSE_ramp_season','MNE_ramp_season','NME_ramp_season','MAE_ramp_season','MNAE_ramp_season','NMAE_ramp_season','ME_ramp2_season','SE2_ramp2_season','RMSS_ramp_season','MR_ramp_season', ... 
    'ME_camp_season','SE_camp_season','ME2_camp_season','SE2_camp_season','MSE_camp_season','R_camp_season', ...
    'R2_camp_season','RMSE_camp_season','MNE_camp_season','NME_camp_season','MAE_camp_season','MNAE_camp_season','NMAE_camp_season','ME_camp2_season','SE2_camp2_season','RMSS_camp_season','MR_camp_season', ... 
    'ME_con_season','SE_con_season','ME2_con_season','SE2_con_season','MSE_con_season','R_con_season', ...
    'R2_con_season','RMSE_con_season','MNE_con_season','NME_con_season','MAE_con_season','MNAE_con_season','NMAE_con_season','ME_con2_season','SE2_con2_season','RMSS_con_season','MR_con_season', ... 
    'network_str','ME_cmaq_network','SE_cmaq_network','ME2_cmaq_network','SE2_cmaq_network','MSE_cmaq_network','R_cmaq_network', ...
    'R2_cmaq_network','RMSE_cmaq_network','MNE_cmaq_network','NME_cmaq_network','MAE_cmaq_network','MNAE_cmaq_network','NMAE_cmaq_network','ME_cmaq2_network','SE2_cmaq2_network','RMSS_cmaq_network','MR_cmaq_network', ...
    'ME_ramp_network','SE_ramp_network','ME2_ramp_network','SE2_ramp_network','MSE_ramp_network','R_ramp_network', ...
    'R2_ramp_network','RMSE_ramp_network','MNE_ramp_network','NME_ramp_network','MAE_ramp_network','MNAE_ramp_network','NMAE_ramp_network','ME_ramp2_network','SE2_ramp2_network','RMSS_ramp_network','MR_ramp_network', ... 
    'ME_camp_network','SE_camp_network','ME2_camp_network','SE2_camp_network','MSE_camp_network','R_camp_network', ...
    'R2_camp_network','RMSE_camp_network','MNE_camp_network','NME_camp_network','MAE_camp_network','MNAE_camp_network','NMAE_camp_network','ME_camp2_network','SE2_camp2_network','RMSS_camp_network','MR_camp_network', ...
    'ME_con_network','SE_con_network','ME2_con_network','SE2_con_network','MSE_con_network','R_con_network', ...
    'R2_con_network','RMSE_con_network','MNE_con_network','NME_con_network','MAE_con_network','MNAE_con_network','NMAE_con_network','ME_con2_network','SE2_con2_network','RMSS_con_network','MR_con_network', ...
    'frm_str','ME_cmaq_frm','SE_cmaq_frm','ME2_cmaq_frm','SE2_cmaq_frm','MSE_cmaq_frm','R_cmaq_frm', ...
    'R2_cmaq_frm','RMSE_cmaq_frm','MNE_cmaq_frm','NME_cmaq_frm','MAE_cmaq_frm','MNAE_cmaq_frm','NMAE_cmaq_frm','ME_cmaq2_frm','SE2_cmaq2_frm','RMSS_cmaq_frm','MR_cmaq_frm', ...
    'ME_ramp_frm','SE_ramp_frm','ME2_ramp_frm','SE2_ramp_frm','MSE_ramp_frm','R_ramp_frm', ...
    'R2_ramp_frm','RMSE_ramp_frm','MNE_ramp_frm','NME_ramp_frm','MAE_ramp_frm','MNAE_ramp_frm','NMAE_ramp_frm','ME_ramp2_frm','SE2_ramp2_frm','RMSS_ramp_frm','MR_ramp_frm', ... 
    'ME_camp_frm','SE_camp_frm','ME2_camp_frm','SE2_camp_frm','MSE_camp_frm','R_camp_frm', ...
    'R2_camp_frm','RMSE_camp_frm','MNE_camp_frm','NME_camp_frm','MAE_camp_frm','MNAE_camp_frm','NMAE_camp_frm','ME_camp2_frm','SE2_camp2_frm','RMSS_camp_frm','MR_camp_frm', ...
    'ME_con_frm','SE_con_frm','ME2_con_frm','SE2_con_frm','MSE_con_frm','R_con_frm', ...
    'R2_con_frm','RMSE_con_frm','MNE_con_frm','NME_con_frm','MAE_con_frm','MNAE_con_frm','NMAE_con_frm','ME_con2_frm','SE2_con2_frm','RMSS_con_frm','MR_con_frm', ... 
    'west_str','ME_cmaq_west','SE_cmaq_west','ME2_cmaq_west','SE2_cmaq_west','MSE_cmaq_west','R_cmaq_west', ...
    'R2_cmaq_west','RMSE_cmaq_west','MNE_cmaq_west','NME_cmaq_west','MAE_cmaq_west','MNAE_cmaq_west','NMAE_cmaq_west','ME_cmaq2_west','SE2_cmaq2_west','RMSS_cmaq_west','MR_cmaq_west', ...
    'ME_ramp_west','SE_ramp_west','ME2_ramp_west','SE2_ramp_west','MSE_ramp_west','R_ramp_west', ...
    'R2_ramp_west','RMSE_ramp_west','MNE_ramp_west','NME_ramp_west','MAE_ramp_west','MNAE_ramp_west','NMAE_ramp_west','ME_ramp2_west','SE2_ramp2_west','RMSS_ramp_west','MR_ramp_west', ... 
    'ME_camp_west','SE_camp_west','ME2_camp_west','SE2_camp_west','MSE_camp_west','R_camp_west', ...
    'R2_camp_west','RMSE_camp_west','MNE_camp_west','NME_camp_west','MAE_camp_west','MNAE_camp_west','NMAE_camp_west','ME_camp2_west','SE2_camp2_west','RMSS_camp_west','MR_camp_west', ...
    'ME_con_west','SE_con_west','ME2_con_west','SE2_con_west','MSE_con_west','R_con_west', ...
    'R2_con_west','RMSE_con_west','MNE_con_west','NME_con_west','MAE_con_west','MNAE_con_west','NMAE_con_west','ME_con2_west','SE2_con2_west','RMSS_con_west','MR_con_west', ... 
    'urban_str','ME_cmaq_urban','SE_cmaq_urban','ME2_cmaq_urban','SE2_cmaq_urban','MSE_cmaq_urban','R_cmaq_urban', ...
    'R2_cmaq_urban','RMSE_cmaq_urban','MNE_cmaq_urban','NME_cmaq_urban','MAE_cmaq_urban','MNAE_cmaq_urban','NMAE_cmaq_urban','ME_cmaq2_urban','SE2_cmaq2_urban','RMSS_cmaq_urban','MR_cmaq_urban', ...
    'ME_ramp_urban','SE_ramp_urban','ME2_ramp_urban','SE2_ramp_urban','MSE_ramp_urban','R_ramp_urban', ...
    'R2_ramp_urban','RMSE_ramp_urban','MNE_ramp_urban','NME_ramp_urban','MAE_ramp_urban','MNAE_ramp_urban','NMAE_ramp_urban','ME_ramp2_urban','SE2_ramp2_urban','RMSS_ramp_urban','MR_ramp_urban', ... 
    'ME_camp_urban','SE_camp_urban','ME2_camp_urban','SE2_camp_urban','MSE_camp_urban','R_camp_urban', ...
    'R2_camp_urban','RMSE_camp_urban','MNE_camp_urban','NME_camp_urban','MAE_camp_urban','MNAE_camp_urban','NMAE_camp_urban','ME_camp2_urban','SE2_camp2_urban','RMSS_camp_urban','MR_camp_urban', ...
    'ME_con_urban','SE_con_urban','ME2_con_urban','SE2_con_urban','MSE_con_urban','R_con_urban', ...
    'R2_con_urban','RMSE_con_urban','MNE_con_urban','NME_con_urban','MAE_con_urban','MNAE_con_urban','NMAE_con_urban','ME_con2_urban','SE2_con2_urban','RMSS_con_urban','MR_con_urban', ... 
    'region_str','ME_cmaq_region','SE_cmaq_region','ME2_cmaq_region','SE2_cmaq_region','MSE_cmaq_region','R_cmaq_region', ...
    'R2_cmaq_region','RMSE_cmaq_region','MNE_cmaq_region','NME_cmaq_region','MAE_cmaq_region','MNAE_cmaq_region','NMAE_cmaq_region','ME_cmaq2_region','SE2_cmaq2_region','RMSS_cmaq_region','MR_cmaq_region', ...
    'ME_ramp_region','SE_ramp_region','ME2_ramp_region','SE2_ramp_region','MSE_ramp_region','R_ramp_region', ...
    'R2_ramp_region','RMSE_ramp_region','MNE_ramp_region','NME_ramp_region','MAE_ramp_region','MNAE_ramp_region','NMAE_ramp_region','ME_ramp2_region','SE2_ramp2_region','RMSS_ramp_region','MR_ramp_region', ... 
    'ME_camp_region','SE_camp_region','ME2_camp_region','SE2_camp_region','MSE_camp_region','R_camp_region', ...
    'R2_camp_region','RMSE_camp_region','MNE_camp_region','NME_camp_region','MAE_camp_region','MNAE_camp_region','NMAE_camp_region','ME_camp2_region','SE2_camp2_region','RMSS_camp_region','MR_camp_region', ...
    'ME_con_region','SE_con_region','ME2_con_region','SE2_con_region','MSE_con_region','R_con_region', ...
    'R2_con_region','RMSE_con_region','MNE_con_region','NME_con_region','MAE_con_region','MNAE_con_region','NMAE_con_region','ME_con2_region','SE2_con2_region','RMSS_con_region','MR_con_region', ... 
    'lendist','ME_cmaq_dist','SE_cmaq_dist','ME2_cmaq_dist','SE2_cmaq_dist','MSE_cmaq_dist','R_cmaq_dist', ...
    'R2_cmaq_dist','RMSE_cmaq_dist','MNE_cmaq_dist','NME_cmaq_dist','MAE_cmaq_dist','MNAE_cmaq_dist','NMAE_cmaq_dist','ME_cmaq2_dist','SE2_cmaq2_dist','RMSS_cmaq_dist','MR_cmaq_dist', ...
    'ME_ramp_dist','SE_ramp_dist','ME2_ramp_dist','SE2_ramp_dist','MSE_ramp_dist','R_ramp_dist', ...
    'R2_ramp_dist','RMSE_ramp_dist','MNE_ramp_dist','NME_ramp_dist','MAE_ramp_dist','MNAE_ramp_dist','NMAE_ramp_dist','ME_ramp2_dist','SE2_ramp2_dist','RMSS_ramp_dist','MR_ramp_dist', ... 
    'ME_camp_dist','SE_camp_dist','ME2_camp_dist','SE2_camp_dist','MSE_camp_dist','R_camp_dist', ...
    'R2_camp_dist','RMSE_camp_dist','MNE_camp_dist','NME_camp_dist','MAE_camp_dist','MNAE_camp_dist','NMAE_camp_dist','ME_camp2_dist','SE2_camp2_dist','RMSS_camp_dist','MR_camp_dist', ...
    'ME_con_dist','SE_con_dist','ME2_con_dist','SE2_con_dist','MSE_con_dist','R_con_dist', ...
    'R2_con_dist','RMSE_con_dist','MNE_con_dist','NME_con_dist','MAE_con_dist','MNAE_con_dist','NMAE_con_dist','ME_con2_dist','SE2_con2_dist','RMSS_con_dist','MR_con_dist', ... 
    'lenobs','ME_cmaq_obs','SE_cmaq_obs','ME2_cmaq_obs','SE2_cmaq_obs','MSE_cmaq_obs','R_cmaq_obs', ...
    'R2_cmaq_obs','RMSE_cmaq_obs','MNE_cmaq_obs','NME_cmaq_obs','MAE_cmaq_obs','MNAE_cmaq_obs','NMAE_cmaq_obs','ME_cmaq2_obs','SE2_cmaq2_obs','RMSS_cmaq_obs','MR_cmaq_obs', ...
    'ME_ramp_obs','SE_ramp_obs','ME2_ramp_obs','SE2_ramp_obs','MSE_ramp_obs','R_ramp_obs', ...
    'R2_ramp_obs','RMSE_ramp_obs','MNE_ramp_obs','NME_ramp_obs','MAE_ramp_obs','MNAE_ramp_obs','NMAE_ramp_obs','ME_ramp2_obs','SE2_ramp2_obs','RMSS_ramp_obs','MR_ramp_obs', ... 
    'ME_camp_obs','SE_camp_obs','ME2_camp_obs','SE2_camp_obs','MSE_camp_obs','R_camp_obs', ...
    'R2_camp_obs','RMSE_camp_obs','MNE_camp_obs','NME_camp_obs','MAE_camp_obs','MNAE_camp_obs','NMAE_camp_obs','ME_camp2_obs','SE2_camp2_obs','RMSS_camp_obs','MR_camp_obs', ...
    'ME_con_obs','SE_con_obs','ME2_con_obs','SE2_con_obs','MSE_con_obs','R_con_obs', ...
    'R2_con_obs','RMSE_con_obs','MNE_con_obs','NME_con_obs','MAE_con_obs','MNAE_con_obs','NMAE_con_obs','ME_con2_obs','SE2_con2_obs','RMSS_con_obs','MR_con_obs', ... 
    'lenmod','ME_cmaq_mod','SE_cmaq_mod','ME2_cmaq_mod','SE2_cmaq_mod','MSE_cmaq_mod','R_cmaq_mod', ...
    'R2_cmaq_mod','RMSE_cmaq_mod','MNE_cmaq_mod','NME_cmaq_mod','MAE_cmaq_mod','MNAE_cmaq_mod','NMAE_cmaq_mod','ME_cmaq2_mod','SE2_cmaq2_mod','RMSS_cmaq_mod','MR_cmaq_mod', ...
    'ME_ramp_mod','SE_ramp_mod','ME2_ramp_mod','SE2_ramp_mod','MSE_ramp_mod','R_ramp_mod', ...
    'R2_ramp_mod','RMSE_ramp_mod','MNE_ramp_mod','NME_ramp_mod','MAE_ramp_mod','MNAE_ramp_mod','NMAE_ramp_mod','ME_ramp2_mod','SE2_ramp2_mod','RMSS_ramp_mod','MR_ramp_mod', ... 
    'ME_camp_mod','SE_camp_mod','ME2_camp_mod','SE2_camp_mod','MSE_camp_mod','R_camp_mod', ...
    'R2_camp_mod','RMSE_camp_mod','MNE_camp_mod','NME_camp_mod','MAE_camp_mod','MNAE_camp_mod','NMAE_camp_mod','ME_camp2_mod','SE2_camp2_mod','RMSS_camp_mod','MR_camp_mod', ...
    'ME_con_mod','SE_con_mod','ME2_con_mod','SE2_con_mod','MSE_con_mod','R_con_mod', ...
    'R2_con_mod','RMSE_con_mod','MNE_con_mod','NME_con_mod','MAE_con_mod','MNAE_con_mod','NMAE_con_mod','ME_con2_mod','SE2_con2_mod','RMSS_con_mod','MR_con_mod', ... 
    'lenlambda2','ME_cmaq_lambda2','SE_cmaq_lambda2','ME2_cmaq_lambda2','SE2_cmaq_lambda2','MSE_cmaq_lambda2','R_cmaq_lambda2', ...
    'R2_cmaq_lambda2','RMSE_cmaq_lambda2','MNE_cmaq_lambda2','NME_cmaq_lambda2','MAE_cmaq_lambda2','MNAE_cmaq_lambda2','NMAE_cmaq_lambda2','ME_cmaq2_lambda2','SE2_cmaq2_lambda2','RMSS_cmaq_lambda2','MR_cmaq_lambda2', ...
    'ME_ramp_lambda2','SE_ramp_lambda2','ME2_ramp_lambda2','SE2_ramp_lambda2','MSE_ramp_lambda2','R_ramp_lambda2', ...
    'R2_ramp_lambda2','RMSE_ramp_lambda2','MNE_ramp_lambda2','NME_ramp_lambda2','MAE_ramp_lambda2','MNAE_ramp_lambda2','NMAE_ramp_lambda2','ME_ramp2_lambda2','SE2_ramp2_lambda2','RMSS_ramp_lambda2','MR_ramp_lambda2', ...
    'ME_camp_lambda2','SE_camp_lambda2','ME2_camp_lambda2','SE2_camp_lambda2','MSE_camp_lambda2','R_camp_lambda2', ...
    'R2_camp_lambda2','RMSE_camp_lambda2','MNE_camp_lambda2','NME_camp_lambda2','MAE_camp_lambda2','MNAE_camp_lambda2','NMAE_camp_lambda2','ME_camp2_lambda2','SE2_camp2_lambda2','RMSS_camp_lambda2','MR_camp_lambda2', ...
    'ME_con_lambda2','SE_con_lambda2','ME2_con_lambda2','SE2_con_lambda2','MSE_con_lambda2','R_con_lambda2', ...
    'R2_con_lambda2','RMSE_con_lambda2','MNE_con_lambda2','NME_con_lambda2','MAE_con_lambda2','MNAE_con_lambda2','NMAE_con_lambda2','ME_con2_lambda2','SE2_con2_lambda2','RMSS_con_lambda2','MR_con_lambda2', ...
    'uniday','ME_cmaq_day','SE_cmaq_day','ME2_cmaq_day','SE2_cmaq_day','MSE_cmaq_day','R_cmaq_day', ...
    'R2_cmaq_day','RMSE_cmaq_day','MNE_cmaq_day','NME_cmaq_day','MAE_cmaq_day','MNAE_cmaq_day','NMAE_cmaq_day','ME_cmaq2_day','SE2_cmaq2_day','RMSS_cmaq_day','MR_cmaq_day', ...
    'ME_ramp_day','SE_ramp_day','ME2_ramp_day','SE2_ramp_day','MSE_ramp_day','R_ramp_day', ...
    'R2_ramp_day','RMSE_ramp_day','MNE_ramp_day','NME_ramp_day','MAE_ramp_day','MNAE_ramp_day','NMAE_ramp_day','ME_ramp2_day','SE2_ramp2_day','RMSS_ramp_day','MR_ramp_day', ...
    'ME_camp_day','SE_camp_day','ME2_camp_day','SE2_camp_day','MSE_camp_day','R_camp_day', ...
    'R2_camp_day','RMSE_camp_day','MNE_camp_day','NME_camp_day','MAE_camp_day','MNAE_camp_day','NMAE_camp_day','ME_camp2_day','SE2_camp2_day','RMSS_camp_day','MR_camp_day', ...
    'ME_con_day','SE_con_day','ME2_con_day','SE2_con_day','MSE_con_day','R_con_day', ...
    'R2_con_day','RMSE_con_day','MNE_con_day','NME_con_day','MAE_con_day','MNAE_con_day','NMAE_con_day','ME_con2_day','SE2_con2_day','RMSS_con_day','MR_con_day', ...
    'ME_campS_season','SE_campS_season','ME2_campS_season','SE2_campS_season','MSE_campS_season','R_campS_season', ...
    'R2_campS_season','RMSE_campS_season','MNE_campS_season','NME_campS_season','MAE_campS_season','MNAE_campS_season','NMAE_campS_season','ME_campS2_season','SE2_campS2_season','RMSS_campS_season','MR_campS_season', ...
    'ME_conS_season','SE_conS_season','ME2_conS_season','SE2_conS_season','MSE_conS_season','R_conS_season', ...
    'R2_conS_season','RMSE_conS_season','MNE_conS_season','NME_conS_season','MAE_conS_season','MNAE_conS_season','NMAE_conS_season','ME_conS2_season','SE2_conS2_season','RMSS_conS_season','MR_conS_season', ...
    'region6_str','ME_cmaq_region6','SE_cmaq_region6','ME2_cmaq_region6','SE2_cmaq_region6','MSE_cmaq_region6','R_cmaq_region6', ...
    'R2_cmaq_region6','RMSE_cmaq_region6','MNE_cmaq_region6','NME_cmaq_region6','MAE_cmaq_region6','MNAE_cmaq_region6','NMAE_cmaq_region6','ME_cmaq2_region6','SE2_cmaq2_region6','RMSS_cmaq_region6','MR_cmaq_region6', ...
    'ME_ramp_region6','SE_ramp_region6','ME2_ramp_region6','SE2_ramp_region6','MSE_ramp_region6','R_ramp_region6', ...
    'R2_ramp_region6','RMSE_ramp_region6','MNE_ramp_region6','NME_ramp_region6','MAE_ramp_region6','MNAE_ramp_region6','NMAE_ramp_region6','ME_ramp2_region6','SE2_ramp2_region6','RMSS_ramp_region6','MR_ramp_region6', ... 
    'ME_camp_region6','SE_camp_region6','ME2_camp_region6','SE2_camp_region6','MSE_camp_region6','R_camp_region6', ...
    'R2_camp_region6','RMSE_camp_region6','MNE_camp_region6','NME_camp_region6','MAE_camp_region6','MNAE_camp_region6','NMAE_camp_region6','ME_camp2_region6','SE2_camp2_region6','RMSS_camp_region6','MR_camp_region6', ...
    'ME_con_region6','SE_con_region6','ME2_con_region6','SE2_con_region6','MSE_con_region6','R_con_region6', ...
    'R2_con_region6','RMSE_con_region6','MNE_con_region6','NME_con_region6','MAE_con_region6','MNAE_con_region6','NMAE_con_region6','ME_con2_region6','SE2_con2_region6','RMSS_con_region6','MR_con_region6', ... 
    'ME_camp6_region6','SE_camp6_region6','ME2_camp6_region6','SE2_camp6_region6','MSE_camp6_region6','R_camp6_region6', ...
    'R2_camp6_region6','RMSE_camp6_region6','MNE_camp6_region6','NME_camp6_region6','MAE_camp6_region6','MNAE_camp6_region6','NMAE_camp6_region6','ME_camp62_region6','SE2_camp62_region6','RMSS_camp6_region6','MR_camp6_region6', ...
    'ME_con6_region6','SE_con6_region6','ME2_con6_region6','SE2_con6_region6','MSE_con6_region6','R_con6_region6', ...
    'R2_con6_region6','RMSE_con6_region6','MNE_con6_region6','NME_con6_region6','MAE_con6_region6','MNAE_con6_region6','NMAE_con6_region6','ME_con62_region6','SE2_con62_region6','RMSS_con6_region6','MR_con6_region6'); 

end