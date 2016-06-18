function [] = exclusionCriteria()
% this function will list the possible exclusion criteria:
% by MP metric, by obs cut offs, by mod cut offs, by lambda1/lambda2 cut
% offs, geographic location, time

% metrics to visualize with % change from kriging to BME: R2, MSE

% load info
load('matfiles/allInfo.mat');

% zh: works best for higher values of pm2.5
% CMAQ                            
% Dist2NMon                        
% Is105,Is110,Is115,Is120,Is125,Is90,Is95                             
% IsFRM,IsIMPROVE,IsSTN,IsTEOM              
% IsRural,IsSuburban,IsUrban                                                    
% IsSpring, IsSummer,IsWinter,IsFall                                                                                                                                            
% IsitWest                         
% R2_2                             
% R_2                               
% beta1_2                                                
% fBias_2                         
% fErr_2                         
% in                                      
% m2DmsBias_2                      
% mBias_2                            
% mDsBias_2                         
% mErr_2                           
% mMod_2                            
% mObs_2                            
% msBias_2                          
% nBias_2                         
% nErr_2                   
% nmBias_2                          
% nmErr_2                          
% nrmsBias_2                       
% rmsBias_2                        
% s2DmsBias_2                       
% sBias_2                            
% vMod_2                            
% vObs_2                                                                      
             
% cut offs
cutoffs = [30 35 40];
cutoffs2 = [35 40 45];
cutoffs = [prctile((mBias_2).^2,80) prctile((mBias_2).^2,85) prctile((mBias_2).^2,90)];
colorz = {'ro-','bo-','go-'};

% % change in MSE from kriging to BME
% each cell is a cut off, first column krig, 2nd column BME, row is radii
MSE_krig2BME = cell(4,1);
for i = 1:length(cutoffs)
    
    MSE_krig2BME{i,1} = NaN*ones(10,2);
    idxcut = zh>=cutoffs(i)&zh<cutoffs2(i);
    idxcut = CMAQ>=cutoffs(i)&CMAQ<cutoffs2(i);
    idxcut = (mBias_2).^2>=cutoffs(i);
    disp(sum(idxcut));
        
    for j = 1:10
        if j == 1
            idxnan = ~isnan(zk_hard_0km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 2
            idxnan = ~isnan(zk_hard_100km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 3
            idxnan = ~isnan(zk_hard_200km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 4
            idxnan = ~isnan(zk_hard_300km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 5
            idxnan = ~isnan(zk_hard_400km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 6
            idxnan = ~isnan(zk_hard_500km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 7
            idxnan = ~isnan(zk_hard_600km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 8
            idxnan = ~isnan(zk_hard_700km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 9
            idxnan = ~isnan(zk_hard_800km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        elseif j == 10
            idxnan = ~isnan(zk_hard_900km);
            MSE_krig2BME{i,1}(j,1) = mean((zk_hard_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
            MSE_krig2BME{i,1}(j,2) = mean((zk_soft_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
        end        
    end   
    
end

figure; hold on;
for i = 1:length(cutoffs)
    plot(0:100:900,100.*(MSE_krig2BME{i}(:,2)-MSE_krig2BME{i}(:,1))./(MSE_krig2BME{i}(:,1)),colorz{i});
    strz{i} = sprintf('>=%0.1f',cutoffs(i));
end
title('Percent change in MSE from kriging to BME by: cutoffs');
xlabel('km');
legend(strz,'location','best');

blah = 5;

% perform cut offs by percentile of values

% loop through all each exclusion criteria individually and then jointly

% another to calculate all the statistics for each method, exclusion,
% increasing radius for each individual exclusions and joint exclusions


% another to make all the plots






end