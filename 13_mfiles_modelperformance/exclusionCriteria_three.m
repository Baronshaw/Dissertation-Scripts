function [] = exclusionCriteria_three()
% this will loop through each continuous exclusion criteria jointly. From
% the previous scripts, a few stand out variables make it through a first
% pass to be considered to look at variables jointly. Another lesson
% learned from the previous scrips is that cutoffs do better than intevals
% and MSE makes more sense than R.

% list of variables I can loop through jointly:
% lambda1_2
% m2DmsBias_2                      
                           
% load info
load('matfiles/allInfo.mat');
colorz1 = {'ro-','bo-','go-','ko-'};
colorz2 = {'rs-','bs-','gs-','ks-'};
colorz3 = {'r*-','b*-','g*-','k*-'};
colorz4 = {'rx-','bx-','gx-','kx-'};

%%% change in MSE from kriging to BME (and other methods)
% each cell is a variables, each subcel is a cut off, first column krig, 
% 2nd column BME, 3rd column CAMP, 4th column DS, 5th column CMAQ, row is radii

% all variables
allvariables = {lambda1_2,m2DmsBias_2,CMAQ,zh,lambda2_2,lambda1_2./lambda2_2};
allvariablesstr = {'lambda1','ME2dMSE','CMAQ','zh','lambda2','lambda1Dlambda2'};
MSE_krig2BME = cell(length(allvariables),length(allvariables));
R_krig2BME = cell(length(allvariables),length(allvariables));
ME_krig2BME = cell(length(allvariables),length(allvariables));
pntsNcalc = cell(length(allvariables),length(allvariables));

% get cut offs
for i = 1:length(allvariables)
    cutoffs{i} = [prctile(allvariables{i},0) prctile(allvariables{i},10) ...
            prctile(allvariables{i},20) prctile(allvariables{i},30) ...
            prctile(allvariables{i},40) prctile(allvariables{i},50) ...
            prctile(allvariables{i},60) prctile(allvariables{i},70) ...
            prctile(allvariables{i},80) prctile(allvariables{i},90) ...
            prctile(allvariables{i},100)];
end

% loop through each variable
for i = 1:length(allvariables)-1
    
    for j = i+1:length(allvariables)

        disp([i j]);
        MSE_krig2BME{i,j} = cell(length(cutoffs{i}),length(cutoffs{j}));
        R_krig2BME{i,j} = cell(length(cutoffs{i}),length(cutoffs{j}));
        ME_krig2BME{i,j} = cell(length(cutoffs{i}),length(cutoffs{j}));
        pntsNcalc{i,j} = cell(length(cutoffs{i}),length(cutoffs{j}));

        % loop through each cutoff
        for k = 1:length(cutoffs{i})-1
            for l = 1:length(cutoffs{j})-1

                MSE_krig2BME{i,j}{k,l} = NaN*ones(10,5);
                R_krig2BME{i,j}{k,l} = NaN*ones(10,5);
                ME_krig2BME{i,j}{k,l} = NaN*ones(10,5);
                pntsNcalc{i,j}{k,l} = NaN*ones(10,5);
                
                if i == 5
                    idxcut = allvariables{i}<cutoffs{i}(k+1) & allvariables{j}>=cutoffs{j}(l+1);
                elseif j == 5
                    idxcut = allvariables{i}>=cutoffs{i}(k) & allvariables{j}<cutoffs{j}(l+1);
                else
                    idxcut = allvariables{i}>=cutoffs{i}(k) & allvariables{j}>=cutoffs{j}(l);
                end

                % loop through each increasing radii
                for m = 1:10
                    if m == 1
                        idxnan = ~isnan(zk_hard_0km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_0km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_0km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_0km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_0km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_0km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_0km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_0km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_0km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_0km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 2
                        idxnan = ~isnan(zk_hard_100km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_100km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_100km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_100km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_100km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_100km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_100km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_100km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_100km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_100km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 3
                        idxnan = ~isnan(zk_hard_200km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_200km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_200km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_200km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_200km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_200km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_200km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_200km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_200km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_200km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 4
                        idxnan = ~isnan(zk_hard_300km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_300km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_300km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_300km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_300km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_300km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_300km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_300km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_300km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_300km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 5
                        idxnan = ~isnan(zk_hard_400km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_400km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_400km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_400km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_400km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_400km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_400km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_400km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_400km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_400km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 6
                        idxnan = ~isnan(zk_hard_500km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_500km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_500km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_500km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_500km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_500km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_500km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_500km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_500km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_500km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 7
                        idxnan = ~isnan(zk_hard_600km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_600km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_600km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_600km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_600km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_600km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_600km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_600km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_600km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_600km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 8
                        idxnan = ~isnan(zk_hard_700km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_700km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_700km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_700km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_700km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_700km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_700km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_700km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_700km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_700km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 9
                        idxnan = ~isnan(zk_hard_800km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_800km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_800km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_800km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_800km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_800km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_800km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_800km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_800km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_800km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    elseif m == 10
                        idxnan = ~isnan(zk_hard_900km);
                        idxnan2 = ~isnan(zk_STDSI_0km) & ~isnan(CMAQ);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_900km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)).^2);
                        
                        R_krig2BME{i,j}{k,l}(m,1) = corr(zk_hard_900km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,2) = corr(zk_soft_900km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,3) = corr(zk_camp_900km(idxnan&idxcut),zh(idxnan&idxcut));
                        R_krig2BME{i,j}{k,l}(m,4) = corr(zk_STDSI_900km(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        R_krig2BME{i,j}{k,l}(m,5) = corr(CMAQ(idxnan&idxcut&idxnan2),zh(idxnan&idxcut&idxnan2));
                        
                        ME_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_900km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_900km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,3) = mean((zk_camp_900km(idxnan&idxcut)-zh(idxnan&idxcut)));
                        ME_krig2BME{i,j}{k,l}(m,4) = mean((zk_STDSI_900km(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        ME_krig2BME{i,j}{k,l}(m,5) = mean((CMAQ(idxnan&idxcut&idxnan2)-zh(idxnan&idxcut&idxnan2)));
                        
                        pntsNcalc{i,j}{k,l}(m,1:3) = sum(idxnan&idxcut); 
                        pntsNcalc{i,j}{k,l}(m,4:5) = sum(idxnan&idxcut&idxnan2);
                    end        
                end 
                
            end
        end

    end
end

% save info
save('matfiles/exclus_three.mat','allvariables','allvariablesstr', ...
    'MSE_krig2BME','R_krig2BME','ME_krig2BME','cutoffs','pntsNcalc');

% % display each variable: MSE
% for i = 1:length(allvariables)-1
%     for j = i+1:length(allvariables)
%     
%         % plot figure
%         figure; hold on;
%         n = 1;
%         strz = [];
%         for k = 1:length(cutoffs{i})
%             
%             if k == 1, colorz = colorz1;
%             elseif k == 2, colorz = colorz2;
%             elseif k == 3, colorz = colorz3; 
%             elseif k == 4, colorz = colorz4; end
%             
%             for l = 1:length(cutoffs{j})
%                 plot(0:100:900,100.*(MSE_krig2BME{i,j}{k,l}(:,2)-MSE_krig2BME{i,j}{k,l}(:,1))./(MSE_krig2BME{i,j}{k,l}(:,1)), ...
%                     colorz{l});
%                 if k==1&l==1, strz{n} = 'baseline';
%                 else strz{n} = sprintf('>=%0.3f, >=%0.3f',cutoffs{i}(k),cutoffs{j}(l)); end
%                 n = n + 1;
%             end
%             
%         end
%         
%         title(sprintf('%% change MSE from kriging to BME by cutoffs: %s, %s', ...
%             allvariablesstr{i},allvariablesstr{j}));
%         xlabel('km');
%         legend(strz,'location','best');
% 
%         % save figure
%         set(gcf,'Position',[0 0 800 600]);
%         print(gcf,'-painters','-dpng','-r600',sprintf('figures/MSE_exclus_cont_double_cutoffs_%s_%s.png', ...
%             allvariablesstr{i},allvariablesstr{j}));
%         
%     end
% end
% close all;

end