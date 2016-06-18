function [] = exclusionCriteria_continuous_singular_intervals()
% this will loop through each continuous exclusion criteria singularly

% list of variables I can loop through:
% zh
% CMAQ 
% Dist2NMon
% lambda1_2
% lambda2_2
% R2_2                             
% R_2                               
% beta1_2                                                
% fBias_2                         
% fErr_2                                                              
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

% load info
load('matfiles/allInfo.mat');
mBias2_2 = (mBias_2).^2;
colorz = {'ro-','bo-','go-','ko-'};

%%% change in MSE from kriging to BME
% each cell is a variables, each subcel is a cut off, first column krig, 
% 2nd column BME, row is radii

% all variables
allvariables = {zh,CMAQ,Dist2NMon,lambda1_2,lambda2_2,R2_2,R_2,beta1_2,fBias_2,fErr_2,m2DmsBias_2, ...
    mBias2_2,mDsBias_2,mErr_2,mMod_2,mObs_2,msBias_2,nBias_2,nErr_2, ...                   
    nmBias_2,nmErr_2,nrmsBias_2,rmsBias_2,s2DmsBias_2,sBias_2,vMod_2,vObs_2};
allvariablesstr = {'zh','CMAQ','Dist2NMon','lambda1_2','lambda2_2','R2_2','R_2','beta1_2','fBias_2','fErr_2','m2DmsBias_2', ...
    'mBias2_2','mDsBias_2','mErr_2','mMod_2','mObs_2','msBias_2','nBias_2','nErr_2', ...                   
    'nmBias_2','nmErr_2','nrmsBias_2','rmsBias_2','s2DmsBias_2','sBias_2','vMod_2','vObs_2'};
MSE_krig2BME = cell(length(allvariables),1);
R_krig2BME = cell(length(allvariables),1);
cutoffs = cell(length(allvariables),1);
cutoffs2 = cell(length(allvariables),1);

% loop through each variable
for i = 1:length(allvariables)
    
    cutoffs{i} = [prctile(allvariables{i},80) prctile(allvariables{i},85) prctile(allvariables{i},90)];
    cutoffs2{i} = [prctile(allvariables{i},85) prctile(allvariables{i},90) prctile(allvariables{i},95)];
    MSE_krig2BME{i} = cell(length(cutoffs{i}),1);
    R_krig2BME{i} = cell(length(cutoffs{i}),1);
    
    % loop through each cutoff
    for j = 1:length(cutoffs{i})

        MSE_krig2BME{i}{j,1} = NaN*ones(10,2);
        R_krig2BME{i}{j,1} = NaN*ones(10,2);
        idxcut = allvariables{i}>=cutoffs{i}(j)&allvariables{i}<cutoffs2{i}(j);
        
        % loop through each increasing radii
        for k = 1:10
            if k == 1
                idxnan = ~isnan(zk_hard_0km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_0km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_0km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 2
                idxnan = ~isnan(zk_hard_100km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_100km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_100km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 3
                idxnan = ~isnan(zk_hard_200km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_200km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_200km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 4
                idxnan = ~isnan(zk_hard_300km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_300km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_300km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 5
                idxnan = ~isnan(zk_hard_400km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_400km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_400km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 6
                idxnan = ~isnan(zk_hard_500km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_500km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_500km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 7
                idxnan = ~isnan(zk_hard_600km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_600km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_600km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 8
                idxnan = ~isnan(zk_hard_700km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_700km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_700km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 9
                idxnan = ~isnan(zk_hard_800km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_800km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_800km(idxnan&idxcut),zh(idxnan&idxcut));
            elseif k == 10
                idxnan = ~isnan(zk_hard_900km);
                MSE_krig2BME{i}{j,1}(k,1) = mean((zk_hard_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                MSE_krig2BME{i}{j,1}(k,2) = mean((zk_soft_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                R_krig2BME{i}{j,1}(k,1) = corr(zk_hard_900km(idxnan&idxcut),zh(idxnan&idxcut));
                R_krig2BME{i}{j,1}(k,2) = corr(zk_soft_900km(idxnan&idxcut),zh(idxnan&idxcut));
            end        
        end   

    end

end

% save info
save('matfiles/exclus_cont_sing_intervals.mat','allvariables','allvariablesstr', ...
    'MSE_krig2BME','R_krig2BME');

% display each variable: MSE
for i = 1:length(allvariables)
    
    % plot figure
    figure; hold on;
    for j = 1:length(cutoffs{i})
        plot(0:100:900,100.*(MSE_krig2BME{i}{j}(:,2)-MSE_krig2BME{i}{j}(:,1))./(MSE_krig2BME{i}{j}(:,1)),colorz{j});
        if j == 1, strz{j} = 'blaseline';
        else, strz{j} = sprintf('>=%0.1f',cutoffs{i}(j)); end
    end
    title(sprintf('%% change MSE from kriging to BME by intervals: %s',allvariablesstr{i}));
    xlabel('km');
    legend(strz,'location','best');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/MSE_exclus_cont_sing_intervals_%s.png',allvariablesstr{i}));

end
close all;

% display each variable: R
for i = 1:length(allvariables)
    
    % plot figure
    figure; hold on;
    for j = 1:length(cutoffs{i})
        plot(0:100:900,100.*(R_krig2BME{i}{j}(:,2)-R_krig2BME{i}{j}(:,1))./(R_krig2BME{i}{j}(:,1)),colorz{j});
        if j == 1, strz{j} = 'blaseline';
        else, strz{j} = sprintf('>=%0.1f',cutoffs{i}(j)); end
    end
    title(sprintf('%% change R from kriging to BME by intervals: %s',allvariablesstr{i}));
    xlabel('km');
    legend(strz,'location','best');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/R_exclus_cont_sing_intervals_%s.png',allvariablesstr{i}));

end
close all;

end