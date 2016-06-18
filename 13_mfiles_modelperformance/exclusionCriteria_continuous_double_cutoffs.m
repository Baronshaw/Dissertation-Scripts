function [] = exclusionCriteria_continuous_double_cutoffs()
% this will loop through each continuous exclusion criteria jointly. From
% the previous scripts, a few stand out variables make it through a first
% pass to be considered to look at variables jointly. Another lesson
% learned from the previous scrips is that cutoffs do better than intevals
% and MSE makes more sense than R.

% list of variables I can loop through jointly:
% CMAQ 
% lambda1_2
% m2DmsBias_2                      
% mBias_2                                                    
% s2DmsBias_2                       
% sBias_2                            

% load info
load('matfiles/allInfo.mat');
mBias2_2 = (mBias_2).^2;
colorz1 = {'ro-','bo-','go-','ko-'};
colorz2 = {'rs-','bs-','gs-','ks-'};
colorz3 = {'r*-','b*-','g*-','k*-'};
colorz4 = {'rx-','bx-','gx-','kx-'};

%%% change in MSE from kriging to BME
% each cell is a variables, each subcel is a cut off, first column krig, 
% 2nd column BME, row is radii

% all variables
allvariables = {CMAQ,lambda1_2,m2DmsBias_2,mBias2_2,s2DmsBias_2,sBias_2};
allvariablesstr = {'CMAQ','lambda1','ME2dMSE','ME2','SE2dMSE','SE'};
MSE_krig2BME = cell(length(allvariables),length(allvariables));

% get cut offs
for i = 1:length(allvariables)
    cutoffs{i} = [prctile(allvariables{i},0) prctile(allvariables{i},80) ...
            prctile(allvariables{i},85) prctile(allvariables{i},90)];
end

% loop through each variable
for i = 1:length(allvariables)-1
    
    for j = i+1:length(allvariables)

        
        MSE_krig2BME{i,j} = cell(length(cutoffs{i}),length(cutoffs{j}));

        % loop through each cutoff
        for k = 1:length(cutoffs{i})
            for l = 1:length(cutoffs{j})

                MSE_krig2BME{i,j}{k,l} = NaN*ones(10,2);
                idxcut = allvariables{i}>=cutoffs{i}(k) & allvariables{j}>=cutoffs{j}(l);

                % loop through each increasing radii
                for m = 1:10
                    if m == 1
                        idxnan = ~isnan(zk_hard_0km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_0km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 2
                        idxnan = ~isnan(zk_hard_100km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_100km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 3
                        idxnan = ~isnan(zk_hard_200km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_200km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 4
                        idxnan = ~isnan(zk_hard_300km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_300km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 5
                        idxnan = ~isnan(zk_hard_400km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_400km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 6
                        idxnan = ~isnan(zk_hard_500km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_500km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 7
                        idxnan = ~isnan(zk_hard_600km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_600km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 8
                        idxnan = ~isnan(zk_hard_700km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_700km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 9
                        idxnan = ~isnan(zk_hard_800km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_800km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    elseif m == 10
                        idxnan = ~isnan(zk_hard_900km);
                        MSE_krig2BME{i,j}{k,l}(m,1) = mean((zk_hard_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                        MSE_krig2BME{i,j}{k,l}(m,2) = mean((zk_soft_900km(idxnan&idxcut)-zh(idxnan&idxcut)).^2);
                    end        
                end 
                
            end
        end

    end
end

% save info
save('matfiles/exclus_cont_double_cutoffs.mat','allvariables','allvariablesstr', ...
    'MSE_krig2BME','cutoffs');

% display each variable: MSE
for i = 1:length(allvariables)-1
    for j = i+1:length(allvariables)
    
        % plot figure
        figure; hold on;
        n = 1;
        strz = [];
        for k = 1:length(cutoffs{i})
            
            if k == 1, colorz = colorz1;
            elseif k == 2, colorz = colorz2;
            elseif k == 3, colorz = colorz3; 
            elseif k == 4, colorz = colorz4; end
            
            for l = 1:length(cutoffs{j})
                plot(0:100:900,100.*(MSE_krig2BME{i,j}{k,l}(:,2)-MSE_krig2BME{i,j}{k,l}(:,1))./(MSE_krig2BME{i,j}{k,l}(:,1)), ...
                    colorz{l});
                if k==1&l==1, strz{n} = 'baseline';
                else strz{n} = sprintf('>=%0.3f, >=%0.3f',cutoffs{i}(k),cutoffs{j}(l)); end
                n = n + 1;
            end
            
        end
        
        title(sprintf('%% change MSE from kriging to BME by cutoffs: %s, %s', ...
            allvariablesstr{i},allvariablesstr{j}));
        xlabel('km');
        legend(strz,'location','best');

        % save figure
        set(gcf,'Position',[0 0 800 600]);
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/MSE_exclus_cont_double_cutoffs_%s_%s.png', ...
            allvariablesstr{i},allvariablesstr{j}));
        
    end
end
close all;

end