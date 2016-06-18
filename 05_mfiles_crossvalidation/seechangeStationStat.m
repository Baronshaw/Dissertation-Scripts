function [] = seechangeStationStat()
% this function created 4/29/2014 will create histograms of the statistics
% by station for both methods (obs and obs + CTM) and look at % change
% going from one method to the other

%%% OBS %%%

soft = 0;  
constant = 0; 
gauss = 1;
if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end
load(sprintf('Xval_10fold%s%s%s_results_test1.mat',softstr,constr,gaussstr));

% variables to work with
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = ckall(:,1) == unisID(i,1) & ckall(:,2) == unisID(i,2);
    obs(i) = mean(zhall(idx));
end
vars_all{1} = { obs ; RMSE_sID ; MAE_sID ; ME_sID ; r2_sID ; MS_sID ; RMSS_sID ; ...
    MR_sID ; std_est_sID ; std_obs_sID };
str_all = { 'obs' ; 'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'std est' ; 'std ob' };

%%% OBS + CTM %%%

soft = 1;  
constant = 0; 
gauss = 1;
if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end
load(sprintf('Xval_10fold%s%s%s_results_test1.mat',softstr,constr,gaussstr));

% variables to work with
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = ckall(:,1) == unisID(i,1) & ckall(:,2) == unisID(i,2);
    obs(i) = mean(zhall(idx));
end
vars_all{2} = { obs ; RMSE_sID ; MAE_sID ; ME_sID ; r2_sID ; MS_sID ; RMSS_sID ; ...
    MR_sID ; std_est_sID ; std_obs_sID };
str_all = { 'obs' ; 'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'std est' ; 'std ob' };

% loop through histograms
for i = 1:length(vars_all{1})
    for j = 1:2
    
        figure; hold on;

        % colorplot
        hist(vars_all{j}{i},100);
        if j == 1
            title(sprintf('%s %s %s %s','nosoft',constr(2:end),gaussstr(2:end),str_all{i}));
        else
            title(sprintf('%s %s %s %s','soft',constr(2:end),gaussstr(2:end),str_all{i})); 
        end
        [y1 x1] = hist(vars_all{1}{i},100);
        [y2 x2] = hist(vars_all{2}{i},100);
        xlim([min([x1,x2]) max([x1,x2])]);
        ylim([min([y1,y2]) max([y1,y2])]);

        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        if j == 1
            print(gcf,'-painters','-dpdf','-r600',sprintf('hist%s_%s%s%s_test1.pdf', ...
                str_all{i},'_nosoft',constr,gaussstr));
        else
            print(gcf,'-painters','-dpdf','-r600',sprintf('hist%s_%s%s%s_test1.pdf', ...
                str_all{i},'_soft',constr,gaussstr));
        end
    
    end
end

end