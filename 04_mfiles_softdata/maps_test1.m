function [] = maps_test1(soft,constant,gauss)
% created 4/15/2014 this function will display cross-validation statistics
% results across space in maps for 2001

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% load data
load(sprintf('../05_mfiles_crossvalidation/Xval_LOOCV%s%s%s_results_test1.mat', ...
    '_nosoft',constr,gaussstr));
% variables to work with
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = ckall(:,1) == unisID(i,1) & ckall(:,2) == unisID(i,2);
    obs(i) = mean(zhall(idx));
end
vars_all{1} = { obs ; RMSE_sID ; MAE_sID ; ME_sID ; r2_sID ; MS_sID ; RMSS_sID ; ...
    MR_sID ; std_est_sID ; std_obs_sID };
str_all = { 'obs' ; 'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'std est' ; 'std ob' };


% load data
load(sprintf('../05_mfiles_crossvalidation/Xval_LOOCV%s%s%s_results_test1.mat', ...
    '_soft',constr,gaussstr));
% variables to work with
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = ckall(:,1) == unisID(i,1) & ckall(:,2) == unisID(i,2);
    obs(i) = mean(zhall(idx));
end
vars_all{2} = { obs ; RMSE_sID ; MAE_sID ; ME_sID ; r2_sID ; MS_sID ; RMSS_sID ; ...
    MR_sID ; std_est_sID ; std_obs_sID };
str_all = { 'obs' ; 'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'std est' ; 'std ob' };

for i = 1:length(vars_all{1})
    for j = 1:2
        
        figure; hold on;

        % country outline
        cd ../09_mfiles_projections
        load('USAcontiguous.mat');
        plotax = ell2lambertcc([x,y],'whiproj2001');
        cd ../04_mfiles_softdata

        % setting axis
        xlabel('km');
        ylabel('km');
        axis([ -3000000 3000000 -2000000 1500000 ]);

        % overlaying the states
        load('../09_mfiles_projections/USAstates5.mat');
        for k= 1:length(X)
            cd ../09_mfiles_projections
            states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
            cd ../04_mfiles_softdata
            plot(states(:,1),states(:,2),'k-');
        end

        % colorplot
        temp = [vars_all{1}{i} vars_all{2}{i}];
        cax = [prctile(temp,0) prctile(temp,100)];
        deciles = prctile(temp,0:10:100);
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',5,[0 0 0]};
        colorplot2(unisID,vars_all{j}{i}',flipud(hot(10)),Property,Value,cax,deciles);
                
        % colorbar
        cmap = flipud(hot(10));  
        colormap(cmap);         
        cbh = colorbar;         
        caxis(cax);
        set(cbh,'ytick',linspace(prctile(temp,0),prctile(temp,100),11));
        set(cbh,'yticklabel',arrayfun(@num2str,deciles,'uni',false));

        if j == 1
            title(sprintf('%s %s %s %s','nosoft',constr(2:end),gaussstr(2:end),str_all{i}));
        else 
            title(sprintf('%s %s %s %s','soft',constr(2:end),gaussstr(2:end),str_all{i}));
        end
         
        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        set(gca,'XTickLabel',get(gca,'XTick')/1000);
        set(gca,'YTickLabel',get(gca,'YTick')/1000);
        if j == 1
            print(gcf,'-painters','-dpdf','-r600',sprintf('LOOCV_%s_%s%s%s_test1.pdf', ...
            str_all{i},'_nosoft',constr,gaussstr));
        else
            print(gcf,'-painters','-dpdf','-r600',sprintf('LOOCV_%s_%s%s%s_test1.pdf', ...
            str_all{i},'soft',constr,gaussstr));
        end        
    
    end

end

end