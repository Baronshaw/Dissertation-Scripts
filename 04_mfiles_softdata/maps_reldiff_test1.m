function [] = maps_reldiff_test1(soft,constant,gauss)
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
vars_all{1} = { RMSE_sID ; MAE_sID ; ME_sID ; r2_sID ; MS_sID ; RMSS_sID ; ...
    MR_sID ; std_est_sID };
str_all = { 'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'std est' };


% load data
load(sprintf('../05_mfiles_crossvalidation/Xval_LOOCV%s%s%s_results_test1.mat', ...
    '_soft',constr,gaussstr));
% variables to work with
unisID = unique(ckall(:,1:2),'rows');
vars_all{2} = { RMSE_sID ; MAE_sID ; ME_sID ; r2_sID ; MS_sID ; RMSS_sID ; ...
    MR_sID ; std_est_sID };
str_all = { 'RMSE' ; 'MAE' ; 'ME' ; 'r2' ; 'MS' ; 'RMSS' ; 'MR' ; 'std est' };

for i = 1:length(vars_all{1})
           
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
    % new - old / old
    temp = 100.*((vars_all{2}{i}' - vars_all{1}{i}')./vars_all{1}{i}');
    cax = [prctile(temp,0) prctile(temp,100)];
    deciles = prctile(temp,0:10:100);
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    colorplot2(unisID,temp,flipud(hot(10)),Property,Value,cax,deciles);

    % colorbar
    cmap = flipud(hot(10));  
    colormap(cmap);         
    cbh = colorbar;         
    caxis(cax);
    set(cbh,'ytick',linspace(prctile(temp,0),prctile(temp,100),11));
    set(cbh,'yticklabel',arrayfun(@num2str,deciles,'uni',false));

    title(sprintf('LOOCV reldiff %s',str_all{i}));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);

    print(gcf,'-painters','-dpdf','-r600',sprintf('LOOCV_reldiff_%s_test1.pdf',str_all{i}));

end

end