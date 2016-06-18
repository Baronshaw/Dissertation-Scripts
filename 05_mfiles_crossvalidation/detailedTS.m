function [] = detailedTS(soft,constant,gauss,dispnum)
% this function created 4/29/2014 will display detailed time series of Xval
% results. This will include nhmax, nsmax, pred, and obs

if nargin < 1, soft = 0; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not
if nargin < 4, dispnum = 5; end % number of time series to randomly display

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% gathering all the data
allhardall = [];
alllambda1all = [];
alllambda2all = [];
for i = 1:10   
    % loading results
    load(sprintf('../matfiles/Xval_10fold_fold%d%s%s%s_test1.mat',i,softstr,constr,gaussstr));
    zkall{i,1} = zk_madd;
    zhall{i,1} = zh_Xval;
    ckall{i,1} = ck;
    vkall{i,1} = vk;  
    foldall{i,1} = i*ones(length(zkall{i,1}),1);
    allhardall = [ allhardall ; allhard ];
    alllambda1all = [ alllambda1all ; alllambda1 ];
    alllambda2all = [ alllambda2all ; alllambda2 ];
end
zkall = cell2mat(zkall);
idx = ~isnan(zkall); zkall = zkall(idx);
zhall = cell2mat(zhall); zhall = zhall(idx);
ckall = cell2mat(ckall); ckall = ckall(idx,:);
vkall = cell2mat(vkall); vkall = vkall(idx);
foldall = cell2mat(foldall); foldall = foldall(idx);
% statistics results
load(sprintf('Xval_10fold%s%s%s_results_test1.mat',softstr,constr,gaussstr));

% variables to work with
unisID = unique(ckall(:,1:2),'rows');
% looping through time series
unisID = unique(ckall(:,1:2),'rows');
rand('seed',0);
picksID = randsample(length(unisID),dispnum);

% map of locations
% country outline
figure; hold on;
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
for j = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
    cd ../04_mfiles_softdata
    plot(states(:,1),states(:,2),'k-');
end
% plot
plot(unisID(picksID(1:dispnum),1),unisID(picksID(1:dispnum),2),'bo');
title('time series locations'); 
% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
set(gca,'XTickLabel',get(gca,'XTick')/1000);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
print(gcf,'-painters','-dpdf','-r600',sprintf('map%d_%s%s%s_test1.pdf', ...
    dispnum,softstr,constr,gaussstr));
    
for i = 1:dispnum   
    
    figure; hold on;
    idx = unisID(picksID(i),1) == ckall(:,1) & unisID(picksID(i),2) == ckall(:,2);
    plot(ckall(idx,3),zkall(idx),'b.'); % predicted
    plot(ckall(idx,3),zhall(idx),'r.'); % observed
    
    cksub = ckall(idx,3);
    for j = 1:length(cksub)
        idx2 = idx & cksub(j) == ckall(:,3);
        len = length(allhardall{idx2});
        plot(repmat(cksub(j),len,1),allhardall{idx2},'c.');
    end
    
    if soft == 1 
        % display soft data
        legend('pred','obs','hard','soft');
    else
        legend('pred','obs','hard');
    end

    ylabel('PM2.5 concentration');
    title(sprintf('detailed time series station num %d (%0.0f km,%0.0f km)\nRMSE=%0.2f, MAE=%0.2f, ME=%0.2f, r2=%0.2f\nMS=%0.2f, RMSS=%0.2f, MR=%0.2f, stdest=%0.2f, stdobs=%0.2f', ...
        i,unisID(picksID(i),1)/1000,unisID(picksID(i),2)/1000,RMSE_sID(picksID(i)),MAE_sID(picksID(i)), ...
        ME_sID(picksID(i)),r2_sID(picksID(i)),MS_sID(picksID(i)),RMSS_sID(picksID(i)), ...
        MR_sID(picksID(i)),std_est_sID(picksID(i)),std_obs_sID(picksID(i))));
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',datestr(get(gca,'XTick')));
    print(gcf,'-painters','-dpdf','-r600',sprintf('detailedTS_num%d_%s%s%s_test1.pdf', ...
        i,softstr,constr,gaussstr));
    
end

end