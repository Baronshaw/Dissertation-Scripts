function [] = testValidation()
% this function will test the validation results obtained

if nargin < 1, soft = 0; end % soft data or not
if nargin < 2, constant = 1; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 0, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

load(sprintf('Xval_10fold%s%s%s_results.mat',softstr,constr,gaussstr));

% %%% Why is yearly r2 < daily r2?
% 
% % figure of MSE by station across the US
% uniloc = unique(ckall(:,1:2),'rows');
% for i = 1:length(uniloc)
%     disp(i);
%     idx = uniloc(i,1) == ckall(:,1) & uniloc(i,2) == ckall(:,2);
%     msetest(i) = mean( (zkall(idx)-zhall(idx)).^2 );
%     
% end
% figure; hold on;
% Property={'Marker','MarkerSize','MarkerEdgeColor'};
% Value ={'o',5,[0 0 0]};
% cax = prctile(msetest,[5 95]);
% colorplot(uniloc,msetest','hot',Property,Value,cax);
% caxis(cax);
% colorbar;
% title('daily MSE by location');
% print(gcf,'-painters','-dpdf','-r600','daily_MSE_by_location.pdf');
% 
% % figure of r2 by station across the US
% uniloc = unique(ckall(:,1:2),'rows');
% for i = 1:length(uniloc)
%     disp(i);
%     idx = uniloc(i,1) == ckall(:,1) & uniloc(i,2) == ckall(:,2);
%     r2test(i) = (corr(zkall(idx),zhall(idx),'type','Pearson')).^2;
%     
% end
% figure; hold on;
% Property={'Marker','MarkerSize','MarkerEdgeColor'};
% Value ={'o',5,[0 0 0]};
% cax = [0 1];
% colorplot(uniloc,r2test','hot',Property,Value,cax);
% caxis(cax);
% colorbar;
% title('daily r2 by location');
% print(gcf,'-painters','-dpdf','-r600','daily_r2_by_location.pdf');
% 
% % figure of yearly MSE by station across the US
% uniloc2 = unique(ckyrs(:,1:2),'rows');
% for i = 1:length(uniloc2)
%     disp(i);
%     idx = uniloc2(i,1) == ckyrs(:,1) & uniloc2(i,2) == ckyrs(:,2);
%     msetest2(i) = mean( (zkyrs(idx)-zhyrs(idx)).^2 );
%     
% end
% figure; hold on;
% Property={'Marker','MarkerSize','MarkerEdgeColor'};
% Value ={'o',5,[0 0 0]};
% cax = prctile(msetest2,[5 95]);
% colorplot(uniloc2,msetest2','hot',Property,Value,cax);
% caxis(cax);
% colorbar;
% title('yearly MSE by location');
% print(gcf,'-painters','-dpdf','-r600','yearly_MSE_by_location.pdf');
% 
% % figure of yearly r2 by station across the US
% uniloc2 = unique(ckyrs(:,1:2),'rows');
% for i = 1:length(uniloc)
%     disp(i);
%     idx = uniloc2(i,1) == ckyrs(:,1) & uniloc2(i,2) == ckyrs(:,2);
%     r2test2(i) = (corr(zkyrs(idx),zhyrs(idx),'type','Pearson')).^2;
%     
% end
% figure; hold on;
% Property={'Marker','MarkerSize','MarkerEdgeColor'};
% Value ={'o',5,[0 0 0]};
% cax = [0 1];
% colorplot(uniloc,r2test2','hot',Property,Value,cax);
% caxis(cax);
% colorbar;
% title('yearly r2 by location');
% print(gcf,'-painters','-dpdf','-r600','yearly_r2_by_location.pdf');
% 
% %%% Why does r2 change by year for daily?
% ckallyrs = datevec(ckall(:,3)); ckallyrs = ckallyrs(:,1);
% uniyrs = 2006:2008;
% for i = 1:length(uniyrs)    
%     figure; hold on;
%     idx = uniyrs(i) == ckallyrs;
%     plot(ckall(idx,3),zhall(idx),'b.');
%     plot(ckall(idx,3),zkall(idx),'r.');
%     legend('obs','est');
%     title(sprintf('daily obs and est for %d with R2=%0.3f',uniyrs(i),r2_year(i+7)));  
% end

% figure of MSE by station/year across the US
ckallyrs = datevec(ckall(:,3)); ckallyrs = ckallyrs(:,1);
uniyrs = 2006:2008;
for i = 1:length(uniyrs)
    
    idx = ckallyrs == uniyrs(i);
    cksub = ckall(idx,:); zhsub = zhall(idx); zksub = zkall(idx);
    uniloc = unique(cksub(:,1:2),'rows');
    msebyyeartest = [];
    for j = 1:length(uniloc)
        disp(j);
        idx = uniloc(j,1) == cksub(:,1) & uniloc(j,2) == cksub(:,2);
        msebyyeartest(j) = mean( (zksub(idx)-zhsub(idx)).^2 );

    end
    figure; hold on;
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = prctile(msebyyeartest,[5 95]);
    colorplot(uniloc,msebyyeartest','hot',Property,Value,cax);
    colorbar;
    caxis(cax);
    title('daily MSE by location');
    %print(gcf,'-painters','-dpdf','-r600','daily_MSE_by_location.pdf');
    
end

uniyrs = 2006:2008;
for i = 1:length(uniyrs)
    figure; hold on;
    idx = uniyrs(i) == ckyrs(:,3);
    plot(ckyrs(idx,3),zhyrs(idx),'b.');
    plot(ckyrs(idx,3),zkyrs(idx),'r.');
    legend('obs','est');
    title(sprintf('yearly obs and est for %d with R2=%0.3f',uniyrs(i),r2_yearyrs(i+7)));
end

%%% Why does r2 change by region?
uniEPA = unique(EPAall);
for i = 1:length(uniEPA)-1
    figure; hold on;
    idx = uniEPA(i) == EPAall;
    plot(ckall(idx,3),zhall(idx),'b.');
    plot(ckall(idx,3),zkall(idx),'r.');
    legend('obs','est');
    title(sprintf('daily obs and est for region %d with R2=%0.3f',uniEPA(i),r2_EPA(i)));
end

uniEPA = unique(EPAyrs);
for i = 1:length(uniEPA)-1
    figure; hold on;
    idx = uniEPA(i) == EPAyrs;
    plot(ckyrs(idx,3),zhyrs(idx),'b.');
    plot(ckyrs(idx,3),zkyrs(idx),'r.');
    legend('obs','est');
    title(sprintf('yearly obs and est for region %d with R2=%0.3f',uniEPA(i),r2_EPA(i)));
end

% look at individual stations
EPA_test = [1 2 3 8 10];
for i = 1:length(EPA_test)
    
    rand('seed',0);
    idx = EPAall == EPA_test(i);
    unicksub = unique(ckall(idx,1:2),'rows');
    temp = randsample(length(unicksub),5);
    cksub = ckall(idx,:); zksub = zkall(idx); zhsub = zhall(idx);
    
    idx = EPAyrs == EPA_test(i);
    cksubyrs = ckyrs(idx,:); zksubyrs = zkyrs(idx); zhsubyrs = zhyrs(idx);
    
    % daily
    for j = 1:5
        idx = unicksub(temp(j),1) == cksub(:,1) & unicksub(temp(j),2) == cksub(:,2);
        figure; hold on;
        plot(cksub(idx,3),zhsub(idx),'b.-');
        plot(cksub(idx,3),zksub(idx),'r.-');
        legend('obs','est');
        title(sprintf('daily loc= %0.3f %03.f region %d',unicksub(temp(j),1),unicksub(temp(j),2),EPA_test(i)));
    end
    
    % yearly
    for j = 1:5
        idx = unicksub(temp(j),1) == cksubyrs(:,1) & unicksub(temp(j),2) == cksubyrs(:,2);
        figure; hold on;
        plot(cksubyrs(idx,3),zhsubyrs(idx),'b.-');
        plot(cksubyrs(idx,3),zksubyrs(idx),'r.-');
        legend('obs','est');
        title(sprintf('yearly loc= %0.3f %03.f region %d',unicksub(temp(j),1),unicksub(temp(j),2),EPA_test(i)));
    end
    
end

save(sprintf('testing_Xval_10fold%s%s%s_results.mat',softstr,constr,gaussstr),...
    'msetest','msetest2','skew','skew2');

end