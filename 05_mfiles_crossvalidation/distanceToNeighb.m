function [] = distanceToNeighb()
% this function will plot the difference in R2 by distance to closest
% neighbor

% does this make sense to have this on a map, or through a histogram?

% load data
load(sprintf('../matfiles/Xvalforcediso_LOOCV__soft_long_gauss_foriso%dkm.mat',0));
zkallh = zk_madd; zhallh = zh_Xval; ckallh = ck; vkallh = vk; 
load(sprintf('../matfiles/Xvalforcediso_LOOCV__nosoft_long_gauss_foriso%dkm.mat',0));
zkalls = zk_madd; zhalls = zh_Xval; ckalls = ck; vkalls = vk;  
zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
zhallh = cell2mat(zhallh); zhalls = cell2mat(zhalls);
idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx); zhallh = zhallh(idx); 
zkalls = zkalls(idx); zhalls = zhalls(idx); 
% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckalls = cell(length(cMS),1);
for j = 1:length(ckalls)
    idx2 = cMS(j,1) == ch(:,1) & cMS(j,2) == ch(:,2);
    ckalls{j,1} = ch(idx2,:);
end
ckalls = cell2mat(ckalls);
ckalls = ckalls(idx,:);

% get distances to closest neighbor
X = cMS'; Y = X;
D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
[b ix] = sort(D,2);
closestNeib = cMS(ix(:,2),:);
distNeib = b(:,2);

% filling in the distances with all observations
[ai bi] = ismember(ckalls(:,1:2),cMS,'rows');
chdists = distNeib(bi);

% sort distances by 10 percentile
% neighb_tile = prctile(distNeib,0:10:100); % I could do it by unique station
for i = [1 2 4 5 10] % 100 50 25 20 10
    
    neighb_tile = prctile(chdists,0:i:100);

    % calculate diffference in R2 for each percentile
    absDiffR2 = NaN*ones(length(neighb_tile)-1,1);
    for j = 1:length(neighb_tile)-1
        idx = chdists >= neighb_tile(j) & chdists < neighb_tile(j+1);
        a = corrcoef(zkalls(idx),zhalls(idx));
        b = corrcoef(zkallh(idx),zhallh(idx));
        absDiffR2(j) = a(2)^2 - b(2)^2;
    end    
    
    % display results
    % should I put confidence intervals here?
    figure; hold on;
    plot(1:length(absDiffR2),absDiffR2,'bo-');
    title(sprintf('Diff in R2 by %d perctiles',100/i));
    xlabel('percentile');
    ylabel('Difference in R2');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('dist2NeighR2_%dtile.pdf',100/i));
    print(gcf,'-painters','-dpng','-r600',sprintf('dist2NeighR2_%dtile.png',100/i));

end

end