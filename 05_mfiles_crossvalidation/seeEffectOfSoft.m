function [] = seeEffectOfSoft(forceddist)
% this function looks at years that contain soft data and see if the
% incorporation of soft data reduced or increased bias and by how much

if nargin < 1, forceddist = 300000; end
forced = floor(forceddist./1000);

soft_years = [2001:2002 2005:2007]; 

% gathering all the data
load(sprintf('../matfiles/Xvalforcediso_LOOCV__soft_long_gauss_foriso%dkm.mat',floor(forceddist./1000)));
zkalls = zk_madd; zhalls = zh_Xval; ckalls = ck; vkalls = vk; 

load(sprintf('../matfiles/Xvalforcediso_LOOCV__nosoft_long_gauss_foriso%dkm.mat',floor(forceddist./1000)));
zkallh = zk_madd; zhallh = zh_Xval; ckallh = ck; vkallh = vk;  

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
for i = 1:length(ckalls)
    idx2 = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckalls{i,1} = ch(idx2,:);
end
ckalls = cell2mat(ckalls);
ckalls = ckalls(idx,:);

temp = datevec(ckalls(:,3));
yrs = temp(:,1);

% looping through each soft year
for i = 1:length(soft_years)
    
    idx = soft_years(i) == yrs;
    
    % looking at absoluate error, count # of times it increased/decreased
    allAEs = abs(zkalls(idx) - zhalls(idx));
    allAEh = abs(zkallh(idx) - zhallh(idx));
    
    idxbetter = allAEs<allAEh;
    numbetter = sum(idxbetter)/length(allAEs);
    
    % loop through each station and count when soft was better
    cksub = ckalls(idx,:);
    zksubs = zkalls(idx); zksubh = zkallh(idx);
    zhsubs = zhalls(idx); zhsubh = zhallh(idx);
    unisub = unique(cksub(:,1:2),'rows');
    unicount = [];
    unirdiff = [];
    unirreldiff = [];
    unir2diff = [];
    unir2reldiff = [];
    for j = 1:size(unisub)
        uniidx = unisub(j,1) == cksub(:,1) & unisub(j,2) == cksub(:,2);
        unicount(j) = sum(idxbetter & uniidx)/sum(uniidx);
        unisum(j) = sum(uniidx);
        a = corrcoef(zksubs(uniidx),zhsubs(uniidx)) - ...
            corrcoef(zksubh(uniidx),zhsubh(uniidx));        
        b = ( corrcoef(zksubs(uniidx),zhsubs(uniidx)) - ...
            corrcoef(zksubh(uniidx),zhsubh(uniidx)) )./corrcoef(zksubs(uniidx),zhsubs(uniidx)) ;        
        c = (corrcoef(zksubs(uniidx),zhsubs(uniidx))).^2 - ...
            (corrcoef(zksubh(uniidx),zhsubh(uniidx))).^2;        
        d = ( (corrcoef(zksubs(uniidx),zhsubs(uniidx))).^2 - ...
            (corrcoef(zksubh(uniidx),zhsubh(uniidx))).^2 )./(corrcoef(zksubs(uniidx),zhsubs(uniidx))).^2 ;
        if sum(uniidx) > 1
            unirdiff(j) = a(2);
            unirreldiff(j) = b(2);
            unir2diff(j) = c(2);
            unir2reldiff(j) = d(2);
        else
            unirdiff(j) = a(1);
            unirreldiff(j) = b(1);
            unir2diff(j) = c(1);
            unir2reldiff(j) = d(1);
        end
    end
    % I had to normalize by 'sum(uniidx)' because of the natural bimodality
    % that comes from sampling frequency
    
%     % make a map of normalized counts 
%     figure; hold on;
%     % country outline
%     cd ../09_mfiles_projections
%     load('USAcontiguous.mat');
%     plotax = ell2lambertcc([x,y],'whiproj2001');
%     cd ../05_mfiles_crossvalidation
%     % setting axis
%     xlabel('km');
%     ylabel('km');
%     axis([ -3000000 3000000 -2000000 1500000 ]);
%     % overlaying the states
%     load('../09_mfiles_projections/USAstates5.mat');
%     for j = 1:length(X)
%         cd ../09_mfiles_projections
%         states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
%         cd ../05_mfiles_crossvalidation
%         plot(states(:,1),states(:,2),'k-');
%     end
%     % colorplot
%     Property={'Marker','MarkerSize','MarkerEdgeColor'};
%     Value ={'o',5,[0 0 0]};
%     cax = [0 1];
%     colorplot(unisub,unicount','hot',Property,Value,cax);    
%     caxis(cax);
%     colorbar;
%     title(sprintf('%% for ea loc in %d where soft data reduced abs bias',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('per_softLowBias_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('per_softLowBias_%d_%dkm.png',soft_years(i),forced));
    
    % make a map of bias r  
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    
    cax = [prctile(unirdiff,5) prctile(unirdiff,95)];
    %a = floor(cax(1)*100)/100; b = ceil(cax(2)*100)/100;
    c = linspace(cax(1),cax(2),15);
%     if sign(cax(1)*cax(2)) == -1 % insert 0
%         temp = find(c<0); temp = temp(end);
%         c = [c(1:temp) 0 c(temp+1:end)];
%     else
    
    %c(1) = cax(1); c(end) = cax(2);
    temp1 = cool(length(c));
    temp2 = hot(length(c));
    temp3 = [temp1(c<0,:) ; temp2(c>=0,:)];
    temp3 = temp3(1:end-1,:);
    %temp3 = [temp1(1:find(c==0)-1,:);temp2(find(c==0):end,:)];
    %if sum(c==0)==0, temp3 = temp2; end
    
    colorplot(unisub,unirdiff',temp3,Property,Value,cax);
    caxis(cax);
    %colorbar;
    cbh = colorbar('YGrid','on');
    set(cbh,'ytick',c);
    
    title(sprintf('bias of R %d %d km',soft_years(i),forced));
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('mapRdiff_%d_%dkm.pdf',soft_years(i),forced));
    print(gcf,'-painters','-dpng','-r600',sprintf('mapRdiff_%d_%dkm.png',soft_years(i),forced));
    
%     % make a map of rel bias r  
%     figure; hold on;
%     % country outline
%     cd ../09_mfiles_projections
%     load('USAcontiguous.mat');
%     plotax = ell2lambertcc([x,y],'whiproj2001');
%     cd ../05_mfiles_crossvalidation
%     % setting axis
%     xlabel('km');
%     ylabel('km');
%     axis([ -3000000 3000000 -2000000 1500000 ]);
%     % overlaying the states
%     load('../09_mfiles_projections/USAstates5.mat');
%     for j = 1:length(X)
%         cd ../09_mfiles_projections
%         states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
%         cd ../05_mfiles_crossvalidation
%         plot(states(:,1),states(:,2),'k-');
%     end
%     % colorplot
%     Property={'Marker','MarkerSize','MarkerEdgeColor'};
%     Value ={'o',5,[0 0 0]};
%     cax = [prctile(unirreldiff,5) prctile(unirreldiff,95)];
%     colorplot(unisub,unirreldiff','hot',Property,Value,cax);    
%     caxis(cax);
%     colorbar;
%     title(sprintf('rel bias of R %d',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('mapRreldiff_%d_%d.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('mapRreldiff_%d_%d.png',soft_years(i),forced));
    
    % make a map of bias r2 
    figure; hold on;
    % country outline
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    % setting axis
    xlabel('km');
    ylabel('km');
    axis([ -3000000 3000000 -2000000 1500000 ]);
    % overlaying the states
    load('../09_mfiles_projections/USAstates5.mat');
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        plot(states(:,1),states(:,2),'k-');
    end
    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    
    cax = [prctile(unir2diff,5) prctile(unir2diff,95)];
    %a = floor(cax(1)*100)/100; b = ceil(cax(2)*100)/100;
    c = linspace(cax(1),cax(2),15);
%     if sign(cax(1)*cax(2)) == -1 % insert 0
%         temp = find(c<0); temp = temp(end);
%         c = [c(1:temp) 0 c(temp+1:end)];
%     else
    
    %c(1) = cax(1); c(end) = cax(2);
    temp1 = cool(length(c));
    temp2 = hot(length(c));
    temp3 = [temp1(c<0,:) ; temp2(c>=0,:)];
    temp3 = temp3(1:end-1,:);
    %temp3 = [temp1(1:find(c==0)-1,:);temp2(find(c==0):end,:)];
    %if sum(c==0)==0, temp3 = temp2; end
    
    colorplot(unisub,unir2diff',temp3,Property,Value,cax);
    caxis(cax);
    %colorbar;
    cbh = colorbar('YGrid','on');
    set(cbh,'ytick',c);

    title(sprintf('bias of R2 %d %d km',soft_years(i),forced));
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('mapR2diff_%d_%dkm.pdf',soft_years(i),forced));
    print(gcf,'-painters','-dpng','-r600',sprintf('mapR2diff_%d_%dkm.png',soft_years(i),forced));
    
%     % make a map of rel bias r2  
%     figure; hold on;
%     % country outline
%     cd ../09_mfiles_projections
%     load('USAcontiguous.mat');
%     plotax = ell2lambertcc([x,y],'whiproj2001');
%     cd ../05_mfiles_crossvalidation
%     % setting axis
%     xlabel('km');
%     ylabel('km');
%     axis([ -3000000 3000000 -2000000 1500000 ]);
%     % overlaying the states
%     load('../09_mfiles_projections/USAstates5.mat');
%     for j = 1:length(X)
%         cd ../09_mfiles_projections
%         states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
%         cd ../05_mfiles_crossvalidation
%         plot(states(:,1),states(:,2),'k-');
%     end
%     % colorplot
%     Property={'Marker','MarkerSize','MarkerEdgeColor'};
%     Value ={'o',5,[0 0 0]};
%     cax = [prctile(unir2reldiff,5) prctile(unir2reldiff,95)];
%     colorplot(unisub,unir2reldiff','hot',Property,Value,cax);    
%     caxis(cax);
%     colorbar;
%     title(sprintf('rel bias of R2 %d',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('mapR2reldiff_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('mapR2reldiff_%d_%dkm.png',soft_years(i),forced));
%     
%     % plot station counts
%     figure; hold on;
%     hist(unicount,100);
%     title(sprintf('# for ea loc in %d where soft data reduced abs bias',soft_years(i)));
%     
%     % create the same plots for each day
%     uniday = unique(cksub(:,3)) - datenum(soft_years(i)-1,12,31);
%     cksubtemp = cksub - datenum(soft_years(i)-1,12,31);
%     unicountday = [];
%     unirdiffday = [];
%     unirreldiffday = [];
%     unir2diffday = [];
%     unir2reldiffday = [];
%     for j = 1:size(uniday,1)
%         uniidx = unisub(j,1) == cksub(:,1) & unisub(j,2) == cksub(:,2);
%         unicountday(j) = sum(idxbetter & uniidx)/sum(uniidx);
%         a = corrcoef(zksubs(uniidx),zhsubs(uniidx)) - ...
%             corrcoef(zksubh(uniidx),zhsubh(uniidx));        
%         b = ( corrcoef(zksubs(uniidx),zhsubs(uniidx)) - ...
%             corrcoef(zksubh(uniidx),zhsubh(uniidx)) )./corrcoef(zksubs(uniidx),zhsubs(uniidx)) ;        
%         c = (corrcoef(zksubs(uniidx),zhsubs(uniidx))).^2 - ...
%             (corrcoef(zksubh(uniidx),zhsubh(uniidx))).^2;        
%         d = ( (corrcoef(zksubs(uniidx),zhsubs(uniidx))).^2 - ...
%             (corrcoef(zksubh(uniidx),zhsubh(uniidx))).^2 )./(corrcoef(zksubs(uniidx),zhsubs(uniidx))).^2 ;
%         if sum(uniidx) > 1
%             unirdiffday(j) = a(2);
%             unirreldiffday(j) = b(2);
%             unir2diffday(j) = c(2);
%             unir2reldiffday(j) = d(2);
%         else
%             unirdiffday(j) = a(1);
%             unirreldiffday(j) = b(1);
%             unir2diffday(j) = c(1);
%             unir2reldiffday(j) = d(1);
%         end
%     end
%     
%     % plotting
%     figure; hold on;
%     hist(unicountday,100);
%     title(sprintf('# of ea day in %d where soft data reduced abs bias',soft_years(i)));
%     
%     % make a time series of normalized counts
%     figure; hold on;
%     plot(uniday,unicountday,'b--');
%     title(sprintf('%% of ea day in %d where soft data reduced abs bias',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('per_softLowBiasTS_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('per_softLowBiasTS_%d_%dkm.png',soft_years(i),forced));
%     
%     % make a time series of bias of R
%     figure; hold on;
%     plot(uniday,unirdiffday,'b--');
%     title(sprintf('bias of R %d',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('absRTS_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('absRTS_%d_%dkm.png',soft_years(i),forced));
%     
%     % make a time series of rel bias of R
%     figure; hold on;
%     plot(uniday,unirreldiffday,'b--');
%     title(sprintf('rel bias of R %d',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('relbiasRTS_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('relbiasRTS_%d_%dkm.png',soft_years(i),forced));
%     
%     % make a time series of bias of R2
%     figure; hold on;
%     plot(uniday,unir2diffday,'b--');
%     title(sprintf('bias of R2 %d',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('absR2TS_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('absR2TS_%d_%dkm.png',soft_years(i),forced));
%     
%     % make a time series of rel bias of R2
%     figure; hold on;
%     plot(uniday,unir2reldiffday,'b--');
%     title(sprintf('rel bias of R2 %d',soft_years(i)));
%     % save figure
%     set(gcf,'Position',[0 0 800 600]);
%     set(gcf,'PaperUnits','inches');    
%     set(gcf,'PaperPosition',[0 0 800 600]./100);
%     set(gcf,'PaperPositionMode','manual');
%     print(gcf,'-painters','-dpdf','-r600',sprintf('relbiasR2TS_%d_%dkm.pdf',soft_years(i),forced));
%     print(gcf,'-painters','-dpng','-r600',sprintf('relbiasR2TS_%d_%dkm.png',soft_years(i),forced));
%     
%     % I also need to look at where the abs bias is the AMOUNT of abs bias that's reduced
%     allAEdiff = (allAEs - allAEh)./allAEs;
%     figure; hold on;
%     hist(allAEdiff,100);
       
end

% Conclusion: There is really no pattern to the counts of reduction in abs
% bias - it seems normally distributed

% there could be a pattern in WHERE the counts are high or low, I need to
% make maps and time series to test this

end