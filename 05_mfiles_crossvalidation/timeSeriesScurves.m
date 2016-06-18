function [] = timeSeriesScurves(soft,constant,gauss)
% this function created 4/16/2014 will create time series of stations and
% then corresponding S curves for each day

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end
deltaT = 365; n = 3; minpnts = 150; numbins = 10; negval = 0;

% comparing the non soft data
load(sprintf('Xval_10fold%s%s%s_results_test1.mat','_nosoft',constr,gaussstr));
zkall_ns = zkall; zhall_ns = zhall; errall_ns = zkall_ns - zhall_ns;

% loading Xval data
load(sprintf('Xval_10fold%s%s%s_results_test1.mat',softstr,constr,gaussstr));
errall = zkall - zhall;
lb = prctile(errall,5);
ub = prctile(errall,95);

% count station that has the highest amount of large errors or with the
% highest proportion of large errors (quanity 'large')
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = unisID(i,1) == ckall(:,1) & unisID(i,2) == ckall(:,2);
    sumoutb(i) = sum( errall(idx) < lb | errall(idx) > ub );    
end

% finding the particular station of interest
idxloc = unisID(:,1) >= -2.5*10^6 & unisID(:,1) <= 2.5*10^6 & unisID(:,2) >= -2*10^6 & unisID(:,2) <= 1.5*10^6;
investsID = unisID(sumoutb'>100 & idxloc,:); % stations that are below 5th or above 
                                   % 95th %tile at least 100 times in 2001
                                   
% map of where trouble stations are located
figure; hold on;
plot(unisID(:,1),unisID(:,2),'bo');
plot(investsID(:,1),investsID(:,2),'ro');

% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpdf','-r600',sprintf('badlocs_%s%s%s_test1.pdf',softstr,constr,gaussstr));
print(gcf,'-painters','-dpng','-r600',sprintf('badlocs_%s%s%s_test1.png',softstr,constr,gaussstr));

% for each trouble station without soft data
for i = 1:size(investsID,1)
    
    % time series
    idx = investsID(i,1) == ckall(:,1) & investsID(i,2) == ckall(:,2);
    figure; hold on;
    plot(ckall(idx,3),zkall_ns(idx),'b.');
    plot(ckall(idx,3),zhall_ns(idx),'r.');
    plot(ckall(idx,3),errall_ns(idx),'k.');
    legend('pred','obs','err');
    title(sprintf('No soft - Time Series for trouble station %0.0f km %0.0f km', ...
        investsID(i,1)/1000,investsID(i,2)/1000));
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('badTSnum%d_%s%s%s_test1.pdf',i,'_nosoft',constr,gaussstr));
    print(gcf,'-painters','-dpng','-r600',sprintf('badTSnum%d_%s%s%s_test1.png',i,'_nosoft',constr,gaussstr));

end

% for each trouble station
for i = 1:size(investsID,1)
    
    % time series
    idx = investsID(i,1) == ckall(:,1) & investsID(i,2) == ckall(:,2);
    figure; hold on;
    plot(ckall(idx,3),zkall(idx),'b.');
    plot(ckall(idx,3),zhall(idx),'r.');
    plot(ckall(idx,3),errall(idx),'k.');
    legend('pred','obs','err');
    title(sprintf('Time Series for trouble station %0.0f km %0.0f km', ...
        investsID(i,1)/1000,investsID(i,2)/1000));
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('badTSnum%d_%s%s%s_test1.pdf',i,softstr,constr,gaussstr));
    print(gcf,'-painters','-dpng','-r600',sprintf('badTSnum%d_%s%s%s_test1.png',i,softstr,constr,gaussstr));

    
    cksub = ckall(idx,3);
    dayzsub = datevec(ckall(idx,3));
    dayzsubdisp = dayzsub(:,1)*10000 + dayzsub(:,2)*100 + dayzsub(:,3);

    % for each day, create S curve in the grid where the station is located
    for j = 1:sum(idx)

        % loading soft data: soft data NOT from test1
        load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
            dayzsubdisp(j),deltaT,n,minpnts,numbins,negval));

        distz = sqrt( (CTMlocs(:,1)-investsID(i,1)).^2 + (CTMlocs(:,2)-investsID(i,2)).^2 );
        [sorteddistz idxMod] = sort(distz);
 
        figure; hold on;

        idxshow = investsID(i,1) == ckall(:,1) & investsID(i,2) == ckall(:,2) & ...
            datenum(dayzsub(j,:)) == ckall(:,3);
        idxMod = distz == sorteddistz(1);
        plot(valMod(idxMod,:),valObs(idxMod,:),'b.');
        plot([get(gca,'XLim')],[get(gca,'XLim')],'k--');
        xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
        ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
        title(sprintf('PM_{2.5} (\\mug/m^3) on %d, \\DeltaT=%d days, n=%d\nx=%0.0f km, y=%0.0f km\npred=%0.2f, obs=%0.2f, err=%0.2f',...        
            dayzsubdisp(j),deltaT,n,CTMlocs(idxMod,1)./1000,CTMlocs(idxMod,2)./1000,...
            zkall(idxshow),zhall(idxshow),errall(idxshow)));

        % add bins
        for k = 1:numbins+1
            plot([perctile_data(idxMod,k) perctile_data(idxMod,k)],get(gca,'YLim'),'r-');
        end

        % calculate means in each bin
        plot(mean_Mod(idxMod,:),mean_Obs(idxMod,:),'ro','MarkerFaceColor','r');

        % linear interpolation in each bin
        xlim = get(gca,'XLim');
        ylim = get(gca,'YLim');
        interpx = linspace(xlim(1),xlim(2));
        interpy = interp1(mean_Mod(idxMod,:),mean_Obs(idxMod,:),interpx,'linear','extrap'); 
        plot(interpx,interpy,'r-');

        % add modeled value
        plot([dailyCTMg(idxMod) dailyCTMg(idxMod)],[ylim(1) meanGivMod(idxMod,1)], ...
            'c--','LineWidth',2);
        plot([xlim(1) dailyCTMg(idxMod)],[meanGivMod(idxMod,1) meanGivMod(idxMod,1)], ...
            'c--','LineWidth',2);
        
        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpdf','-r600',sprintf('TSnum%d_Scurve%d_%s%s%s_test1.pdf',i,dayzsubdisp(j),softstr,constr,gaussstr));
        print(gcf,'-painters','-dpng','-r600',sprintf('TSnum%d_Scurve%d_%s%s%s_test1.png',i,dayzsubdisp(j),softstr,constr,gaussstr));
  
        if length(findall(0,'type','figure')) >= 20, close all; end
        
    end

end

end