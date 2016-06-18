function [] = timeSeriesScurvesMean(soft,constant,gauss)
% this function created 4/17/2014 will create time series of stations and
% then corresponding mean bins of the S curves for each day

% how do these compare with "good stations"?

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end
deltaT = 365; n = 3; minpnts = 150; numbins = 10; negval = 0;

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
 
meansModAll = cell(size(investsID,1),1);
meansObsAll = cell(size(investsID,1),1);

% for each trouble station
for i = 1:size(investsID,1)
 
    idx = ckall(:,1) == investsID(i,1) & ckall(:,2) == investsID(i,2);
    cksub = ckall(idx,3);
    dayzsub = datevec(ckall(idx,3));
    dayzsubdisp = dayzsub(:,1)*10000 + dayzsub(:,2)*100 + dayzsub(:,3);
    meansModAll{i} = NaN*ones(size(dayzsubdisp,1),10);
    meansObsAll{i} = NaN*ones(size(dayzsubdisp,1),10);

    % for each day, create S curve in the grid where the station is located
    for j = 1:sum(idx)

        % loading soft data: soft data NOT from test1
        load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
            dayzsubdisp(j),deltaT,n,minpnts,numbins,negval));

        distz = sqrt( (CTMlocs(:,1)-investsID(i,1)).^2 + (CTMlocs(:,2)-investsID(i,2)).^2 );
        [sorteddistz idxMod] = sort(distz);
 
        % calculate means in each bin
        meansModAll{i}(j,:) = mean_Mod(idxMod(1),:);
        meansObsAll{i}(j,:) = mean_Obs(idxMod(1),:);
       
    end
    
    % plotting
    figure; hold on;
    plot(datenum(dayzsub),meansModAll{i},'.');
    legend('1','2','3','4','5','6','7','8','9','10');
    title(sprintf('Bad Station %d Mean Mod %s %s %s',i,softstr(2:end),constr(2:end),gaussstr(2:end)));
    xlabel('day');
    ylabel('Mod conc');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('TSnum%d_ScurveModMeans_%s%s%s_test1.pdf',i,dayzsubdisp(j),softstr,constr,gaussstr));
    print(gcf,'-painters','-dpng','-r600',sprintf('TSnum%d_ScurveModMeans%d_%s%s%s_test1.png',i,dayzsubdisp(j),softstr,constr,gaussstr));

     % plotting
    figure; hold on;
    plot(datenum(dayzsub),meansObsAll{i},'.');
    legend('1','2','3','4','5','6','7','8','9','10');
    title(sprintf('Bad Station %d Mean Obs %s %s %s',i,softstr(2:end),constr(2:end),gaussstr(2:end)));
    xlabel('day');
    ylabel('Obs conc');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('TSnum%d_ScurveObsMeans_%s%s%s_test1.pdf',i,dayzsubdisp(j),softstr,constr,gaussstr));
    print(gcf,'-painters','-dpng','-r600',sprintf('TSnum%d_ScurveObsMeans%d_%s%s%s_test1.png',i,dayzsubdisp(j),softstr,constr,gaussstr));

end

% loading Xval data
load(sprintf('../matfiles/Xval_10fold%s%s%s_results.mat','_nosoft',constr,gaussstr));
errall_nosoft = zkall - zhall;


end