function [diagnostics] = plotScurve_diagnostics2(showspace,dayWHI,daysWHIdisp,location,correct,n,deltaT,numbins,minpnts,modplots,...
    negval)
% this function will create plots for given mean and variances for a given
% location/day

% parameters subject to change
if nargin < 1, showspace = 1; end % showing closest stations or seasons
if nargin < 2, dayWHI = [ 2001 06 15 ]; end
if nargin < 3, daysWHIdisp = 20010615; end  
if nargin < 4, location = rand(1,2); end
if nargin < 5, correct = rand; end
if nargin < 6, n = 3; end % number of closest stations
if nargin < 7, deltaT = 365; end % number of days in the interval
if nargin < 8, numbins = 10; end % number of bins for each plot
if nargin < 9, minpnts = 150; end % min number of point to each plot
if nargin < 10, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 11, negval = 0; end % 0 = there are no negative predicted values

% added later
n = 6;
deltaT = 180;

diagnostics = NaN*ones(3,1);

% loading data
load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat',daysWHIdisp,deltaT,n,minpnts,numbins,negval));
CTMlocs = round(CTMlocs);

figure; hold on;

idxMod = location(1) == CTMlocs(:,1) & location(2) == CTMlocs(:,2);

% adding colors by space or season
colors = {'r','b','c','k','m','g','y',[0.5 0.5 0.5]};
if showspace == 1
    len = length(cMSObsx(idxMod,:));
    npick = sum(~isnan(cMSObsx(idxMod,:)));
    daypick = length(tMEObs)./len;    
    for i = 1:npick
        plot(valMod(idxMod,1+daypick*(i-1):daypick*i), ...
            valObs(idxMod,1+daypick*(i-1):daypick*i),'.','Color',colors{i});
        legendstr{i} = sprintf('station %d',i);
    end
else
    yrz = floor(tMEObs./10^4);
    moz = floor((tMEObs - yrz*10^4)./10^2);
    daz = tMEObs - yrz*10^4 - moz*10^2;
    for i = 1:4
        if i == 1, idxTime = moz == 1 | moz ==2 | moz == 12; if sum(idxTime)>0, legendstr{i} = 'winter'; end
        elseif i == 2, idxTime = moz == 3 | moz == 4 | moz == 5; if sum(idxTime)>0, legendstr{i} = 'spring'; end
        elseif i == 3, idxTime = moz == 6 | moz == 7 | moz == 8; if sum(idxTime)>0, legendstr{i} = 'summer'; end
        else idxTime = moz == 9 | moz == 10 | moz == 11; if sum(idxTime)>0, legendstr{i} = 'fall'; end 
        end
        plot(valMod(idxMod,idxTime),valObs(idxMod,idxTime),'.','Color',colors{i});
    end
    legendstr = legendstr(~cellfun('isempty',legendstr));
end

legend(legendstr);
plot([get(gca,'XLim')],[get(gca,'XLim')],'k--');

xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
title(sprintf('PM_{2.5} (\\mug/m^3) on %d, \\DeltaT=%d days, n=%d\nx=%0.0f km, y=%0.0f km',...        
    daysWHIdisp,deltaT,n,CTMlocs(idxMod,1)./1000,CTMlocs(idxMod,2)./1000));

diagnostics(2) = sum(~isnan(valObs(idxMod,:)));

% add bins
for j = 1:numbins+1
    plot([perctile_data(idxMod,j) perctile_data(idxMod,j)],get(gca,'YLim'),'r-');
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

diagnostics(1) = meanGivMod(idxMod,1);
temp = sort([perctile_data(idxMod,:) meanGivMod(idxMod,1)]);
diagnostics(3) = find(temp==meanGivMod(idxMod,1)) - 1;

% add correct value
plot(dailyCTMg(idxMod),correct,'rs','MarkerSize',7);

end