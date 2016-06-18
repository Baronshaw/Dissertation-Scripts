function [] = compareTradModPer2Xval()
% this function will compare the traditional model performance with the
% performance of the crossvalidation statistics by station

% load traditional model performance results
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
coordObs = round(coordObs);
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
yrmodanum = datenum(yr,mo,da);
unilocs = unique(coordObs(:,1:2),'rows');
load('../matfiles/dispTradModPerform.mat');

% load cross validation results
load('../matfiles/Xvalresultsbystation.mat');

% pairing cMS with unilocs
indexunilocs = NaN*ones(length(unilocs),1);
for i = 1:length(unilocs)
    dists = sqrt( (unilocs(i,1)-cMS(:,1)).^2 + (unilocs(i,2)-cMS(:,2)).^2 );
    a = find(dists==min(dists));
    if min(dists) < 100 
        indexunilocs(i) = a; 
    end
end

% loop through each statistic
statstr = { 'ME' ; 'R2' ; 'MAE' ; 'RMSE' };
Ttoplot = { TlME ; TlR2 ; TlMAE ; TlRMSE };
Xtoplot = { MEs ; r2s ; MAEs ; RMSEs };
for i = 1:length(statstr)

    figure; hold on;
    plot(Ttoplot{i},Xtoplot{i}(indexunilocs),'b.');
    title(sprintf('Model vs Xval stats by station 2001 for %s',statstr{i}));
    xlabel('Model Stats');
    ylabel('Xval Stat');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('ByStation2001_ModVXval_%s.png',statstr{i})); 
    
end

end