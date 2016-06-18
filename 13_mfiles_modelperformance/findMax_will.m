function [] = findMax_will(yrz)
% this function will find the day with the max observed data and create
% maps for fixed modeled values up to that max obs value

if nargin < 1, yrz = 2001; end

% load CMAQ paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yrz));
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
dayz = datenum([yr mo da]);
unidayz = unique(dayz);

% getting min/max/mean/median obs for each day
len = length(unidayz);
Min = NaN*ones(len,1); Max = NaN*ones(len,1); 
Mean = NaN*ones(len,1); Median = NaN*ones(len,1);
for i = 1:len
    disp(i);
    idx = unidayz(i) == dayz;
    Min(i) = nanmin(Obs(idx)); Max(i) = nanmax(Obs(idx));
    Mean(i) = nanmean(Obs(idx)); Median(i) = nanmedian(Obs(idx));
end
figure; plot(unidayz,Min,unidayz,Max,unidayz,Mean,unidayz,Median);
disp(datestr(unidayz(Max==max(Max))));
unicoord = unique(coordObs,'rows');
figure; plot(unicoord(:,1),unicoord(:,2),'bo',coordObs(Obs==max(Max),1),coordObs(Obs==max(Max),2),'ro')

% save figure 
figure; hold on;
plot(unidayz,Mean,'b.-');
title('Mean PM for each day in 2001');
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd mmm'));
print(gcf,'-painters','-dpng','-r600','figures/MeanObs_2001.png');
           
% save results
save('matfiles/statDayObs.mat','unidayz','Min','Max','Mean','Median');

end