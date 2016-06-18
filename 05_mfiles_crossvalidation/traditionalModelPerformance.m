function [] = traditionalModelPerformance()
% this function will create statistics, maps, and time sereies of 
% traditional model performance measures

% load paired Obs and CTM data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
coordObs = round(coordObs);
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
yrmodanum = datenum(yr,mo,da);

% loop through each station and find R2
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
cMS = round(cMS);
for i = 1:length(cMS)
    idx1 = sqrt( (cMS(i,1)-coordObs(:,1)).^2 + (cMS(i,2)-coordObs(:,2)).^2 );
    idx = idx1 < 36; % inside the grid
    if sum(idx) > 0 
        r2Trad(i,1) = (corr(Obs(idx),Mod(idx),'type','Pearson')).^2;
    else
        r2Trad(i,1) = NaN;
    end
end

% plot and save figure
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
for k = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    plot(states(:,1),states(:,2),'k-');
end
% colorplot
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'o',5,[0 0 0]};       
cax = [prctile(r2Trad,5) prctile(r2Trad,95)];
colorplot(cMS,r2Trad,hot(10),Property,Value,cax);       
caxis(cax);
colorbar;
title('map of R2 Obs v CTM in 2001');
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','map_TradR2_ObsvCTM.png'); 

% condencing lambda1 and lambda2

% load each soft data day
summLam1 = NaN*ones(length(Obs),1);
summBias = NaN*ones(length(Obs),1);
n = 3; 
deltaT = 365; 
numbins = 10; 
minpnts = 150; 
negval = 0;
temp = datevec(datenum(2001,1,1):datenum(2001,12,31));
dayWHI = temp(:,1:3);
dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
for i = 1:length(dayWHIdisps)
    disp(i);
    load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
        dayWHIdisps(i),deltaT,n,minpnts,numbins,negval));
    
    % calculate midpoint
    midpnts = arrayfun(@(j) (perctile_data(j,2:end)-perctile_data(j,1:end-1))./2 ...
        + perctile_data(j,1:end-1), 1:length(perctile_data),'UniformOutput',false);
    midpnts = cell2mat(midpnts');
    
    % summarize bias/lambda1
    biassumize = arrayfun(@(j) (1./length(midpnts(j,1:end))).*sum(midpnts(j,1:end)-mean_Obs(j,1:end)), ...
        1:length(midpnts),'UniformOutput',false);
    biassumize = cell2mat(biassumize');
    lam1sumize = arrayfun(@(j) mean(mean_Obs(j,1:end)),1:length(midpnts),'UniformOutput',false);
    lam1sumize = cell2mat(lam1sumize');
    
    % pairing summarized lambda1's with observed value
    idx = yrmodanum == datenum(dayWHI(i,:));
    X = coordObs(idx,:)';
    Y = CTMlocs';
    D = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
    [D Didx] = sort(D,2);
    lam12Get = lam1sumize(Didx(:,1));
    bias2Get = biassumize(Didx(:,1));
    
    % putting summiarized lambda1's in order
    [ai bi] = ismember([coordObs(idx,:) yrmodanum(idx)],[coordObs yrmodanum],'rows'); 
    summLam1(bi) = lam12Get;
    summBias(bi) = bias2Get;

end

save('../matfiles/TradMethod.mat','coordObs','yrmodanum','summLam1','summBias');

blah = 5;

end