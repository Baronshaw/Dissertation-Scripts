function [] = displayTradModPerform()
% this function will display the measures of traditional model performance
% including all the statistics, maps, and time series. This will be
% compared to the RAMP method at the locations where both exists.

% load paired Obs and CTM data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
coordObs = round(coordObs);
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
yrmodanum = datenum(yr,mo,da);

% statistics overall
TME = mean(Mod-Obs); 
TNME = sum(Mod-Obs)/sum(Obs); 
TSE = std(Mod-Obs); 
TR = corr(Obs,Mod,'type','Pearson'); 
TR2 = TR.^2;
TMAE = (1./length(Mod)).*sum(abs(Mod-Obs));
TNMAE = sum(abs(Mod-Obs))./sum(Mod); 
TRMSE = sqrt( mean((Mod-Obs).^2) );
TNRMSE = sqrt( sum((Mod-Obs).^2) )./sqrt(sum(Obs));
allstats = [ TME ; TNME ; TSE ; TR ; TR2 ; TMAE ; TNMAE ; TRMSE ; TNRMSE ];

% average across days
unidays = unique(yrmodanum);
TdME = NaN*ones(size(unidays,1),1); TdNME = NaN*ones(size(unidays,1),1); 
TdSE = NaN*ones(size(unidays,1),1); TdR = NaN*ones(size(unidays,1),1);
TdR2 = NaN*ones(size(unidays,1),1); TdMAE = NaN*ones(size(unidays,1),1); 
TdNMAE = NaN*ones(size(unidays,1),1); TdRMSE = NaN*ones(size(unidays,1),1); 
TdNRMSE = NaN*ones(size(unidays,1),1); 
for i = 1:length(unidays)
    disp(i);
    idx = unidays(i) == yrmodanum;
    TdME(i,1) = mean(Mod(idx)-Obs(idx)); 
    TdNME(i,1) = sum(Mod(idx)-Obs(idx))/sum(Obs(idx)); 
    TdSE(i,1) = std(Mod(idx)-Obs(idx)); 
    TdR(i,1) = corr(Obs(idx),Mod(idx),'type','Pearson'); 
    TdR2(i,1) = TdR(i).^2;
    TdMAE(i,1) = (1./length(Mod(idx))).*sum(abs(Mod(idx)-Obs(idx)));
    TdNMAE(i,1) = sum(abs(Mod(idx)-Obs(idx)))./sum(Mod(idx)); 
    TdRMSE(i,1) = sqrt( mean((Mod(idx)-Obs(idx)).^2) );
    TdNRMSE(i,1) = sqrt( sum((Mod(idx)-Obs(idx)).^2) )./sqrt(sum(Obs(idx)));
    
end

% average across locations
unilocs = unique(coordObs(:,1:2),'rows');
TlME = NaN*ones(size(unilocs,1),1); TlNME = NaN*ones(size(unilocs,1),1); 
TlSE = NaN*ones(size(unilocs,1),1); TlR = NaN*ones(size(unilocs,1),1);
TlR2 = NaN*ones(size(unilocs,1),1); TlMAE = NaN*ones(size(unilocs,1),1); 
TlNMAE = NaN*ones(size(unilocs,1),1); TlRMSE = NaN*ones(size(unilocs,1),1); 
TlNRMSE = NaN*ones(size(unilocs,1),1);
for i = 1:length(unilocs)
    disp(i);
    idx = unilocs(i,1) == coordObs(:,1) & unilocs(i,2) == coordObs(:,2);
    TlME(i,1) = mean(Mod(idx)-Obs(idx)); 
    TlNME(i,1) = sum(Mod(idx)-Obs(idx))/sum(Obs(idx)); 
    TlSE(i,1) = std(Mod(idx)-Obs(idx)); 
    TlR(i,1) = corr(Obs(idx),Mod(idx),'type','Pearson'); 
    TlR2(i,1) = TlR(i).^2;
    TlMAE(i,1) = (1./length(Mod(idx))).*sum(abs(Mod(idx)-Obs(idx)));
    TlNMAE(i,1) = sum(abs(Mod(idx)-Obs(idx)))./sum(Mod(idx)); 
    TlRMSE(i,1) = sqrt( mean((Mod(idx)-Obs(idx)).^2) );
    TlNRMSE(i,1) = sqrt( sum((Mod(idx)-Obs(idx)).^2) )./sqrt(sum(Obs(idx)));
    
end

% saving all variables
save('../matfiles/dispTradModPerform.mat', ...    
    'TME','TNME','TSE','TR','TR2','TMAE','TNMAE','TRMSE','TNRMSE', ... 
    'TdME','TdNME','TdSE','TdR','TdR2','TdMAE','TdNMAE','TdRMSE','TdNRMSE', ...
    'TlME','TlNME','TlSE','TlR','TlR2','TlMAE','TlNMAE','TlRMSE','TlNRMSE');

% display maps
cd ../BMELIB2.0b
startup
cd ../05_mfiles_crossvalidation
statstr = { 'ME' ; 'NME' ; 'SE' ; 'R' ; 'R2' ; 'MAE' ; 'NMAE' ; 'RMSE' ; 'NRMSE' };
toplot = { TlME ; TlNME ; TlSE ; TlR ; TlR2 ; TlMAE ; TlNMAE ; TlRMSE ; TlNRMSE };
for i = 1:length(statstr)   
    
    figure; hold on;
    toshow = toplot{i}; 
    
    % country outline   
    cd ../09_mfiles_projections
    load('USAcontiguous.mat');
    plotax = ell2lambertcc([x,y],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};       
    cax = [prctile(toshow,5) prctile(toshow,95)];
    colorplot(unilocs,toshow,hot(10),Property,Value,cax);       
    caxis(cax);
    colorbar;
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
    title(sprintf('map of %s in 2001',statstr{i}));

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('map_Trad_%s.png',statstr{i})); 
  
end

% load RAMP results
load('../matfiles/dispL1L2Mod.mat');

% find the modeling locations that correspond to the observed locations
subsumMod = NaN*ones(length(Obs),1);
subsumLam1 = NaN*ones(length(Obs),1);
subMESCurve = NaN*ones(length(Obs),1);
subNMESCurve = NaN*ones(length(Obs),1);
subSESCurve = NaN*ones(length(Obs),1);
subRSCurve = NaN*ones(length(Obs),1);
subR2SCurve = NaN*ones(length(Obs),1);
subMAESCurve = NaN*ones(length(Obs),1);
subNMAESCurve = NaN*ones(length(Obs),1);
subRMSESCurve = NaN*ones(length(Obs),1);
subNRMSESCurve = NaN*ones(length(Obs),1);
for i = 1:length(unilocs)
    disp(i);
    % find the locations of the RAMP
    distsidx = sqrt( (unilocs(i,1)-allLocs(:,1)).^2 + (unilocs(i,2)-allLocs(:,2)).^2 );
    [sorted sortidx] = sort(distsidx); 
    tempsumMod = sumMod(sortidx(1:365)); tempsumLam1 = sumLam1(sortidx(1:365)); 
    tempallLocs = allLocs(sortidx(1:365),:); 
    tempallLocs(:,1:2) = repmat(unilocs(i,:),length(tempallLocs(:,3)),1);
    [lia lib] = ismember(tempallLocs,[coordObs yrmodanum],'rows');
    lib(lib==0) = [];
    subsumMod(lib) = tempsumMod(lia);
    subsumLam1(lib) = tempsumLam1(lia);
    
    tempMESCurve = MESCurve(sortidx(1:365)); tempNMESCurve = NMESCurve(sortidx(1:365));
    tempSESCurve = SESCurve(sortidx(1:365)); tempRSCurve = RSCurve(sortidx(1:365));
    tempR2SCurve = R2SCurve(sortidx(1:365)); tempMAESCurve = MAESCurve(sortidx(1:365));
    tempNMAESCurve = NMAESCurve(sortidx(1:365)); tempRMSESCurve = RMSESCurve(sortidx(1:365));
    tempNRMSESCurve = NRMSESCurve(sortidx(1:365));
    subMESCurve(lib) = tempMESCurve(lia); subNMESCurve(lib) = tempNMESCurve(lia);
    subSESCurve(lib) = tempSESCurve(lia); subRSCurve(lib) = tempRSCurve(lia);
    subR2SCurve(lib) = tempR2SCurve(lia); subMAESCurve(lib) = tempMAESCurve(lia);
    subNMAESCurve(lib) = tempNMAESCurve(lia); subRMSESCurve(lib) = tempRMSESCurve(lia);
    subNRMSESCurve(lib) = tempNRMSESCurve(lia);
    
end

% recalculate overall statistics for RAMP and SRAMP
ME = mean(subsumMod-subsumLam1); 
NME = sum(subsumMod-subsumLam1)/sum(subsumLam1); 
SE = std(subsumMod-subsumLam1); 
R = corr(subsumLam1,subsumMod,'type','Pearson'); 
R2 = R.^2;
MAE = (1./length(subsumMod)).*sum(abs(subsumMod-subsumLam1));
NMAE = sum(abs(subsumMod-subsumLam1))./sum(subsumMod); 
RMSE = sqrt( mean((subsumMod-subsumLam1).^2) );
NRMSE = sqrt( sum((subsumMod-subsumLam1).^2) )./sqrt(sum(subsumLam1));
SME = mean(subMESCurve);
SNME = mean(subNMESCurve);
SSE = mean(subSESCurve);
SR = mean(subRSCurve);
SR2 = mean(subR2SCurve);
SMAE = mean(subMAESCurve);
SNMAE = mean(subNMAESCurve);
SRMSE = mean(subRMSESCurve);
SNRMSE = mean(subNRMSESCurve);
allstats = [ ME SME TME ; NME SNME TNME ; SE SSE TSE ; R SR TR ; R2 SR2 TR2 ; MAE SMAE TMAE ; ...
    NMAE SNMAE TNMAE ; RMSE SRMSE TRMSE ; NRMSE SNRMSE TNRMSE ];

% average across days
unidays = unique(yrmodanum);
dME = NaN*ones(size(unidays,1),1); dNME = NaN*ones(size(unidays,1),1); 
dSE = NaN*ones(size(unidays,1),1); dR = NaN*ones(size(unidays,1),1);
dR2 = NaN*ones(size(unidays,1),1); dMAE = NaN*ones(size(unidays,1),1); 
dNMAE = NaN*ones(size(unidays,1),1); dRMSE = NaN*ones(size(unidays,1),1); 
dNRMSE = NaN*ones(size(unidays,1),1); 
SdME = NaN*ones(size(unidays,1),1); SdNME = NaN*ones(size(unidays,1),1); 
SdSE = NaN*ones(size(unidays,1),1); SdR = NaN*ones(size(unidays,1),1);
SdR2 = NaN*ones(size(unidays,1),1); SdMAE = NaN*ones(size(unidays,1),1); 
SdNMAE = NaN*ones(size(unidays,1),1); SdRMSE = NaN*ones(size(unidays,1),1); 
SdNRMSE = NaN*ones(size(unidays,1),1); 
for i = 1:length(unidays)
    disp(i);
    idx = unidays(i) == yrmodanum;
    dME(i,1) = mean(subsumMod(idx)-subsumLam1(idx)); 
    dNME(i,1) = sum(subsumMod(idx)-subsumLam1(idx))/sum(subsumLam1(idx)); 
    dSE(i,1) = std(subsumMod(idx)-subsumLam1(idx)); 
    dR(i,1) = corr(subsumLam1(idx),subsumMod(idx),'type','Pearson'); 
    dR2(i,1) = dR(i).^2;
    dMAE(i,1) = (1./length(subsumMod(idx))).*sum(abs(subsumMod(idx)-subsumLam1(idx)));
    dNMAE(i,1) = sum(abs(subsumMod(idx)-subsumLam1(idx)))./sum(subsumMod(idx)); 
    dRMSE(i,1) = sqrt( mean((subsumMod(idx)-subsumLam1(idx)).^2) );
    dNRMSE(i,1) = sqrt( sum((subsumMod(idx)-subsumLam1(idx)).^2) )./sqrt(sum(subsumLam1(idx)));
    SdME(i,1) = mean(subMESCurve(idx));
    SdNME(i,1) = mean(subNMESCurve(idx));
    SdSE(i,1) = mean(subSESCurve(idx));
    SdR(i,1) = mean(subRSCurve(idx));
    SdR2(i,1) = mean(subR2SCurve(idx));
    SdMAE(i,1) = mean(subMAESCurve(idx));
    SdNMAE(i,1) = mean(subNMAESCurve(idx));
    SdRMSE(i,1) = mean(subRMSESCurve(idx));
    SdNRMSE(i,1) = mean(subNRMSESCurve(idx)); 
    
end

% display time series (only at collocated locations)
statstr = { 'ME' ; 'NME' ; 'SE' ; 'R' ; 'R2' ; 'MAE' ; 'NMAE' ; 'RMSE' ; 'NRMSE' };
toplot = { dME ; dNME ; dSE ; dR ; dR2 ; dMAE ; dNMAE ; dRMSE ; dNRMSE };
Stoplot = { SdME ; SdNME ; SdSE ; SdR ; SdR2 ; SdMAE ; SdNMAE ; SdRMSE ; SdNRMSE };
Ttoplot = { TdME ; TdNME ; TdSE ; TdR ; TdR2 ; TdMAE ; TdNMAE ; TdRMSE ; TdNRMSE };
for i = 1:length(statstr)
    
    % plot figure
    figure; hold on;
    plot(unidays,toplot{i},'r.',unidays,Stoplot{i},'b.',unidays,Ttoplot{i},'g.');
    legend('stats','Scurvestat','Tradstats','Location','Best');
    title(sprintf('%s across days in 2001',statstr{i}));
    xlabel('days in 2001');
    ylabel(sprintf('%s',statstr{i}));
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('TS_TradRAMP_%s.png',statstr{i})); 
       
end

% saving all variables
save('../matfiles/dispTradModPerform.mat','subsumLam1','subsumMod', ...
    'subMESCurve','subNMESCurve','subSESCurve','subRSCurve','subR2SCurve','subMAESCurve', ...
    'subNMAESCurve','subRMSESCurve','subNRMSESCurve','allstats', ...
    'ME','NME','SE','R','R2','MAE','NMAE','RMSE','NRMSE', ... 
    'SME','SNME','SSE','SR','SR2','SMAE','SNMAE','SRMSE','SNRMSE','allstats', ...
    'dME','dNME','dSE','dR','dR2','dMAE','dNMAE','dRMSE','dNRMSE', ...
    'SdME','SdNME','SdSE','SdR','SdR2','SdMAE','SdNMAE','SdRMSE','SdNRMSE', ...
    'TME','TNME','TSE','TR','TR2','TMAE','TNMAE','TRMSE','TNRMSE', ... 
    'TdME','TdNME','TdSE','TdR','TdR2','TdMAE','TdNMAE','TdRMSE','TdNRMSE', ...
    'TlME','TlNME','TlSE','TlR','TlR2','TlMAE','TlNMAE','TlRMSE','TlNRMSE');

end