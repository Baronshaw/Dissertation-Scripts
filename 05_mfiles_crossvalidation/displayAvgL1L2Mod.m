function [] = displayAvgL1L2Mod()
% this function will display overall staistics, maps, and time series for
% the more traditional model performance statistics for the RAMP method

% load data
load('../matfiles/avgL1L2Mod.mat');

% calculate statistics overall
ME = mean(sumMod-sumLam1); 
NME = sum(sumMod-sumLam1)/sum(sumLam1); 
SE = std(sumMod-sumLam1); 
R = corr(sumLam1,sumMod,'type','Pearson'); 
R2 = R.^2;
MAE = (1./length(sumMod)).*sum(abs(sumMod-sumLam1));
NMAE = sum(abs(sumMod-sumLam1))./sum(sumMod); 
RMSE = sqrt( mean((sumMod-sumLam1).^2) );
NRMSE = sqrt( sum((sumMod-sumLam1).^2) )./sqrt(sum(sumLam1));
L2 = mean(sumLam2); % not sure about this one
SME = mean(MESCurve);
SNME = mean(NMESCurve);
SSE = mean(SESCurve);
SR = mean(RSCurve);
SR2 = mean(R2SCurve);
SMAE = mean(MAESCurve);
SNMAE = mean(NMAESCurve);
SRMSE = mean(RMSESCurve);
SNRMSE = mean(NRMSESCurve);
SL2 = mean(L2SCurve);
allstats = [ ME SME ; NME SNME ; SE SSE ; R SR ; R2 SR2 ; MAE SMAE ; ...
    NMAE SNMAE ; RMSE SRMSE ; NRMSE SNRMSE ; L2 SL2 ];

% average across days
unidays = unique(allLocs(:,3));
dME = NaN*ones(size(unidays,1),1); dNME = NaN*ones(size(unidays,1),1); 
dSE = NaN*ones(size(unidays,1),1); dR = NaN*ones(size(unidays,1),1);
dR2 = NaN*ones(size(unidays,1),1); dMAE = NaN*ones(size(unidays,1),1); 
dNMAE = NaN*ones(size(unidays,1),1); dRMSE = NaN*ones(size(unidays,1),1); 
dNRMSE = NaN*ones(size(unidays,1),1); dL2 = NaN*ones(size(unidays,1),1);
SdME = NaN*ones(size(unidays,1),1); SdNME = NaN*ones(size(unidays,1),1); 
SdSE = NaN*ones(size(unidays,1),1); SdR = NaN*ones(size(unidays,1),1);
SdR2 = NaN*ones(size(unidays,1),1); SdMAE = NaN*ones(size(unidays,1),1); 
SdNMAE = NaN*ones(size(unidays,1),1); SdRMSE = NaN*ones(size(unidays,1),1); 
SdNRMSE = NaN*ones(size(unidays,1),1); SdL2 = NaN*ones(size(unidays,1),1);
for i = 1:length(unidays)
    disp(i);
    idx = unidays(i) == allLocs(:,3);
    dME(i,1) = mean(sumMod(idx)-sumLam1(idx)); 
    dNME(i,1) = sum(sumMod(idx)-sumLam1(idx))/sum(sumLam1(idx)); 
    dSE(i,1) = std(sumMod(idx)-sumLam1(idx)); 
    dR(i,1) = corr(sumLam1(idx),sumMod(idx),'type','Pearson'); 
    dR2(i,1) = dR(i).^2;
    dMAE(i,1) = (1./length(sumMod(idx))).*sum(abs(sumMod(idx)-sumLam1(idx)));
    dNMAE(i,1) = sum(abs(sumMod(idx)-sumLam1(idx)))./sum(sumMod(idx)); 
    dRMSE(i,1) = sqrt( mean((sumMod(idx)-sumLam1(idx)).^2) );
    dNRMSE(i,1) = sqrt( sum((sumMod(idx)-sumLam1(idx)).^2) )./sqrt(sum(sumLam1(idx)));
    dL2(i,1) = mean(sumLam2(idx)); % not sure about this one
    SdME(i,1) = mean(MESCurve(idx));
    SdNME(i,1) = mean(NMESCurve(idx));
    SdSE(i,1) = mean(SESCurve(idx));
    SdR(i,1) = mean(RSCurve(idx));
    SdR2(i,1) = mean(R2SCurve(idx));
    SdMAE(i,1) = mean(MAESCurve(idx));
    SdNMAE(i,1) = mean(NMAESCurve(idx));
    SdRMSE(i,1) = mean(RMSESCurve(idx));
    SdNRMSE(i,1) = mean(NRMSESCurve(idx));
    SdL2(i,1) = mean(L2SCurve(idx)); 
    
end

% average across locations
unilocs = unique(allLocs(:,1:2),'rows');
lME = NaN*ones(size(unilocs,1),1); lNME = NaN*ones(size(unilocs,1),1); 
lSE = NaN*ones(size(unilocs,1),1); lR = NaN*ones(size(unilocs,1),1);
lR2 = NaN*ones(size(unilocs,1),1); lMAE = NaN*ones(size(unilocs,1),1); 
lNMAE = NaN*ones(size(unilocs,1),1); lRMSE = NaN*ones(size(unilocs,1),1); 
lNRMSE = NaN*ones(size(unilocs,1),1); lL2 = NaN*ones(size(unilocs,1),1);
SlME = NaN*ones(size(unilocs,1),1); SlNME = NaN*ones(size(unilocs,1),1); 
SlSE = NaN*ones(size(unilocs,1),1); SlR = NaN*ones(size(unilocs,1),1);
SlR2 = NaN*ones(size(unilocs,1),1); SlMAE = NaN*ones(size(unilocs,1),1); 
SlNMAE = NaN*ones(size(unilocs,1),1); SlRMSE = NaN*ones(size(unilocs,1),1); 
SlNRMSE = NaN*ones(size(unilocs,1),1); SlL2 = NaN*ones(size(unilocs,1),1);
for i = 1:length(unilocs)
    disp(i);
    idx = unilocs(i,1) == allLocs(:,1) & unilocs(i,2) == allLocs(:,2);
    lME(i,1) = mean(sumMod(idx)-sumLam1(idx)); 
    lNME(i,1) = sum(sumMod(idx)-sumLam1(idx))/sum(sumLam1(idx)); 
    lSE(i,1) = std(sumMod(idx)-sumLam1(idx)); 
    lR(i,1) = corr(sumLam1(idx),sumMod(idx),'type','Pearson'); 
    lR2(i,1) = lR(i).^2;
    lMAE(i,1) = (1./length(sumMod(idx))).*sum(abs(sumMod(idx)-sumLam1(idx)));
    lNMAE(i,1) = sum(abs(sumMod(idx)-sumLam1(idx)))./sum(sumMod(idx)); 
    lRMSE(i,1) = sqrt( mean((sumMod(idx)-sumLam1(idx)).^2) );
    lNRMSE(i,1) = sqrt( sum((sumMod(idx)-sumLam1(idx)).^2) )./sqrt(sum(sumLam1(idx)));
    lL2(i,1) = mean(sumLam2(idx)); % not sure about this one
    SlME(i,1) = mean(MESCurve(idx));
    SlNME(i,1) = mean(NMESCurve(idx));
    SlSE(i,1) = mean(SESCurve(idx));
    SlR(i,1) = mean(RSCurve(idx));
    SlR2(i,1) = mean(R2SCurve(idx));
    SlMAE(i,1) = mean(MAESCurve(idx));
    SlNMAE(i,1) = mean(NMAESCurve(idx));
    SlRMSE(i,1) = mean(RMSESCurve(idx));
    SlNRMSE(i,1) = mean(NRMSESCurve(idx));
    SlL2(i,1) = mean(L2SCurve(idx));
    
end

% saving all variables
save('../matfiles/dispL1L2Mod.mat','sumLam1','sumMod','sumLam2','allLocs', ...
    'MESCurve','NMESCurve','SESCurve','RSCurve','R2SCurve','MAESCurve', ...
    'NMAESCurve','RMSESCurve','NRMSESCurve','L2SCurve', ...
    'ME','NME','SE','R','R2','MAE','NMAE','RMSE','NRMSE','L2', ... 
    'SME','SNME','SSE','SR','SR2','SMAE','SNMAE','SRMSE','SNRMSE','SL2','allstats', ...
    'dME','dNME','dSE','dR','dR2','dMAE','dNMAE','dRMSE','dNRMSE','dL2', ...
    'SdME','SdNME','SdSE','SdR','SdR2','SdMAE','SdNMAE','SdRMSE','SdNRMSE','SdL2', ... 
    'lME','lNME','lSE','lR','lR2','lMAE','lNMAE','lRMSE','lNRMSE','lL2', ...
    'SlME','SlNME','SlSE','SlR','SlR2','SlMAE','SlNMAE','SlRMSE','SlNRMSE','SlL2');

% display time series
statstr = { 'ME' ; 'NME' ; 'SE' ; 'R' ; 'R2' ; 'MAE' ; 'NMAE' ; 'RMSE' ; 'NRMSE' ; 'L2' };
Sstatstr = { 'SME' ; 'SNME' ; 'SSE' ; 'SR' ; 'SR2' ; 'SMAE' ; 'SNMAE' ; 'SRMSE' ; 'SNRMSE' ; 'SL2' };
toplot = { dME ; dNME ; dSE ; dR ; dR2 ; dMAE ; dNMAE ; dRMSE ; dNRMSE ; dL2 };
Stoplot = { SdME ; SdNME ; SdSE ; SdR ; SdR2 ; SdMAE ; SdNMAE ; SdRMSE ; SdNRMSE ; SdL2 };
for i = 1:length(statstr)
    
    % plot figure
    figure; hold on;
    plot(unidays,toplot{i},'r.',unidays,Stoplot{i},'b.');
    legend('stats','Scurvestat','Location','Best');
    title(sprintf('%s across days in 2001',statstr{i}));
    xlabel('days in 2001');
    ylabel(sprintf('%s',statstr{i}));
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('TS_RAMP_%s.png',statstr{i})); 
       
end

% display maps
cd ../BMELIB2.0b
startup
cd ../05_mfiles_crossvalidation
toplot = { lME ; lNME ; lSE ; lR ; lR2 ; lMAE ; lNMAE ; lRMSE ; lNRMSE ; lL2 };
Stoplot = { SlME ; SlNME ; SlSE ; SlR ; SlR2 ; SlMAE ; SlNMAE ; SlRMSE ; SlNRMSE ; SlL2 };
for i = 1:length(statstr)
    
    cax = [prctile(toplot{i},5) prctile(toplot{i},95)];
    
    for j = 1:2
        
        if j == 1 
            toshow = toplot{i}; tostr = statstr{i};
        else
            toshow = Stoplot{i}; tostr = Sstatstr{i};
        end
        
        % plot figure
        % figure; hold on;
        % country outline
        cd ../09_mfiles_projections
        load('USAcontiguous.mat');
        plotax = ell2lambertcc([x,y],'whiproj2001');
        cd ../05_mfiles_crossvalidation
        [xg yg Zg] = plotField(unilocs,toshow,[],plotax);      
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
        title(sprintf('map of %s in 2001',tostr));
        
        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpng','-r600',sprintf('map_RAMP_%s.png',tostr)); 
        
    end
    
end

end