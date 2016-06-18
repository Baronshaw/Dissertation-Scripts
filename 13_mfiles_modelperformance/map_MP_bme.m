function [] = map_MP_bme()
% create maps at the monitoring location of these different measures for
% bme performance

strwhat = { '_nosoft' ; '_soft' ; '_diff' ; '_reldiff' };
strwhat2 = { 'Kriging' ; 'RAMP' ; 'Difference' ; 'Relative Difference' };

% getting color bar range
load(sprintf('matfiles/traditional_performance_bme%s.mat','_soft'));
valplot = { numM' ; mObsM' ; mModM' ; mBiasM' ; nBiasM' ; nmBiasM' ; fBiasM' ; ...
    mErrM' ; nErrM' ; nmErrM' ; fErrM' ; RM' ; R2M' ; sBiasM' ; msBiasM' ; ...
    rmsBiasM' ; nrmsBiasM' ; mDsBiasM' ; m2DmsBiasM' ; s2DmsBiasM' ; beta1M' ; vObsM' ; vModM'; ...
    mBiasM.^2' ; sBiasM.^2' ; ...
    mBias.^2./mObs.^2 ; sBias.^2./mObs.^2 ; msBias.^2./mObs.^2 };
uniPMplot = uniPMM;
for i = 1:length(valplot)
    caxs{i} = [prctile(valplot{i},5) prctile(valplot{i},90)];
end

for l = 1:length(strwhat) 
    
    % load file
    load(sprintf('matfiles/traditional_performance_bme%s.mat',strwhat{l}));

    % measure names/values
    strnm = {'Number of Paired Values';'Mean Observed';'Mean Modeled';'Mean Error'; ...
        'Mean Normalized Error';'Normalized Mean Error';'Fractional Error';'Mean Absolute Error'; ...
        'Mean Normalized Absolute Error';'Normalized Mean Absolute Error';'Fractional Absolute Error'; ...
        'Correlation';'Correlation Squared';'Standard Error';'Mean Squared Error'; ...
        'Root Mean Squared Error';'Normalized Root Mean Squared Error'; ...
        'Mean Error DIV Standard Error';'Mean Error Squared DIV Mean Squared Error'; ...
        'Standard Error Squared DIV Mean Squared Error';'beta1';'Variance of Observed';'Variance of Modeled'; ...
        'Mean Error Squared' ; 'Standard Error Squared' ; ...
        'Mean Error Squared DIV Mean Observed Squared' ; 'Standard Error Squared DIV Mean Observed Squared' ; 'Mean Squared Error DIV Mean Observed Squared'};
    valplot = { numM' ; mObsM' ; mModM' ; mBiasM' ; nBiasM' ; nmBiasM' ; fBiasM' ; ...
        mErrM' ; nErrM' ; nmErrM' ; fErrM' ; RM' ; R2M' ; sBiasM' ; msBiasM' ; ...
        rmsBiasM' ; nrmsBiasM' ; mDsBiasM' ; m2DmsBiasM' ; s2DmsBiasM' ; beta1M' ; vObsM' ; vModM' ; ...
        mBiasM.^2' ; sBiasM.^2' ; ...
        mBiasM.^2'./mObsM.^2' ; sBiasM.^2'./mObsM.^2' ; msBiasM.^2'./mObsM.^2' };
    valplotoverall = { num ; mObs ; mMod ; mBias ; nBias ; nmBias ; fBias ; ...
        mErr ; nErr ; nmErr ; fErr ; R ; R2 ; sBias ; msBias ; ...
        rmsBias ; nrmsBias ; mDsBias ; m2DmsBias ; s2DmsBias ; beta1 ; vObs ; vMod ; ...
        mBias.^2 ; sBias.^2 ; ...
        mBias.^2./mObs.^2 ; sBias.^2./mObs.^2 ; msBias.^2./mObs.^2 } ;

%     %%% histogram
%     for i = 1:length(valplot)
% 
%         figure; hold on;
%         hist(valplot{i},100);
%         title(sprintf('%s histogram %s for 2001',strwhat2{l},strnm{i})); 
% 
%         % save figure
%         set(gcf,'Position',[0 0 800 600]);
%         set(gcf,'PaperUnits','inches');    
%         set(gcf,'PaperPosition',[0 0 800 600]./100);
%         set(gcf,'PaperPositionMode','manual');
%         print(gcf,'-painters','-dpng','-r600',sprintf('figures/hist_%s_bme%s.png',strnm{i},strwhat{l}));
% 
%     end
%     close all;

    %%% maps overall
    for i = 1:length(valplot)

        figure; hold on;

        % country outline
        cd ../09_mfiles_projections
        load('USAcontiguous.mat');
        plotax = ell2lambertcc([x,y],'whiproj2001');
        cd ../13_mfiles_modelperformance

        % setting axis
        xlabel('km');
        ylabel('km');
        axis([ -3000000 3000000 -2000000 1500000 ]);

        % overlaying the states
        load('../09_mfiles_projections/USAstates5.mat');
        allstates = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
            {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
        for k = 1:length(allstates)
            cd ../09_mfiles_projections
            states = ell2lambertcc([allstates(k).Lon',allstates(k).Lat'],'whiproj2001');
            cd ../13_mfiles_modelperformance
            plot(states(:,1),states(:,2),'k-');
        end

        % colorplot
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',5,[0 0 0]};   
        if length(unique(valplot{i})) == 1            
            cax = [unique(valplot{i})-1 unique(valplot{i})+1];
            valplot{i}(1) = valplot{i}(1) + 0.01;
        elseif sum(isfinite(valplot{i})) == 0
            cax = [0 1]; 
        elseif caxs{i}(1) == caxs{i}(2)
            cax = [prctile(valplot{i},5) prctile(valplot{i},90)];
        elseif l == 0 | l == 1
            cax = caxs{i};        
        else
            cax = [prctile(valplot{i},5) prctile(valplot{i},90)];
        end
        
        % getting right color scheme
        colorz = redpink;
        if strcmp(strnm{i},'Correlation')==1 | strcmp(strnm{i},'Correlation Squared')==1 
            colorz = flipud(colorz); 
        end
            
        colorplot(uniPMplot,valplot{i},colorz,Property,Value,cax);
        caxis(cax);
        colorbar;
        title(sprintf('%s %s for 2001\noverall 2001 = %0.3f',strwhat2{l},strnm{i},valplotoverall{i})); 

        % save figure
        set(gcf,'Position',[0 0 800 600]);       
        set(gca,'YTickLabel',get(gca,'YTick')/1000);
        set(gca,'XTickLabel',get(gca,'XTick')/1000);
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_bme%s.png',strnm{i},strwhat{l}));

    end
    close all;

end

end