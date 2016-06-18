function [] = map_MP_cmaq()
% create maps at the monitoring location of these different measures
% overall for the year and by season

% load file
load('matfiles/traditional_performance.mat');

% measure names/values
strnm = { 'number of paired modeled and obs'; 'mean obs value' ; 'mean modeled value' ; ...
    'mean bias' ; 'normalized bias' ; 'normalized mean bias' ; 'fractional bias' ; 'mean error' ; ...
    'normalized error' ; 'normalized mean error' ; 'fractional error' ; 'correlation' ; 'correlation squared' ; ...
    'standard bias' ; 'mean squared bias' ; 'root mean squared bias' ; 'normalized root mean squared bias' ; ...
    'mean bias DIV standard bias' ; 'mean bias squared DIV mean squared bias' ; ...
    'variance of bias DIV mean squared bias' ; 'beta1' ; 'variance of obs' ; 'variance of mod' };
valplot = { numM' ; mObsM' ; mModM' ; mBiasM' ; nBiasM' ; nmBiasM' ; fBiasM' ; ...
    mErrM' ; nErrM' ; nmErrM' ; fErrM' ; RM' ; R2M' ; sBiasM' ; msBiasM' ; ...
    rmsBiasM' ; nrmsBiasM' ; mDsBiasM' ; m2DmsBiasM' ; s2DmsBiasM' ; beta1M' ; vObsM' ; vModM' };

%%% histogram
for i = 1:length(valplot)
    
    figure; hold on;
    hist(valplot{i},100);
    title(sprintf('histogram %s for 2001',strnm{i})); 

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/hist_%s.png',strnm{i}));
    
end
close all;

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
    for j = 1:length(X)
        cd ../09_mfiles_projections
        states = ell2lambertcc([X{j},Y{j}],'whiproj2001');
        cd ../13_mfiles_modelperformance
        plot(states(:,1),states(:,2),'k-');
    end

    % colorplot
    Property={'Marker','MarkerSize','MarkerEdgeColor'};
    Value ={'o',5,[0 0 0]};
    cax = [prctile(valplot{i},5) prctile(valplot{i},90)];
    colorplot(uniPMM,valplot{i},'hot',Property,Value,cax);
    caxis(cax);
    colorbar;
    title(sprintf('%s for 2001',strnm{i})); 

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    set(gca,'YTickLabel',get(gca,'YTick')/1000);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s.png',strnm{i}));

end
close all;

%%% maps by season
% load file
load('matfiles/traditional_performanceII.mat');

% measure names/values
strnmII = { 'winter' ; 'spring' ; 'summer' ; 'fall' };
valplotII = { numMS ; mObsMS ; mModMS ; mBiasMS ; nBiasMS ; nmBiasMS ; fBiasMS ; ...
    mErrMS ; nErrMS ; nmErrMS ; fErrMS ; RMS ; R2MS ; sBiasMS ; msBiasMS ; ...
    rmsBiasMS ; nrmsBiasMS ; mDsBiasMS ; m2DmsBiasMS ; s2DmsBiasMS ; beta1MS ; vObsMS ; vModMS }; 

for i = 1:length(valplot)
    for j = 1:length(strnmII)
      
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
        for k = 1:length(X)
            cd ../09_mfiles_projections
            states = ell2lambertcc([X{k},Y{k}],'whiproj2001');
            cd ../13_mfiles_modelperformance
            plot(states(:,1),states(:,2),'k-');
        end

        % colorplot
        Property={'Marker','MarkerSize','MarkerEdgeColor'};
        Value ={'o',5,[0 0 0]};
        idx = isfinite(valplotII{i}(:,j)); 
        cax = [prctile(valplot{i},5) prctile(valplot{i},90)];
        colorplot(uniPMM(idx,:),valplotII{i}(idx,j),'hot',Property,Value,cax);       
        caxis(cax);
        colorbar;
        title(sprintf('%s in %s for 2001',strnm{i},strnmII{j})); 

        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        set(gca,'XTickLabel',get(gca,'XTick')/1000);
        set(gca,'YTickLabel',get(gca,'YTick')/1000);
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_%s.png',strnm{i},strnmII{j}));
        
        drawnow;
        frameI = getframe(gcf);
        imI = frame2im(frameI);
        [imindI,cmI] = rgb2ind(imI,256);
        outfileI = sprintf('figures/%s_season.gif',strnm{i});
        % On the first loop, create the file. In subsequent loops, append.
        if j == 1
            imwrite(imindI,cmI,outfileI,'gif','DelayTime',1,'loopcount',inf);
        else
            imwrite(imindI,cmI,outfileI,'gif','DelayTime',1,'writemode','append');
        end
    
    end
    close all;
end

end