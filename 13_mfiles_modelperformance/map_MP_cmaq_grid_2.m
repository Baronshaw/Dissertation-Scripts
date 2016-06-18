function [] = map_MP_cmaq_grid_2()
% create maps of traditional model performance on a grid

% load data
load('matfiles/traditional_performance_grid_2.mat')
yr = floor(yrNday./10^4); uniyr = unique(yr);
mo = floor((yrNday - yr*10^4)./10^2);
da = yrNday - yr*10^4 - mo.*10^2;
yrznum = datevec(datenum(uniyr,1:12,1));

% measure names/values
% % CMAQ naming
% strnm = { 'number of paired modeled and obs'; 'mean obs value' ; 'mean modeled value' ; ...
%     'mean bias' ; 'normalized bias' ; 'normalized mean bias' ; 'fractional bias' ; 'mean error' ; ...
%     'normalized error' ; 'normalized mean error' ; 'fractional error' ; 'correlation' ; 'correlation squared' ; ...
%     'standard bias' ; 'mean squared bias' ; 'root mean squared bias' ; 'normalized root mean squared bias' ; ...
%     'mean bias DIV standard bias' ; 'mean bias squared DIV mean squared bias' ; ...
%     'variance of bias DIV mean squared bias' ; 'beta1' ; 'variance of obs' ; 'variance of mod' ; ...
%     'mean bias squared' ; 'standard bias squared' ; ...
%     'mean bias squared DIV mean obs squared' ; 'standard bias squared DIV mean obs squared' ; 'mean squared bias DIV mean obs squared' };
strnm = {'Number of Paired Values';'Mean Observed';'Mean Modeled';'Mean Error'; ...
        'Mean Normalized Error';'Normalized Mean Error';'Fractional Error';'Mean Absolute Error'; ...
        'Mean Normalized Absolute Error';'Normalized Mean Absolute Error';'Fractional Absolute Error'; ...
        'Correlation';'Correlation Squared';'Standard Error';'Mean Squared Error'; ...
        'Root Mean Squared Error';'Normalized Root Mean Squared Error'; ...
        'Mean Error DIV Standard Error';'Mean Error Squared DIV Mean Squared Error'; ...
        'Standard Error Squared DIV Mean Squared Error';'beta1';'Variance of Observed';'Variance of Modeled'; ...
        'Mean Error Squared' ; 'Standard Error Squared' ; ...
        'Mean Error Squared DIV Mean Observed Squared' ; 'Standard Error Squared DIV Mean Observed Squared' ; 'Mean Squared Error DIV Mean Observed Squared'};
valplot = { num ; mObs ; mMod ; mBias ; nBias ; nmBias ; fBias ; ...
    mErr ; nErr ; nmErr ; fErr ; R ; R2 ; sBias ; msBias ; ...
    rmsBias ; nrmsBias ; mDsBias ; m2DmsBias ; s2DmsBias ; beta1 ; vObs ; vMod ; ...
    mBias.^2 ; sBias.^2 ; ...
    mBias.^2./mObs.^2 ; sBias.^2./mObs.^2 ; msBias.^2./mObs.^2 }; 

% loop through each day on the month
for i = 1:length(valplot)
    for j = 1:length(yrznum)
        disp([i j]);
        idx = yrznum(j,1) == yr & yrznum(j,2) == mo & yrznum(j,3) == da;
        figure; hold on;

        load('../09_mfiles_projections/USAcontiguous.mat');
        cd ../09_mfiles_projections
        plotax = ell2lambertcc([x,y],'whiproj2001');
        cd ../13_mfiles_modelperformance

        % plotting 
        lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
        %colorz = colormap(jet);
        colorz = redpink;
        if strcmp(strnm{i},'Correlation')==1 | strcmp(strnm{i},'Correlation Squared')==1, colorz = flipud(colorz); end
        [xg yg Zg] = plotField(CTMlocs,valplot{i}(:,idx),lax,[plotax(:,1) plotax(:,2)],colorz);
        caxis([prctile(valplot{i}(:),5) prctile(valplot{i}(:),90)]);   
        colorbar;
        axis([lax(1)+500 lax(2)-10000 lax(3)+500 lax(4)-10000]);

        % setting axis        
        xlabel('km');
        ylabel('km');

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

        % title 
        % title(sprintf('%s of PM_{2.5} (\\mug/m^3) on %s',strnm{i},datestr(yrznum(j,:))));
        title(sprintf('%s of PM_{2.5} on %s',strnm{i},datestr(yrznum(j,:))));    

        % save figure 
        set(gcf,'Position',[0 0 800 600]);       
        set(gca,'YTickLabel',get(gca,'YTick')/1000);
        set(gca,'XTickLabel',get(gca,'XTick')/1000);
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/Two_%s_%s_grid.png',strnm{i},datestr(yrznum(j,:))));

        drawnow;
        frameI = getframe(gcf);
        imI = frame2im(frameI);
        [imindI,cmI] = rgb2ind(imI,256);
        outfileI = sprintf('figures/Two_%s_grid.gif',strnm{i});
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