function [] = map_MP_cmaq_gridlamper()
% create maps of traditional model performance on a grid and  look at a few
% different dates by percentile and also by increasing modeled value

% loading BME function
cd ../BMELIB2.0b
startup
cd ../13_mfiles_modelperformance

% load data
load('matfiles/traditional_performance_grid.mat');
valplot = { num ; mObs ; mMod ; mBias ; nBias ; nmBias ; fBias ; ...
    mErr ; nErr ; nmErr ; fErr ; R ; R2 ; sBias ; msBias ; ...
    rmsBias ; nrmsBias ; mDsBias ; m2DmsBias ; s2DmsBias ; beta1 ; vObs ; vMod ; ...
    mBias.^2 ; sBias.^2 ; ...
    mBias.^2./mObs.^2 ; sBias.^2./mObs.^2 ; msBias.^2./mObs.^2};  
for i = 1:length(valplot)   
    cax{i} = [prctile(valplot{i}(:),5) prctile(valplot{i}(:),90)]; 
end

yr = floor(yrNday./10^4); uniyr = unique(yr);
mo = floor((yrNday - yr*10^4)./10^2);
da = yrNday - yr*10^4 - mo.*10^2;
yrznum = datevec(datenum(uniyr,[1 7],1));
[a b] = size(yrznum);

%%% look at increasing modeled values
tic
load('matfiles/par_traditional_performance_gridbin.mat');
toc

strnm = {'Number of Paired Values';'Mean Observed';'Mean Modeled';'Mean Error'; ...
        'Mean Normalized Error';'Normalized Mean Error';'Fractional Error';'Mean Absolute Error'; ...
        'Mean Normalized Absolute Error';'Normalized Mean Absolute Error';'Fractional Absolute Error'; ...
        'Correlation';'Correlation Squared';'Standard Error';'Mean Squared Error'; ...
        'Root Mean Squared Error';'Normalized Root Mean Squared Error'; ...
        'Mean Error DIV Standard Error';'Mean Error Squared DIV Mean Squared Error'; ...
        'Standard Error Squared DIV Mean Squared Error';'beta1';'Variance of Observed';'Variance of Modeled'; ...
        'Mean Error Squared' ; 'Standard Error Squared' ; ...
        'Mean Error Squared DIV Mean Observed Squared' ; 'Standard Error Squared DIV Mean Observed Squared' ; 'Mean Squared Error DIV Mean Observed Squared'};
len = length(mBias);
mBias2 = cell(len,1); sBias2 = cell(len,1); mBias2DmObs2 = cell(len,1);
sBias2DmObs2 = cell(len,1); msBias2DmObs2 = cell(len,1);
for i = 1:len
    mBias2{i} = mBias{i}.^2; sBias2{i} = sBias{i}.^2; 
    mBias2DmObs2{i} = mBias{i}.^2./mObs{i}.^2; 
    sBias2DmObs2{i} = sBias{i}.^2./mObs{i}.^2; 
    msBias2DmObs2{i} = msBias{i}.^2./mObs{i}.^2; 
end
valplot = { num ; mObs ; mMod ; mBias ; nBias ; nmBias ; fBias ; ...
    mErr ; nErr ; nmErr ; fErr ; R ; R2 ; sBias ; msBias ; ...
    rmsBias ; nrmsBias ; mDsBias ; m2DmsBias ; s2DmsBias ; beta1 ; vObs ; vMod ; ...
    mBias2 ; sBias2 ; ...
    mBias2DmObs2 ; sBias2DmObs2 ; msBias2DmObs2 }; 

%%% look at increasing percentile
% loop through a few days 
for l = 1:length(strnm)
    for i = 1:a
        disp(i);
        idx = yrNday == yrznum(i,1)*10000+yrznum(i,2)*100+yrznum(i,3);

        for j = 1:10

            figure; hold on;

            load('../09_mfiles_projections/USAcontiguous.mat');
            cd ../09_mfiles_projections
            plotax = ell2lambertcc([x,y],'whiproj2001');
            cd ../13_mfiles_modelperformance

            % plotting 
            lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
            [xg yg Zg] = plotField(CTMlocs,valplot{l}{idx}(:,j),lax,[plotax(:,1) plotax(:,2)]);
            caxis(cax{l});   
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
            title(sprintf('%s of PM_{2.5} on %s, %% %d',strnm{l},datestr(yrznum(i,:)),j*10));

            % save figure 
            set(gcf,'Position',[0 0 800 600]);       
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_%s_per%d_grid.png',strnm{l},datestr(yrznum(i,:)),j*10));

            drawnow;
            frameI = getframe(gcf);
            imI = frame2im(frameI);
            [imindI,cmI] = rgb2ind(imI,256);
            outfileI = sprintf('figures/%s_per_%s_grid.gif',strnm{l},datestr(yrznum(i,:)));
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

%%% look at increasing modeled values

% load data
load('matfiles/traditional_performance_grid.mat');
valplot = { mBias ; mBias.^2 ; sBias.^2 ; mBias.^2+sBias.^2 };
cax = [];
for i = 1:length(valplot)   
    cax{i} = [prctile(valplot{i}(:),5) prctile(valplot{i}(:),90)]; 
end

tic
load('matfiles/traditional_performance_extendbin.mat');
toc

% measure names/values
temp = repmat([NaN modplots],length(CTMlocs),1);
modplotsrep = cell(length(numGivMod),1);
for i = 1:length(modplotsrep), modplotsrep{i} = temp; end
valplot1 = arrayfun( @(x) mObsGivMod{x}-modplotsrep{x}, 1:length(mObsGivMod), 'UniformOutput', false)';
valplot2 = arrayfun( @(x) (mObsGivMod{x}-modplotsrep{x}).^2, 1:length(mObsGivMod), 'UniformOutput', false)';
valplot3 = arrayfun( @(x) (mObsGivMod{x}-modplotsrep{x}).^2+vObsGivMod{x}, 1:length(mObsGivMod), 'UniformOutput', false)';
valplot = { valplot1 ; valplot2 ; vObsGivMod ; valplot3 };
strnm = { 'Mean Error' ; 'Mean Error Squared' ; 'Standard Error Squared' ; 'Mean Squared Error' };

% loop through a few days 
for l = 1:length(strnm) % measure
    for i = 1:a % day
        disp(i);
        idx = yrNday == yrznum(i,1)*10000+yrznum(i,2)*100+yrznum(i,3);
        for j = 1:length(modplots) % modeled values

            figure; hold on;

            load('../09_mfiles_projections/USAcontiguous.mat');
            cd ../09_mfiles_projections
            plotax = ell2lambertcc([x,y],'whiproj2001');
            cd ../13_mfiles_modelperformance

            % plotting 
            lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
            [xg yg Zg] = plotField(CTMlocs,valplot{l}{idx}(:,j+1),lax,[plotax(:,1) plotax(:,2)]);
            caxis(cax{l});   
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
            title(sprintf('%s of PM_{2.5} on %s, giv Mod = %d',strnm{l},datestr(yrznum(i,:)),modplots(j)));

            % save figure 
            set(gcf,'Position',[0 0 800 600]);       
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_%s_Mod%d_grid.png',strnm{l},datestr(yrznum(i,:)),modplots(j)));

            drawnow;
            frameI = getframe(gcf);
            imI = frame2im(frameI);
            [imindI,cmI] = rgb2ind(imI,256);
            outfileI = sprintf('figures/%s_Mod_%s_grid.gif',strnm{l},datestr(yrznum(i,:)));
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

end