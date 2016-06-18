function []  = map_MP_cmaq_gridPatch()
% this function will show where patching occured and the results of said
% patching

% loading BME function
cd ../BMELIB2.0b
startup
cd ../13_mfiles_modelperformance

% load data
load('matfiles/traditional_performance_grid.mat');
valplot = { mBias ; mBias.^2 ; sBias.^2 ; mBias.^2+sBias.^2 };
cax = [];
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
load('matfiles/traditional_performance_extendbinFixed.mat');
toc

% measure names/values
temp = repmat([NaN modplots],length(CTMlocs),1);
modplotsrep = cell(length(numGivModF),1);
for i = 1:length(modplotsrep), modplotsrep{i} = temp; end
valplot1 = arrayfun( @(x) modplotsrep{x}-mObsGivModF{x}, 1:length(mObsGivModF), 'UniformOutput', false)';
valplot2 = arrayfun( @(x) (modplotsrep{x}-mObsGivModF{x}).^2, 1:length(mObsGivModF), 'UniformOutput', false)';
valplot3 = arrayfun( @(x) (modplotsrep{x}-mObsGivModF{x}).^2+vObsGivModF{x}, 1:length(mObsGivModF), 'UniformOutput', false)';
valplot = { valplot1 ; valplot2 ; vObsGivModF ; valplot3 };
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
            toplot = valplot{l}{idx}(:,j+1);
            toplot(flagGrid{j}(:,idx)) = NaN;
            lax = [min(plotax(:,1))-100000 max(plotax(:,1))+100000 min(plotax(:,2))-100000 max(plotax(:,2))+100000];
            colorz = redpink;
            if strcmp(strnm{i},'Correlation')==1 | strcmp(strnm{i},'Correlation Squared')==1, colorz = flipud(colorz); end        
            [xg yg Zg] = plotField(CTMlocs,toplot,lax,[plotax(:,1) plotax(:,2)],colorz);
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
            title(sprintf('Patched %s of PM_{2.5} on %s, giv Mod = %d',strnm{l},datestr(yrznum(i,:)),modplots(j)));

            % save figure 
            set(gcf,'Position',[0 0 800 600]);       
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            print(gcf,'-painters','-dpng','-r600',sprintf('figures/%s_%s_Mod%d_gridF.png',strnm{l},datestr(yrznum(i,:)),modplots(j)));

            drawnow;
            frameI = getframe(gcf);
            imI = frame2im(frameI);
            [imindI,cmI] = rgb2ind(imI,256);
            outfileI = sprintf('figures/%s_Mod_%s_gridF.gif',strnm{l},datestr(yrznum(i,:)));
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

% look at the patches 
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
        [xg yg Zg] = plotField(CTMlocs,double(flagGrid{j}(:,idx)),lax,[plotax(:,1) plotax(:,2)]);
        caxis([-1 2]);   
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
        title(sprintf('Location of Patches on %s, giv Mod = %d',datestr(yrznum(i,:)),modplots(j)));

        % save figure 
        set(gcf,'Position',[0 0 800 600]);       
        set(gca,'YTickLabel',get(gca,'YTick')/1000);
        set(gca,'XTickLabel',get(gca,'XTick')/1000);
        print(gcf,'-painters','-dpng','-r600',sprintf('figures/PatchesLocs_%s_Mod%d.png',datestr(yrznum(i,:)),modplots(j)));

        drawnow;
        frameI = getframe(gcf);
        imI = frame2im(frameI);
        [imindI,cmI] = rgb2ind(imI,256);
        outfileI = sprintf('figures/PatchesLocs_Mod_%s_grid.gif',datestr(yrznum(i,:)));
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