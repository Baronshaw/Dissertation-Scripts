function [] = gifTraditionalVSSCurve()
% this function will create the gifs of displayTraditionalVSSCruve.m

% load results
cd ../BMELIB2.0b
startup
cd ../05_mfiles_crossvalidation
n = 3; 
deltaT = 365; 
numbins = 10; 
minpnts = 150; 
negval = 0;
deltaT = 90; exstr = '_T90';
load(sprintf('../matfiles/correctTradVSSCurve%s.mat',exstr));
unilocs = unique(allLocs(:,1:2),'rows');
unidays = unique(allLocs(:,3));
dayz = datevec(datenum(2001,1:12,1));
% dayz = datevec(datenum(2001,1,1):datenum(2001,12,31));
dayWHI = dayz(:,1:3);
dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);

% loop through each model performance statistics
statstr = { 'ME' ; 'NME' ; 'SE' ; 'R' ; 'R2' ; 'MAE' ; 'NMAE' ; 'RMSE' ; 'NRMSE' };
methstr = { 'Traditional' ; 'SCurve' };
Ttoplot = { TME ; TNME ; TSE ; TR ; TR2 ; TMAE ; TNMAE ; TRMSE ; TNRMSE };
Stoplot = { SME ; SNME ; SSE ; SR ; SR2 ; SMAE ; SNMAE ; SRMSE ; SNRMSE };

% go through all the maps
for i = 1:length(statstr)
    
    cax = [prctile(Ttoplot{i},5) prctile(Ttoplot{i},95)];
    
    for j = 1:2
        
        if j == 1, toshow = Ttoplot{i}; else toshow = Stoplot{i}; end
        
        for k = 1:length(dayz)
            
            idx = datenum(dayz(k,:)) == allLocs(:,3);
        
            % plot figure
            cd ../09_mfiles_projections
            load('USAcontiguous.mat');
            plotax = ell2lambertcc([x,y],'whiproj2001');
            cd ../05_mfiles_crossvalidation
            [xg yg Zg] = plotField(unilocs,toshow(idx),[],plotax);      
            caxis(cax);
            colorbar;
            % setting axis
            xlabel('km');
            ylabel('km');
            axis([ -3000000 3000000 -2000000 1500000 ]);
            % overlaying the states
            load('../09_mfiles_projections/USAstates5.mat');
            for l = 1:length(X)
                cd ../09_mfiles_projections
                states = ell2lambertcc([X{l},Y{l}],'whiproj2001');
                cd ../05_mfiles_crossvalidation
                plot(states(:,1),states(:,2),'k-');
            end
            title(sprintf('%s map of %s on %s',methstr{j},statstr{i},datestr(datenum(dayz(k,:)),'mm/dd/yyyy')));

            % save figure
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            drawnow;
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            outfile = sprintf('../modperplots/%s_map_%s%s.gif',methstr{j},statstr{i},exstr);
            % On the first loop, create the file. In subsequent loops, append.
            if k == 1
                imwrite(imind,cm,outfile,'gif','DelayTime',1,'loopcount',inf);
            else
                imwrite(imind,cm,outfile,'gif','DelayTime',1,'writemode','append');
            end
            close all;
        end
       
    end
    
end

% go through all lambda2s
temp = datevec(datenum(2001,1:12,1));
dayWHI = temp(:,1:3);
dayWHIdisps = dayWHI(:,1)*10^4 + dayWHI(:,2)*10^2 + dayWHI(:,3);
numdisp = 5;
cd ../09_mfiles_projections
load('USAcontiguous.mat');
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../05_mfiles_crossvalidation
in = inpolygon(unilocs(:,1),unilocs(:,2),plotax(:,1),plotax(:,2));
rand('seed',0);
idxlocs = randsample(sum(in),numdisp);
unilocssub = unilocs(in,:);
for j = 1:numdisp

    for k = 1:length(dayz)
    
        % load soft data info
        load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
        dayWHIdisps(k),deltaT,n,minpnts,numbins,negval));

        idx = unilocssub(idxlocs(j),1) == CTMlocs(:,1) & unilocssub(idxlocs(j),2) == CTMlocs(:,2);

        % plot figure
        figure; hold on;
        plot(mean_Mod(idx,:),var_Mod(idx,:),'b.-');
        xlim([0 20]);
        ylim([0 20]);
        xlabel('Mean Modeled');
        ylabel('Lambda2');
        title(sprintf('Lambda2 VS Modeled at (%.0fkm,%.0fkm) Location #%d on %s', ...
            floor(unilocssub(idxlocs(j),1)./1000),floor(unilocssub(idxlocs(j),2)./1000),j,...
            datestr(datenum(dayz(k,:)),'mm/dd/yyyy')));

        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        outfile = sprintf('../modperplots/L2VSMod_Loc%d%s.gif',j,exstr);
        % On the first loop, create the file. In subsequent loops, append.
        if k == 1
            imwrite(imind,cm,outfile,'gif','DelayTime',1,'loopcount',inf);
        else
            imwrite(imind,cm,outfile,'gif','DelayTime',1,'writemode','append');
        end
        close all;
    end

end

end