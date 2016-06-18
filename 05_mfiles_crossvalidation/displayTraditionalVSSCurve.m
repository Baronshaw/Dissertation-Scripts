function [] = displayTraditionalVSSCurve(numdisp)
% this function will take the results of correctTraditionalVSSCurve.m and
% create maps for each day and time series for select locations for each
% model performance statistics

if nargin < 1, numdisp = 5; end

% load results
cd ../BMELIB2.0b
startup
cd ../05_mfiles_crossvalidation
exstr = '_T90';
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
            print(gcf,'-painters','-dpng','-r600',sprintf('../modperplots/%s_map_%s_%d%s.png',methstr{j},statstr{i},dayWHIdisps(k),exstr)); 
        
            if mod(k,20) == 0, close all; end
            
        end
        
        close all
        
    end
    
end

% go through some locations
cd ../09_mfiles_projections
load('USAcontiguous.mat');
plotax = ell2lambertcc([x,y],'whiproj2001');
cd ../05_mfiles_crossvalidation
in = inpolygon(unilocs(:,1),unilocs(:,2),plotax(:,1),plotax(:,2));
rand('seed',0);
idxlocs = randsample(sum(in),numdisp);
unilocssub = unilocs(in,:);

% creat map of selected locations
figure; hold on;
load('../09_mfiles_projections/USAstates5.mat');
for i = 1:length(X)
    cd ../09_mfiles_projections
    states = ell2lambertcc([X{i},Y{i}],'whiproj2001');
    cd ../05_mfiles_crossvalidation
    plot(states(:,1),states(:,2),'k-');
end
for i = 1:numdisp
    text(unilocssub(idxlocs(i),1),unilocssub(idxlocs(i),2),num2str(i),'fontname','Arial','color','blue')
end
% save figure
set(gcf,'Position',[0 0 800 600]);
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperPosition',[0 0 800 600]./100);
set(gcf,'PaperPositionMode','manual');
print(gcf,'-painters','-dpng','-r600','../modperplots/TS_locations.png'); 

for i = 1:length(statstr)

    for j = 1:numdisp

        idx = unilocssub(idxlocs(j),1) == allLocs(:,1) & unilocssub(idxlocs(j),2) == allLocs(:,2);

        % plot figure
        figure; hold on;
        plot(allLocs(idx,3),Ttoplot{i}(idx),'b.',allLocs(idx,3),Stoplot{i}(idx),'r.');
        legend('Traditional','SCurve','Location','Best');
        xlabel('month/year in 2001');
        ylabel(sprintf('%s',statstr{i}));
        set(gca,'XTickLabel',datestr(get(gca,'XTick'),'mm/yy'));
        title(sprintf('Time Series of %s at (%.0fkm,%.0fkm) Location #%d', ...
            statstr{i},floor(unilocssub(idxlocs(j),1)./1000),floor(unilocssub(idxlocs(j),2)./1000),j));

        % save figure
        set(gcf,'Position',[0 0 800 600]);
        set(gcf,'PaperUnits','inches');    
        set(gcf,'PaperPosition',[0 0 800 600]./100);
        set(gcf,'PaperPositionMode','manual');
        print(gcf,'-painters','-dpng','-r600',sprintf('../modperplots/TS_%s_%d%s.png',statstr{i},j,exstr)); 

    end 
    close all
end

end