function [] = QAplots(monthsdate,testrun,strfolder,soft,constant,gauss)
% this function will create time series for numpick locations showing 1)
% residual plots and 2) plots with 5 closest stations

if nargin < 1, monthsdate = 102:104; end % months to show
if nargin < 2, testrun = 0; end % file to load
if nargin < 3, strfolder = 'sum06'; end
if nargin < 4, soft = 1; end % soft data or not
if nargin < 5, constant = 0; end % constant offset or not
if nargin < 6, gauss = 1; end % gaussian soft data or not

if testrun == 1 
    teststr = 'test_'; 
    pathstr = 'matfiles_QAQC';
elseif testrun == 0
    teststr = ''; 
    pathstr = 'matfiles_est';
elseif testrun == 2
    teststr = 'BMEest_'; 
    pathstr = 'matfiles_QAQC';
end

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 0, gaussstr = '_nongauss'; else gaussstr = '_gauss'; end

% load data and combining all months of interest
for i = 1:length(monthsdate)
    
    monthsplot(i) = datenum(1998,monthsdate(i),1);
    inputdates1 = datevec(datenum(1998,monthsdate(i),1));
    inputdates2 = datevec(datenum(1998,monthsdate(i)+1,1)-1);
    startdate = inputdates1(:,1:3);
    enddate = inputdates2(:,1:3);
    daydisp1 = startdate(1)*10000 + startdate(2)*100 + startdate(3);
    daydisp2 = enddate(1)*10000 + enddate(2)*100 + enddate(3);
    
    if soft == 0 
        load(sprintf('../%s/%sWHIMS_%d_%d.mat', ...
            pathstr,teststr,daydisp1,daydisp2)); 
    else
        load(sprintf('../%s/%sWHIMS_%d_%d%s%s%s.mat', ...
            pathstr,teststr,daydisp1,daydisp2,softstr,constr,gaussstr));  
    end
    load(sprintf('../%s/%smI_%d_%d.mat',pathstr,teststr,daydisp1,daydisp2));
    
    if soft == 1
        cks{i,1} = ck; zks{i,1} = zk; vks{i,1} = vk;
        temp1s{i,1} = temp1; temp1as{i,1} = temp1a;
        temp1bs{i,1} = temp1b; mItests{i,1} = mItest; 
    else
        cks{i,1} = ck; zks{i,1} = zk; vks{i,1} = vk;
        temp1s{i,1} = temp1'; temp1as{i,1} = temp1a';
        temp1bs{i,1} = temp1b'; mItests{i,1} = mItest; 
    end
    
    if soft == 1
        allhards{i,1} = allhard; alllambda1s{i,1} = alllambda1;
        alllambda2s{i,1} = alllambda2; allchs{i,1} = allch; allcss{i,1} = allcs;
        temp2s{i,1} = temp2; 
    end
          
end

cks = cell2mat(cks);
zks = cell2mat(zks);
vks = cell2mat(vks);
temp1s = cell2mat(temp1s);
temp1as = cell2mat(temp1as);
temp1bs = cell2mat(temp1bs);
mItests = cell2mat(mItests);

if soft == 1
    temp2s = cell2mat(temp2s); 
end

if soft == 1
    % combining all the cells properly
    a = cellfun(@length,allhards); a = [0;a];
    b = cumsum(a);
    c = sum(cellfun(@length,allhards));
    allhard = cell(c,1);
    alllambda1 = cell(c,1);
    alllambda2 = cell(c,1);
    allch = cell(c,1);
    allcs = cell(c,1);
    for i = 1:length(b)-1
        allhard(1+b(i):b(i+1)) = allhards{i,1};
        alllambda1(1+b(i):b(i+1),:) = alllambda1s{i,1};
        alllambda2(1+b(i):b(i+1),:) = alllambda2s{i,1};
        allch(1+b(i):b(i+1),:) = allchs{i,1};
        allcs(1+b(i):b(i+1),:) = allcss{i,1};
    end
    allhards = allhard; alllambda1s = alllambda1; alllambda2s = alllambda2;
    allchs = allch; allcss = allcs;
end

% selecting stations to plot
load('../matfiles_QAQC/forQA.mat');
if testrun == 1
    plotx = simuxQA;
    ploty = simuyQA;
    plotidx = simuQA;
    if strcmp(strfolder,'whole')
        plotx = plotx(1:5);
        ploty = ploty(1:5);
        plotidx = plotidx(1:5);
    end
elseif testrun == 0
    plotx = partxQA;
    ploty = partyQA;
    plotidx = partQA;
    if strcmp(strfolder,'whole')
        plotx = plotx(1:5);
        ploty = ploty(1:5);
        plotidx = plotidx(1:5);
    end
end

% getting all the hard data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
zh = zd;
ch = pd;
cMSObs = unique(ch(:,1:2),'rows');

if soft == 1
    % getting all the soft data
    CTMyears = [2001 2002 2005:2007]; 
    for i = 1:length(CTMyears) 
        disp(CTMyears(i));
        load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',CTMyears(i)));
        cssA{i,1} = css;
        lambda1A{i,1} = lambda1;
    end
    css = cell2mat(cssA);
    lambda1 = cell2mat(lambda1A);
end

% looping through all the figures
for i = 1:length(plotidx)
    
    %%% figure 1, mean with min and max
    idx = plotx(i) == cks(:,1) & ploty(i) == cks(:,2);
    figure; hold on; box on;
    plot(cks(idx,3),zks(idx),'-b',cks(idx,3),temp1bs(idx),'--r',cks(idx,3),temp1as(idx),'--g');
    if soft == 1, plot(cks(idx,3),temp2s(idx),'--k'); end
    
    % plot displays
    set(gca,'XTick',monthsplot);
    set(gca,'XTickLabel',datestr(get(gca,'XTick'),'mmmyyyy'));
    ylabel('residual concentration');
    if testrun == 1
        title(sprintf('residual (fake ID # %d)',plotidx(i)));
    else
        title(sprintf('residual (ID # %d)',plotidx(i)));
    end
    if soft == 1 & sum(~isnan(temp2s(idx))) > 0 % soft data used
        legend('mean','min hard','max hard','mean soft');
    else
        legend('mean','min hard','max hard');
    end
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    
    if strcmp(strfolder,'whole')
        print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_QAQC/%s/%smeanminmax%s_%0.5d_%d.pdf',strfolder,teststr,softstr,plotidx(i),monthsdate(1)));
        print(gcf,'-painters','-dpng','-r600',sprintf('../plots_QAQC/%s/%smeanminmax%s_%0.5d_%d.png',strfolder,teststr,softstr,plotidx(i),monthsdate(1)));
    else
        print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_QAQC/%s/%smeanminmax%s_%0.5d.pdf',strfolder,teststr,softstr,plotidx(i)));
        print(gcf,'-painters','-dpng','-r600',sprintf('../plots_QAQC/%s/%smeanminmax%s_%0.5d.png',strfolder,teststr,softstr,plotidx(i)));
    end
    
    %%% figure 2, mean with 5 closest stations
    zkdisplay = zks + mItests;
    figure; hold on; box on;
    plot(cks(idx,3),zkdisplay(idx),'-b','LineWidth',2);
    
    % closest stations
    alldists = sqrt( (plotx(i)-cMSObs(:,1)).^2 + (ploty(i)-cMSObs(:,2)).^2 );
    [sorted sortidx] = sort(alldists);
    
    % getting closest stations
    j = 0;
    n = 0;
    colors = {'r--.','g--.','c--.','m--.','y--.'};
    while j < 5
        idxclose = ch(:,1) == cMSObs(sortidx(n+1),1) & ch(:,2) == cMSObs(sortidx(n+1),2);
        subchplot = ch(idxclose,:);
        subzhplot = zh(idxclose);
        [a b] = ismember(cks(idx,3),subchplot(:,3));
        b(b==0) = [];        
        if sum(a) > 0
            j = j + 1; % increases if finds closest station
            plot(subchplot(b,3),subzhplot(b),colors{j});
        end
        n = n + 1; % increases to get next closest station
    end
    
    if soft == 1
        % soft data
        [a b] = ismember(css(:,3),cks(idx,3));
        cMSCTM = unique(css(a,1:2),'rows');
        alldists = sqrt( (plotx(i)-cMSCTM(:,1)).^2 + (ploty(i)-cMSCTM(:,2)).^2 );
        [sorted sortidx] = sort(alldists);
        if sum(~isnan(temp2s(idx))) > 0 & size(cMSCTM,1) >= 2 % soft data used
            j = 0;
            n = 0;
            colors = {'bs--','rs--'};
            while j < 2
                idxclose = css(:,1) == cMSCTM(sortidx(n+1),1) & css(:,2) == cMSCTM(sortidx(n+1),2);
                subcsplot = css(idxclose,:);
                sublamplot = lambda1(idxclose);
                [a b] = ismember(cks(idx,3),subcsplot(:,3));
                b(b==0) = [];
                if sum(a) > 0
                    j = j + 1; % increases if finds closest station
                    plot(subcsplot(b,3),sublamplot(b),colors{j});
                end
                n = n + 1; % increases to get next closest station
            end
        end
    end
    
    % plot displays
    set(gca,'XTick',monthsplot);
    set(gca,'XTickLabel',datestr(get(gca,'XTick'),'mmmyyyy'));
    ylabel('PM_{2.5} concentration');
    if testrun == 1
        title(sprintf('PM2.5 (fake ID # %d)',plotidx(i)));
    else
        title(sprintf('PM2.5 (ID # %d)',plotidx(i)));
    end
    if soft == 1 & sum(~isnan(temp2s(idx))) > 0 % soft data used
        legend('mean','1st','2nd','3rd','4th','5th','1st s','2nd s');
    else
        legend('mean','1st','2nd','3rd','4th','5th');
    end
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    if strcmp(strfolder,'whole')
        print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_QAQC/%s/%smeanclosest%s_%0.5d_%d.pdf',strfolder,teststr,softstr,plotidx(i),monthsdate(1)));
        print(gcf,'-painters','-dpng','-r600',sprintf('../plots_QAQC/%s/%smeanclosest%s_%0.5d_%d.png',strfolder,teststr,softstr,plotidx(i),monthsdate(1)));
    else
        print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_QAQC/%s/%smeanclosest%s_%0.5d_%d.pdf',strfolder,teststr,softstr,plotidx(i)));
        print(gcf,'-painters','-dpng','-r600',sprintf('../plots_QAQC/%s/%smeanclosest%s_%0.5d_%d.png',strfolder,teststr,softstr,plotidx(i)));
    end
    
end

end