function [] = plotObsWSoft(best)
% this function will create time series that plots observed with estimates
% (LOOCV) using obs and obs+ctm along with lambda1's with lambda2's error
% bars

if nargin < 1, best = 0; end % getting the best performances

if best == 1, beststr = 'best'; else beststr = 'worst'; end

forceddist = 0; 
soft_years = 2001; 
numplots = 10;

% gathering all the data
load(sprintf('../matfiles/Xvalforcediso_LOOCV__nosoft_long_gauss_foriso%dkm.mat',floor(forceddist./1000)));
zkallh = zk_madd; zhallh = zh_Xval; ckallh = ck; vkallh = vk;  

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckallh = cell(length(cMS),1);
for i = 1:length(ckallh)
    disp(i);
    idx2 = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckallh{i,1} = ch(idx2,:);
end
ckallh = cell2mat(ckallh);

load(sprintf('../matfiles/Xvalforcediso_LOOCV__soft_long_gauss_foriso%dkm_timezone_n6.mat',floor(forceddist./1000)));
zkalls = zk_madd; zhalls = zh_Xval; ckalls = ck; vkalls = vk; 

% loading data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
temp2001 = datevec(pd(:,3)); temp2001 = temp2001(:,1);
zh = zd;
ch = pd;
cMS = unique(ch(:,1:2),'rows');
ckalls = cell(length(cMS),1);
for i = 1:length(ckalls)
    disp(i);
    idx2 = cMS(i,1) == ch(:,1) & cMS(i,2) == ch(:,2);
    ckalls{i,1} = ch(idx2 & temp2001==2001,:);
end
ckalls = cell2mat(ckalls);

zkallh = cell2mat(zkallh); zkalls = cell2mat(zkalls);
zhallh = cell2mat(zhallh); zhalls = cell2mat(zhalls);

[lia lib] = ismember(ckallh,ckalls,'rows');
ckallh = ckallh(lia,:);
zkallh = zkallh(lia); zhallh = zhallh(lia);

idx = ~isnan(zkallh) & ~isnan(zkalls); 
zkallh = zkallh(idx); zhallh = zhallh(idx); 
zkalls = zkalls(idx); zhalls = zhalls(idx); 
ckallh = ckallh(idx,:); ckalls = ckalls(idx,:);

% looking at performance by station/year
temp = datevec(ckalls(:,3));
yrs = temp(:,1);
diffR2 = NaN*ones(length(cMS),1);
for i = 1:length(cMS)
    
    idx = cMS(i,1) == ckalls(:,1) & cMS(i,2) == ckalls(:,2) & yrs == soft_years;
    if sum(idx) > 0
        corrH = (corr(zhallh(idx),zkallh(idx))).^2;
        corrS = (corr(zhalls(idx),zkalls(idx))).^2;
        diffR2(i) = corrS - corrH;
    end
    
end

% loading soft mean trend
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 300000 100 50],soft_years));

% starting writing a text file of results
fileID = fopen(sprintf('%s_n6.txt',beststr),'w');
fprintf(fileID,'%s stations\n',beststr);

% plotting time series, plot obs, est hard, est soft, lam1's with lam2's
if best == 1
    [aorig borig] = sort(diffR2,'descend');
    a = aorig(~isnan(aorig)); b = borig(~isnan(aorig)); 
else
    [a b] = sort(diffR2);
end
for i = 1:numplots
    
    fprintf(fileID,'station %d\tdiffR2\t%0.4f\n',i,a(i));
    fprintf(fileID,'\tlocations\t%0.4f\t%0.4f\n',cMS(b(i),1),cMS(b(i),2));
    
    figure; hold on;
    
    idx = cMS(b(i),1) == ckalls(:,1) & cMS(b(i),2) == ckalls(:,2) & yrs == soft_years;
    plot(ckalls(idx,3),zhallh(idx),'ks');
    plot(ckalls(idx,3),zkallh(idx),'bo-');
    plot(ckalls(idx,3),zkalls(idx),'ro-');
    
    % gather all the lambda1's
    tempval = cell2mat(alllambda1{b(i)});
    templ2 = cell2mat(alllambda2{b(i)});
    tempcs = cell2mat(allcs{b(i)});
    s_yrs = datevec(tempcs(:,3)); s_yrs = s_yrs(:,1);
    tempcs = tempcs(s_yrs==2001,:);
    tempval = tempval(s_yrs==2001);
    templ2 = templ2(s_yrs==2001);
    [lia lib] = ismember(round(tempcs),round(pI),'rows');
    lib(lib==0)=[];
    lam1 = tempval + mI(lib);
    lam2 = templ2;

    % within each plot, plot some S-curves (with the worst performance)
    % find the worst performances for monitoring station
    ckallsub = ckalls(idx,:);
    zhallsub = zhallh(idx);
    zkallsubh = zkallh(idx);
    zkallsubs = zkalls(idx);
    stationErr = abs(zhallh(idx)-zkallh(idx)) - abs(zhalls(idx)-zkalls(idx));
    if best == 1
        [a2orig b2orig] = sort(stationErr,'descend');
        a2 = a2orig(~isnan(a2orig)); b2 = b2orig(~isnan(b2orig));
    else
        [a2 b2] = sort(stationErr);
    end
    
    % add the lam1's with lam2's
    errorbar(tempcs(:,3),lam1,sqrt(lam2),'rx');
    
    % the ones I'll be viewing
    idxB = [];
    for j = 1:10
        idxB = [ idxB ; find( tempcs(:,3) == ckallsub(b2(j),3) ) ];
    end
    errorbar(tempcs(idxB,3),lam1(idxB),sqrt(lam2(idxB)),'cx');
    
    legend('obs','est h','est s','lam1');
    title(sprintf('obs, est, and soft with diffR2=%f',diffR2(b(i))));
    set(gca,'XTickLabel',datestr(datevec(get(gca,'XTick'))));
    xlabel('time');
    ylabel('PM2.5 conc');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperPosition',[0 0 800 600]./100);
    set(gcf,'PaperPositionMode','manual');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/diag_%s_station_%d_n6.pdf',beststr,i));
    print(gcf,'-painters','-dpng','-r600',sprintf('../plots/diag_%s_station_%d_n6.png',beststr,i));
    
    for j = 1:numplots
        
        fprintf(fileID,'loc %d\tday\t%s\n',j,datestr(ckallsub(b2(j),3),'mm/dd/yyyy'));
        fprintf(fileID,'\tdiff soft est - obs\t%0.4f\n',zkallsubs(b2(j))-zhallsub(b2(j)));
        fprintf(fileID,'\tdiff hard est - obs\t%0.4f\n',zkallsubh(b2(j))-zhallsub(b2(j)));
                
        cd ../04_mfiles_softdata
        dayWHI = datevec(ckallsub(b2(j),3)); dayWHI = dayWHI(:,1:3);
        daysWHIdisp = dayWHI(1)*10^4 + dayWHI(2)*10^2 + dayWHI(3);
        
        % finding all the correct soft locations
        idxsoft = tempcs(:,3) == ckallsub(b2(j),3);
        locations = tempcs(idxsoft,:);
        correct = zhallsub(b2(j));
        
        fprintf(fileID,'\tlam1''s similar?\t%0.4f\n',var(lam1(idxsoft)));
        texta = zkallsubs(b2(j))-zhallsub(b2(j));
        textb = zkallsubh(b2(j))-zhallsub(b2(j));
        fprintf(fileID,'\tlam1''s influencial?\t%d\n',(texta-textb)./textb>0.1);
        
        for k = 1:3

            location = locations(k,1:2);             
            [diagnostics] = plotScurve_diagnostics(dayWHI,daysWHIdisp,location,correct);
            
            % save figure
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/diag_%s_station_%d_%d_Scurve_%d_n6.pdf',...
                beststr,i,daysWHIdisp,k));
            print(gcf,'-painters','-dpng','-r600',sprintf('../plots/diag_%s_station_%d_%d_Scurve_%d_n6.png',...
                beststr,i,daysWHIdisp,k));
            
            [diagnostics] = plotScurve_diagnostics2(1,dayWHI,daysWHIdisp,location,correct);
            
            % save figure
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/diag_%s_station_%d_%d_Scurve_%d_space_n6.pdf',...
                beststr,i,daysWHIdisp,k));
            print(gcf,'-painters','-dpng','-r600',sprintf('../plots/diag_%s_station_%d_%d_Scurve_%d_space_n6.png',...
                beststr,i,daysWHIdisp,k));
            
            [diagnostics] = plotScurve_diagnostics2(0,dayWHI,daysWHIdisp,location,correct);
            
            % save figure
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            print(gcf,'-painters','-dpdf','-r600',sprintf('../plots/diag_%s_station_%d_%d_Scurve_%d_time_n6.pdf',...
                beststr,i,daysWHIdisp,k));
            print(gcf,'-painters','-dpng','-r600',sprintf('../plots/diag_%s_station_%d_%d_Scurve_%d_time_n6.png',...
                beststr,i,daysWHIdisp,k));            
            
            if k == 1, fprintf(fileID,'s-curves\tlam1\t%0.4f\tsparse?\t%d\tquantile\t%d\n',...
                    diagnostics(1),diagnostics(2),diagnostics(3));
            else fprintf(fileID,'\tlam1\t%0.4f\tsparse?\t%d\tquantile\t%d\n',...
                    diagnostics(1),diagnostics(2),diagnostics(3)); 
            end
            
        end
        cd ../05_mfiles_crossvalidation
    end
    
    close all;
       
end

fclose(fileID);

end