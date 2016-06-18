function [] = exclusionCriteria_plots2()
% new plots from the 02/17/2016

%%% one exclusion criteria

% % change in MSE from kriging to METHOD with CMAQ>70,80,90th %tile
% load data
load('matfiles/exclus_three.mat');
colorz = {'ro-','bo-','go-','rs-','bs-','gs-','r*-','b*-','g*-'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
figure; hold on;
n = 1;
for i = [1 3 4]
    for j = 1:3  
        temp1 = MSE_krig2BME{2,3}{1,6+j}(:,1);
        temp2 = MSE_krig2BME{2,3}{1,6+j}(:,i+1);
        plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
        strz{n} = sprintf('%s, CMAQ >=%0.3f - %d pnts', ...
            methodz{i+1},cutoffs{1}(6+j),pntsNcalc{1,2}{1,6+j}(1,1));
        n = n + 1;
    end
end

legend(strz);
title('% change in MSE from krig to METHOD for ME2/MSE>0th%tile across cross validation radii');
xlabel('km');
ylabel('% change in MSE');  
ylim([-50 100]);
plot([0 900],[0 0],'k-');
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures/CMAQ_percentile_ME2dMSE0th_krig2METHOD.png');
close all;

% % change in MSE from kriging to METHOD with lambda1>70,80,90th %tile
% load data
load('matfiles/exclus_three.mat');
colorz = {'ro-','bo-','go-','rs-','bs-','gs-','r*-','b*-','g*-'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
figure; hold on;
n = 1;
for i = [1 3 4]
    for j = 1:3   
        temp1 = MSE_krig2BME{1,2}{6+j,1}(:,1);
        temp2 = MSE_krig2BME{1,2}{6+j,1}(:,i+1);
        plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
        strz{n} = sprintf('%s, lambda1 >=%0.3f - %d pnts', ...
            methodz{i+1},cutoffs{1}(6+j),pntsNcalc{1,2}{1,6+j}(1,1));
        n = n + 1;
    end
end

legend(strz);
title('% change in MSE from krig to METHOD for ME2/MSE>0th%tile across cross validation radii');
xlabel('km');
ylabel('% change in MSE');  
ylim([-50 100]);
plot([0 900],[0 0],'k-');
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures/lambda1_percentile_ME2dMSE0th_krig2METHOD.png');
close all;

%%% two exclusion critieria

% with ME2/MSE>70th %tile, changee in MSE from kriging to METHOD with
% CMAQ>70,80,90th %tile
% % change in MSE from kriging to METHOD with CMAQ>70,80,90th %tile
% load data
load('matfiles/exclus_three.mat');
colorz = {'ro-','bo-','go-','rs-','bs-','gs-','r*-','b*-','g*-'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
figure; hold on;
n = 1;
for i = [1 3 4]
    for j = 1:3   
        temp1 = MSE_krig2BME{2,3}{8,6+j}(:,1);
        temp2 = MSE_krig2BME{2,3}{8,6+j}(:,i+1);
        plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
        strz{n} = sprintf('%s, CMAQ >=%0.3f - %d pnts', ...
            methodz{i+1},cutoffs{1}(6+j),pntsNcalc{1,2}{1,6+j}(1,1));
        n = n + 1;
    end
end

legend(strz);
title('% change in MSE from krig to METHOD for ME2/MSE>70th%tile across cross validation radii');
xlabel('km');
ylabel('% change in MSE');  
ylim([-50 100]);
plot([0 900],[0 0],'k-');
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures/CMAQ_percentile_ME2dMSE70th_krig2METHOD.png');
close all;

% with ME2/MSE>70th %tile, change in MSE from kriging to METHOD with
% lambda1>70,80,90th %tile

load('matfiles/exclus_three.mat');
colorz = {'ro-','bo-','go-','rs-','bs-','gs-','r*-','b*-','g*-'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
figure; hold on;
n = 1;
for i = [1 3 4]
    for j = 1:3   
        temp1 = MSE_krig2BME{1,2}{6+j,8}(:,1);
        temp2 = MSE_krig2BME{1,2}{6+j,8}(:,i+1);
        plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
        strz{n} = sprintf('%s, lambda1 >=%0.3f - %d pnts', ...
            methodz{i+1},cutoffs{1}(6+j),pntsNcalc{1,2}{1,6+j}(1,1));
        n = n + 1;
    end
end

legend(strz);
title('% change in MSE from krig to METHOD for ME2/MSE>70th%tile across cross validation radii');
xlabel('km');
ylabel('% change in MSE');  
ylim([-50 100]);
plot([0 900],[0 0],'k-');
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures/lambda1_percentile_ME2dMSE70th_krig2METHOD.png');
close all;

%%% relative to DS

% plots of % change in MSE/R for different Xval radius 200/300/400 from DS
% to METHOD
% load data
load('matfiles/exclus_cont_sing_cutoffs.mat');

colorz = {'ro-','bs-','c*-'};
Xvalrad = {'200','300','400'};
varz = {'ME2dMSE'};

% for MSE
for i = 1:3 % each of the Xval radius
    
    figure; hold on;
    n = 1;
    for j = [1 2 5] % each of the methods compared to DS

        % loop through each cross validation radius
        for k = 1:10        
            % quickly getting all values for plot
            temp1 = []; temp2 = [];
            for l = 1:10
                temp1(l) = MSE_krig2BME{11}{l}(i+2,4); 
                temp2(l) = MSE_krig2BME{11}{l}(i+2,j); 
            end
        end    
        plot(0:10:90,100.*(temp2-temp1)./temp1,colorz{n});
        n = n + 1;
    end
    plot([0 90],[0 0],'k--');

    title(sprintf('%% change in MSE from DS to METHOD across percentiles for %s for %skm cross validation radius',varz{1},Xvalrad{i}));
    xlabel('>= percentile');
    ylabel('% change in MSE from DS to METHOD');
    legend('krig','BME','CMAQ','location','best');
    ylim([-40 60]);

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/MSE_MEdMSE_percentile_%skm_DS2METHOD.png',Xvalrad{i}));

end
close all;

% for R
for i = 1:3 % each of the Xval radius
    
    figure; hold on;
    n = 1;
    for j = [1 2 5] % each of the methods compared to kriging

        % loop through each cross validation radius
        for k = 1:10        
            % quickly getting all values for plot
            temp1 = []; temp2 = [];
            for l = 1:10
                temp1(l) = R_krig2BME{11}{l}(i+2,4); 
                temp2(l) = R_krig2BME{11}{l}(i+2,j); 
            end
        end    
        plot(0:10:90,100.*(temp2-temp1)./temp1,colorz{n});
        n = n + 1;
    end
    plot([0 90],[0 0],'k--');

    title(sprintf('%% change in R from DS to METHOD across percentiles for %s for %skm cross validation radius',varz{1},Xvalrad{i}));
    xlabel('>= percentile');
    ylabel('% change in R from DS to METHOD');
    legend('krig','BME','CMAQ','location','best');
    ylim([-30 30]);

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/R_MEdMSE_percentile_%skm_DS2METHOD.png',Xvalrad{i}));

end
close all;

end