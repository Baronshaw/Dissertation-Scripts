function [] = exclusionCriteria_plots3()
% new plots from the 02/17/2016 afternoon

%%% one exclusion criteria

% % change in MSE from kriging to METHOD with lambda1>70,80,90th %tile
% load data
load('matfiles/exclus_three.mat');
colorz = {'ro-','bs-','cx-','ro-.','bs-.','cx-.'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
figure; hold on;
n = 1;
for i = [1 3]
    for j = [0 1 3]   
        if j == 0 % baseline
            temp1 = MSE_krig2BME{1,2}{1,1}(:,1);
            temp2 = MSE_krig2BME{1,2}{1,1}(:,i+1);
            plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
            strz{n} = sprintf('%s, baseline - %d pnts', ...
                methodz{i+1},pntsNcalc{1,2}{1,1}(1,1));
        else
            temp1 = MSE_krig2BME{1,2}{6+j,1}(:,1);
            temp2 = MSE_krig2BME{1,2}{6+j,1}(:,i+1);
            plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
            strz{n} = sprintf('%s, lambda1 >=%0.3f - %d pnts', ...
                methodz{i+1},cutoffs{1}(6+j),pntsNcalc{1,2}{6+j,1}(1,1));
        end
        n = n + 1;
    end
end

legend(strz);
title('% change in MSE from krig to METHOD for ME2/MSE>0th%tile across cross validation radii');
xlabel('leave-one-out cross-validation exclusion radius (km)');
ylabel('% change in MSE');  
ylim([-40 55]);
plot([0 900],[0 0],'k-');
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures2/lambda1_2nd_percentile_ME2dMSE0th_krig2METHOD.png');
close all;

%%% two exclusion critieria

% with ME2/MSE>70th %tile, change in MSE from kriging to METHOD with
% lambda1>70,80,90th %tile

load('matfiles/exclus_three.mat');
colorz = {'ro-','bs-','cx-','ro-.','bs-.','cx-.'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
figure; hold on;
n = 1;
for i = [1 3]
    for j = [0 1 3]  
        if j == 0
            temp1 = MSE_krig2BME{1,2}{1,8}(:,1);
            temp2 = MSE_krig2BME{1,2}{1,8}(:,i+1);
            plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
            strz{n} = sprintf('%s, baseline - %d pnts', ...
                methodz{i+1},pntsNcalc{1,2}{1,8}(1,1));
        else
            temp1 = MSE_krig2BME{1,2}{6+j,8}(:,1);
            temp2 = MSE_krig2BME{1,2}{6+j,8}(:,i+1);
            plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{n}); 
            strz{n} = sprintf('%s, lambda1 >=%0.3f - %d pnts', ...
                methodz{i+1},cutoffs{1}(6+j),pntsNcalc{1,2}{6+j,8}(1,1));
        end
        n = n + 1;
    end
end

legend(strz);
title('% change in MSE from krig to METHOD for ME2/MSE>70th%tile across cross validation radii');
xlabel('leave-one-out cross-validation exclusion radius (km)');
ylabel('% change in MSE');  
ylim([-40 55]);
plot([0 900],[0 0],'k-');
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures2/lambda1_2nd_percentile_ME2dMSE70th_krig2METHOD.png');
close all;

%%% relative to DS

% plots of % change in MSE/R for different Xval radius 200 from DS
% to METHOD
% load data
load('matfiles/exclus_cont_sing_cutoffs.mat');

colorz = {'ro-','bs-','c*-'};
Xvalrad = {'200','300','400'};
varz = {'ME2dMSE'};

% for MSE
for i = 1:1 % each of the Xval radius
    
    figure; hold on;
    n = 1;
    for j = [2 5] % each of the methods compared to DS

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
    xlabel('>= percentile exclusion given ME2dMSE');
    ylabel('% change in MSE from DS to METHOD');
    legend('BME','CMAQ','location','best');
    ylim([-40 60]);

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures2/MSE_2nd_MEdMSE_percentile_%skm_DS2METHOD.png',Xvalrad{i}));

end
close all;

end