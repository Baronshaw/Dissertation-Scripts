function [] = exclusionCriteria_plots()
% this function will create some exclusion criteria plots

% plot idea #1: five plots (one for each method), of ME2/MSE as a function
% of percentile. Each plot has 10 lines (one for each cross-validation
% radius)

% load data
load('matfiles/exclus_cont_sing_cutoffs.mat');

colorz = {'ro-','bo-','go-','ko-','rx-','bx-','gx-','kx-','r*-','b*-'};
methodz = {'krig','BME','CAMP','DS','CMAQ'};
varz = {'ME2dMSE'};

% change in MSE
for i = 1:5
    
    figure; hold on;

    % loop through each cross validation radius
    for j = 1:10
        
        % quickly getting all values for plot
        temp = [];
        for k = 1:10
            temp(k) = MSE_krig2BME{11}{k}(j,i); 
        end
        
        plot(0:10:90,temp,colorz{j});
        
    end
    
    ylim([10 100]);
    title(sprintf('MSE for %s across percentiles for %s',methodz{i},varz{1}));
    xlabel('>= percentile');
    ylabel('MSE');
    legend('0km','100km','200km','300km','400km','500km','600km', ...
        '700km','800km','900km','location','best');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/MSE_MEdMSE_percentile_%s.png',methodz{i}));
    
end
close all;

% change in R
for i = 1:5
    
    figure; hold on;

    % loop through each cross validation radius
    for j = 1:10
        
        % quickly getting all values for plot
        temp = [];
        for k = 1:10
            temp(k) = R_krig2BME{11}{k}(j,i); 
        end
        
        plot(0:10:90,temp,colorz{j});
        
    end
    
    ylim([0.3 1]);
    title(sprintf('R for %s across percentiles for %s',methodz{i},varz{1}));
    xlabel('>= percentile');
    ylabel('R');
    legend('0km','100km','200km','300km','400km','500km','600km', ...
        '700km','800km','900km','location','best');
    
    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/R_MEdMSE_percentile_%s.png',methodz{i}));
    
end
close all;

% plot idea #2: plot the % change in MSE from kriging to BME as a function
% of increasing radius ME2/MSE>=80th%ile with increasing cutoffs for
% lambda1. The plot has three lines (1 for each lambda1 cutoff)

% load data
load('matfiles/exclus_three.mat');

figure; hold on;
for i = 1:3   
    temp1 = MSE_krig2BME{1,2}{i+7,9}(:,1);
    temp2 = MSE_krig2BME{1,2}{i+7,9}(:,2);
    plot(0:100:900,100.*(temp2-temp1)./temp1,colorz{i}); 
    strz{i} = sprintf('lambda1 >=%0.3f - %d pnts',cutoffs{1}(i+7),pntsNcalc{1,2}{i+7,9}(1,1));
end

legend(strz);
title('% change in MSE from krig to BME for ME2/MSE>70th%ile across cross validation radii');
xlabel('km');
ylabel('% change in MSE');   
    
% save figure
set(gcf,'Position',[0 0 800 600]);
print(gcf,'-painters','-dpng','-r600','figures/lambda1_percentile_ME2dMSE70th_krig2BME.png');
close all;

% plot idea #3: % change in MSE from krig in METHOD for the remaining four
% methods as a function of percentile for ME2/MSE for 200/300/400 km cross 
% validation radius

% load data
load('matfiles/exclus_cont_sing_cutoffs.mat');

colorz = {'ro-','bo-','go-','ko-','rx-','bx-','gx-','kx-','r*-','b*-'};
varz = {'ME2dMSE'};
Xvalrad = {'200','300','400'};

% for MSE
for i = 1:3 % each of the Xval radius
    
    figure; hold on;

    for j = 1:4 % each of the methods compared to kriging

        % loop through each cross validation radius
        for k = 1:10        
            % quickly getting all values for plot
            temp1 = []; temp2 = [];
            for l = 1:10
                temp1(l) = MSE_krig2BME{11}{l}(i+2,1); 
                temp2(l) = MSE_krig2BME{11}{l}(i+2,j+1); 
            end
        end    
        plot(0:10:90,100.*(temp2-temp1)./temp1,colorz{j});

    end
    plot([0 90],[0 0],'k--');

    title(sprintf('%% change in MSE from kriging to METHOD across percentiles for %s for %skm cross validation radius',varz{1},Xvalrad{i}));
    xlabel('>= percentile');
    ylabel('% change in MSE from kriging to METHOD');
    legend('BME','CAMP','DS','CMAQ','location','best');
    ylim([-40 100]);

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/MSE_MEdMSE_percentile_%skm_krig2METHOD.png',Xvalrad{i}));

end
close all;

% for R
for i = 1:3 % each of the Xval radius
    
    figure; hold on;

    for j = 1:4 % each of the methods compared to kriging

        % loop through each cross validation radius
        for k = 1:10        
            % quickly getting all values for plot
            temp1 = []; temp2 = [];
            for l = 1:10
                temp1(l) = R_krig2BME{11}{l}(i+2,1); 
                temp2(l) = R_krig2BME{11}{l}(i+2,j+1); 
            end
        end    
        plot(0:10:90,100.*(temp2-temp1)./temp1,colorz{j});

    end
    plot([0 90],[0 0],'k--');

    title(sprintf('%% change in R from kriging to METHOD across percentiles for %s for %skm cross validation radius',varz{1},Xvalrad{i}));
    xlabel('>= percentile');
    ylabel('% change in MSE from kriging to METHOD');
    legend('BME','CAMP','DS','CMAQ','location','best');
    ylim([-30 30]);

    % save figure
    set(gcf,'Position',[0 0 800 600]);
    print(gcf,'-painters','-dpng','-r600',sprintf('figures/R_MEdMSE_percentile_%skm_krig2METHOD.png',Xvalrad{i}));

end
close all;

end