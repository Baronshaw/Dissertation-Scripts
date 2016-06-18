function [] = visualize_covmodel_DS()
% this function will visualize the covariance model for the DS method

% loading data
for i = 2001:2002
    load(sprintf('../matfiles/prepCTMandObs_%d.mat',i));
    Modall{i-2000,1} = Mod; Obsall{i-2000,1} = Obs;
    coordObsall{i-2000,1} = coordObs; cht{i-2000,1} = yrmodaObs;
end
zhm = cell2mat(Modall); zho = cell2mat(Obsall);
ch = cell2mat(coordObsall); cht = cell2mat(cht);

% modifying dates
yrall = floor(cht./10000);
moall = floor((cht - yrall*10000)/100);
daall = cht  - yrall*10000 - moall*100;
cht = datenum(yrall,moall,daall);
unidates = unique(cht);
len = length(unique(cht));
unidays = datevec(unidates);
unidays = unidays(:,1).*10000 + unidays(:,2).*100 + unidays(:,3);

% loading experimental covariance
load('matdata/DS_results_v1.mat');

for i = 1:30:len
    disp(i);
    % covariance parameters
    A11 = betac{i}(1); A12 = betac{i}(2); A22 = betac{i}(3); phi0 = betac{i}(4); phi1 = betac{i}(5);
    
    for j = 1:10
        for k = j:10
            
            % more covariance parameters
            xB1 = problems{i}(j); xB2 = problems{i}(k); tau2 = problems{i}(end);
            % modeled covariance
            x = 0:max(rLags{i}{j,k});            
            y = (A11-A12.*xB1).*(A11-A12.*xB2).*exp(-1./phi0.*x) + ...
                (A12-A22.*xB1).*(A12-A22.*xB2).*exp(-1./phi1.*x) + tau2.*exp(-0.5.*x);

            % visualize covariance
            figure; hold on;
            plot(rLags{i}{j,k},Crtests{i}{j,k},'b.'); % experimental covariance
            plot(x,y,'r-');
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            xlabel('km');
            ylabel('Covariance');
            title(sprintf('DS Covariance for PM2.5 on %d with %dth %%tile and %dth %%tile',unidays(i),j*10,k*10));

            % save figure
            set(gcf,'Position',[0 0 800 600]);
            set(gcf,'PaperUnits','inches');    
            set(gcf,'PaperPosition',[0 0 800 600]./100);
            set(gcf,'PaperPositionMode','manual');
            set(gca,'XTickLabel',get(gca,'XTick')/1000);
            set(gca,'YTickLabel',get(gca,'YTick')/1000);
            print(gcf,'-painters','-dpng','-r600',sprintf('covariance_models/Cov_%d_%0.2d_%0.2d.png',unidays(i),j,k));

        end
        close all;
    end
    
end
close all;
end