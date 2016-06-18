function [] = plotEffectsOfLambdas(numplots)
% this function will make plots of the effect that the soft data has on the
% estimate 

if nargin < 1, numplots = 5; end % number of plots for each scenario

soft_years = [2001:2002 2005 2006:2007];
for i = 1:length(soft_years)
    disp(soft_years(i));
    % load results
    load(sprintf('help_lambda_%d.mat',soft_years(i)));
    load(sprintf('dataneighb_%d.mat',soft_years(i)));
    actuallambda_help2 = ~actuallambda_help;
    test1 = nansum(lambda_help2,2); 
    % I = improved, E = estimate, W = worsened, S = soft
    
    % soft improved, but shouldn't have
    idxIEWS = actuallambda_help == 1 & test1 == 0;
    a = find(idxIEWS==1);
    plotIEWS = randsample(length(a),numplots);
    
    % soft improved, and should have
    idxIEIS = actuallambda_help == 1 & test1 == 3;
    b = find(idxIEIS==1);
    plotIEIS = randsample(length(b),numplots);
    
    % soft worsened, and should have
    idxWEWS = actuallambda_help == 0 & test1 == 0;
    c = find(idxWEWS==1);
    plotWEWS = randsample(length(c),numplots);
    
    % soft worsened, but shouldn't have
    idxWEIS = actuallambda_help == 0 & test1 == 3;    
    d = find(idxWEIS==1);
    plotWEIS = randsample(length(d),numplots);
    
    figure; hold on;
    for j = 1:numplots        
        plot(zh_lambda(a(plotIEWS(j))),j,'bx');
        plot(zkh_lambda(a(plotIEWS(j))),j,'bs');
        plot(zks_lambda(a(plotIEWS(j))),j,'bo');
        for k = 1:3
            plot(zslocal{a(plotIEWS(j))}(k),j,'b.');           
        end        
    end
    title(sprintf('examples where good soft est, bad lambda1s, %d',soft_years(i)));
    ylim([0 numplots+1]);
    legend('obs','hard est','soft est','lam1','lam1','lam1');
    print(gcf,'-painters','-dpng','-r600',sprintf('catIEWS_%d.png',soft_years(i)));
    
    figure; hold on;
    for j = 1:numplots        
        plot(zh_lambda(b(plotIEIS(j))),j,'bx');
        plot(zkh_lambda(b(plotIEIS(j))),j,'bs');
        plot(zks_lambda(b(plotIEIS(j))),j,'bo');
        for k = 1:3
            plot(zslocal{b(plotIEIS(j))}(k),j,'b.');           
        end        
    end
    title(sprintf('examples where good soft est, good lambda1s, %d',soft_years(i)));
    ylim([0 numplots+1]);
    legend('obs','hard est','soft est','lam1','lam1','lam1');
    print(gcf,'-painters','-dpng','-r600',sprintf('catIEIS_%d.png',soft_years(i)));
    
    figure; hold on;
    for j = 1:numplots        
        plot(zh_lambda(c(plotWEWS(j))),j,'bx');
        plot(zkh_lambda(c(plotWEWS(j))),j,'bs');
        plot(zks_lambda(c(plotWEWS(j))),j,'bo');
        for k = 1:3
            plot(zslocal{c(plotWEWS(j))}(k),j,'b.');           
        end        
    end
    title(sprintf('examples where bad soft est, bad lambda1s, %d',soft_years(i)));
    ylim([0 numplots+1]);
    legend('obs','hard est','soft est','lam1','lam1','lam1');
    print(gcf,'-painters','-dpng','-r600',sprintf('catWEWS_%d.png',soft_years(i)));
    
    figure; hold on;
    for j = 1:numplots        
        plot(zh_lambda(d(plotWEIS(j))),j,'bx');
        plot(zkh_lambda(d(plotWEIS(j))),j,'bs');
        plot(zks_lambda(d(plotWEIS(j))),j,'bo');
        for k = 1:3
            plot(zslocal{d(plotWEIS(j))}(k),j,'b.');           
        end        
    end
    title(sprintf('examples where bad soft est, good lambda1s, %d',soft_years(i)));
    ylim([0 numplots+1]);
    legend('obs','hard est','soft est','lam1','lam1','lam1');
    print(gcf,'-painters','-dpng','-r600',sprintf('catWEIS_%d.png',soft_years(i)));
    
end

end