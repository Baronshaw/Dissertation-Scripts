function [] = plotCovModel()
% this function will plot all the covariance models

% all the mean trend models
meanNow = {'short','intermediate','long','very_long'}; 
markerOrder = {'o' ; 's' ; 'x' ; '*' };
plotcolor = cool(6); plotcolor([1 end],:) = [];

% first, display of all the spatial experimental covariances
figure; hold on;
markerOrder = {'o' ; 's' ; 'x' ; '*' };
plotcolor = cool(6);
for i = 1:length(meanNow)
    load(sprintf('../matfiles/expcov_%s.mat',meanNow{i}));
    plot(rLag,Crtest,markerOrder{i},'Color',plotcolor(i+1,:),...
        'LineWidth',2);
end
legend(meanNow);
set(gca,'XTickLabel',get(gca,'XTick')/1000);
xlabel('spatial lag (km)');
ylabel('spatial covariance');
title(sprintf('spatial covariance (r,\\tau=0)'));
h = gcf;
print(h,'-painters','-dpdf','-r600',sprintf('../plots_covModel/expcov_r.pdf'));

% second, display of all the temporal experimental covariances
figure; hold on;
markerOrder = {'o' ; 's' ; 'x' ; '*' };
plotcolor = cool(6);
for i = 1:length(meanNow)
    load(sprintf('../matfiles/expcov_%s.mat',meanNow{i}));
    plot(tLag,Cttest,markerOrder{i},'Color',plotcolor(i+1,:),...
        'LineWidth',2);
end
legend(meanNow);
xlabel('temporal lag (days)');
ylabel('temporal covariance');
title(sprintf('temporal covariance (r=0,\\tau)'));
h = gcf;
print(h,'-painters','-dpdf','-r600',sprintf('../plots_covModel/expcov_t.pdf'));

load('../matfiles/covmod_modinfo.mat');
clear s g p

% plotting all mean trends for a given model on one plot
for i = 1:length(modstr)-1 % '-1' because can't disp jointexpexp
    
    % in space
    figure; hold on;    
    for j = 1:length(meanNow)
        
        load(sprintf('../matfiles/expcov_%s.mat',meanNow{j}));
        Cr1 = Crtest(1);
        len = length(Crtest);
        
        cd covariance_models
        % powered exponential family 
        modstr{j,1} = 'powered exponential';
        s{j,1} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
            'Upper',[2,Inf],'Startpoint',[1,300000]);
        g{j,1} = fittype( @(alp1,ar1,x) Cr1*exp(-(3*x./ar1).^alp1),'options',s{j,1} );
        p{j,1} = 2;

        % gaussian 
        modstr{j,2} = 'gaussian';
        s{j,2} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[5000]);
        g{j,2}= fittype( @(ar1,x) Cr1*exp(-(sqrt(3)*x/ar1).^2),'options',s{j,2} );
        p{j,2} = 1;

        % matern 
        modstr{j,3} = 'matern';
        s{j,3} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
            'Upper',[Inf,Inf],'Startpoint',[5,50000]);
        g{j,3} = fittype( @(nu1,ar1,x) ...
            Cr1 .* 2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1),'options',s{j,3} );
        p{j,3} = 2;

        % exponential 
        modstr{j,4} = 'exponential';
        s{j,4} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[500000]);
        g{j,4} = fittype( @(ar1,x) Cr1*exp(-3*x./ar1),'options',s{j,4} );
        p{j,4} = 1;

        % nugget 
        modstr{j,5} = 'nugget';
        s{j,5} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[1],'Startpoint',[0.5]);
        g{j,5} = fittype( @(ar1,x) Cr1*exp(-3*x./ar1),'options',s{j,5} );
        p{j,5} = 1;

        % cubic 
        modstr{j,6} = 'cubic';
        s{j,6} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[50000]);
        g{j,6} = fittype( 'cubic(ar1,Cr1,x)','options',s{j,6},'problem','Cr1');
        p{j,6} = 1;

        % spherical 
        modstr{j,7} = 'spherical';
        s{j,7} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[500000]);
        g{j,7} = fittype( 'spherical(ar1,Cr1,x)','options',s{j,7},'problem','Cr1');
        p{j,7} = 1;

        % pentaspherical 
        modstr{j,8} = 'pentaspherical';
        s{j,8} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[500000]);
        g{j,8} = fittype( 'pentaspherical(ar1,Cr1,x)','options',s{j,8},'problem','Cr1');
        p{j,8} = 1;

        % cosine hole 
        modstr{j,9} = 'cosine hole';
        s{j,9} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[50000]);
        g{j,9} = fittype( @(ar1,x) Cr1*cos(pi*x./ar1),'options',s{j,9} );
        p{j,9} = 1;

        % sine hole 
        modstr{j,10} = 'sine hole';
        s{j,10} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[5000]);
        g{j,10} = fittype( 'sinehole(ar1,Cr1,x)','options',s{j,10},'problem','Cr1');
        p{j,10} = 1;

        % exponential exponential 
        modstr{j,11} = 'exponential exponential';
        s{j,11} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
            'Upper',[1 Inf Inf],'Startpoint',[0.5 100000 300000]);
        g{j,11} = fittype( @(alp,ar1,ar2,x) Cr1.*((alp).*exp(-3*x./ar1)+(1-alp).*exp(-3*x./ar2)),'options',s{j,11} );
        p{j,11} = 3;

        % powered exponential powered exponential 
        modstr{j,12} = 'powered exponential powered exponential';
        s{j,12} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,2,Inf,2,Inf],'Startpoint',[0.5,1,5000,1,5000]);
        g{j,12} = fittype( @(alp,alp1,ar1,alp2,ar2,x) Cr1*( (alp).*exp(-(3*x./ar1)).^alp1 + (1-alp).*exp(-(3*x./ar2).^alp2 )),...
            'options',s{j,12} );
        p{j,12} = 5;

        % gaussian gaussian 
        modstr{j,13} = 'gaussian gaussian';
        s{j,13} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
            'Upper',[1 Inf Inf],'Startpoint',[0.5 50000 50000]);
        g{j,13}= fittype( @(alp,ar1,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + (1-alp)*exp(-sqrt(3)*x/ar2).^2),'options',s{j,13} );
        p{j,13} = 3;

        % matern matern
        modstr{j,14} = 'matern matern';
        s{j,14} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,50,50000,50,50000]);
        g{j,14} = fittype( @(alp,nu1,ar1,nu2,ar2,x) ...
            Cr1 .*( (alp).*(2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1))  + ...
            (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{j,14} );
        p{j,14} = 5;

        % exponential powered exponential 
        modstr{j,15} = 'exponential powered exponential';
        s{j,15} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
            'Upper',[1,Inf,2,Inf],'Startpoint',[0.5,500000,1,500000]);
        g{j,15} = fittype( @(alp,ar1,alp2,ar2,x) Cr1* ( (alp).*exp(-3*x./ar1) + ...
            (1-alp).*exp(-(3*x./ar2).^alp2) ),'options',s{j,15} );
        p{j,15} = 4;

        % exponential gaussian 
        modstr{j,16} = 'exponential gaussian';
        s{j,16} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
            'Upper',[1 Inf Inf],'Startpoint',[0.5 100000 300000]);
        g{j,16} = fittype( @(alp,ar1,ar2,x) Cr1.*((alp).*exp(-3*x./ar1)+(1-alp).*exp(-3*x./ar2).^2),'options',s{j,16} );
        p{j,16} = 3;

        % exponential matern 
        modstr{j,17} = 'exponential matern';
        s{j,17} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf],'Startpoint',[0.5,50000,2,50000]);
        g{j,17} = fittype( @(alp,ar1,nu2,ar2,x) Cr1* ( (alp).*exp(-3*x./ar1) + ...
            (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{j,17} );
        p{j,17} = 4;

        % exponential nugget 
        modstr{j,18} = 'exponential nugget';
        s{j,18} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
            'Upper',[1 Inf 1],'Startpoint',[0.5 500000 0.5]);
        g{j,18} = fittype( @(alp,ar1,ar2,x) Cr1.*((alp).*exp(-3*x./ar1)+(1-alp).*exp(-3*x./ar2)),'options',s{j,18} );
        p{j,18} = 3;

        % powered exponential gaussian 
        modstr{j,19} = 'powered exponential gaussian';
        s{j,19} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
            'Upper',[1,2,Inf,Inf],'Startpoint',[0.5,1,5000,5000]);
        g{j,19} = fittype( @(alp,alp1,ar1,ar2,x) Cr1*( (alp).*exp(-(3*x./ar1)).^alp1 + (1-alp).*exp(-(sqrt(3)*x/ar2).^2)),...
            'options',s{j,19} );
        p{j,19} = 4;

        % powered exponential nugget 
        modstr{j,20} = 'powered exponential nugget';
        s{j,20} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
            'Upper',[1,2,Inf,1],'Startpoint',[0.5,1,500000,0.5]);
        g{j,20} = fittype( @(alp,alp1,ar1,ar2,x) Cr1*( (alp).*exp(-(3*x./ar1)).^alp1 + (1-alp).*exp(-3*x./ar2)),...
            'options',s{j,20} );
        p{j,20} = 4;

        % gaussian matern 
        modstr{j,21} = 'gaussian matern';
        s{j,21} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],...
            'Upper',[1 Inf Inf Inf],'Startpoint',[0.5 50000 20 50000]);
        g{j,21}= fittype( @(alp,ar1,nu2,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + ...
            (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{j,21} );
        p{j,21} = 4;

        % gaussian nugget  
        modstr{j,22} = 'gaussian nugget';
        s{j,22} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
            'Upper',[1 Inf 1],'Startpoint',[0.5 500000 0.5]);
        g{j,22}= fittype( @(alp,ar1,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + (1-alp)*exp(-sqrt(3)*x/ar2)),'options',s{j,22} );
        p{j,22} = 3;    

        % matern nugget 
        modstr{j,23} = 'matern nugget';
        s{j,23} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
            'Upper',[1,Inf,Inf,1],'Startpoint',[0.5,5,50000,0.5]);
        g{j,23} = fittype( @(alp,nu1,ar1,ar2,x) ...
            Cr1 .*( (alp).*(2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1))  + ...
            (1-alp)*exp(-sqrt(3)*x/ar2) ),'options',s{j,23} );
        p{j,23} = 4;

        % joint exponential 
        modstr{j,24} = 'joint exponential';
        s{j,24} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,500000,500000,500,500]);
        g{j,24} = fittype( 'jointexpexp(apl,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,24},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'},'coefficients',{'apl','ar1','ar2','at1','at2'} );
        p{j,24} = 3;

        cd ..
        
        if i == 3 | i==14 | i==17 | i==21 | i==23, 
            [f gof output] = fit( rLag(2:end)', Crtest(2:end), g{j,i} );
        elseif i==6 | i==7 | i==8 | i==10, 
            [f gof output] = fit( rLag', Crtest, g{j,i}, 'problem', Cr1 );
        else
            [f gof output] = fit( rLag', Crtest, g{j,i} );
        end
        
        aicval(j) = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{j,i};
        r2(j) = gof.rsquare;
                
        hl(j) = plot(rLag,Crtest,markerOrder{j},'Color',plotcolor(j,:),...
            'LineWidth',2);        
        h(j) = plot(f);
        set(h(j),{'Color'},num2cell(plotcolor(j,:),2));

    end

    ylabel(sprintf('Covariance C(r,\\tau=0)'));
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    xlabel('spatial lag (km)');
    title(sprintf('Covariance for r, %s\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f',...
        modstr{1,i},meanNow{1},aicval(1),r2(1),meanNow{2},aicval(2),r2(2),...
        meanNow{3},aicval(3),r2(3),meanNow{4},aicval(4),r2(4)));        
    legend(hl,meanNow,'Location','Best');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_covModel/covmod_r_%s.pdf',...
        modstr{1,i}));
        
    % in time 
    figure; hold on;
    
    for j = 1:length(meanNow)
        
        load(sprintf('../matfiles/expcov_%s.mat',meanNow{j}));
        Cr1 = Cttest(1); 
        len = length(Cttest);
    
        % update some models for time

        % gaussian 
        modstr{j,2} = 'gaussian';
        s{j,2} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
            'Upper',[Inf],'Startpoint',[5000]);
        g{j,2}= fittype( @(ar1,x) Cr1*exp(-(sqrt(3)*x/ar1).^2),'options',s{j,2} );
        p{j,2} = 1;

        % matern matern
        modstr{j,14} = 'matern matern';
        s{j,14} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,50,500,50,500]);
        g{j,14} = fittype( @(alp,nu1,ar1,nu2,ar2,x) ...
            Cr1 .*( (alp).*(2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1))  + ...
            (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{j,14} );
        p{j,14} = 5;

        % gaussian matern 
        modstr{j,21} = 'gaussian matern';
        s{j,21} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],...
            'Upper',[1 Inf Inf Inf],'Startpoint',[0.5 5000 20 5000]);
        g{j,21}= fittype( @(alp,ar1,nu2,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + ...
            (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{j,21} );
        p{j,21} = 4;

        % joint exponential 
        modstr{j,24} = 'joint exponential';
        s{j,24} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,5000,5000,500,500]);
        g{j,24} = fittype( 'jointexpexp(apl,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,24},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'},'coefficients',{'apl','ar1','ar2','at1','at2'} );
        p{j,24} = 3;

        if i == 3 | i==14 | i==17 | i==21 | i==23, [f gof output] = fit( tLag(2:end)', Cttest(2:end)', g{j,i} );
        elseif i==6 | i==7 | i==8 | i==10, 
            [f gof output] = fit( tLag', Cttest', g{j,i}, 'problem', Cr1 );
        else
            [f gof output] = fit( tLag', Cttest', g{j,i} );
        end
        
        aicval(j) = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{j,i};
        r2(j) = gof.rsquare;
                
        hl(j) = plot(tLag,Cttest,markerOrder{j},'Color',plotcolor(j,:),...
            'LineWidth',2);        
        h(j) = plot(f);
        set(h(j),{'Color'},num2cell(plotcolor(j,:),2));

    end

    ylabel(sprintf('Covariance C(r=0,\\tau)'));
    xlabel('temporal lag (days)');
    title(sprintf('Covariance for t, %s\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f',...
        modstr{1,i},meanNow{1},aicval(1),r2(1),meanNow{2},aicval(2),r2(2),...
        meanNow{3},aicval(3),r2(3),meanNow{4},aicval(4),r2(4)));        
    legend(hl,meanNow,'Location','Best');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_covModel/covmod_t_%s.pdf',...
        modstr{1,i}));
        
end

close all;

end