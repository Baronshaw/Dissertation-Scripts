function [] = getCovModel()
% this function will see the fit of many different models

% all the mean trend models
meanNow = {'short','intermediate','long','very_long'}; 

% fitting some models to the data, one mean trend model at a time
for i = 1:length(meanNow)
    
    % loading data
    load(sprintf('../matfiles/expcov_%s.mat',meanNow{i})); 
    len = length(Crtest);
    Cr1 = Crtest(1);
    
    cd covariance_models
    % powered exponential family 
    modstr{i,1} = 'powered exponential';
    s{i,1} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
    	'Upper',[2,Inf],'Startpoint',[1,300000]);
    g{i,1} = fittype( @(alp1,ar1,x) Cr1*exp(-(3*x./ar1).^alp1),'options',s{i,1} );
    p{i,1} = 2;
    
    % gaussian 
    modstr{i,2} = 'gaussian';
    s{i,2} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[50000]);
    g{i,2}= fittype( @(ar1,x) Cr1*exp(-(sqrt(3)*x/ar1).^2),'options',s{i,2} );
    p{i,2} = 1;
    
    % matern 
    modstr{i,3} = 'matern';
    s{i,3} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
    	'Upper',[Inf,Inf],'Startpoint',[5,50000]);
    g{i,3} = fittype( @(nu1,ar1,x) ...
        Cr1 .* 2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1),'options',s{i,3} );
    p{i,3} = 2;
    
    % exponential 
    modstr{i,4} = 'exponential';
    s{i,4} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[500000]);
    g{i,4} = fittype( @(ar1,x) Cr1*exp(-3*x./ar1),'options',s{i,4} );
    p{i,4} = 1;
    
    % nugget 
    modstr{i,5} = 'nugget';
    s{i,5} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[1],'Startpoint',[0.5]);
    g{i,5} = fittype( @(ar1,x) Cr1*exp(-3*x./ar1),'options',s{i,5} );
    p{i,5} = 1;
    
    % cubic 
    modstr{i,6} = 'cubic';
    s{i,6} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[50000]);
    g{i,6} = fittype( 'cubic(ar1,Cr1,x)','options',s{i,6},'problem','Cr1');
    p{i,6} = 1;
    
    % spherical 
    modstr{i,7} = 'spherical';
    s{i,7} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[500000]);
    g{i,7} = fittype( 'spherical(ar1,Cr1,x)','options',s{i,7},'problem','Cr1');
    p{i,7} = 1;
    
    % pentaspherical 
    modstr{i,8} = 'pentaspherical';
    s{i,8} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[500000]);
    g{i,8} = fittype( 'pentaspherical(ar1,Cr1,x)','options',s{i,8},'problem','Cr1');
    p{i,8} = 1;
    
    % cosine hole 
    modstr{i,9} = 'cosine hole';
    s{i,9} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[50000]);
    g{i,9} = fittype( @(ar1,x) Cr1*cos(pi*x./ar1),'options',s{i,9} );
    p{i,9} = 1;
    
    % sine hole 
    modstr{i,10} = 'sine hole';
    s{i,10} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[5000]);
    g{i,10} = fittype( 'sinehole(ar1,Cr1,x)','options',s{i,10},'problem','Cr1');
    p{i,10} = 1;
    
    % exponential exponential 
    modstr{i,11} = 'exponential exponential';
    s{i,11} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
    	'Upper',[1 Inf Inf],'Startpoint',[0.5 100000 300000]);
    g{i,11} = fittype( @(alp,ar1,ar2,x) Cr1.*((alp).*exp(-3*x./ar1)+(1-alp).*exp(-3*x./ar2)),'options',s{i,11} );
    p{i,11} = 3;
    
    % powered exponential powered exponential 
    modstr{i,12} = 'powered exponential powered exponential';
    s{i,12} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,2,Inf,2,Inf],'Startpoint',[0.5,1,5000,1,5000]);
    g{i,12} = fittype( @(alp,alp1,ar1,alp2,ar2,x) Cr1*( (alp).*exp(-(3*x./ar1)).^alp1 + (1-alp).*exp(-(3*x./ar2).^alp2 )),...
        'options',s{i,12} );
    p{i,12} = 5;
    
    % gaussian gaussian 
    modstr{i,13} = 'gaussian gaussian';
    s{i,13} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
    	'Upper',[1 Inf Inf],'Startpoint',[0.5 500000 500000]);
    g{i,13}= fittype( @(alp,ar1,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + (1-alp)*exp(-sqrt(3)*x/ar2).^2),'options',s{i,13} );
    p{i,13} = 3;
    
    % matern matern
    modstr{i,14} = 'matern matern';
    s{i,14} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,50,50000,50,50000]);
    g{i,14} = fittype( @(alp,nu1,ar1,nu2,ar2,x) ...
        Cr1 .*( (alp).*(2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1))  + ...
        (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{i,14} );
    p{i,14} = 5;
    
    % exponential powered exponential 
    modstr{i,15} = 'exponential powered exponential';
    s{i,15} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
    	'Upper',[1,Inf,2,Inf],'Startpoint',[0.5,500000,1,500000]);
    g{i,15} = fittype( @(alp,ar1,alp2,ar2,x) Cr1* ( (alp).*exp(-3*x./ar1) + ...
        (1-alp).*exp(-(3*x./ar2).^alp2) ),'options',s{i,15} );
    p{i,15} = 4;
    
    % exponential gaussian 
    modstr{i,16} = 'exponential gaussian';
    s{i,16} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
    	'Upper',[1 Inf Inf],'Startpoint',[0.5 100000 300000]);
    g{i,16} = fittype( @(alp,ar1,ar2,x) Cr1.*((alp).*exp(-3*x./ar1)+(1-alp).*exp(-3*x./ar2).^2),'options',s{i,16} );
    p{i,16} = 3;
    
    % exponential matern 
    modstr{i,17} = 'exponential matern';
    s{i,17} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf],'Startpoint',[0.5,50000,2,50000]);
    g{i,17} = fittype( @(alp,ar1,nu2,ar2,x) Cr1* ( (alp).*exp(-3*x./ar1) + ...
        (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{i,17} );
    p{i,17} = 4;
    
    % exponential nugget 
    modstr{i,18} = 'exponential nugget';
    s{i,18} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
    	'Upper',[1 Inf 1],'Startpoint',[0.5 500000 0.5]);
    g{i,18} = fittype( @(alp,ar1,ar2,x) Cr1.*((alp).*exp(-3*x./ar1)+(1-alp).*exp(-3*x./ar2)),'options',s{i,18} );
    p{i,18} = 3;
    
    % powered exponential gaussian 
    modstr{i,19} = 'powered exponential gaussian';
    s{i,19} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
    	'Upper',[1,2,Inf,Inf],'Startpoint',[0.5,1,50000,50000]);
    g{i,19} = fittype( @(alp,alp1,ar1,ar2,x) Cr1*( (alp).*exp(-(3*x./ar1)).^alp1 + (1-alp).*exp(-(sqrt(3)*x/ar2).^2)),...
        'options',s{i,19} );
    p{i,19} = 4;
    
    % powered exponential nugget 
    modstr{i,20} = 'powered exponential nugget';
    s{i,20} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
    	'Upper',[1,2,Inf,1],'Startpoint',[0.5,1,500000,0.5]);
    g{i,20} = fittype( @(alp,alp1,ar1,ar2,x) Cr1*( (alp).*exp(-(3*x./ar1)).^alp1 + (1-alp).*exp(-3*x./ar2)),...
        'options',s{i,20} );
    p{i,20} = 4;
    
    % gaussian matern 
    modstr{i,21} = 'gaussian matern';
    s{i,21} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],...
    	'Upper',[1 Inf Inf Inf],'Startpoint',[0.5 50000 20 50000]);
    g{i,21}= fittype( @(alp,ar1,nu2,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + ...
        (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{i,21} );
    p{i,21} = 4;
    
    % gaussian nugget  
    modstr{i,22} = 'gaussian nugget';
    s{i,22} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],...
    	'Upper',[1 Inf 1],'Startpoint',[0.5 500000 0.5]);
    g{i,22}= fittype( @(alp,ar1,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + (1-alp)*exp(-sqrt(3)*x/ar2)),'options',s{i,22} );
    p{i,22} = 3;    
    
    % matern nugget 
    modstr{i,23} = 'matern nugget';
    s{i,23} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
    	'Upper',[1,Inf,Inf,1],'Startpoint',[0.5,5,50000,0.5]);
    g{i,23} = fittype( @(alp,nu1,ar1,ar2,x) ...
        Cr1 .*( (alp).*(2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1))  + ...
        (1-alp)*exp(-sqrt(3)*x/ar2) ),'options',s{i,23} );
    p{i,23} = 4;
    
    % joint exponential 
    modstr{i,24} = 'joint exponential';
    s{i,24} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,50000,50000,500,500]);
    g{i,24} = fittype( 'jointexpexp(apl,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,24},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'},'coefficients',{'apl','ar1','ar2','at1','at2'} );
    p{i,24} = 3;
    
    cd ..
    
    save('../matfiles/covmod_modinfo.mat','modstr','s','g','p');
    
    % in space
    for j = 1:length(modstr) 
                                
        disp(j);
        if j == 3 | j==14 | j==17 | j==21 | j==23, [f gof output] = fit( rLag(2:end)', Crtest(2:end), g{i,j} );
        elseif j==6 | j==7 | j==8 | j==10 
            [f gof output] = fit( rLag', Crtest, g{i,j}, 'problem', Cr1 );
        elseif j == 24
            x = [rLag' ; zeros(length(tLag),1)];
            y = [zeros(length(rLag),1) ; tLag'];
            z = [Crtest ; Cttest'];
            [f gof output] = fit([x,y],z,g{i,j},'problem',Cr1);
        else
            [f gof output] = fit( rLag', Crtest, g{i,j} );
        end
        aicval = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{i,j};
        
        if j < length(modstr) % because you can't disp jointexpexp
            figure; hold on;
            plot(rLag,Crtest,'bo');
            plot(f,'r');
            ylabel(sprintf('Covariance C(r,\\tau=0)'));
            xlabel('spatial lag (meters)');
            title(sprintf('Covariance for r, %s, %s, r^{2}=%0.3f',meanNow{i},modstr{i,j},gof.rsquare));
            h = gcf;
            print(h,'-painters','-dpdf','-r600',sprintf('plots_covModel/covmod_r_%s_%s.pdf',...
                meanNow{i},modstr{i,j}));
            sr = s{i,j}; gr = g{i,j}; pr = p{i,j};
        end
        
        save(sprintf('../matfiles/covmod_test_r_%s_%s.mat',meanNow{i},modstr{i,j}), ...
            'f','gof','rLag','Crtest','output','aicval','sr','gr','pr');
        disp(modstr{i,j});
        disp(f);
        disp([gof.rsquare aicval]);
        
    end
    
    % update some models for time
    
    % gaussian 
    modstr{i,2} = 'gaussian';
    s{i,2} = fitoptions('Method','NonlinearLeastSquares','Lower',[0],...
    	'Upper',[Inf],'Startpoint',[5000]);
    g{i,2}= fittype( @(ar1,x) Cr1*exp(-(sqrt(3)*x/ar1).^2),'options',s{i,2} );
    p{i,2} = 1;
    
    % matern matern
    modstr{i,14} = 'matern matern';
    s{i,14} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,50,500,50,500]);
    g{i,14} = fittype( @(alp,nu1,ar1,nu2,ar2,x) ...
        Cr1 .*( (alp).*(2.^(1-nu1)./gamma(nu1) .* (x./ar1).^nu1 .* besselk(nu1,x./ar1))  + ...
        (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{i,14} );
    p{i,14} = 5;
    
    % gaussian matern 
    modstr{i,21} = 'gaussian matern';
    s{i,21} = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],...
    	'Upper',[1 Inf Inf Inf],'Startpoint',[0.5 5000 20 5000]);
    g{i,21}= fittype( @(alp,ar1,nu2,ar2,x) Cr1*( alp*exp(-sqrt(3)*x/ar1).^2 + ...
        (1-alp).*(2.^(1-nu2)./gamma(nu2) .* (x./ar2).^nu2 .* besselk(nu2,x./ar2) ) ),'options',s{i,21} );
    p{i,21} = 4;
    
    % joint exponential 
    modstr{i,24} = 'joint exponential';
    s{i,24} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,5000,5000,500,500]);
    g{i,24} = fittype( 'jointexpexp(apl,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,24},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'},'coefficients',{'apl','ar1','ar2','at1','at2'} );
    p{i,24} = 3;
    
    % in time
    for j = 1:length(modstr)-1 % '-1' because can't disp jointexpexp
        disp(j);
        if j == 3 | j==14 | j==17 | j==21 | j==23, [f gof output] = fit( tLag(2:end)', Cttest(2:end)', g{i,j} );
        elseif j==6 | j==7 | j==8 | j==10, 
            [f gof output] = fit( tLag', Cttest', g{i,j}, 'problem', Cr1 );
        elseif j == 24
            x = [rLag' ; zeros(length(tLag),1)];
            y = [zeros(length(rLag),1) ; tLag'];
            z = [Crtest ; Cttest'];
            [f gof output] = fit([x,y],z,g{i,j},'problem',Cr1);
        else
            [f gof output] = fit( tLag', Cttest', g{i,j} );
        end
        aicval = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{i,j};
        
        if j < length(modstr) % because you can't disp jointexpexp
            figure; hold on;
            plot(tLag,Cttest,'bo');
            plot(f,'r');
            ylabel(sprintf('Covariance C(r=0,\\tau)'));
            xlabel('temporal lag (days)');
            title(sprintf('Covariance for t, %s, %s, r^{2}=%0.3f',meanNow{i},modstr{i,j},gof.rsquare));
            h = gcf;
            print(h,'-painters','-dpdf','-r600',sprintf('plots_covModel/covmod_t_%s_%s.pdf',...
                meanNow{i},modstr{i,j}));
        end
        
        save(sprintf('../matfiles/covmod_test_t_%s_%s.mat',meanNow{i},modstr{i,j}), ...
            'f','gof','tLag','Cttest','output','aicval');
        disp(modstr{i,j});
        disp(f);
        disp([gof.rsquare aicval]);
        
    end
    
    close all;
    
end

end