function [] = plotCovModelJoint()
% this function will plot the joint covariance models

% all the mean trend models
meanNow = {'short','intermediate','long','very_long'};
markerOrder = {'o' ; 's' ; 'x' ; '*' };
plotcolor = cool(6); plotcolor([1 end],:) = [];
load('../matfiles/covmod_modinfo_joint.mat');
    
% plotting all mean trends for a given model on one plot
for i = 1:length(modstr)
    
    for j = 1:length(meanNow)
   
        % loading data
        load(sprintf('../matfiles/expcov_%s.mat',meanNow{j})); 
        len = length(Crtest);
        Cr1 = Crtest(1);

        cd covariance_models
        % joint powered exponential 
        modstr{j,1} = 'joint powered exponential';
        s{j,1} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
            'Upper',[2,2,Inf,Inf],'Startpoint',[1,1,1500000,5]); 
        g{j,1} = fittype( 'jointpowexp(alpr1,alpt1,ar1,at1,Cr1,x,y)','options',s{j,1},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,1} = 4;

        % joint gaussian 
        modstr{j,2} = 'joint gaussian';
        s{j,2} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
            'Upper',[Inf,Inf],'Startpoint',[300000,5000]);
        g{j,2} = fittype( 'jointgau(ar1,at1,Cr1,x,y)','options',s{j,2},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,2} = 2;    

        % joint exponential
        modstr{j,3} = 'joint exponential';
        s{j,3} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
            'Upper',[Inf,Inf],'Startpoint',[300000,5000]);
        g{j,3} = fittype( 'jointexp(ar1,at1,Cr1,x,y)','options',s{j,3},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,3} = 2;

        % joint powered exponential powered exponential 
        modstr{j,4} = 'joint powered exponential powered exponential';
        s{j,4} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0,0,0,0,0],...
            'Upper',[1,2,2,2,2,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1,1,1,1,100000,100000,5000,5000]);
        g{j,4} = fittype( 'jointpowexppowexp(alp,alpr1,alpt1,alpt2,alpr2,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,4},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,4} = 9;

        % joint gaussian gaussian 
        modstr{j,5} = 'joint gaussian gaussian';
        s{j,5} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1000000,1000000,500,500]);
        g{j,5} = fittype( 'jointgaugau(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,5}, ...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,5} = 5;

        % joint exponential exponential 
        modstr{j,6} = 'joint exponential exponential';
        s{j,6} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1500000,1500000,500,500]);
        g{j,6} = fittype( 'jointexpexp(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,6},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'} );
        p{j,6} = 5;

        % joint exponential powered exponential 
        modstr{j,7} = 'joint exponential powered exponential';
        s{j,7} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0,0,0],...
            'Upper',[1,Inf,2,Inf,Inf,2,Inf],'Startpoint',[0.5,300000,1,300000,5000,1,5000]);
        g{j,7} = fittype('jointexppowexp(alp,ar1,aplr2,ar2,at1,alpt2,at2,Cr1,x,y)','options',s{j,7},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,7} = 7;

        % joint exponential gaussian 
        modstr{j,8} = 'joint exponential gaussian';
        s{j,8} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
            'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,100000,300000,5000,5000]);
        g{j,8} = fittype( 'jointexpgau(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,8}, ...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,8} = 5;   

        % joint powered exponential gaussian space
        modstr{j,9} = 'joint powered exponential gaussian';
        s{j,9} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0,0,0],...
            'Upper',[1,2,2,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1,1,300000,300000,5000,5000]);
        g{j,9} = fittype( 'jointpowexpgau(alp,alpr1,alpt1,ar1,ar2,at1,at2,Cr1,x,y)','options',s{j,9},...
            'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
        p{j,9} = 7;
        
        cd ..

        x = [rLag' ; zeros(length(tLag),1)];
        y = [zeros(length(rLag),1) ; tLag'];
        z = [Crtest ; Cttest'];
        [f gof output] = fit([x,y],z,g{j,i},'problem',Cr1);    

        aicval(j) = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{j,i};
        r2(j) = gof.rsquare;

        figure(1); hold on;
        hlr(j) = plot(rLag,Crtest,markerOrder{j},'Color',plotcolor(j,:),...
            'LineWidth',2); 
        xplot = linspace(rLag(1),rLag(end));
        yplot = feval(f,[xplot' zeros(100,1)]);
        h(j) = plot(xplot,yplot);
        set(h(j),{'Color'},num2cell(plotcolor(j,:),2));

        figure(2); hold on;
        hlt(j) = plot(tLag,Cttest,markerOrder{j},'Color',plotcolor(j,:),...
            'LineWidth',2); 
        xplot = linspace(tLag(1),tLag(end));
        yplot = feval(f,[zeros(100,1) xplot']);
        h(j) = plot(xplot,yplot);
        set(h(j),{'Color'},num2cell(plotcolor(j,:),2));
            
    end

    figure(1);
    ylabel(sprintf('Covariance C(r,\\tau=0)'));
    set(gca,'XTickLabel',get(gca,'XTick')/1000);
    xlabel('spatial lag (km)');
    title(sprintf('Covariance for r, %s\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f',...
        modstr{1,i},meanNow{1},aicval(1),r2(1),meanNow{2},aicval(2),r2(2),...
        meanNow{3},aicval(3),r2(3),meanNow{4},aicval(4),r2(4)));        
    legend(hlr,meanNow,'Location','Best');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_covModel/covmod_r_%s_joint.pdf',...
        modstr{1,i}));
    
    figure(2);
    ylabel(sprintf('Covariance C(r=0,\\tau)'));
    xlabel('temporal lag (days)');
    title(sprintf('Covariance for t, %s\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f\n%s: aic=%0.3f, r^{2}=%0.3f',...
        modstr{1,i},meanNow{1},aicval(1),r2(1),meanNow{2},aicval(2),r2(2),...
        meanNow{3},aicval(3),r2(3),meanNow{4},aicval(4),r2(4)));        
    legend(hlt,meanNow,'Location','Best');
    print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_covModel/covmod_t_%s_joint.pdf',...
        modstr{1,i}));
     
    close all;
          
end

end