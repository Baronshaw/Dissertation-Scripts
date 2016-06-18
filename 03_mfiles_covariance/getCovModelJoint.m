function [] = getCovModelJoint()
% this function will display the experimental covariance models in 
% space/time jointly 

% all the mean trend models
meanNow = {'short','intermediate','long','very_long'}; 

% fitting some models to the data, one mean trend model at a time
for i = 1:length(meanNow)
    
    % loading data
    load(sprintf('../matfiles/expcov_%s.mat',meanNow{i})); 
    len = length(Crtest);
    Cr1 = Crtest(1);
    
    cd covariance_models
    % joint powered exponential 
    modstr{i,1} = 'joint powered exponential';
    s{i,1} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],...
    	'Upper',[2,2,Inf,Inf],'Startpoint',[1,1,1500000,5]); 
    g{i,1} = fittype( 'jointpowexp(alpr1,alpt1,ar1,at1,Cr1,x,y)','options',s{i,1},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,1} = 4;
    
    % joint gaussian 
    modstr{i,2} = 'joint gaussian';
    s{i,2} = fitoptions('Method','NonlinearLeastSquares','Lower',[00],...
    	'Upper',[Inf,Inf],'Startpoint',[300000,5000]);
    g{i,2} = fittype( 'jointgau(ar1,at1,Cr1,x,y)','options',s{i,2},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,2} = 2;    
    
    % joint exponential
    modstr{i,3} = 'joint exponential';
    s{i,3} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],...
    	'Upper',[Inf,Inf],'Startpoint',[300000,5000]);
    g{i,3} = fittype( 'jointexp(ar1,at1,Cr1,x,y)','options',s{i,3},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,3} = 2;
    
    % joint powered exponential powered exponential 
    modstr{i,4} = 'joint powered exponential powered exponential';
    s{i,4} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0,0,0,0,0],...
    	'Upper',[1,2,2,2,2,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1,1,1,1,100000,100000,5000,5000]);
    g{i,4} = fittype( 'jointpowexppowexp(alp,alpr1,alpt1,alpt2,alpr2,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,4},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,4} = 9;
    
    % joint gaussian gaussian 
    modstr{i,5} = 'joint gaussian gaussian';
    s{i,5} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1000000,1000000,500,500]);
    g{i,5} = fittype( 'jointgaugau(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,5}, ...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,5} = 5;
    
    % joint exponential exponential 
    modstr{i,6} = 'joint exponential exponential';
    s{i,6} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1500000,1500000,500,500]);
    g{i,6} = fittype( 'jointexpexp(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,6},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'} );
    p{i,6} = 5;
    
    % joint exponential powered exponential 
    modstr{i,7} = 'joint exponential powered exponential';
    s{i,7} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0,0,0],...
    	'Upper',[1,Inf,2,Inf,Inf,2,Inf],'Startpoint',[0.5,300000,1,300000,5000,1,5000]);
    g{i,7} = fittype('jointexppowexp(alp,ar1,aplr2,ar2,at1,alpt2,at2,Cr1,x,y)','options',s{i,7},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,7} = 7;
    
    % joint exponential gaussian 
    modstr{i,8} = 'joint exponential gaussian';
    s{i,8} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    	'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,100000,300000,5000,5000]);
    g{i,8} = fittype( 'jointexpgau(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,8}, ...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,8} = 5;   
    
    % joint powered exponential gaussian space
    modstr{i,9} = 'joint powered exponential gaussian';
    s{i,9} = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0,0,0],...
    	'Upper',[1,2,2,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1,1,300000,300000,5000,5000]);
    g{i,9} = fittype( 'jointpowexpgau(alp,alpr1,alpt1,ar1,ar2,at1,at2,Cr1,x,y)','options',s{i,9},...
        'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'});
    p{i,9} = 7;

    cd ..
    
    save('../matfiles/covmod_modinfo_joint.mat','modstr','s','g','p');
    
    % calculate the covariance jointly, show space first
    for j = 1:length(modstr)
                                
        disp(j);
        x = [rLag' ; zeros(length(tLag),1)];
        y = [zeros(length(rLag),1) ; tLag'];
        z = [Crtest ; Cttest'];
        [f gof output] = fit([x,y],z,g{i,j},'problem',Cr1);

        aicval = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{i,j};
        save(sprintf('../matfiles/covmod_r_%s_%s_joint.mat',meanNow{i},modstr{i,j}), ...
            'f','gof','rLag','Crtest','output','aicval','s','g','p');
        disp(modstr{i,j});
        disp(f);
        disp([gof.rsquare aicval]);
        
    end
    
    % calculate the covariance jointly, show time second
    for j = 1:length(modstr)
                                
        disp(j);
        x = [rLag' ; zeros(length(tLag),1)];
        y = [zeros(length(rLag),1) ; tLag'];
        z = [Crtest ; Cttest'];
        [f gof output] = fit([x,y],z,g{i,j},'problem',Cr1);

        aicval = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p{i,j};
        save(sprintf('../matfiles/covmod_t_%s_%s_joint.mat',meanNow{i},modstr{i,j}), ...
            'f','gof','rLag','Crtest','output','aicval','s','g','p');
        disp(modstr{i,j});
        disp(f);
        disp([gof.rsquare aicval]);
        
    end
    
end

end