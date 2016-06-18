function [] = CovModelConstant()
% this function calculate the experimental covariance of the PM2.5 data
% after removing a constant mean trend. Much of the code here is similar to
% the code in tempCovMod.m.

% loading all the data
for i = 1999:2010
    load(sprintf('../matfiles/prepObs_%d.mat',i)); 
    ch{i-1998,1} = coordObs;
    yrs = floor(yrmodaObs./10000);
    mos = floor( (yrmodaObs-10000*yrs) ./ 100 );
    das = yrmodaObs - 10000*yrs - 100*mos;
    cht{i-1998,1} = datenum(yrs,mos,das);
    zh{i-1998,1} = Obs;
end
ch = [cell2mat(ch) cell2mat(cht)];
zh = cell2mat(zh);
mO = mean(zh);
  
zh_mrmvd = zh - mO;
% convert from stv2stg
[Z_mrmvd,cMS,tME,nanratio]=valstv2stg(ch,zh_mrmvd);

% in order to come up with rLag, I need to know what the spatial
% distances are and to equal spacing in a log space
DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));
DME = abs(bsxfun(@minus,tME,tME'))';

% spatial covariance
rLag = [0 prctile(unique(DMS(:)),[0.25 0.5 0.75 1 1.5 2:10 12.5 15:5:50])];
rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];

% temporal covariance
tLag = [0:10 15:5:40 50:25:150];
tTol = [0 repmat(0.5,1,10) repmat(2.5,1,6) repmat(25,1,5)];

% trying it a new way that's mine to compare to stcov
tic
[Crtest nprtest]=stcov(Z_mrmvd,cMS,tME,Z_mrmvd,cMS,tME,rLag,rTol,0,0);
toc
tic
[Cttest npttest]=stcov(Z_mrmvd,cMS,tME,Z_mrmvd,cMS,tME,0,0,tLag,tTol); 
toc

% diplaying results
figure; hold on;
plot(rLag,Crtest,'bo');
figure; hold on;
plot(tLag,Cttest,'bo');

% saving results of experimental covariance
save('../matfiles/expcov_constant.mat', ...
    'Z_mrmvd','cMS','tME','rLag','rTol','tLag','tTol','Crtest','nprtest',...
    'Cttest','npttest','mO');

% instead of doing the exhaustive search of the best covariance fittings,
% we can just assuming that the covariance model has already been chosen

% joint exponential exponential 
Cr1 = Crtest(1);
modstr = 'joint exponential exponential';
s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0,0],...
    'Upper',[1,Inf,Inf,Inf,Inf],'Startpoint',[0.5,1500000,1500000,500,500]);
g = fittype( 'jointexpexp(alp,ar1,ar2,at1,at2,Cr1,x,y)','options',s,...
    'problem',{'Cr1'},'independent',{'x','y'},'dependent',{'z'} );
p = 5;
        
x = [rLag' ; zeros(length(tLag),1)];
y = [zeros(length(rLag),1) ; tLag'];
z = [Crtest ; Cttest'];
[f gof output] = fit([x,y],z,g,'problem',Cr1); 

len = length(Crtest);
aicval = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p;
r2 = gof.rsquare;

save(sprintf('../matfiles/covmod_%s_%s_joint.mat','constant',modstr), ...
	'f','gof','rLag','Crtest','tLag','Cttest','output','aicval','s','g','p');

figure; hold on;
plot(rLag,Crtest,'s','Color','b','LineWidth',2); 
xplot = linspace(rLag(1),rLag(end));
yplot = feval(f,[xplot' zeros(100,1)]);
plot(xplot,yplot);
ylabel(sprintf('Covariance C(r,\\tau=0)'));
set(gca,'XTickLabel',get(gca,'XTick')/1000);
xlabel('spatial lag (km)');
title(sprintf('Covariance for r, %s\n%s: aic=%0.3f, r^{2}=%0.3f',modstr,'constant',aicval,r2));        
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_covModel/covmod_r_%s_constant_joint.pdf',modstr));
    
figure; hold on;
plot(tLag,Cttest,'s','Color','b','LineWidth',2); 
xplot = linspace(tLag(1),tLag(end));
yplot = feval(f,[zeros(100,1) xplot']);
plot(xplot,yplot);
ylabel(sprintf('Covariance C(r=0,\\tau)'));
xlabel('temporal lag (days)');
title(sprintf('Covariance for t, %s\n%s: aic=%0.3f, r^{2}=%0.3f',modstr,'constant',aicval,r2));        
print(gcf,'-painters','-dpdf','-r600',sprintf('../plots_covModel/covmod_t_%s_constant_joint.pdf',modstr));

end