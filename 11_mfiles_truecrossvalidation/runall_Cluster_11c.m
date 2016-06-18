function [] = runall_Cluster_11c(FOLDIDX)
% this function will run the covariance model for the true cross validation

if nargin < 1, FOLDIDX = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../11_mfiles_truecrossvalidation

%%% experimental covariance %%%

% experimental covariance
smoothingParam = [900000 300000 100 50];     
load(sprintf('../matfiles/meanTrend_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam), ...
    'ch_nonXval','zh_nonXval','ch_Xval','zh_Xval','mI_nonXval');
zh_mtr = zh_nonXval - mI_nonXval; % mtr = mean trend removed
[Z_mtr,cMS,tME,nanratio]=valstv2stg(ch_nonXval,zh_mtr);

% in order to come up with rLag, I need to know what the spatial
% distances are and to equal spacing in a log space
DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));
DME = abs(bsxfun(@minus,tME,tME'))';

rLag = [0 prctile(unique(DMS(:)),[0.25 0.5 0.75 1 1.5 2:10 12.5 15:5:50])];
rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];
tLag = [0:10 15:5:40 50:25:150];
tTol = [0 repmat(0.5,1,10) repmat(2.5,1,6) repmat(25,1,5)];

[Crtest nprtest]=stcov(Z_mtr,cMS,tME,Z_mtr,cMS,tME,rLag,rTol,0,0);
[Cttest npttest]=stcov(Z_mtr,cMS,tME,Z_mtr,cMS,tME,0,0,tLag,tTol);

% diplaying results
figure; hold on;
plot(rLag,Crtest,'bo');
figure; hold on;
plot(tLag,Cttest,'bo');

%%% model fitting %%%

% automatic fit for covariance model, joint exp/exp
Cr1 = Crtest(1);
len = length(Crtest);
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
aicval = len*log(2*pi/len) + len + 2 + len*log(mean(output.residuals.^2)) + 2*p;

% display results
disp(f);
disp([gof.rsquare aicval]);

%%% save results %%%

save(sprintf('../matfiles/covmod_true10fold_%d.mat',FOLDIDX), ...
    'Z_mtr','cMS','tME','rLag','rTol','tLag','tTol','Crtest','nprtest','Cttest','npttest', ...
    'f','gof','rLag','Crtest','output','aicval','s','g','p');

end