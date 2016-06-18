function [] = runall_Cluster_11a(FOLDIDX)
% the true 10-fold cross validation mean trend

if nargin < 1, FOLDIDX = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../04_mfiles_softdata

load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
smoothingParam = [900000 300000 100 50];
zh = zd;
ch = pd;

% picking the corss-validation locations in space 
rand('seed',0);
cMS = unique(ch(:,1:2),'rows');
len = length(cMS); 
len2 = 1:floor(len/10):len; 
randidx = randperm(len);

% removing the X-val locations/getting the X-val locations
if FOLDIDX == 10
    [locA locB] = ismember(ch(:,1:2),cMS(randidx(len2(FOLDIDX):len2(FOLDIDX+1)),:),'rows');
else
    [locA locB] = ismember(ch(:,1:2),cMS(randidx(len2(FOLDIDX):len2(FOLDIDX+1)-1),:),'rows');
end
ch_nonXval = ch(~locA,:);
zh_nonXval = zh(~locA);
ch_Xval = ch(locA,:);
zh_Xval = zh(locA);

% mean trend calculation
cd ../10_mfiles_newmeantrend
tic
[mI_nonXval]=expKernelSmooth_stv_parallel(ch_nonXval,zh_nonXval,smoothingParam,ch_nonXval);
toc
cd ../11_mfiles_truecrossvalidation 

% saving results
save(sprintf('../matfiles/meanTrend_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam), ...
    'ch_nonXval','zh_nonXval','ch_Xval','zh_Xval','mI_nonXval');

end