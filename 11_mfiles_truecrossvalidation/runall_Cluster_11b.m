function [] = runall_Cluster_11b(FOLDIDX)
% the true 10-fold cross validation mean trend for soft data

if nargin < 1, FOLDIDX = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../04_mfiles_softdata

% loading all hard data
load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,300000,100,50));
smoothingParam = [2*900000 300000 100 50]; % '2' added later
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

% load soft data locations
CTMyears = [2001 2002 2005 2006 2006 2007];
pI = cell(length(CTMyears),1);
for i = 1:length(CTMyears);    
    if CTMyears(i) == 2006 && i == 5
        load(sprintf('../matfiles/prepCTM_%dsum.mat',CTMyears(i)));  
    else
        load(sprintf('../matfiles/prepCTM_%d.mat',CTMyears(i)));  
    end    
    pI{i,1} = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];   
end
pI = cell2mat(pI);

% mean trend calculation
cd ../10_mfiles_newmeantrend
tic
[mIsoft_nonXval]=expKernelSmooth_stv_parallel(ch_nonXval,zh_nonXval,smoothingParam,pI);
toc
cd ../11_mfiles_truecrossvalidation 

% saving results
save(sprintf('../matfiles/meanTrendSoft_true10fold_%d_%d_%d_%d_%d.mat',FOLDIDX,smoothingParam), ...
    'ch_nonXval','zh_nonXval','ch_Xval','zh_Xval','pI','mIsoft_nonXval');

end