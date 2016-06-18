function [] = patch_MP_cmaq_grid(yrz)
% this function will find the location/day outside the range of the SCurve
% for fixed modeled values and performance an interpolation
% Notes: interpolation was only done in 2D
%        'flagged' mean there was a problematic extrapolation

% bsub -x -q day -n 9 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "patch_MP_cmaq_grid" -logfile "patch_MP_cmaq_grid.out"

if nargin < 1, yrz = 2001; end

% loop through each day of the year
yrNday = datevec(datenum(yrz,1,1):datenum(yrz,12,31));
yrNday = yrNday(:,1).*10^4 + yrNday(:,2).*10^2 + yrNday(:,3);

% variables 
load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(1)));
r = length(CTMlocs); c = length(yrNday);
maxGrid = NaN*ones(r,c); minGrid = NaN*ones(r,c); 

% load soft data files
for i = 1:length(yrNday)
    disp(i);
    load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));
    Mod = valMod;    
    maxGrid(:,i) = nanmax(Mod,[],2);
    minGrid(:,i) = nanmin(Mod,[],2);    
end

% flag values outside of min/max range
[r c] = size(maxGrid);
modplots = 0:5:50;
flagGrid = cell(length(modplots),1);
for i = 1:length(modplots)
    disp(modplots(i));
    flagGrid{i} = NaN*ones(r,c);
    flagGrid{i} = maxGrid < modplots(i);
    disp(sum(flagGrid{i}(:))/(r*c));
end

% reshape data
Xpre = unique(round(CTMlocs(:,1)));
lenx = length(Xpre);
Ypre = unique(round(CTMlocs(:,2)));
leny = length(Ypre);
X = repmat(Xpre',leny,1);
Y = repmat(Ypre,1,lenx);
[dummy sortidx] = ismember([X(:) Y(:)],round(CTMlocs),'rows');
% checking: isequal(round(CTMlocs(sortidx,:)),[X(:) Y(:)]) == 1
[dummy origsortidx] = ismember(round(CTMlocs),[X(:) Y(:)],'rows');
% checking: a = [X(:) Y(:)]; isequal(a(origsortidx,:),round(CTMlocs)) == 1
V = reshape(maxGrid(sortidx,1),[leny lenx]);
% checking: isequal([X(:) Y(:) V(:)],[round(CTMlocs(sortidx,:)) maxGrid(sortidx,1)]) == 1;

% looking at the metrics
load('matfiles/traditional_performance_extendbin.mat');

% measure names/values
valplot = { numGivMod ; mObsGivMod ; mModGivMod ; mBiasGivMod ; nBiasGivMod ; nmBiasGivMod ; fBiasGivMod ; ...
    mErrGivMod ; nErrGivMod ; nmErrGivMod ; fErrGivMod ; RGivMod ; R2GivMod ; sBiasGivMod ; msBiasGivMod ; ...
    rmsBiasGivMod ; nrmsBiasGivMod ; mDsBiasGivMod ; m2DmsBiasGivMod ; s2DmsBiasGivMod ; beta1GivMod ; vObsGivMod ; vModGivMod };
[r c] = size(valplot{1});
[r2 c2] = size(valplot{1}{1});
valplotfixed = cell(length(valplot),1);

for i = 1:length(valplot) % loop through the metrics
    tic
    valplotfixed{i} = cell(r,1);
    
    for j = 1:r % looping through the days
        
        valplotfixed{i}{j} = NaN*ones(r2,c2);
        
        for k = 2:c2 % looping through the fixed modeled values
            disp([i j k]);
            idxflag = flagGrid{k-1}(:,j) == 1;% locating flag
            idxV = reshape(idxflag(sortidx),[leny lenx]); % reshape flag index
            V = reshape(valplot{i}{j}(sortidx,k),[leny lenx]); % reshape metric
            Vflag = V; Vflag(idxV) = NaN; % filling in flagged vals with NaN
            
            % performance the interpolation
            Vfixed = inpaint_nans(Vflag); % this fn is from the file exchange
            
            % reshaping original metric
            dummy = Vfixed(:);
            valplotfixed{i}{j}(:,k) = dummy(origsortidx);
        
        end        
    end
    toc
end

% rename metrics
numGivModF = valplotfixed{1}; 
mObsGivModF = valplotfixed{2};
mModGivModF = valplotfixed{3};
mBiasGivModF = valplotfixed{4}; 
nBiasGivModF = valplotfixed{5}; 
nmBiasGivModF = valplotfixed{6}; 
fBiasGivModF = valplotfixed{7};
mErrGivModF = valplotfixed{8};
nErrGivModF = valplotfixed{9};
nmErrGivModF = valplotfixed{10};
fErrGivModF = valplotfixed{11}; 
RGivModF = valplotfixed{12}; 
R2GivModF = valplotfixed{13};
sBiasGivModF = valplotfixed{14};
msBiasGivModF = valplotfixed{15};
rmsBiasGivModF = valplotfixed{16};
nrmsBiasGivModF = valplotfixed{17};
mDsBiasGivModF = valplotfixed{18};
m2DmsBiasGivModF = valplotfixed{19};
s2DmsBiasGivModF = valplotfixed{20};
beta1GivModF = valplotfixed{21}; 
vObsGivModF = valplotfixed{22}; 
vModGivModF = valplotfixed{23};

% save results
save('matfiles/traditional_performance_extendbinFixed.mat', ...
    'numGivModF','mObsGivModF','mModGivModF','mBiasGivModF','nBiasGivModF','nmBiasGivModF', ...
    'fBiasGivModF','mErrGivModF','nErrGivModF','nmErrGivModF','fErrGivModF','RGivModF', ... 
    'R2GivModF','sBiasGivModF','msBiasGivModF','rmsBiasGivModF','nrmsBiasGivModF','mDsBiasGivModF', ...
    'm2DmsBiasGivModF','s2DmsBiasGivModF','beta1GivModF','vObsGivModF','vModGivModF', ...
    'CTMlocs','yrNday','modplots','maxGrid','minGrid','flagGrid'); 

end