function [] = stat_MP_cmaq_extendbin(yrz)
% this function extends the work of stat_MP_cmaq_gridbin.m by calculating
% all the measures for fixed modeled values

% bsub -x -q week -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "stat_MP_cmaq_extendbin" -logfile "stat_MP_cmaq_extendbin.out"

if nargin < 1, yrz = 2001; end

% loop through each day of the year
yrNday = datevec(datenum(yrz,1,1):datenum(yrz,12,31));
yrNday = yrNday(:,1).*10^4 + yrNday(:,2).*10^2 + yrNday(:,3);
len = length(yrNday);

% load saved measures
load('matfiles/par_traditional_performance_gridbin.mat');
load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(1)));
modplots = modplots;

% variables 
len2 = size(num{1},1);
numGivMod = cell(len,1); mObsGivMod = cell(len,1); mModGivMod = cell(len,1);
mBiasGivMod = cell(len,1); nBiasGivMod = cell(len,1); nmBiasGivMod = cell(len,1);
fBiasGivMod = cell(len,1); mErrGivMod = cell(len,1); nErrGivMod = cell(len,1);
nmErrGivMod = cell(len,1); fErrGivMod = cell(len,1); RGivMod = cell(len,1);
R2GivMod = cell(len,1); sBiasGivMod = cell(len,1); msBiasGivMod = cell(len,1);
rmsBiasGivMod = cell(len,1); nrmsBiasGivMod = cell(len,1); mDsBiasGivMod = cell(len,1);
m2DmsBiasGivMod = cell(len,1); s2DmsBiasGivMod = cell(len,1); beta1GivMod = cell(len,1);
vObsGivMod = cell(len,1); vModGivMod = cell(len,1);

% go through each saved measure and find the measure for the modeled values
% in the grid and for fixed modeled values
matlabpool open 12
parfor i = 1:len
    disp(i);    
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),num{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    numGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),mObs{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    mObsGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),mMod{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    mModGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),mBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    mBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),nBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    nBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),nmBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    nmBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),fBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    fBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),mErr{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    mErrGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),nErr{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    nErrGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),nmErr{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    nmErrGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),fErr{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    fErrGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),R{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    RGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),R2{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    R2GivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),sBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    sBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),msBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    msBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),rmsBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    rmsBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),nrmsBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    nrmsBiasGivMod{i} = cell2mat(temp'); 
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),mDsBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    mDsBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),m2DmsBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    m2DmsBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),s2DmsBias{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    s2DmsBiasGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),beta1{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    beta1GivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),vObs{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    vObsGivMod{i} = cell2mat(temp');
end
parfor i = 1:len
    disp(i);
    temp2 = load(sprintf('../matfiles/PM2p5_%d_365_3_150_10_neg0.mat',yrNday(i)));    
    temp = arrayfun(@(x) interp1(temp2.mean_Mod(x,:),vMod{i}(x,:), ...
            [temp2.dailyCTMg(x) temp2.modplots],'linear','extrap'),1:len2,'UniformOutput',false); 
    vModGivMod{i} = cell2mat(temp');
end
matlabpool close

% save results
save('matfiles/traditional_performance_extendbin.mat', ...
    'numGivMod','mObsGivMod','mModGivMod','mBiasGivMod','nBiasGivMod','nmBiasGivMod', ...
    'fBiasGivMod','mErrGivMod','nErrGivMod','nmErrGivMod','fErrGivMod','RGivMod', ... 
    'R2GivMod','sBiasGivMod','msBiasGivMod','rmsBiasGivMod','nrmsBiasGivMod','mDsBiasGivMod', ...
    'm2DmsBiasGivMod','s2DmsBiasGivMod','beta1GivMod','vObsGivMod','vModGivMod', ...
    'CTMlocs','yrNday','modplots'); 

end