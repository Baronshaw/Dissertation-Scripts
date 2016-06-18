function [] = performanceInfo()
% this function will get all the information needed to run all the scripts
% in this folder

% load longitude, latitude, day, observed value, CTM value for 2001
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
xproj = coordObs(:,1);
yproj = coordObs(:,2);
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;

% get lambda1 and lambda2 value
yrmoda = datenum(yr,mo,da);
uniyrmoda = unique(yrmoda); 
load('../matfiles/PM2p5_soft_yr2001.mat');
clear limi probdens
L1 = NaN*ones(length(Mod),1);
L2 = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == css(:,3);
    cssSub = css(idx,1:2); lambda1sub = lambda1(idx); lambda2sub = lambda2(idx);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(cssSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1sub(idx);
    L2ordered = lambda2sub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    L1(bidx) = L1ordered;
    L2(bidx) = L2ordered;
    
end

% rename variables
lambda1 = L1; lambda2 = L2;

% East/West (using the 100 degree meridian)
cd ../09_mfiles_projections
M100 = ell2lambertcc([-100 40],'whiproj2001');
cd ../17_mfiles_validatelambda1
IsWest = NaN*ones(length(Obs),1);
IsWest(coordObs(:,1)<M100(1)) = 1;
IsWest(coordObs(:,1)>M100(1)) = 0;

% distance to closest monitor
unilocs = unique(coordObs,'rows');
X = unilocs';
Y = unilocs';
pairdist = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
[distz dummy] = sort(pairdist,2);
distz = distz(:,2);
Dist2NMon = NaN*ones(length(Obs),1);
for i = 1:length(unilocs)
    disp(i);
    idx = unilocs(i,1) == coordObs(:,1) & unilocs(i,2) == coordObs(:,2);
    Dist2NMon(idx) = distz(i);
end

% getting original information
load(sprintf('../datafiles/Observed_PM2p5/MasterDaily_PM2p5_%d.mat',2001));

% convert to projection
cd ../09_mfiles_projections
ogProj = ell2lambertcc([longitude latitude],'whiproj2001');
cd ../17_mfiles_validatelambda1
unidatenum = unique(yrmodaObs);
ranks_2 = NaN*ones(length(coordObs),1);
locations_2 = NaN*ones(length(coordObs),1);
for i = 1:length(unidatenum)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = unidatenum(i) == dates;
    ogProjSub = ogProj(idx,1:2); rankssub = ranks(idx); locationssub = locations(idx);
    idx = unidatenum(i) == yrmodaObs;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(ogProjSub,cMSsub);
    if max(dist) > 10, disp(max(dist)); end
    
    % putting all the values in place for the final variables
    ranksordered = rankssub(idx);
    locationsordered = locationssub(idx);
    [aidx bidx] = ismember([cMSsub repmat(unidatenum(i),length(cMSsub),1)],[coordObs yrmodaObs],'rows');
    ranks_2(bidx) = ranksordered;
    locations_2(bidx) = locationsordered;
    
end

% network
IsIMPROVE = NaN*ones(length(Obs),1); IsIMPROVE(ranks_2==5) = 1;
IsSTN = NaN*ones(length(Obs),1); IsSTN(ranks_2==4) = 1;

% TEOM/FRM
IsFRM = NaN*ones(length(Obs),1); IsFRM(ranks>=1&ranks<=3) = 1;
IsTEOM = NaN*ones(length(Obs),1); IsTEOM(ranks_2>=4&ranks_2<=6) = 1;

% urban/rural/suburban
IsUrban = NaN*ones(length(Obs),1); IsUrban(locations_2==5) = 1;
IsRural = NaN*ones(length(Obs),1); IsRural(locations_2==2) = 1;
IsSuburban = NaN*ones(length(Obs),1); IsSuburban(locations_2==3) = 1;

% get region
if ~exist('regions.mat')
    allregions = shaperead('../13_mfiles_modelperformance/FWS_LCC/FWS_LCC.shp');
    for k = 1:length(allregions)
        tic
        cd ../09_mfiles_projections
        allregions_p{k,1} = ell2lambertcc([allregions(k).X',allregions(k).Y'],'whiproj2001');
        cd ../17_mfiles_validatelambda1
        inregions{k,1} = inpolygon(coordObs(:,1),coordObs(:,2),allregions_p{k,1}(:,1),allregions_p{k,1}(:,2));
        toc
    end
    save('regions.mat','allregions_p','inregions');
else
    load('regions.mat')
end

% get region string
allregions = shaperead('../13_mfiles_modelperformance/FWS_LCC/FWS_LCC.shp');
for k = 1:length(allregions)
    region_str{k,1} = allregions(k).area_names;
end

% get season
IsWinter = NaN*ones(length(Obs),1); IsWinter(mo==1|mo==2|mo==12) = 1;
IsSpring = NaN*ones(length(Obs),1); IsSpring(mo==3|mo==4|mo==5) = 1;
IsSummer = NaN*ones(length(Obs),1); IsSummer(mo==6|mo==7|mo==8) = 1;
IsFall = NaN*ones(length(Obs),1); IsFall(mo==9|mo==10|mo==11) = 1;

% get lambda1 and lambda2 from the CAMP method
yrmoda = datenum(yr,mo,da);
uniyrmoda = unique(yrmoda); 
load('CAMPmethod.mat');
lambda1CAMP = NaN*ones(length(Mod),1);
lambda2CAMP = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3));
    distCTMvSub = distCTMv(idx,:); lambda1CAMPsub = meanGivMod(idx,1); lambda2CAMPsub = varGivMod(idx,1);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(distCTMvSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1CAMPsub(idx);
    L2ordered = lambda2CAMPsub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    lambda1CAMP(bidx) = L1ordered;
    lambda2CAMP(bidx) = L2ordered;
    
end

% get lambda1 and lambda2 from the regionalized CAMP method
load('CAMPmethod_regional.mat');
idx = ~isnan(inregion6v); yrmodaCTMv = yrmodaCTMv(idx,:); distCTMv = distCTMv(idx,:);
meanGivMod = meanGivMod(idx,:); varGivMod = varGivMod(idx,:);
lambda1CAMP6 = NaN*ones(length(Mod),1);
lambda2CAMP6 = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3));
    distCTMvSub = distCTMv(idx,:); lambda1CAMP6sub = meanGivMod(idx,1); lambda2CAMP6sub = varGivMod(idx,1);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(distCTMvSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1CAMP6sub(idx);
    L2ordered = lambda2CAMP6sub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    lambda1CAMP6(bidx) = L1ordered;
    lambda2CAMP6(bidx) = L2ordered;
    
end

% get lambda1 and lambda2 from the seasonal CAMP method
load('CAMPmethod_seasonal.mat');
lambda1CAMPS = NaN*ones(length(Mod),1);
lambda2CAMPS = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3));
    distCTMvSub = distCTMv(idx,:); lambda1CAMPSsub = meanGivMod(idx,1); lambda2CAMPSsub = meanGivMod(idx,1);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(distCTMvSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1CAMPSsub(idx);
    L2ordered = lambda2CAMPSsub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    lambda1CAMPS(bidx) = L1ordered;
    lambda2CAMPS(bidx) = L2ordered;
    
end

% get lambda1 and lambda2 from the Constant method
yrmoda = datenum(yr,mo,da);
uniyrmoda = unique(yrmoda); 
load('Constantmethod.mat');
lambda1Constant = NaN*ones(length(Mod),1);
lambda2Constant = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3));
    distCTMvSub = distCTMv(idx,:); lambda1Constantsub = meanGivMod(idx,1); lambda2Constantsub = varGivMod(idx,1);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(distCTMvSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1Constantsub(idx);
    L2ordered = lambda2Constantsub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    lambda1Constant(bidx) = L1ordered;
    lambda2Constant(bidx) = L2ordered;
    
end

% get lambda1 and lambda2 from the regionalized Constant method
load('ConstantRegionmethod.mat');
idx = ~isnan(inregion6v); yrmodaCTMv = yrmodaCTMv(idx,:); distCTMv = distCTMv(idx,:);
meanGivMod = meanGivMod(idx,:); varGivMod = varGivMod(idx,:);
lambda1Constant6 = NaN*ones(length(Mod),1);
lambda2Constant6 = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3));
    distCTMvSub = distCTMv(idx,:); lambda1Constant6sub = meanGivMod(idx,1); lambda2Constant6sub = varGivMod(idx,1);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(distCTMvSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1Constant6sub(idx);
    L2ordered = lambda2Constant6sub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    lambda1Constant6(bidx) = L1ordered;
    lambda2Constant6(bidx) = L2ordered;
    
end

% get lambda1 and lambda2 from the seasonal CAMP method
load('ConstantSeasonmethod.mat');
lambda1ConstantS = NaN*ones(length(Mod),1);
lambda2ConstantS = NaN*ones(length(Mod),1);
for i = 1:length(uniyrmoda)
    disp(i);
    % finding the nearest lambda1 locations for each observed value
    idx = uniyrmoda(i) == datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3));
    distCTMvSub = distCTMv(idx,:); lambda1ConstantSsub = meanGivMod(idx,1); lambda2ConstantSsub = meanGivMod(idx,1);
    idx = uniyrmoda(i) == yrmoda;
    cMSsub = [xproj(idx) yproj(idx)];
    [idx, dist] = knnsearch(distCTMvSub,cMSsub);
    
    % putting all the values in place for the final variables
    L1ordered = lambda1ConstantSsub(idx);
    L2ordered = lambda2ConstantSsub(idx);
    [aidx bidx] = ismember([cMSsub repmat(uniyrmoda(i),length(cMSsub),1)],[coordObs yrmoda],'rows');
    lambda1ConstantS(bidx) = L1ordered;
    lambda2ConstantS(bidx) = L2ordered;
    
end

% save all variables
save('performanceInfo.mat','Obs','Mod','coordObs','xproj','yproj', ...
    'yrmodaObs','yr','mo','da','yrmoda','lambda1','lambda2', ...
    'IsWest','Dist2NMon','IsIMPROVE','IsSTN','IsFRM','IsTEOM','IsUrban', ...
    'IsRural','IsSuburban','inregions','allregions_p','region_str', ...
    'IsWinter','IsSpring','IsSummer','IsFall', ...
    'lambda1CAMP','lambda2CAMP','lambda1CAMP6','lambda2CAMP6','inregion6', ...
    'lambda1CAMPS','lambda2CAMPS','lambda1Constant','lambda2Constant', ...
    'lambda1Constant6','lambda2Constant6','lambda1ConstantS','lambda2ConstantS');

end