function [] = antiPerform()
% this function will try to anticipate performance: soft/hard, DS, CMAQ, 
% model performance, location information at different radii

% bsub -x -q day -n 3 -R /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "antiPerform(0)" -logfile "stat_MP_cmaq_extendbin0.out"

%%% getting observed data and locations

load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone__20150408.mat', ...
        '_soft','_long','_gauss',0));
zh = cell2mat(zh_Xval);
ck_base = cell2mat(ckXval);

% subset to 2001
datez = datevec(ck_base(:,3));
idx = datez(:,1) == 2001;
zh = zh(idx); ck_base = ck_base(idx,:); 

%%% getting hard/soft mean/variance for increasing radius

for i = 0:100000:900000
    disp(i);
    % load hard
    if i == 0
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone__20150408.mat', ...
            '_nosoft','_long','_gauss',floor(i./1000)));
    else
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_20150408.mat', ...
        '_nosoft','_long','_gauss',floor(i./1000)));
    end
    
    % make sure locations match up/subset to 2001
    ck = cell2mat(ckXval);
    datez = datevec(ck(:,3));
    idx = datez(:,1) == 2001;
    if ~isequal(ck_base,ck(idx,:)), disp('mismatch locations');  end % this should be true
    zk = cell2mat(zk_madd);
    vk = cell2mat(vk);
        
    if i == 0
        zk_hard_0km = zk(idx);
        vk_hard_0km = vk(idx);
    elseif i == 100000
        zk_hard_100km = zk(idx);
        vk_hard_100km = vk(idx);
    elseif i == 200000
        zk_hard_200km = zk(idx);
        vk_hard_200km = vk(idx);
    elseif i == 300000
        zk_hard_300km = zk(idx);
        vk_hard_300km = vk(idx);
    elseif i == 400000
        zk_hard_400km = zk(idx);
        vk_hard_400km = vk(idx);
    elseif i == 500000
        zk_hard_500km = zk(idx);
        vk_hard_500km = vk(idx);
    elseif i == 600000
        zk_hard_600km = zk(idx);
        vk_hard_600km = vk(idx);
    elseif i == 700000
        zk_hard_700km = zk(idx);
        vk_hard_700km = vk(idx);
    elseif i == 800000
        zk_hard_800km = zk(idx);
        vk_hard_800km = vk(idx);
    elseif i == 900000
        zk_hard_900km = zk(idx);
        vk_hard_900km = vk(idx);
    end

    % load soft
    if i == 0
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone__20150408.mat', ...
            '_soft','_long','_gauss',floor(i./1000)));
    else
        load(sprintf('../matfiles/Xvalforcediso_LOOCV_%s%s%s_foriso%dkm_timezone_20150408.mat', ...
        '_soft','_long','_gauss',floor(i./1000)));
    end
    
    % make sure locations match up/subset to 2001
    ck = cell2mat(ckXval);
    datez = datevec(ck(:,3));
    idx = datez(:,1) == 2001;
    if ~isequal(ck_base,ck(idx,:)), disp('mismatch locations');  end % this should be true
    zk = cell2mat(zk_madd);
    vk = cell2mat(vk);
        
    if i == 0
        zk_soft_0km = zk(idx);
        vk_soft_0km = vk(idx);
    elseif i == 100000
        zk_soft_100km = zk(idx);
        vk_soft_100km = vk(idx);
    elseif i == 200000
        zk_soft_200km = zk(idx);
        vk_soft_200km = vk(idx);
    elseif i == 300000
        zk_soft_300km = zk(idx);
        vk_soft_300km = vk(idx);
    elseif i == 400000
        zk_soft_400km = zk(idx);
        vk_soft_400km = vk(idx);
    elseif i == 500000
        zk_soft_500km = zk(idx);
        vk_soft_500km = vk(idx);
    elseif i == 600000
        zk_soft_600km = zk(idx);
        vk_soft_600km = vk(idx);
    elseif i == 700000
        zk_soft_700km = zk(idx);
        vk_soft_700km = vk(idx);
    elseif i == 800000
        zk_soft_800km = zk(idx);
        vk_soft_800km = vk(idx);
    elseif i == 900000
        zk_soft_900km = zk(idx);
        vk_soft_900km = vk(idx);
    end
    
end

%%% getting mean/variance for CAMP for increasing radius

for i = 0:100000:900000
    disp(i);
           
    % load CAMP
    load(sprintf('../12_mfiles_othermethods/matdata/Xvalforcediso_LOOCV_CAMP_foriso%dkm_timezone_20150408.mat',floor(i./1000)));
        
    % make sure locations match up/subset to 2001
    ck = cell2mat(ckXval);
    datez = datevec(ck(:,3));
    idx = datez(:,1) == 2001;
    if ~isequal(ck_base,ck(idx,:)), disp('mismatch locations');  end % this should be true
    zk = cell2mat(zk_madd);
    vk = cell2mat(vk);
        
    if i == 0
        zk_camp_0km = zk(idx);
        vk_camp_0km = vk(idx);
    elseif i == 100000
        zk_camp_100km = zk(idx);
        vk_camp_100km = vk(idx);
    elseif i == 200000
        zk_camp_200km = zk(idx);
        vk_camp_200km = vk(idx);
    elseif i == 300000
        zk_camp_300km = zk(idx);
        vk_camp_300km = vk(idx);
    elseif i == 400000
        zk_camp_400km = zk(idx);
        vk_camp_400km = vk(idx);
    elseif i == 500000
        zk_camp_500km = zk(idx);
        vk_camp_500km = vk(idx);
    elseif i == 600000
        zk_camp_600km = zk(idx);
        vk_camp_600km = vk(idx);
    elseif i == 700000
        zk_camp_700km = zk(idx);
        vk_camp_700km = vk(idx);
    elseif i == 800000
        zk_camp_800km = zk(idx);
        vk_camp_800km = vk(idx);
    elseif i == 900000
        zk_camp_900km = zk(idx);
        vk_camp_900km = vk(idx);
    end
    
end

%%% getting mean/variance for three downscaler methods for increasing radius

for i = 0:100000:900000
    disp(i);
    %%% load static downscaler
    load(sprintf('../12_mfiles_othermethods/matdata/LOOCE_DSforcedisolation_%dkm.mat',floor(i./1000)));
    
    % cell2mat
    temp1 = cell2mat(cellfun(@length,zk,'UniformOutput',false));
    temp2 = cell2mat(arrayfun(@(x,y) repmat(x,y,1),unidates,temp1,'UniformOutput',false));
    zk = cell2mat(zk);
    vk = cell2mat(vk);
    ck = [cell2mat(ckall) temp2];    
    
    % finding closest observed data location
    zk_SDS = NaN*ones(length(ck_base),1); 
    vk_SDS = NaN*ones(length(ck_base),1);
    
    unitime = unique(ck_base(:,3));
    disp('static');
    for j = 1:length(unitime)
        idx = unitime(j) == ck(:,3);
        idx_base = unitime(j) == ck_base(:,3);
        zkidx = zk(idx); vkidx = vk(idx);
        [idxmin d] = knnsearch(ck(idx,1:2),ck_base(idx_base,1:2));
        zksub = zkidx(idxmin); zksub(d>36000) = NaN; 
        vksub = vkidx(idxmin); vksub(d>36000) = NaN;
        
        [aidx bidx] = ismember(ck_base,ck_base(idx_base,:),'rows');
        bidx(bidx==0) = [];
        if ~isempty(zksub), zk_SDS(aidx) = zksub(bidx); end
        if ~isempty(vksub), vk_SDS(aidx) = vksub(bidx); end
    end
    
    %%% loading space/time downscaler number 1
    load(sprintf('../12_mfiles_othermethods/matdata/LOOCE_DSforcedisolation_%dkm_add_ind_muli_ind.mat',floor(i./1000)));
    zk = cell2mat(zk);
    vk = cell2mat(vk);
    ck = cell2mat(ck);

    % finding closest observed data location
    zk_STDSI = NaN*ones(length(ck_base),1); 
    vk_STDSI = NaN*ones(length(ck_base),1);
    
    unitime = unique(ck_base(:,3));
    disp('s/t ds 1');
    for j = 1:length(unitime)
        idx = unitime(j) == ck(:,3);
        idx_base = unitime(j) == ck_base(:,3);
        zkidx = zk(idx); vkidx = vk(idx);
        [idxmin d] = knnsearch(ck(idx,1:2),ck_base(idx_base,1:2));
        zksub = zkidx(idxmin); zksub(d>36000) = NaN; 
        vksub = vkidx(idxmin); vksub(d>36000) = NaN;
        
        [aidx bidx] = ismember(ck_base,ck_base(idx_base,:),'rows');
        bidx(bidx==0) = [];
        if ~isempty(zksub), zk_STDSI(aidx) = zksub(bidx); end
        if ~isempty(vksub), vk_STDSI(aidx) = vksub(bidx); end
    end
    
    %%% loading space/time downscaler number 2
    load(sprintf('../12_mfiles_othermethods/matdata/LOOCE_DSforcedisolation_%dkm_add_dyn_muli_ind.mat',floor(i./1000)));
    zk = cell2mat(zk);
    vk = cell2mat(vk);
    ck = cell2mat(ck);

    % finding closest observed data location
    zk_STDSII = NaN*ones(length(ck_base),1); 
    vk_STDSII = NaN*ones(length(ck_base),1);
    
    unitime = unique(ck_base(:,3));
    disp('s/t ds 2');
    for j = 1:length(unitime)
        idx = unitime(j) == ck(:,3);
        idx_base = unitime(j) == ck_base(:,3);
        zkidx = zk(idx); vkidx = vk(idx);
        [idxmin d] = knnsearch(ck(idx,1:2),ck_base(idx_base,1:2));
        zksub = zkidx(idxmin); zksub(d>36000) = NaN; 
        vksub = vkidx(idxmin); vksub(d>36000) = NaN;
        
        [aidx bidx] = ismember(ck_base,ck_base(idx_base,:),'rows');
        bidx(bidx==0) = [];
        if ~isempty(zksub), zk_STDSII(aidx) = zksub(bidx); end
        if ~isempty(vksub), vk_STDSII(aidx) = vksub(bidx); end
    end
        
    if i == 0
        zk_SDS_0km = zk_SDS;
        vk_SDS_0km = vk_SDS;
        zk_STDSI_0km = zk_STDSI;
        vk_STDSI_0km = vk_STDSI;
        zk_STDSII_0km = zk_STDSII;
        vk_STDSII_0km = vk_STDSII;
    elseif i == 100000
        zk_SDS_100km = zk_SDS;
        vk_SDS_100km = vk_SDS;
        zk_STDSI_100km = zk_STDSI;
        vk_STDSI_100km = vk_STDSI;
        zk_STDSII_100km = zk_STDSII;
        vk_STDSII_100km = vk_STDSII;
    elseif i == 200000
        zk_SDS_200km = zk_SDS;
        vk_SDS_200km = vk_SDS;
        zk_STDSI_200km = zk_STDSI;
        vk_STDSI_200km = vk_STDSI;
        zk_STDSII_200km = zk_STDSII;
        vk_STDSII_200km = vk_STDSII;
    elseif i == 300000
        zk_SDS_300km = zk_SDS;
        vk_SDS_300km = vk_SDS;
        zk_STDSI_300km = zk_STDSI;
        vk_STDSI_300km = vk_STDSI;
        zk_STDSII_300km = zk_STDSII;
        vk_STDSII_300km = vk_STDSII;
    elseif i == 400000
        zk_SDS_400km = zk_SDS;
        vk_SDS_400km = vk_SDS;
        zk_STDSI_400km = zk_STDSI;
        vk_STDSI_400km = vk_STDSI;
        zk_STDSII_400km = zk_STDSII;
        vk_STDSII_400km = vk_STDSII;
    elseif i == 500000
        zk_SDS_500km = zk_SDS;
        vk_SDS_500km = vk_SDS;
        zk_STDSI_500km = zk_STDSI;
        vk_STDSI_500km = vk_STDSI;
        zk_STDSII_500km = zk_STDSII;
        vk_STDSII_500km = vk_STDSII;
    elseif i == 600000
        zk_SDS_600km = zk_SDS;
        vk_SDS_600km = vk_SDS;
        zk_STDSI_600km = zk_STDSI;
        vk_STDSI_600km = vk_STDSI;
        zk_STDSII_600km = zk_STDSII;
        vk_STDSII_600km = vk_STDSII;
    elseif i == 700000
        zk_SDS_700km = zk_SDS;
        vk_SDS_700km = vk_SDS;
        zk_STDSI_700km = zk_STDSI;
        vk_STDSI_700km = vk_STDSI;
        zk_STDSII_700km = zk_STDSII;
        vk_STDSII_700km = vk_STDSII;
    elseif i == 800000
        zk_SDS_800km = zk_SDS;
        vk_SDS_800km = vk_SDS;
        zk_STDSI_800km = zk_STDSI;
        vk_STDSI_800km = vk_STDSI;
        zk_STDSII_800km = zk_STDSII;
        vk_STDSII_800km = vk_STDSII;
    elseif i == 900000
        zk_SDS_900km = zk_SDS;
        vk_SDS_900km = vk_SDS;
        zk_STDSI_900km = zk_STDSI;
        vk_STDSI_900km = vk_STDSI;
        zk_STDSII_900km = zk_STDSII;
        vk_STDSII_900km = vk_STDSII;
    end
    
end

%%% load CMAQ data

load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
ckt = datenum(yr,mo,da);
ck = [coordObs ckt];

% finding closest observed data location
CMAQ = NaN*ones(length(ck_base),1);
zk = Mod;

unitime = unique(ck_base(:,3));
for i = 1:length(unitime)
    idx = unitime(i) == ck(:,3);
    idx_base = unitime(i) == ck_base(:,3);
    zkidx = zk(idx); 
    [idxmin d] = knnsearch(ck(idx,1:2),ck_base(idx_base,1:2));
    zksub = zkidx(idxmin); zksub(d>36000) = NaN; 

    [aidx bidx] = ismember(ck_base,ck_base(idx_base,:),'rows');
    bidx(bidx==0) = [];
    if ~isempty(zksub), CMAQ(aidx) = zksub(bidx); end
end

%%% load all CMAQ model performance metrics

load('matfiles/traditional_performance_grid.mat');
yr = floor(yrNday./10000);
mo = floor((yrNday - yr.*10000)./100);
da = yrNday - yr.*10000 - mo.*100;
ctCTM = datenum(yr,mo,da)'; ckt = repmat(ctCTM,length(CTMlocs),1);
ck = [repmat(CTMlocs,length(yrNday),1) ckt(:)];

% finding closest observed data location
zk = 1:length(ck);
zkfinal = NaN*ones(length(ck_base),1);

unitime = unique(ck_base(:,3));
for i = 1:length(unitime)
    idx = unitime(i) == ck(:,3);
    idx_base = unitime(i) == ck_base(:,3);
    zkidx = zk(idx); 
    [idxmin d] = knnsearch(ck(idx,1:2),ck_base(idx_base,1:2));
    zksub = zkidx(idxmin); zksub(d>36000) = NaN; 

    [aidx bidx] = ismember(ck_base,ck_base(idx_base,:),'rows');
    bidx(bidx==0) = [];
    if ~isempty(zksub), zkfinal(aidx) = zksub(bidx); end
end
idxfinal = ~isnan(zkfinal);
zkfinal(isnan(zkfinal)) = [];

fErr_2=NaN*ones(length(ck_base),1); prefErr=fErr(:); fErr_2(idxfinal)=prefErr(zkfinal);  
mMod_2=NaN*ones(length(ck_base),1); premMod=mMod(:); mMod_2(idxfinal)=premMod(zkfinal);
nmBias_2=NaN*ones(length(ck_base),1); prenmBias=nmBias(:); nmBias_2(idxfinal)=prenmBias(zkfinal);       
s2DmsBias_2=NaN*ones(length(ck_base),1); pres2DmsBias=s2DmsBias(:); s2DmsBias_2(idxfinal)=pres2DmsBias(zkfinal);    
R_2=NaN*ones(length(ck_base),1); preR=R(:); R_2(idxfinal)=preR(zkfinal);
m2DmsBias_2=NaN*ones(length(ck_base),1); prem2DmsBias=m2DmsBias(:); m2DmsBias_2(idxfinal)=prem2DmsBias(zkfinal);        
mObs_2=NaN*ones(length(ck_base),1); premObs=mObs(:); mObs_2(idxfinal)=premObs(zkfinal);
nmErr_2=NaN*ones(length(ck_base),1); prenmErr=nmErr(:); nmErr_2(idxfinal)=prenmErr(zkfinal);      
sBias_2=NaN*ones(length(ck_base),1); presBias=sBias(:); sBias_2(idxfinal)=presBias(zkfinal);     
R2_2=NaN*ones(length(ck_base),1); preR2=R2(:); R2_2(idxfinal)=preR2(zkfinal);      
mBias_2=NaN*ones(length(ck_base),1); premBias=mBias(:); mBias_2(idxfinal)=premBias(zkfinal);        
msBias_2=NaN*ones(length(ck_base),1); premsBias=msBias(:); msBias_2(idxfinal)=premsBias(zkfinal);     
nrmsBias_2=NaN*ones(length(ck_base),1); prenrmsBias=nrmsBias(:); nrmsBias_2(idxfinal)=prenrmsBias(zkfinal);     
vMod_2=NaN*ones(length(ck_base),1); prevMod=vMod(:); vMod_2(idxfinal)=prevMod(zkfinal);   
beta1_2=NaN*ones(length(ck_base),1); prebeta1=beta1(:); beta1_2(idxfinal)=prebeta1(zkfinal);      
mDsBias_2=NaN*ones(length(ck_base),1); premDsBias=mDsBias(:); mDsBias_2(idxfinal)=premDsBias(zkfinal);      
nBias_2=NaN*ones(length(ck_base),1); prenBias=nBias(:); nBias_2(idxfinal)=prenBias(zkfinal);    
vObs_2=NaN*ones(length(ck_base),1); prevObs=vObs(:); vObs_2(idxfinal)=prevObs(zkfinal);             
fBias_2=NaN*ones(length(ck_base),1); prefBias=fBias(:); fBias_2(idxfinal)=prefBias(zkfinal);       
mErr_2=NaN*ones(length(ck_base),1); premErr=mErr(:); mErr_2(idxfinal)=premErr(zkfinal);      
nErr_2=NaN*ones(length(ck_base),1); prenErr=nErr(:); nErr_2(idxfinal)=prenErr(zkfinal);       
rmsBias_2=NaN*ones(length(ck_base),1); prermsBias=rmsBias(:); rmsBias_2(idxfinal)=prermsBias(zkfinal); 

%%% get lambda1 and lambda2

load(sprintf('../matfiles/PM2p5_soft_yr%d.mat',2001));
ck = css;

% finding closest observed data location
lambda1_2 = NaN*ones(length(ck_base),1);
lambda2_2 = NaN*ones(length(ck_base),1);
zk = lambda1;
zk2 = lambda2;

unitime = unique(ck_base(:,3));
for i = 1:length(unitime)
    idx = unitime(i) == ck(:,3);
    idx_base = unitime(i) == ck_base(:,3);
    zkidx = zk(idx); zk2idx = zk2(idx);
    [idxmin d] = knnsearch(ck(idx,1:2),ck_base(idx_base,1:2));
    zksub = zkidx(idxmin); zksub(d>36000) = NaN; 
    zk2sub = zk2idx(idxmin); zk2sub(d>36000) = NaN;

    [aidx bidx] = ismember(ck_base,ck_base(idx_base,:),'rows');
    bidx(bidx==0) = [];
    if ~isempty(zksub), lambda1_2(aidx) = zksub(bidx); end
    if ~isempty(zk2sub), lambda2_2(aidx) = zk2sub(bidx); end
end

%%% load all location information (e.g. East/West, network, distance to 
% closest monitor, TEOM/FRM, etc.)

% East/West (using the 100 degree meridian)
cd ../09_mfiles_projections
M100 = ell2lambertcc([-100 40],'whiproj2001');
cd ../13_mfiles_modelperformance
IsitWest = NaN*ones(length(ck_base),1);
IsitWest(ck_base(:,1)<M100(1)) = 1;
IsitWest(ck_base(:,1)>M100(1)) = 0;

% distance to closest monitor
unilocs = unique(ck_base(:,1:2),'rows');
X = unilocs';
Y = unilocs';
pairdist = sqrt( bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y) );
[distz dummy] = sort(pairdist,2);
distz = distz(:,2);
Dist2NMon = NaN*ones(length(ck_base),1);
for i = 1:length(unilocs)
    idx = unilocs(i,1) == ck_base(:,1) & unilocs(i,2) == ck_base(:,2);
    Dist2NMon(idx) = distz(i);
end

% getting original information
load(sprintf('../datafiles/Observed_PM2p5/MasterDaily_PM2p5_%d.mat',2001));
yr = floor(dates./10000);
mo = floor((dates - yr.*10000)./100);
da = dates - yr.*10000 - mo.*100;
ogObs = datenum(yr,mo,da);

% convert to projection
cd ../09_mfiles_projections
ogProj = ell2lambertcc([longitude latitude],'whiproj2001');
cd ../13_mfiles_modelperformance
[lia lib] = ismember([ogProj ogObs],ck_base,'rows');
lib(lib==0) = [];
ranks_2 = ranks(lib);
locations_2 = locations(lib);

% network
IsIMPROVE = NaN*ones(length(ck_base),1); IsIMPROVE(ranks_2==5) = 1;
IsSTN = NaN*ones(length(ck_base),1); IsSTN(ranks_2==4) = 1;

% TEOM/FRM
IsFRM = NaN*ones(length(ck_base),1); IsFRM(ranks>=1&ranks<=3) = 1;
IsTEOM = NaN*ones(length(ck_base),1); IsTEOM(ranks_2>=4&ranks_2<=6) = 1;

% urban/rural/suburban
IsUrban = NaN*ones(length(ck_base),1); IsUrban(locations_2==5) = 1;
IsRural = NaN*ones(length(ck_base),1); IsRural(locations_2==2) = 1;
IsSuburban = NaN*ones(length(ck_base),1); IsSuburban(locations_2==3) = 1;

% get region
if ~exist('matfiles/inregion.mat')
    allregions = shaperead('FWS_LCC/FWS_LCC.shp');
    for k = 1:length(allregions)
        tic
        cd ../09_mfiles_projections
        allregions_p{k,1} = ell2lambertcc([allregions(k).X',allregions(k).Y'],'whiproj2001');
        cd ../13_mfiles_modelperformance
        in{k,1} = inpolygon(ck_base(:,1),ck_base(:,2),allregions_p{k,1}(:,1),allregions_p{k,1}(:,2));
        toc
    end
    save('matfiles/inregion.mat','in','allregions');
else
    load('matfiles/inregion.mat')
end

% get season
[yr mo da] = datevec(ck_base(:,3));
IsWinter = NaN*ones(length(ck_base),1); IsWinter(mo==1|mo==2|mo==12) = 1;
IsSpring = NaN*ones(length(ck_base),1); IsSpring(mo==3|mo==4|mo==5) = 1;
IsSummer = NaN*ones(length(ck_base),1); IsSummer(mo==6|mo==7|mo==8) = 1;
IsFall = NaN*ones(length(ck_base),1); IsFall(mo==9|mo==10|mo==11) = 1;

% get longitude
cd ../09_mfiles_projections
testing = [-125:5:-90]'; 
M = ell2lambertcc([testing 40*ones(length(testing),1)],'whiproj2001');
cd ../13_mfiles_modelperformance
Is125 = NaN*ones(length(ck_base),1); Is125(ck_base(:,1)<M(1,1)) = 1;
Is120 = NaN*ones(length(ck_base),1); Is120(ck_base(:,1)<M(2,1)) = 1;
Is115 = NaN*ones(length(ck_base),1); Is115(ck_base(:,1)<M(3,1)) = 1;
Is110 = NaN*ones(length(ck_base),1); Is110(ck_base(:,1)<M(4,1)) = 1;
Is105 = NaN*ones(length(ck_base),1); Is105(ck_base(:,1)<M(5,1)) = 1;
Is100 = NaN*ones(length(ck_base),1); Is100(ck_base(:,1)<M(6,1)) = 1;
Is95 = NaN*ones(length(ck_base),1); Is95(ck_base(:,1)<M(7,1)) = 1;
Is90 = NaN*ones(length(ck_base),1); Is90(ck_base(:,1)<M(8,1)) = 1;

%%% save results
save('matfiles/allInfo.mat','zh','ck_base','zk_hard_0km','vk_hard_0km', ...
    'zk_hard_100km','vk_hard_100km','zk_hard_200km','vk_hard_200km', ...
    'zk_hard_300km','vk_hard_300km','zk_hard_400km','vk_hard_400km', ...
    'zk_hard_500km','vk_hard_500km','zk_hard_600km','vk_hard_600km', ...
    'zk_hard_700km','vk_hard_700km','zk_hard_800km','vk_hard_800km', ...
    'zk_hard_900km','vk_hard_900km','zk_soft_0km','vk_soft_0km', ...
    'zk_soft_100km','vk_soft_100km','zk_soft_200km','vk_soft_200km', ...
    'zk_soft_300km','vk_soft_300km','zk_soft_400km','vk_soft_400km', ...
    'zk_soft_500km','vk_soft_500km','zk_soft_600km','vk_soft_600km', ...
    'zk_soft_700km','vk_soft_700km','zk_soft_800km','vk_soft_800km', ...
    'zk_soft_900km','vk_soft_900km', ...
    'zk_camp_0km','vk_camp_0km','zk_camp_100km','vk_camp_100km', ...
    'zk_camp_200km','vk_camp_200km','zk_camp_300km','vk_camp_300km', ...
    'zk_camp_400km','vk_camp_400km','zk_camp_500km','vk_camp_500km', ...
    'zk_camp_600km','vk_camp_600km','zk_camp_700km','vk_camp_700km', ...
    'zk_camp_800km','vk_camp_800km','zk_camp_900km','vk_camp_900km', ...
    'zk_SDS_0km','vk_SDS_0km','zk_STDSI_0km','vk_STDSI_0km','zk_STDSII_0km','vk_STDSII_0km', ...
    'zk_SDS_100km','vk_SDS_100km','zk_STDSI_100km','vk_STDSI_100km','zk_STDSII_100km','vk_STDSII_100km', ...
    'zk_SDS_200km','vk_SDS_200km','zk_STDSI_200km','vk_STDSI_200km','zk_STDSII_200km','vk_STDSII_200km', ...
    'zk_SDS_300km','vk_SDS_300km','zk_STDSI_300km','vk_STDSI_300km','zk_STDSII_300km','vk_STDSII_300km', ...
    'zk_SDS_400km','vk_SDS_400km','zk_STDSI_400km','vk_STDSI_400km','zk_STDSII_400km','vk_STDSII_400km', ...
    'zk_SDS_500km','vk_SDS_500km','zk_STDSI_500km','vk_STDSI_500km','zk_STDSII_500km','vk_STDSII_500km', ...
    'zk_SDS_600km','vk_SDS_600km','zk_STDSI_600km','vk_STDSI_600km','zk_STDSII_600km','vk_STDSII_600km', ...
    'zk_SDS_700km','vk_SDS_700km','zk_STDSI_700km','vk_STDSI_700km','zk_STDSII_700km','vk_STDSII_700km', ...
    'zk_SDS_800km','vk_SDS_800km','zk_STDSI_800km','vk_STDSI_800km','zk_STDSII_800km','vk_STDSII_800km', ...
    'zk_SDS_900km','vk_SDS_900km','zk_STDSI_900km','vk_STDSI_900km','zk_STDSII_900km','vk_STDSII_900km', ...
    'CMAQ', ...
    'fErr_2','mMod_2','nmBias_2','s2DmsBias_2','R_2','m2DmsBias_2','mObs_2', ...
    'nmErr_2','sBias_2','R2_2','mBias_2','msBias_2','nrmsBias_2','vMod_2', ...       
    'beta1_2','mDsBias_2','nBias_2','vObs_2','fBias_2','mErr_2','nErr_2','rmsBias_2', ...
    'lambda1_2','lambda2_2', ...
    'IsitWest','Dist2NMon','IsIMPROVE','IsSTN','IsFRM','IsTEOM','IsUrban','IsRural','IsSuburban', ...
    'in','IsWinter','IsSpring','IsSummer','IsFall','Is125','Is120','Is115','Is110','Is105','Is95','Is90');

end