function [] = simulateObs(yr2simu)
% this function will take the results of generatelambda and the existing
% CMAQ modeled data and generate observed concentrations

if nargin < 1, yr2simu = 2002; end

% set seed
rng(0);

% gather all selected lambda values
load(sprintf('PM2p5_meanGivMod_yr%d.mat',yr2simu));
lambda1_select = meanGivMod_all;
lambda2_select = varGivMod_all;
modplots = modplots_all(1,:);

% getting paired modeled data and space/time locations of obs
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yr2simu));
% use coordObs, yrmodaObs, Mod
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
dayrmoda = datenum(yr,mo,da);
uniyrmoda = unique(dayrmoda);

% calculate observed data where z_hat ~ N(lambda1,lambda2)
mObs_calc = NaN*ones(length(Mod),1);
vObs_calc = NaN*ones(length(Mod),1);
Obssimu_calc = NaN*ones(length(Mod),1);
% 1) for every coordObs location see Mod, lambda1_select, lambda2_select
% 2) simulate lambda1, lambda2 through interpolation
for i = 1:length(uniyrmoda)
    disp(i);
    % subsetting all selected fields
    idx1 = uniyrmoda(i) == css(:,3); 
    cssSub = css(idx1,:); 
    lambda1_selectSub = lambda1_select(idx1,:);
    lambda2_selectSub = lambda2_select(idx1,:);

    % subsetting paired data
    idx2 = uniyrmoda(i) == dayrmoda;
    coordObsSub = coordObs(idx2,:);
    ModSub = Mod(idx2);
    
    % find closest points
    [idxa, dist] = knnsearch(cssSub(:,1:2),coordObsSub);

    % calculate lambda1, lambda2 through interpolation
    % unfortunately all the inputs for 'interp1' can't be a matrix, so I
    % had to loop through each row
    mObs_interp = NaN*ones(length(idxa),1); 
    vObs_interp = NaN*ones(length(idxa),1);
    x = modplots; 
    y = lambda1_selectSub(idxa,:); 
    z = lambda2_selectSub(idxa,:);
    for j = 1:length(ModSub)
        mObs_interp(j) = interp1(x,y(j,:),ModSub(j),'linear','extrap');
        vObs_interp(j) = interp1(x,z(j,:),ModSub(j),'linear','extrap');
    end
    %mObs_interp = interp1(meanMod_allSub(idx,:),lambda1_simuSub(idx,:),ModSub,'linear','extrap');
    %vObs_interp = interp1(meanMod_allSub(idx,:),lambda2_simuSub(idx,:),ModSub,'linear','extrap');
    
    % from lambda1 and lambda2 above, do the observed simulation
    len = length(mObs_interp);
    Obssimu_interp = mObs_interp + sqrt(vObs_interp).*randn(len,1);
    
    % putting final values in the correct order
    [aidx bidx] = ismember([coordObsSub repmat(uniyrmoda(i),length(coordObsSub),1)],[coordObs dayrmoda],'rows');
    mObs_calc(bidx) = mObs_interp;
    vObs_calc(bidx) = vObs_interp;
    Obssimu_calc(bidx) = Obssimu_interp;
        
end

% save results
save(sprintf('simulateObs_%d.mat',yr2simu),'coordObs','yrmodaObs','Mod', ...
    'css','lambda1_select','lambda2_select','mObs_calc','vObs_calc','Obssimu_calc','modplots'); 

end