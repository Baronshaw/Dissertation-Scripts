function [] = simulateObs(yr2simu)
% this function will take the results of generatelambda and the existing
% CMAQ modeled data and generate observed concentrations

if nargin < 1, yr2simu = 2002; end

% gather all generated lambda1 values
lambda1_gen = cell(1,10);
for i = 1:10
    load(sprintf('RFlambda1_%d_decile_%0.2d.mat',yr2simu,i));
    % use cGen and FcGen
    lambda1_gen{1,i} = FcGen + 16; % I added 16 so no values are negative
end
lambda1_gen = cell2mat(lambda1_gen);

% get mean modeled space/time locations
load(sprintf('PM2p5_meanMod_yr%d.mat',yr2simu));
% use css, mean_Mod

% getting paired modeled data and space/time locations of obs
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yr2simu));
% use coordObs, yrmodaObs, Mod
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
dnyrmoda = datenum(yr,mo,da);
uniyrmoda = unique(dnyrmoda);

% calculate observed data where z_hat ~ N(lambda1,lambda2)
mObs_calc = NaN*ones(length(Mod),1);
% 1) for every coordObs location see Mod, FcGen
% 2) simulate lambda1 through interpolation
for i = 1:length(uniyrmoda)
    disp(i);
    % subsetting all generated random fields
    idx1 = uniyrmoda(i) == cGen(:,3); 
    cGenSub = cGen(idx1,:); 
    lambda1_genSub = lambda1_gen(idx1,:);
    
    % subsetting all the gathered modeled data
    idx2 = uniyrmoda(i) == css(:,3);
    cssSub = css(idx2,:);
    meanMod_allSub = meanMod_all(idx2,:);

    % subsetting paired data
    idx2 = uniyrmoda(i) == dnyrmoda;
    dnyrmodaSub = dnyrmoda(idx2);
    coordObsSub = coordObs(idx2,:);
    ModSub = Mod(idx2);
    
    % find closest points
    [idxa, dist] = knnsearch(cssSub(:,1:2),coordObsSub);
    [idxb, dist] = knnsearch(cGenSub(:,1:2),coordObsSub);
    
    % calculate lambda1, lambda2 through interpolation
    % unfortunately all the inputs for 'interp1' can't be a matrix, so I
    % had to loop through each row
    mObs_interp = NaN*ones(length(idxa),1); 
    x = meanMod_allSub(idxa,:); 
    y = lambda1_genSub(idxb,:); 
    for j = 1:length(ModSub)
        mObs_interp(j) = interp1(x(j,:),y(j,:),ModSub(j),'linear','extrap');
    end
    %mObs_interp = interp1(meanMod_allSub(idx,:),lambda1_simuSub(idx,:),ModSub,'linear','extrap');
    %vObs_interp = interp1(meanMod_allSub(idx,:),lambda2_simuSub(idx,:),ModSub,'linear','extrap');
    
    % putting final values in the correct order
    [aidx bidx] = ismember([coordObsSub repmat(uniyrmoda(i),length(coordObsSub),1)],[coordObs dnyrmoda],'rows');
    mObs_calc(bidx) = mObs_interp;
        
end

% save results
save(sprintf('simulateObs_%d.mat',yr2simu),'coordObs','yrmodaObs','Mod','mObs_calc', ...
    'css','cGen','lambda1_gen','meanMod_all','perctile_data_all');

end