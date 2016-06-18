function [] = recalculatelambdaGivPer()
% this function will perform a recalculation of RAMP and find lambda1 and
% lambda2 for modeled values across percentiles from generated lambda1s and
% lambda2s for 2001

% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "recalculatelambdaGivPer" -logfile "recalculatelambdaGivPer.out"

% loading BME function
cd ../BMELIB2.0b
startup
cd ../17_mfiles_validatelambda1

% load generated Obs for 2001 and 2002
load('simulateObs_2001.mat')
Mod1 = Mod; coordObs1 = coordObs; yrmodaObs1 = yrmodaObs; 
mObs_calc1 = mObs_calc; Obssimu_calc1 = Obssimu_calc;
css1 = css; 
load('simulateObs_2002.mat')
Mod2 = Mod; coordObs2 = coordObs; yrmodaObs2 = yrmodaObs; 
mObs_calc2 = mObs_calc; Obssimu_calc2 = Obssimu_calc;
Mod = [Mod1;Mod2]; coordObs = [coordObs1;coordObs2]; 
yrmodaObs = [yrmodaObs1;yrmodaObs2]; mObs_calc = [mObs_calc1;mObs_calc2];
Obssimu_calc = [Obssimu_calc1;Obssimu_calc2];
css = css1; 
clear Mod1 Mod2 coordObs1 coordObs2 yrmodaObs1 yrmodaObs2 
clear mObs_calc1 mObs_calc2 Obssimu_calc1 Obssimu_calc2
clear css1 vObs_calc modplots

% convert from long to wide format
if ~exist('wideformat.mat')
    unicoordObs = unique(coordObs,'rows');
    uniyrmoda = unique(yrmodaObs);
    Modw = NaN*ones(length(unicoordObs),length(uniyrmoda));
    mObs_calcw = NaN*ones(length(unicoordObs),length(uniyrmoda));
    Obssimu_calcw = NaN*ones(length(unicoordObs),length(uniyrmoda));
    for i = 1:length(unicoordObs)
        idx = coordObs(:,1) == unicoordObs(i,1) & coordObs(:,2) == unicoordObs(i,2);
        coordObsSub = coordObs(idx,:);
        yrmodaObsSub = yrmodaObs(idx);
        ModSub = Mod(idx); 
        mObs_calcSub = mObs_calc(idx,:); 
        Obssimu_calcSub = Obssimu_calc(idx,:);
        [aidx bidx] = ismember(uniyrmoda,yrmodaObsSub);
        Modw(i,aidx) = ModSub;
        mObs_calcw(i,aidx) = mObs_calcSub;
        Obssimu_calcw(i,aidx) = Obssimu_calcSub;
    end
    save('wideformat.mat','unicoordObs','uniyrmoda','Modw','mObs_calcw','Obssimu_calcw');
else
    load('wideformat.mat')
end

% convert dates
yr = floor(uniyrmoda./10000);
mo = floor((uniyrmoda - yr.*10000)./100);
da = uniyrmoda - yr.*10000 - mo.*100;
dayrmoda = datenum(yr,mo,da);

% initialize final matrix
lambda1_recalc = NaN*ones(length(css),10);
lambda2_recalc = NaN*ones(length(css),10);
lambda1_recalcsimu = NaN*ones(length(css),10);
lambda2_recalcsimu = NaN*ones(length(css),10);
perctile_recalc = NaN*ones(length(css),11);
meanMod_recalc = NaN*ones(length(css),10);

% loop through each css locations for 2001 and calculate lambda1 and
% lambda2 across percentiles
unicss = unique(css(:,3));
for i = 1:length(unicss)
    tic
    disp(i);
    % subsetting all CMAQ grids to an individual day
    idx = unicss(i) == css(:,3);    
    cssSub = css(idx,:);

    % subsetting paired data to +/- 180 days
    idx = dayrmoda >= unicss(i)-180 & dayrmoda <= unicss(i)+180;
    ModSub = Modw(:,idx); 
    mObs_calcSub = mObs_calcw(:,idx);
    Obssimu_calcSub = Obssimu_calcw(:,idx);

    % generated lambda1 and lambda2
    lambda1_recalculated = NaN*ones(length(cssSub),10);
    lambda2_recalculated = NaN*ones(length(cssSub),10);
    meanMod_recalculated = NaN*ones(length(cssSub),10);
    
    % finding number of closest neighbors
    k = 3;
    [idxneib, dist] = knnsearch(unicoordObs,cssSub(:,1:2),'K',k);
    r = length(idxneib);
    pntsum = cell2mat( arrayfun( @(x) sum(sum(~isnan(ModSub(idxneib(x,:),:)))), 1:r,'UniformOutput',false) );
    ksize = k*ones(r,1);
    while sum(pntsum>=150) < r
        k = k + 1;
        idxpnt = pntsum < 150;
        [idxneib, dist] = knnsearch(unicoordObs,cssSub(:,1:2),'K',k);        
        pntsumnew = cell2mat( arrayfun( @(x) sum(sum(~isnan(ModSub(idxneib(x,:),:)))), 1:r,'UniformOutput',false) );
        pntsum(idxpnt) = pntsumnew(idxpnt);
        ksize(idxpnt) = k;
    end
    
    % finding closest neighbors
    for j = 3:k
        idx = j == ksize;
        idxneib(idx,j+1:k) = NaN;
    end

    % getting all the modeled data by grid  
    tempa = arrayfun( @(x) ModSub(idxneib(x,~isnan(idxneib(x,:))),:), 1:length(idxneib), 'UniformOutput', false );
    tempb = cellfun( @(x) x(:), tempa, 'UniformOutput', false );
    Modcell = cellfun( @(x) x(~isnan(x)), tempb, 'UniformOut', false );

    % getting all the calculated lambda1s by grid
    tempa = arrayfun( @(x) mObs_calcSub(idxneib(x,~isnan(idxneib(x,:))),:), 1:length(idxneib), 'UniformOutput', false );
    tempb = cellfun( @(x) x(:), tempa, 'UniformOutput', false );
    mObs_calccell = cellfun( @(x) x(~isnan(x)), tempb, 'UniformOut', false );
    
    % getting all the calculated lambda1s by grid through the simulation
    tempa = arrayfun( @(x) Obssimu_calcSub(idxneib(x,~isnan(idxneib(x,:))),:), 1:length(idxneib), 'UniformOutput', false );
    tempb = cellfun( @(x) x(:), tempa, 'UniformOutput', false );
    Obssimu_calccell = cellfun( @(x) x(~isnan(x)), tempb, 'UniformOut', false );
    
    % coming up with percentiles of modeled data
    perctile_dataSub = cell2mat( cellfun( @(x) prctile(x,0:10:100), Modcell, 'UniformOut', false )' );

    % loop through each of 10 percentiles to have all the lambda1s and lambda2s
    for j = 1:10
        idxper = arrayfun( @(x) Modcell{x}>=perctile_dataSub(x,j) & Modcell{x}<perctile_dataSub(x,j+1), ...
            1:length(idxneib), 'UniformOutput', false );
        meanMod_cell = arrayfun( @(x) mean(Modcell{x}(idxper{x})), ...
            1:length(idxneib), 'UniformOutput', false );
        lambda1_cell = arrayfun( @(x) mean(mObs_calccell{x}(idxper{x})), ...
            1:length(idxneib), 'UniformOutput', false );
        lambda2_cell = arrayfun( @(x) var(mObs_calccell{x}(idxper{x})), ...
            1:length(idxneib), 'UniformOutput', false );
        lambda1_cellsimu = arrayfun( @(x) mean(Obssimu_calccell{x}(idxper{x})), ...
            1:length(idxneib), 'UniformOutput', false );
        lambda2_cellsimu = arrayfun( @(x) var(Obssimu_calccell{x}(idxper{x})), ...
            1:length(idxneib), 'UniformOutput', false );
        meanMod_recalculated(:,j) = cell2mat(meanMod_cell);
        lambda1_recalculated(:,j) = cell2mat(lambda1_cell);
        lambda2_recalculated(:,j) = cell2mat(lambda2_cell);
        lambda1_recalculatedsimu(:,j) = cell2mat(lambda1_cellsimu);
        lambda2_recalculatedsimu(:,j) = cell2mat(lambda2_cellsimu);
    end

    % put recalculated lambda1 and lambda2 in place
    [aidx bidx] = ismember(cssSub,css,'rows');
    lambda1_recalc(bidx,:) = lambda1_recalculated;
    lambda2_recalc(bidx,:) = lambda2_recalculated;
    lambda1_recalcsimu(bidx,:) = lambda1_recalculatedsimu;
    lambda2_recalcsimu(bidx,:) = lambda2_recalculatedsimu;
    perctile_recalc(bidx,:) = perctile_dataSub;
    meanMod_recalc(bidx,:) = meanMod_recalculated;
    
    % temp save
    if mod(i,20)==0 
        save('temp_recalculatelambdaGivPer.mat','css','lambda1_recalc', ...
            'lambda2_recalc','lambda1_recalcsimu','lambda2_recalcsimu', ...
            'perctile_recalc','meanMod_recalc');
    end
    toc
end % each unique day

% save results
save('recalculatelambdaGivPer.mat','css','lambda1_recalc', ...
    'lambda2_recalc','lambda1_recalcsimu','lambda2_recalcsimu', ...
    'perctile_recalc','meanMod_recalc');

end