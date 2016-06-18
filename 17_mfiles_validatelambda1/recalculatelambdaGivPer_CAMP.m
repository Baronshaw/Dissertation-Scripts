function [] = recalculatelambdaGivPer_CAMP()
% this function will perform a recalculation of CAMP and find lambda1 and
% lambda2 for modeled values across percentiles from generated lambda1s and
% lambda2s for 2001

% bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "recalculatelambdaGivPer_CAMP" -logfile "recalculatelambdaGivPer_CAMP.out"

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
perctile_recalc = repmat(prctile(Mod,0:10:100),length(css),1);

% calculate lambda1 and lambda2 across percentiles

%%% recalculation

% mean/variance in each bin
prctbins = prctile(Mod,0:10:100);
for i = 1:length(prctbins)-1
    mean_Mod(i,1) = mean(Mod(Mod>=prctbins(i)&Mod<prctbins(i+1)));
    mean_Obs(i,1) = mean(mObs_calc(Mod>=prctbins(i)&Mod<prctbins(i+1)));
    var_Obs(i,1) = var(mObs_calc(Mod>=prctbins(i)&Mod<prctbins(i+1)));
end
meanMod_recalc = repmat(mean_Mod',length(css),1);
lambda1_recalc = repmat(mean_Obs',length(css),1);
lambda2_recalc = repmat(var_Obs',length(css),1);

%%% simulation

% mean/variance in each bin
prctbins = prctile(Mod,0:10:100);
for i = 1:length(prctbins)-1
    mean_Mod(i,1) = mean(Mod(Mod>=prctbins(i)&Mod<prctbins(i+1)));
    mean_Obs(i,1) = mean(Obssimu_calc(Mod>=prctbins(i)&Mod<prctbins(i+1)));
    var_Obs(i,1) = var(Obssimu_calc(Mod>=prctbins(i)&Mod<prctbins(i+1)));
end
lambda1_recalcsimu = repmat(mean_Obs',length(css),1);
lambda2_recalcsimu = repmat(var_Obs',length(css),1);

% save results
save('recalculatelambdaGivPer_CAMP.mat','css','lambda1_recalc', ...
    'lambda2_recalc','lambda1_recalcsimu','lambda2_recalcsimu', ...
    'perctile_recalc','meanMod_recalc');

end