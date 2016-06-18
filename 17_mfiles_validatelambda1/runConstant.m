function [] = runConstant()
% this function will run a constant error correction for CMAQ

modplots = 0:5:50;

% gather all modeled and obs for 2001
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;

% load all modeled data for 2001
load(sprintf('../matfiles/prepCTM_%d.mat',2001));
yrv = yrmodaCTMv(:,1);
mov = yrmodaCTMv(:,2);
dav = yrmodaCTMv(:,3);

%%% domain-wide constant error correction

meanGivMod = repmat(dailyCTMv - mean(Mod-Obs),1,1+length(modplots));
varGivMod = var(Mod)*ones(length(dailyCTMv),1+length(modplots));

% save results
save('Constantmethod.mat','modplots','meanGivMod','varGivMod', ...
    'dailyCTMv','distCTMv','yrmodaCTMv');

%%% region-wide constant error correction

% load results about region information
load('CAMPmethod_regional.mat');
% use 'inregion6', 'inregion6v'

meanGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
for i = 1:length(unique(inregion6))    
    idx = inregion6 == i;
    idxv = inregion6v == i;
    meanGivMod(idxv,:) = repmat(dailyCTMv(idxv) - mean(Mod(idx)-Obs(idx)),1,1+length(modplots));
    varGivMod(idxv,:) = var(Mod(idx)).*ones(sum(idxv),1+length(modplots));   
end

% save results
save('ConstantRegionmethod.mat','modplots','meanGivMod','varGivMod', ...
    'dailyCTMv','distCTMv','yrmodaCTMv','inregion6','inregion6v');

%%% season-wide constant error correction

meanGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
for i = 1:4
    
    if i == 1
        idx = mo == 1 | mo == 2 | mo == 12;
        idxv = mov == 1 | mov == 2 | mov == 12;
    elseif i == 2
        idx = mo == 3 | mo == 4 | mo == 5;
        idxv = mov == 3 | mov == 4 | mov == 5;
    elseif i == 3
        idx = mo == 6 | mo == 7 | mo == 8;
        idxv = mov == 6 | mov == 7 | mov == 8;
    else 
        idx = mo == 9 | mo == 10 | mo == 11;
        idxv = mov == 9 | mov == 10 | mov == 11;
    end  
    
    meanGivMod(idxv,:) = repmat(dailyCTMv(idxv) - mean(Mod(idx)-Obs(idx)),1,1+length(modplots));
    varGivMod(idxv,:) = var(Mod(idx)).*ones(sum(idxv),1+length(modplots)); 
    
end

% save results
save('ConstantSeasonmethod.mat','modplots','meanGivMod','varGivMod', ...
    'dailyCTMv','distCTMv','yrmodaCTMv');

end