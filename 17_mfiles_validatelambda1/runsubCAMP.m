function [] = runsubCAMP()
% this function will perform a regionalized CAMP method

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

%%% getting location information

% state shapefile from: https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
% see which points are in which states
if ~exist('states.mat')
    allstates = shaperead('gz_2010_us_040_00_500k/gz_2010_us_040_00_500k.shp');
    for k = 1:length(allstates)
        tic
        cd ../09_mfiles_projections
        allstates_p{k,1} = ell2lambertcc([allstates(k).X',allstates(k).Y'],'whiproj2001');
        cd ../17_mfiles_validatelambda1
        instates{k,1} = inpolygon(coordObs(:,1),coordObs(:,2),allstates_p{k,1}(:,1),allstates_p{k,1}(:,2));
        instatesv{k,1} = inpolygon(distCTMv(:,1),distCTMv(:,2),allstates_p{k,1}(:,1),allstates_p{k,1}(:,2));
        toc
        namestates{k,1} = allstates(k).NAME;
    end
    save('states.mat','allstates_p','instates','instatesv','namestates');
else
    load('states.mat')
end

% indexes of regions
idxnortheast = [ find(strcmp(namestates,'Maine')) ; find(strcmp(namestates,'New Hampshire')) ; ...
    find(strcmp(namestates,'Vermont')) ; find(strcmp(namestates,'New York')) ; ...
    find(strcmp(namestates,'Massachusetts')) ; find(strcmp(namestates,'Rhode Island')) ; ...
    find(strcmp(namestates,'Connecticut')) ; find(strcmp(namestates,'New Jersey')) ; ...
    find(strcmp(namestates,'Delaware')) ; find(strcmp(namestates,'Maryland')) ; ...
    find(strcmp(namestates,'Pennsylvania')) ; find(strcmp(namestates,'District of Columbia')) ];
idxsoutheast = [ find(strcmp(namestates,'Virginia')) ; find(strcmp(namestates,'West Virginia')) ; ...
    find(strcmp(namestates,'Tennessee')) ; find(strcmp(namestates,'North Carolina')) ; ...
    find(strcmp(namestates,'South Carolina')) ; find(strcmp(namestates,'Mississippi')) ; ...
    find(strcmp(namestates,'Alabama')) ; find(strcmp(namestates,'Georgia')) ; ...
    find(strcmp(namestates,'Florida')) ];
idxuppermidwest = [ find(strcmp(namestates,'Minnesota')) ; find(strcmp(namestates,'Iowa')) ; ...
    find(strcmp(namestates,'Missouri')) ; find(strcmp(namestates,'Kentucky')) ; ...
    find(strcmp(namestates,'Ohio')) ; find(strcmp(namestates,'Indiana')) ; ...
    find(strcmp(namestates,'Illinois')) ; find(strcmp(namestates,'Michigan')) ; ...
    find(strcmp(namestates,'Wisconsin')) ];
idxlowermidwest = [ find(strcmp(namestates,'Arkansas')) ; find(strcmp(namestates,'Louisiana')) ; ...
    find(strcmp(namestates,'Oklahoma')) ; find(strcmp(namestates,'Texas')) ];
idxrockymountains = [ find(strcmp(namestates,'Montana')) ; find(strcmp(namestates,'North Dakota')) ; ...
    find(strcmp(namestates,'South Dakota')) ; find(strcmp(namestates,'Nebraska')) ; ...
    find(strcmp(namestates,'Kansas')) ; find(strcmp(namestates,'Colorado')) ; ...
    find(strcmp(namestates,'New Mexico')) ; find(strcmp(namestates,'Arizona')) ; ...
    find(strcmp(namestates,'Utah')) ; find(strcmp(namestates,'Wyoming')) ; ...
    find(strcmp(namestates,'Idaho')) ; find(strcmp(namestates,'Nevada')) ];
idxpacificcoast = [ find(strcmp(namestates,'Washington')) ; find(strcmp(namestates,'Oregon')) ; ...
    find(strcmp(namestates,'California')) ];

%%% create final index of regions
inregion6 = NaN*ones(length(Mod),1);
inregion6v = NaN*ones(length(dailyCTMv),1);
strregion6 = {'northeast';'southeast';'uppermidwest';'lowermidwest';'rockymountains';'pacificcoast'};

%%% by regions: northeast, southeast, upper midwest, lower midwest, rocky mountains, pacific coast
% regions chosen from: 'Assessment of bias-adjusted PM2.5 air quality 
% forecasts over the continental United States during 2007'

meanGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
[r c] = size(meanGivMod);

for i = 1:6
    disp(i);
    % subsetting to the correct region
    if i == 1
        len = length(idxnortheast);
        idx = 0; idxv = 0;
        for j = 1:len
            idx = idx + instates{idxnortheast(j)};
            idxv = idxv + instatesv{idxnortheast(j)};
        end
        idx = logical(idx); idxv = logical(idxv);
        inregion6(idx) = 1; inregion6v(idxv) = 1;
    elseif i == 2
        len = length(idxsoutheast);
        idx = 0; idxv = 0;
        for j = 1:len
            idx = idx + instates{idxsoutheast(j)};
            idxv = idxv + instatesv{idxsoutheast(j)};
        end
        idx = logical(idx); idxv = logical(idxv);        
        inregion6(idx) = 2; inregion6v(idxv) = 2;
    elseif i == 3
        len = length(idxuppermidwest);
        idx = 0; idxv = 0;
        for j = 1:len
            idx = idx + instates{idxuppermidwest(j)};
            idxv = idxv + instatesv{idxuppermidwest(j)};
        end
        idx = logical(idx); idxv = logical(idxv);
        inregion6(idx) = 3; inregion6v(idxv) = 3;
    elseif i == 4
        len = length(idxlowermidwest);
        idx = 0; idxv = 0;
        for j = 1:len
            idx = idx + instates{idxlowermidwest(j)};
            idxv = idxv + instatesv{idxlowermidwest(j)};
        end
        idx = logical(idx); idxv = logical(idxv);
        inregion6(idx) = 4; inregion6v(idxv) = 4;
    elseif i == 5
        len = length(idxrockymountains);
        idx = 0; idxv = 0;
        for j = 1:len
            idx = idx + instates{idxrockymountains(j)};
            idxv = idxv + instatesv{idxrockymountains(j)};
        end
        idx = logical(idx); idxv = logical(idxv);
        inregion6(idx) = 5; inregion6v(idxv) = 5;
    else
        len = length(idxpacificcoast);
        idx = 0; idxv = 0;
        for j = 1:len
            idx = idx + instates{idxpacificcoast(j)};
            idxv = idxv + instatesv{idxpacificcoast(j)};
        end
        idx = logical(idx); idxv = logical(idxv);
        inregion6(idx) = 6; inregion6v(idxv) = 6;       
    end  
    
    % 10 percentile bins
    prctbins = prctile(Mod(idx),0:10:100);

    % mean/variance in each bin
    for j = 1:length(prctbins)-1
        mean_Mod{i}(j,1) = mean(Mod(Mod>=prctbins(j)&Mod<prctbins(j+1)&idx));
        mean_Obs{i}(j,1) = mean(Obs(Mod>=prctbins(j)&Mod<prctbins(j+1)&idx));
        var_Obs{i}(j,1) = var(Obs(Mod>=prctbins(j)&Mod<prctbins(j+1)&idx));
    end

    % calculate lambda1 and lambda2 for all modeled values
    for j = 1:c
        if j == 1
            meanGivMod(idxv,j) = interp1(mean_Mod{i},mean_Obs{i},dailyCTMv(idxv),'linear','extrap');
            blah = interp1(mean_Mod{i},mean_Obs{i},dailyCTMv(idxv),'linear','extrap');
            disp(sum(isnan(blah)));
            varGivMod(idxv,j) = interp1(mean_Mod{i},var_Obs{i},dailyCTMv(idxv),'linear','extrap');
        else
            meanGivMod(idxv,j) = interp1(mean_Mod{i},mean_Obs{i},repmat(modplots(j-1),sum(idxv),1),'linear','extrap');
            varGivMod(idxv,j) = interp1(mean_Mod{i},var_Obs{i},repmat(modplots(j-1),sum(idxv),1),'linear','extrap');
        end
    end
    
end

% if values are negative, make them zero, but show S-curves
meanGivMod(meanGivMod<0) = 0;

% save results
save('CAMPmethod_regional.mat','modplots','meanGivMod','varGivMod', ...
    'mean_Mod','mean_Obs','var_Obs','dailyCTMv','distCTMv','yrmodaCTMv', ...
    'inregion6','inregion6v','strregion6');

%%%  CAMP by seasons
meanGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMv),1+length(modplots));
[r c] = size(meanGivMod);
for i = 1:4
    disp(i);
    % subsetting to the correct season
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
    
    % 10 percentile bins
    prctbins = prctile(Mod,0:10:100);

    % mean/variance in each bin
    for j = 1:length(prctbins)-1
        mean_Mod{i}(j,1) = mean(Mod(Mod>=prctbins(j)&Mod<prctbins(j+1)&idx));
        mean_Obs{i}(j,1) = mean(Obs(Mod>=prctbins(j)&Mod<prctbins(j+1)&idx));
        var_Obs{i}(j,1) = var(Obs(Mod>=prctbins(j)&Mod<prctbins(j+1)&idx));
    end

    % calculate lambda1 and lambda2 for all modeled values
    for j = 1:c
        if j == 1
            meanGivMod(idxv,j) = interp1(mean_Mod{i},mean_Obs{i},dailyCTMv(idxv),'linear','extrap');
            varGivMod(idxv,j) = interp1(mean_Mod{i},var_Obs{i},dailyCTMv(idxv),'linear','extrap');
        else
            meanGivMod(idxv,j) = interp1(mean_Mod{i},mean_Obs{i},repmat(modplots(j-1),sum(idxv),1),'linear','extrap');
            varGivMod(idxv,j) = interp1(mean_Mod{i},var_Obs{i},repmat(modplots(j-1),sum(idxv),1),'linear','extrap');
        end
    end
    
end

% if values are negative, make them zero, but show S-curves
meanGivMod(meanGivMod<0) = 0;

% save results
save('CAMPmethod_seasonal.mat','modplots','meanGivMod','varGivMod', ...
    'mean_Mod','mean_Obs','var_Obs','dailyCTMv','distCTMv','yrmodaCTMv');

end