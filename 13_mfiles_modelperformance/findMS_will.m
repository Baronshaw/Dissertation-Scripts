function [] = findMS_will()
% this function will take select grid locations and find these monitoring
% stations used to construct their SCurves

% finding the exact locations to do the comparisons
SCurveLoc = { [-475000,470000],[-650000,480000] ; [1100000,-780000],[975000,-700000] ; ...
    [-170000,-1350000],[-170000,-1200000] ; [-550000,-1040000],[-800000,-940000] ; ...
    [-1425000,-750000],[-1425000,-600000] ; [-1950000,-400000],[-1800000,-400000] };


% load paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',2001));
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
dayz = datenum([yr mo da]);
unidayz = unique(dayz);
uniCoord = unique(coordObs,'rows');

% find mointoring stations closest to each grid
for i = 1:length(SCurveLoc)
    
    % first of the pair
    distz = sqrt( (SCurveLoc{i,1}(1)-uniCoord(:,1)).^2 + (SCurveLoc{i,1}(2)-uniCoord(:,2)).^2 );
    [sorted sortidx] = sort(distz);
    sortloc = uniCoord(sortidx,:);
    len = 0;
    for j = 1:3
        len = len + sum(sortloc(j,1)==coordObs(:,1) & sortloc(j,2)==coordObs(:,2));
    end
    while len <= 150
        j = j + 1;
        len = len + sum(sortloc(j,1)==coordObs(:,1) & sortloc(j,2)==coordObs(:,2));
    end
    closeMS{i,1}{1} = sortloc(1:j,:);
    
    % second of the pair
    distz = sqrt( (SCurveLoc{i,2}(1)-uniCoord(:,1)).^2 + (SCurveLoc{i,2}(2)-uniCoord(:,2)).^2 );
    [sorted sortidx] = sort(distz);
    sortloc = uniCoord(sortidx,:);
    len = 0;
    for j = 1:3
        len = len + sum(sortloc(j,1)==coordObs(:,1) & sortloc(j,2)==coordObs(:,2));
    end
    while len <= 150
        j = j + 1;
        len = len + sum(sortloc(j,1)==coordObs(:,1) & sortloc(j,2)==coordObs(:,2));
    end
    closeMS{i,1}{2} = sortloc(1:j,:);
        
end

save('temp.mat','SCurveLoc','closeMS');

end