function [] = get_info(yrz)
% this function will get basic information about a data file including:
% 1) which percentile (of 10 deciles) that each datum falls into
% 2) season of a given datum
% 3) EPA region of a location
% 4) monitoring network of a given location
% 5) if a location is rural/urban/suburban

if nargin < 1, yrz = 2001; end

% load CMAQ paired data
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yrz));
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;

% percentile by observed data values
PMdec = NaN*ones(length(Obs),1);
PMd = prctile(Obs,0:10:100);
for i = 1:length(PMd)-1
    idx = Obs >= PMd(i) & Obs < PMd(i+1);
    PMdec(idx) = i;
end
idx = Obs == PMd(end); PMdec(idx) = length(PMd) - 1;

% season: 1 = winter, 2 = spring, 3 = summer, 4 = fall
PMs = NaN*ones(length(Obs),1);
str = {'spring','summer','fall','winter'};
for i = 1:length(str)
    for j = yrz-1:yrz
        [YY(i,j-1999) MO(i,j-1999) DD(i,j-1999) dum dum] = getSeasonStartDate(str{i},j);
    end
end
idx = datenum(yr,mo,da) >= datenum(YY(4,1),MO(4,1),DD(4,1)) & datenum(yr,mo,da) < datenum(YY(1,2),MO(1,2),DD(1,2));
PMs(idx) = 1;
idx = datenum(yr,mo,da) >= datenum(YY(4,2),MO(4,2),DD(4,2));
PMs(idx) = 1;
idx = datenum(yr,mo,da) >= datenum(YY(1,2),MO(1,2),DD(1,2)) & datenum(yr,mo,da) < datenum(YY(2,2),MO(2,2),DD(2,2));
PMs(idx) = 2;
idx = datenum(yr,mo,da) >= datenum(YY(2,2),MO(2,2),DD(2,2)) & datenum(yr,mo,da) < datenum(YY(3,2),MO(3,2),DD(3,2));
PMs(idx) = 3;
idx = datenum(yr,mo,da) >= datenum(YY(3,2),MO(3,2),DD(3,2)) & datenum(yr,mo,da) < datenum(YY(4,2),MO(4,2),DD(4,2));
PMs(idx) = 4;

% EPA region
% go through original station locations
load(sprintf('../datafiles/Observed_PM2p5/MasterDaily_PM2p5_%d.mat',yrz));
alllon = longitude;
alllat = latitude;
allID = location;
[uniloc uniidx] = unique([alllon alllat],'rows');
unilon = uniloc(:,1);
unilat = uniloc(:,2);
uniID = allID(uniidx);

% convert them to coordinate system
load('../05_mfiles_crossvalidation/projexample.mat');
cd('../09_mfiles_projections');
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load Ellipsoids.mat
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
projuniID = ell2lambertcc([unilon unilat],'whiproj2001'); 
cd ..
% match estimations locations with station id
sIDall = NaN*ones(size(coordObs,1),1);
for i = 1:length(projuniID)
    dists = sqrt((projuniID(i,1)-coordObs(:,1)).^2+(projuniID(i,2)-coordObs(:,2)).^2);
    idx = dists < 10;
    sIDall(idx) = uniID(i);
end
stateall = floor(sIDall./10^7);
cd('13_mfiles_modelperformance');
% match id's with EPA region
% from https://aqs.epa.gov/aqsweb/codes/data/StateCountyCodes.csv
M = csvread('../05_mfiles_crossvalidation/StateCountyCodes_mod.csv');
temp = unique(M,'rows');
stateIDs = temp(:,1);
EPAregionsIDs = temp(:,2);
EPAall = NaN*ones(size(coordObs,1),1);
for i = 1:length(stateIDs)
    idx = stateIDs(i) == stateall;
    EPAall(idx) = EPAregionsIDs(i);
end

% rural/urban/suburban: under 'Location'
% rural = 2, suburban = 3, urban = 5
unilocation = locations(uniidx);
PMl = NaN*ones(size(coordObs,1),1);
for i = 1:length(projuniID)
    dists = sqrt((projuniID(i,1)-coordObs(:,1)).^2+(projuniID(i,2)-coordObs(:,2)).^2);
    idx = dists < 10;
    PMl(idx) = unilocation(i);
end

% monitoring network: 'Rank' 1-4,6 = AQS, 'Rank' 5 = IMPROVE
uniranks = ranks(uniidx);
PMr = NaN*ones(size(coordObs,1),1);
for i = 1:length(projuniID)
    dists = sqrt((projuniID(i,1)-coordObs(:,1)).^2+(projuniID(i,2)-coordObs(:,2)).^2);
    idx = dists < 10;
    PMr(idx) = uniranks(i);
end

% save results
save('matfiles/info.mat','PMd','PMdec','PMs','EPAall','PMl','PMr','yr','mo','da');

end