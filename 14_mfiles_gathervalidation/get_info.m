function [PMs,EPAall] = get_info(ck)
% this function will get basic information about a data file including:
% 1) season of a given datum
% 2) EPA region of a location
% 3) if a location is rural/urban/suburban

% initialize parameters
if nargin < 1, ck = rand(100,3); end

% load dates
[yr mo da] = datevec(ck(:,3));

% season: 1 = winter, 2 = spring, 3 = summer, 4 = fall
PMs = NaN*ones(size(ck,1),1);
idx = mo == 1 |mo == 2 | mo == 12;
PMs(idx) = 1;
idx = mo == 3 |mo == 4 | mo == 5;
PMs(idx) = 2;
idx = mo == 6 |mo == 7 | mo == 8;
PMs(idx) = 3;
idx = mo == 9 |mo == 10 | mo == 11;
PMs(idx) = 4;

% EPA region
% go through original station locations
for i = 1999:2010
    load(sprintf('../datafiles/Observed_PM2p5/MasterDaily_PM2p5_%d.mat',i));
    alllon{i-1998,1} = longitude;
    alllat{i-1998,1} = latitude;
    allID{i-1998,1} = location;
end
alllon = cell2mat(alllon); alllat = cell2mat(alllat); allID = cell2mat(allID);
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
sIDall = NaN*ones(size(ck,1),1);
for i = 1:length(projuniID)
    dists = sqrt((projuniID(i,1)-ck(:,1)).^2+(projuniID(i,2)-ck(:,2)).^2);
    idx = dists < 10;
    sIDall(idx) = uniID(i);
end
stateall = floor(sIDall./10^7);
cd('14_mfiles_gathervalidation');
% match id's with EPA region
% from https://aqs.epa.gov/aqsweb/codes/data/StateCountyCodes.csv
M = csvread('../05_mfiles_crossvalidation/StateCountyCodes_mod.csv');
temp = unique(M,'rows');
stateIDs = temp(:,1);
EPAregionsIDs = temp(:,2);
EPAall = NaN*ones(size(ck,1),1);
for i = 1:length(stateIDs)
    idx = stateIDs(i) == stateall;
    EPAall(idx) = EPAregionsIDs(i);
end

end