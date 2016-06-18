function [coordObs Obs Mod yrmodaObs] = prepCTMandObs(years)

if nargin < 1, years = 2001; end

% load the paired modeled and observed data
if years == 2006
    fid = fopen(sprintf('../datafiles/Paired_%s/%d/sitecmp_%d_12k_spr06_%s.txt', ...
        'PM2p5',years,years,'PM2p5'));
else
    fid = fopen(sprintf('../datafiles/Paired_%s/%d/sitecmp_%d_36k_%s.txt', ...
        'PM2p5',years,years,'PM2p5'));
end
ObsVMod = textscan(fid,'%f%f%f%f%f%s%s%f%f', 'delimiter', '\t');
fclose(fid);
SiteID = ObsVMod{1};
LatitudePaired = ObsVMod{2};
LongitudePaired = ObsVMod{3};
ColumnPaired = ObsVMod{4};
RowPaired = ObsVMod{5};
TimeOnPaired = ObsVMod{6};
TimeOffPaired = ObsVMod{7};
Obs = ObsVMod{8};
Mod = ObsVMod{9};

TimeOnPaired = datevec(char(TimeOnPaired));
TimeOffPaired = datevec(char(TimeOffPaired));
TimeOnPaired = TimeOnPaired(:,1)*10^4 + TimeOnPaired(:,2)*10^2 + TimeOnPaired(:,3);
TimeOffPaired = TimeOffPaired(:,1)*10^4 + TimeOffPaired(:,2)*10^2 + TimeOffPaired(:,3);
yrmodaObs = TimeOnPaired;

% projection location data
cd ../09_mfiles_projections
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load Ellipsoids.mat
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
coordObs = ell2lambertcc([LongitudePaired LatitudePaired],'whiproj2001'); 
cd ../01_mfiles_prepdata

end