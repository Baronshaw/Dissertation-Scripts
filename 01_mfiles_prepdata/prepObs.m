function [coordObs Obs Mod yrmodaObs] = prepObs(years)
% this function will prepare the observational data for years that we do
% not have corresponding modeled data

if nargin < 1, years = 1999; end

% load the paired modeled and observed data
load(sprintf('../datafiles/Paired_PM2p5/MasterDaily_PM2p5_%d.mat',years));
Obs = values;
Mod = [];
yrmodaObs = dates;

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
coordObs = ell2lambertcc([longitude latitude],'whiproj2001'); 
cd ../01_mfiles_prepdata

end