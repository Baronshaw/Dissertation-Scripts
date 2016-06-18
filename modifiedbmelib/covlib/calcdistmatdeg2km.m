function distmat = calcdistmatdeg2km(lonlat1,lonlat2)

% calcdistmatdeg2km - Calculate distance matrix in km based on Lat/Lon
%                     (Nov 28,2012)
%
% distmat = calcdistmatdeg2km(lonlat1,lonlat2)
%
% INPUT:
%
% lonlat1(n by 2): Spatial coordinates 1
%                  1st column is longitude (x)
%                  2nd column is latitude (y)
% lonlat2(m by 2): Spatial coordinates 2
%                  1st column is longitude (x)
%                  2nd column is latitude (y)
%
% OUTPUT :
%
% distmat(n by m): Distance(km) matrix
%

nloc1 = size(lonlat1,1);
nloc2 = size(lonlat2,1);
distmat = zeros(nloc1,nloc2);

for i = 1:nloc1
    dist = distance([lonlat1(i,2),lonlat1(i,1)],...
                    [lonlat2(:,2),lonlat2(:,1)]);
    distmat(i,:) = deg2km(dist');
end
