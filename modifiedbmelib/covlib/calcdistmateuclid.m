function distmat = calcdistmateuclid(sptlcoord1,sptlcoord2)

% calcdistmateuclid - Calculate euclidian distance matrix (Nov 28,2012)
%
% distmat = calcdistmateuclid(sptlcoord1,sptlcoord2)
%
% INPUT:
%
% sptlcoord1(n by 2): Spatial Coordinate 1
%                     1st column is x
%                     2nd column is y
%
% sptlcoord2(m by 2): Spatial Coordinate 2
%                     1st column is x
%                     2nd column is y
%
% OUTPUT :
%
% distmat(n by m): Distance matrix
%

nloc1 = size(sptlcoord1,1);
nloc2 = size(sptlcoord2,1);
distmat = zeros(nloc1,nloc2);

for i = 1:nloc1
    distmat(i,:) = sqrt((sptlcoord1(i,1)-sptlcoord2(:,1)).^2 + ...
                        (sptlcoord1(i,2)-sptlcoord2(:,2)).^2)';
end
