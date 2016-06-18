function distmat = calcdistmattemp(tempcoord1,tempcoord2)

% calcdistmattemp - Calculate temporal distance matrix (Dec 05,2012)
%
% distmat = calcdistmattemp(tempcoord1,tempcoord2)
%
% INPUT:
%
% tempcoord1(1 by n): Temporal location1 
% tempcoord2(1 by m): Temporal location2 
%
% OUTPUT :
%
% distmat(n by m): Temporal distance matrix
%

nloc1 = size(tempcoord1,2);
nloc2 = size(tempcoord2,2);
distmat = zeros(nloc1,nloc2);

if ~isempty(distmat)
    for i = 1:nloc1
        distmat(i,:) = abs(tempcoord1(i) - tempcoord2)';
    end
end
