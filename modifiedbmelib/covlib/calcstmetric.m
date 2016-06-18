function stmetric = calcstmetric(sill,sptlrange,temprange)

% calcstmetric - Calculate Space/Time metric based on sill, spatial range,
%                and temporal range (Nov 29,2012)
%
% stmetric = calcstmetric(sill,sptlrange,temprange)
%
% INPUT:
%
% sill(1 by n): Partial sill
% sptlrange(1 by n): Spatial range
% temprange(1 by n): Temporal range
%
% OUTPUT :
%
% stmetric(scalar): Space/Time metric
%

scovrange = 0;
tcovrange = 0;

totalsill = sum(sill);

for i = 1:size(sill,2)
    scovrange = scovrange + (sill(i)/totalsill) * sptlrange(i);
    tcovrange = tcovrange + (sill(i)/totalsill) * temprange(i);
end

stmetric = scovrange/tcovrange;
