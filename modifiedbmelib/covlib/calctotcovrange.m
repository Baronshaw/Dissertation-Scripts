function totcovrange = calctotcovrange(sill,covrange)

% calcsptlrange - Calculate total spatial covarance range (Dec 10,2012)
%
% sptlrange = calcsptlrange(sptlsill,sptlrange)
%
% INPUT:
%
% sptlsill(1 by n): Partial sill
% covrange(1 by n): Covariance range
%
% OUTPUT :
%
% totcovrange(scalar): Total covariance range
%

totcovrange = 0;
totalsill = sum(sill);

for i = 1:size(sill,2)
    totcovrange = totcovrange + (sill(i)/totalsill) * covrange(i);
end
