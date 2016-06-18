function [sill,covrange] = separatecovparam(covmodel,covparam)

% separatecovparam - Separate covariance parameter into sill and range
%                    (Dec 05,2012)
%
% [sill,range] = separatecovparam(covmodel,covparam)
%
% INPUT:
%
% covmdl(string/cell): Covariance model
% covparam(1 by 1 or 2/cell): Covariance parameter 
%
% OUTPUT :
%
% sill(1 by n): Sill (n = # of cell)
% covrange(1 by n): Covariance range
%

if iscell(covmodel);
    nummodel = size(covmodel,2);
else
    nummodel = 1;
    covmodel = {covmodel};
    covparam = {covparam};
end

sill = zeros(size(covmodel)) * NaN;
covrange = zeros(size(covmodel)) * NaN;

for i = 1:nummodel
    sill(i) = covparam{i}(1);
    if strcmp(covmodel{i},'nuggetC')
        covrange(i) = 0;
    else
        covrange(i) = covparam{i}(2);
    end
end
