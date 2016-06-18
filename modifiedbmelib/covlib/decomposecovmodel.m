function [sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam)

% decomposecovmodel - Decompose separable covariance model and parameter
%                     into spatial and temporal component (Nov 29,2012)
%
% [sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
%      decomposecovmodel(covmodel,covparam)
%
% INPUT:
%
% covmdl(string/cell): Covariance model
% covparam(1 by 2 or 3/cell): Covariance parameter 
%
% OUTPUT :
%
% sptlcovmodel(cell): Spatial covariance model
% sptlcovparam(cell): Spatial covariance parameter
% tempcovmodel(cell): Temporal covariance model
% tempcovparam(cell): Temporal covariance parameter
%

if iscell(covmodel);
    nummodel = size(covmodel,2);
else
    nummodel = 1;
    covmodel = {covmodel};
    covparam = {covparam};
end

for i = 1:nummodel
    idxslash = findstr(covmodel{i},'/');
    if isempty(idxslash);
        error('Covariance model error');
    end
    sptlcovmodel{i} = covmodel{i}(1:idxslash-1);
    tempcovmodel{i} = covmodel{i}(idxslash+1:end);

    if strcmp(sptlcovmodel{i},'nuggetC')
        sptlcovparam{i} = covparam{i}(1);
        if strcmp(tempcovmodel{i},'nuggetC')
            tempcovparam{i} = covparam{i}(1);
        else
            tempcovparam{i} = [covparam{i}(1),covparam{i}(2)];
        end
    else
        sptlcovparam{i} = [covparam{i}(1),covparam{i}(2)];
        if strcmp(tempcovmodel{i},'nuggetC')
            tempcovparam{i} = covparam{i}(1);
        else
            tempcovparam{i} = [covparam{i}(1),covparam{i}(3)];
        end
    end
end
