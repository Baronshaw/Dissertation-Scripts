function covmat = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam)
    
% calcstcovmat - Calculate space/time covariance matrix based on 
%                distance matrix
%
% covmat = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam)
%
% INPUT:
%
% sptldistmat(n by m): Spatial distance matrix
% tempdistmat(n by m): Temporal distance matrix
% covmdl(string/cell): Covariance model
% covparam(1 by 2 or 3/cell): Covariance parameter 
%
% OUTPUT :
%
% covmat(n by m): Covariance matrix
%

[sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam);

covmat = zeros(size(sptldistmat));

for i = 1:size(sptlcovmodel,2)
    sill = sptlcovparam{i}(1);
    if strcmp(sptlcovmodel{i},'nuggetC')
        scovparam = 1;
    else
        scovparam = [1,sptlcovparam{i}(2)];
    end
    scovmat = calcsinglecovmat(sptldistmat,sptlcovmodel{i},...
                               scovparam);
    if strcmp(tempcovmodel{i},'nuggetC')
        tcovparam = 1;
    else
        tcovparam = [1,tempcovparam{i}(2)];
    end
    tcovmat = calcsinglecovmat(tempdistmat,tempcovmodel{i},...
                               tcovparam);
    
    covmat = covmat + sill.*scovmat.*tcovmat;
end
