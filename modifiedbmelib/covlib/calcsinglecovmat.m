function covmat = calcsinglecovmat(distmat,covmdl,covparam)

% calcsinglecovmat - Calculate sigle covariance model matrix (Nov 29,2012)
%
% covmat = calcsinglecovmat(distmat,covmdl,covparam)
%
% INPUT:
%
% distmat(n by m): Distance matric
% covmdl(string): Covariance model
% covparam(1 by 2 or 3): Covariance parameter 
%
% OUTPUT :
%
% covmat(n by m): Value of objective function
%

cmdstr = ['mdlfun = @',covmdl,';'];
eval(cmdstr);

covmat = mdlfun(distmat,covparam);
