function fval = objfunsttwoexp(covparam,slag,scovval,sweight,...
                                        tlag,tcovval,tweight)

% objfunsttwoexp - Objective function of space/time two exponential model 
%                  for least square fit (Nov 28,2012)
%
% fval = objfunsttwoexp(covparam,lag,covval,weight)
%
% INPUT:
%
% covparam(1 by 6): Covariance parameter
% slag(1 by m): Spatial lag
% scovval(1 by m): Experimental covariance
% sweigh(1 by m): Weight for general least square covarinace model fit
% tlag(1 by m): Temporal lag
% tcovval(1 by m): Experimental covariance
% tweigh(1 by m): Weight for general least square covarinace model fit
%
% OUTPUT :
%
% fval(scalar): Value of objective function
%

scovmdlval = exponentialC(slag,[covparam(1),covparam(2)]) + ...
             exponentialC(slag,[covparam(4),covparam(5)]);
tcovmdlval = exponentialC(tlag,[covparam(1),covparam(3)]) + ...
             exponentialC(tlag,[covparam(4),covparam(6)]);

fval = sweight*((scovmdlval - scovval).^2)' + ...
       tweight*((tcovmdlval - tcovval).^2)';
