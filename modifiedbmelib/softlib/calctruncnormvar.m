function varval = calctruncnormvar(mu,sig)

% calctruncnormvar - Calculate the variance of the truncated normal
%                    distribution based on the mean and standard deviation
%                    of the corresponding normal distirubution
%                    (Dec 05,2012)
%
% expval = calctruncnormvar(mu,sig);
%
% INPUTS :
%  
% mu(scalr): the mean of the normal distribution
% sig(scalr): the standard deviation of the normal distribution
%
% OUTPUT :
%
% expval(scalar): the expected value of the truncated normal distribution
%                 (one side, truncated at 0)
%

varval = sig^2 * (1 - ((normpdf((0-mu)/sig))/...
                       (1 - normcdf((0-mu)/sig))) * ...
                      (((normpdf((0-mu)/sig))/...
                        (1 - normcdf((0-mu)/sig))) - ((0-mu)/sig)));
