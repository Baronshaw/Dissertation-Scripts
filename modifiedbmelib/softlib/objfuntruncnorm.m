function fval = objfuntruncnorm(param,mu_orig,sig_orig)

% objfuntruncnorm - Objective function of the parameters of the truncated
%                   normal distribution for least square fit (Dec 05,2012)
%
% fval = objfuntruncnorm(param,mu_orig,sig_orig)
%
% INPUT:
%
% param(1 by 2): [mu,sigma] parameter of the truncated normal
% mu_orig(scalr): the mean of the normal distribution
% sig_orig(scalr): the standard deviation of the normal distribution
%
% OUTPUT :
%
% fval(scalar): Value of objective function
%

mu = param(1);
sig = param(2);

expval = calctruncnormexp(mu,sig);
varval = calctruncnormvar(mu,sig);

fval = ((expval - mu_orig)/mu_orig)^2 + ...
       ((varval - sig_orig^2)/(sig_orig^2))^2;
