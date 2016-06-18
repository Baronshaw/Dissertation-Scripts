function [mu_new,sig_new] = lsparamfittruncnorm(mu_orig,sig_orig)

% lsparamfittruncnorm - Calculate the mu and sigma parameter for the 
%                       truncated normal distribution with
%                       the mean = mu and the variance = sigma2
%                       (Dec 05,2012)
%
% param = lsparamfittruncnorm(mu_orig,sig_orig);
%
% INPUTS :
%  
% mu_orig(scalr): the mean of the normal distribution
% sig_orig(scalr): the standard deviation of the normal distribution
%
% OUTPUT :
%
% mu_new(scalr): the mu parameter for the truncated normal distribution
% sig_new(scalr): the sigma parameter for the truncated normal distribution
%

initval = [mu_orig,sig_orig];

lb = [realmin,realmin];
ub = [Inf,Inf];

fmops = optimset('Algorithm','interior-point','DerivativeCheck','off',...
                 'Diagnostic','off','Display','off','FunValCheck','on',...
                 'MaxFunEvals',10000,'MaxIter',10000,...
                 'TolX',1e-16,'TolFun',1e-16);

[param,fval,exitflag,output] = ...
    fmincon(@(x) objfuntruncnorm(x,mu_orig,sig_orig),initval,[],[],[],[],...
            lb,ub,[],fmops);

if exitflag < 1
    disp('--- msg ---')
    disp(output.message);
end

mu_new = param(1);
sig_new = param(2);
