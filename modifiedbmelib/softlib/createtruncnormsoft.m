function [softpdftype,nl,limi,probdens] = createtruncnormsoft(mu,sig)

% createtruncnormsoft - Create soft data based on truncated normal
%                       distribution (Dec 05,2012)
%
% [softpdftype,nl,limi,probdens] = createtruncnormsoft(mu,sig)
%
% INPUT:
%
% mu(scalr): the mean of the normal distribution
% sig(scalr): the standard deviation of the normal distribution
%
% OUTPUT:
%
% softpdftype(scalar=2): Indicates the type of soft pdf representing 
%                        the probabilitic soft data.  
% nl(n by 1): Vector of number of limits  nl (see probasyntax)
% limi(n by l3): matrix representing the limits  values (see probasyntax)
% probdens(n by 13): matrix representing the probability densities
%                    (see probasyntax)
%

limi = zeros(size(mu,1),13);
probdens = zeros(size(mu,1),13);
nl = ones(size(mu)) * 13;
softpdftype = 2;

for i = 1:size(mu,1)
    
    if mu - sig < 0
        plist = [0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.95,...
                 0.99,0.999,0.9999];
    else
        plist = [0.0001,0.001,0.01,0.05,0.3,0.4,0.5,0.6,0.7,0.95,...
                 0.99,0.999,0.9999];
    end
    [limi(i,:),probdens(i,:)] = truncnormprob(mu(i),sig(i),0,Inf,plist);
end
