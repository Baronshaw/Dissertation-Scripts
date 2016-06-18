function [px,py] = truncnormprob(mu,sig,lb,ub,plist)

% truncnormprob - Calculate the probability density based on the value of 
%                 the cumulative probability values
%
% y = truncnorm(mu,sig,lb,ub,plist)
%
% INPUT:
%
% mu(scalar): Mean
% sig(scalar): Standard Deviation
% lb(scalar): Lower Bound
% ub(scalar): Upper Bound
% plist(n by 1): Cumulative probability value
%
% OUTPUT:
%
% px(n by 1): Value at Cumulative probability value
% py(n by 1): Probabilty density
%

xmax = mu + 4 * sig;
x = 0:xmax/10000:xmax;
y = truncnormpdf(x,mu,sig,lb,ub);

cumprob = zeros(1,size(x,2)-1);
for i = 1:length(x)-1
    cumprob(i) = trapz(x(1:i+1),y(1:i+1));
end

px = zeros(size(plist)) * NaN;
py = zeros(size(plist)) * NaN;
for i = 1:length(plist)
    idx = find(cumprob>plist(i),1,'first');
    px(i) = x(idx);
    py(i) = y(idx);
end
