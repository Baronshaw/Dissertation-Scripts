function y = truncnormpdf(x,mu,sig,lb,ub)

% truncnormpdf - Truncated normal PDF (Dec 05,2012)
%
% y = truncnormpdf(x,mu,sig,lb,ub)
%
% INPUT:
%
% x(n by 1): values
% mu(scalar): mean of the corresponding normal distribution
% sig(scalar): standard deviation of the corresponding normal distribution
% lb(scalar): lower bound
% ub(scalar): upper bound
%
% OUTPUT:
%
% y(n by 1): Truncated normal pdf
%

y = zeros(size(x));
for i = 1:length(x)
    if x(i) >= lb && x(i) <= ub
        y(i) = normpdf(x(i),mu,sig)/...
               (normcdf(ub,mu,sig) - normcdf(lb,mu,sig));
    else
        y(i) = 0;
    end
end
