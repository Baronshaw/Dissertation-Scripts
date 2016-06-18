function [RMSE] = myErrFunMST(BETA,VAL,MOD)
% This is for the space/time downscaler - mean trend.

betas = BETA(1:2);
rhos = BETA(3:4);
etas = BETA(5:6);

[dummy T] = size(MOD);
lambda = 1/T;

for i = 1:T

    if i == 1
        beta0t(i,1) = rhos(1)*betas(1)+etas(1);
        beta1t(i,1) = rhos(2)*betas(2)+etas(2);
    else
        beta0t(i,1) = betatfun(beta0t(i-1,1),rhos(1),etas(1));
        beta1t(i,1) = betatfun(beta1t(i-1,1),rhos(2),etas(2));
    end
    meantrend = beta0t(i,1) + beta1t(i,1)*MOD(:,i);
    idx = ~isnan(meantrend);
    errors(i,1) = lambda.*sqrt(mean((VAL(idx,i)-meantrend(idx)).^2));
    
end

RMSE = sum(errors);

end