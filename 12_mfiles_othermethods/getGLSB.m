function [beta0,beta1] = getGLSB(betac,problems,chsub,zhosub,zhmsub,zhmp)
% this function will get the GLS estimator of beta0 and beta1

leng = length(zhosub);
X = [ones(leng,1) zhmsub];
xh = zhmsub;
Dxh = sqrt(bsxfun(@plus,dot(chsub,chsub,2),dot(chsub,chsub,2)')-2*(chsub*chsub'));
zhosub1 = zhosub; zhosub2 = zhosub;
Chh = zeros(leng,leng);
for i = 1:10
    for j = 1:10
        idx1 = zhosub1 >= zhmp(i) & zhosub1 < zhmp(i+1); 
        idx2 = zhosub2 >= zhmp(j) & zhosub2 < zhmp(j+1);
        A11 = betac(1); A12 = betac(2); A22 = betac(3); 
        xB1 = problems(i); xB2 = problems(j); tau2 = problems(end);
        phi0 = betac(4); phi1 = betac(5);
        Chhtemp = (A11-A12.*xB1).*(A11-A12.*xB2).*exp(-1./phi0.*Dxh) + ...
                (A12-A22.*xB1).*(A12-A22.*xB2).*exp(-1./phi1.*Dxh) + tau2.*exp(-0.5.*Dxh);
        idx1C = repmat(idx1,1,leng);
        idx2C = repmat(idx2',leng,1);
        idxC = idx1C & idx2C;
        Chhtemp(~idxC) = 0;
        Chh = Chh + Chhtemp;
    end
end
B = inv(X'*pinv(Chh)*X)*X'*pinv(Chh)*xh;

beta0 = B(1);
beta1 = B(2);

end