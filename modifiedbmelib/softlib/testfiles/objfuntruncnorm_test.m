function objfuntruncnorm_test()

mu_orig = 1;
sig_orig = 2;

expval = calctruncnormexp(mu_orig,sig_orig);
varval = calctruncnormvar(mu_orig,sig_orig);

mu = 0:0.1:10;
sig = 0:0.1:10;


fval = ones(size(mu,2) * size(sig,2),1) * NaN;
params = [];

midx = 1;
for i = 1:size(mu,2)
    for j = 1:size(sig,2)
        param = [mu(i),sig(j)];
        fval(midx) = objfuntruncnorm(param,expval,sqrt(varval));
        params = [params;param];
        midx = midx + 1;
    end
end

disp(['index of min(f): ',num2str(transpose(find(fval==min(fval))))]);
disp(['min(f)         : ',num2str(min(fval))]);

disp(['original mu      : ',num2str(mu_orig)]);
disp(['original sigma   : ',num2str(sig_orig)]);
disp(['calculated mu    : ',...
      num2str(transpose(params(fval == min(fval),1)))]);
disp(['calculated sigma : ',...
      num2str(transpose(params(fval == min(fval),2)))]);
