function lsparamfittruncnorm_test()

mu = 3;
sig = 2;

expval = calctruncnormexp(mu,sig);
varval = calctruncnormvar(mu,sig);

[mu_new,sig_new] = lsparamfittruncnorm(expval,sqrt(varval));

disp(['the mean of truncnorm with mu = ',num2str(mu),' and sig = ',...
      num2str(sig),' is ',num2str(expval)]);
disp(['the variance of truncnorm with mu = ',num2str(mu),' and sig = ',...
      num2str(sig),' is ',num2str(varval)]);

disp(['the mu parameter of truncnorm with the mean = ',num2str(expval),...
      ' and variance = ',num2str(varval),' is ',num2str(mu_new)]);
disp(['the sig parameter of truncnorm with the mean = ',num2str(expval),...
      ' and variance = ',num2str(varval),' is ',num2str(sig_new)]);
