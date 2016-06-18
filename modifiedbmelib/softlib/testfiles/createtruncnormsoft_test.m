function createtruncnormsoft_test()

x = -10:0.001:10;

mu = 1;
sigma = 2;

ynorm = normpdf(x,mu,sigma);
ytruncnorm1 = truncnormpdf(x,mu,sigma,0,Inf);
[softpdftype,nl,limi,probdens] = createtruncnormsoft(mu,sigma);

figure;
hold on;
plot(x,ynorm,'-r');
plot(x,ytruncnorm1,'-b');
plot(limi,probdens,'-g','LineWidth',2);
