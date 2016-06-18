function truncnormprob_test()

mu = 1;
sig = 2;
lb = 0;
ub = Inf;
plist = 0.1:0.1:0.9;

[px,py] = truncnormprob(mu,sig,lb,ub,plist);

x = 0:0.01:10;
y = truncnormpdf(x,mu,sig,lb,ub);

figure;
hold on;
plot(x,y,'-r');
plot(px,py,'-go');
