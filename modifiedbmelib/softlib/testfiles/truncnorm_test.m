function truncnorm_test()

x = -10:0.001:30;
mu = 10;
sig = 2;
lb = 0;
ub = Inf;
ytrunc = truncnormpdf(x,mu,sig,lb,ub);
ynorm = normpdf(x,mu,sig);
[softpdftype,nl,limisoft,probdenssoft] = createtruncnormsoft(mu,sig);

auctrunc = trapz(x,ytrunc);
aucnorm = trapz(x,ynorm);
aucsoft = trapz(limisoft,probdenssoft);
disp(['AUC - truncated normal: ',num2str(auctrunc)]);
disp(['AUC - normal distribution: ',num2str(aucnorm)]);
disp(['AUC - soft: ',num2str(aucsoft)]);

cumprob = zeros(1,size(x,2)-1);
for i = 1:length(x)-1
    cumprob(i) = trapz(x(1:i+1),ytrunc(1:i+1));
end

% plist = [0.0001,0.001,0.01,0.05,0.3,0.4,0.5,0.6,0.7,0.95,0.99,0.999,0.9999];
% [px,py] = truncnormprob(mu,sig,lb,ub,plist);

[softpdftype,nl,liminorm,probdensnorm] = probaGaussian(mu,sig^2);

figure;
hold on;
h1 = plot(x,ytrunc,'-r','LineWidth',2);
h2 = plot(x,ynorm,'-b','LineWidth',2);
h3 = plot(x(1:end-1),cumprob,'--r','LineWidth',2);
% h4 = plot(px,plist,'go');
% h5 = plot(px,py,'--co','LineWidth',2);
h6 = plot(limisoft,probdenssoft,'--yo','LineWidth',2);
h7 = plot(liminorm,probdensnorm,'--mo','LineWidth',2);
legend([h1,h2,h3,h6,h7],'Truncated Normal PDF','Normal PDF',...
       'Truncated Normal CDF','Truncated Normal Soft PDF',...
       'Normal Soft PDF',2);
xlabel('X');
ylabel('Y');
title(['Truncated Normal (',num2str(mu),'-',num2str(sig),'-',...
       num2str(lb),'-',num2str(ub),')'],'FontSize',18);
ylim([0,1]);
