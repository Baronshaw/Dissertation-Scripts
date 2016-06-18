function [] = plotdemoCAMP()
% this function will create a figure demonstrating a example of the
% on-the-fly RAMP method

pah = 1;
i = 1;

% load data
load('matfiles/pah_data.mat');

% defining variables
idx = ~isnan(val(:,pah+4));
zPAH = log(val(idx,pah+4));
zPMsub = log(val(idx,4));
zPM = log(val(:,4));
cPM = [ProjectX ProjectY Time];
cPAH = [ProjectX(idx) ProjectY(idx) Time(idx)];
buff = std(zPMsub);
idx = zPMsub >= zPM(i)-buff & zPMsub <= zPM(i)+buff;
mPAHs = mean(zPAH(idx));
vPAHs = var(zPAH(idx));

figure; hold on;

plot(zPM(i),sqrt(vPAHs),'g+','MarkerSize',10,'LineWidth',3);
plot(zPMsub,zPAH,'b.','MarkerSize',20);

a = get(gca,'YLim');
b = get(gca,'XLim');

plot([zPM(i)-buff zPM(i)-buff],[a(1) a(2)],'c--','LineWidth',5);
plot([zPM(i)+buff zPM(i)+buff],[a(1) a(2)],'c--','LineWidth',5);

xlabel(sprintf('log(PM)'));
ylabel(sprintf('log(%s)',valdispname{pah}));

plot([zPM(i) zPM(i)],[a(1) mPAHs],'k--','LineWidth',3);
plot([b(1) zPM(i)],[mPAHs mPAHs],'k--','LineWidth',3);
plot(zPM(i),mPAHs,'rx','MarkerSize',20,'LineWidth',5);

set(gcf,'Position',[0 0 800 500]); 
print(gcf,'-painters','-dpng','-r600','figures/plotdemoCAMP.png');

end