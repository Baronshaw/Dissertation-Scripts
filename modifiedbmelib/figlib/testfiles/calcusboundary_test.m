function calcusboundary_test()

figure;
plotusstate;

[xmin,xmax,ymin,ymax] = calcusboundary();
plot(xmin,ymin,'r.');
plot(xmin,ymax,'r.');
plot(xmax,ymin,'r.');
plot(xmax,ymax,'r.');

marginratio = 1/20;
[xmin,xmax,ymin,ymax] = calcusboundary(marginratio);
plot(xmin,ymin,'b.');
plot(xmin,ymax,'b.');
plot(xmax,ymin,'b.');
plot(xmax,ymax,'b.');
