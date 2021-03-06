function maskcontigusa(maskfillcolor,masklinecolor,stboundarycolor)

marginratio = 1/20;
[xmin,xmax,ymin,ymax] = calcusboundary(marginratio);

load('USAcontiguous.mat');
x(1576:1602) = [];
y(1576:1602) = [];
pVal{1} = [x,y];
plotmask(pVal,[xmin,xmax,ymin,ymax],maskfillcolor,masklinecolor);
plotusstate(stboundarycolor);
