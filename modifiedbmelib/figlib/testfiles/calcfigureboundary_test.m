function calcfigureboundary_test()

sptlcoord = [10,120;20,180];
marginratio = 1/20;
[xmin,xmax,ymin,ymax] = calcfigureboundary(sptlcoord,marginratio);

delval = (sptlcoord(2,2) - sptlcoord(1,2)) * marginratio;

disp(['input xmin:',num2str(sptlcoord(1,1)-delval),...
      ' xmax:',num2str(sptlcoord(2,1)+delval),...
      ' ymin:',num2str(sptlcoord(1,2)-delval),...
      ' ymax:',num2str(sptlcoord(2,2)+delval)]);
disp(['output xmin:',num2str(xmin),' xmax:',num2str(xmax),...
      ' ymin:',num2str(ymin),' ymax:',num2str(ymax)]);
