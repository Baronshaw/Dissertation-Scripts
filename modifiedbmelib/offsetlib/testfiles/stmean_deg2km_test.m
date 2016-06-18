function stmean_deg2km_test()

lonlat = [0,0;0,1;0,2;0,3];
tME = [0,1];
Z = [1,2;2,3;3,4;4,5];

p = [1000,deg2km(1),10,10];
[ms1,mss1,mt1,mts1] = stmean_deg2km(Z,lonlat,tME,p);

p = [10,1,10,10];
[ms2,mss2,mt2,mts2] = stmean2(Z,lonlat,tME,p);

disp('(ms,mss) Use great circle distance');
disp(num2str([ms1,mss1]));

disp('(ms,mss) Use euclid distance');
disp(num2str([ms2,mss2]));
