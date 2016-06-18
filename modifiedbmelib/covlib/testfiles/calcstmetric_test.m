function calcstmetric_test()

sill = [0.15,0.85];
sptlrange = [3,120];
temprange = [2,400];

scovrange = 0.15 * 3 + 0.85 * 120;
tcovrange = 0.15 * 2 + 0.85 * 400;
stmetric1 = scovrange/tcovrange;

stmetric2 = calcstmetric(sill,sptlrange,temprange);

disp(['Manual Calculation: ',num2str(stmetric1)]);
disp(['Output of calcstmetric: ',num2str(stmetric2)]);
