function calcstcovmom_test2()

csvdata = csvread('data04.csv');

[vcoord,vval,dummy] = avedupli(csvdata(:,[1,2,3]),csvdata(:,7));
[Z,cMS,tME,nanratio] = valstv2stg(vcoord,vval);

sdistmat = calcdistmateuclid(cMS,cMS);
tdistmat = calcdistmattemp(tME,tME);

slagbound = calclagbound(sdistmat,10);
tlagbound = calclagbound(tdistmat,10);

[slag1,scovval1,snumpair1,sweight,tlag1,tcovval1,tnumpair1,tweight] = ...
    calcstcovmom(sdistmat,tdistmat,Z,slagbound,tlagbound);

[ds1,dt1,c1,o1] = crosscovarioST2(vcoord,vcoord,...
                                  vval,vval,slagbound,0);

[ds2,dt2,c2,o2] = crosscovarioST2(vcoord,vcoord,...
                                  vval,vval,0,tlagbound);

disp('calcstcovmom_test2')
disp(['Num Pair (mom)            :',num2str(snumpair1)]);
disp(['Num Pair (crosscovarioST2):',num2str([size(vval,1),o1'])]);
disp(['Num Pair (mom)            :',num2str(tnumpair1)]);
disp(['Num Pair (crosscovarioST2):',num2str([size(vval,1),o2])]);
                              
figure;
subplot(2,1,1)
hold on;
h1 = plot(slag1,scovval1,'-ro');
h2 = plot(ds1,c1,'--b.');

subplot(2,1,2)
hold on;
plot(tlag1,tcovval1,'-ro');
plot(dt2,c2,'--b.');
