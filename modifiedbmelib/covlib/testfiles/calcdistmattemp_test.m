function calcdistmattemp_test()

tempcoord1 = [0,1,2,3,4];
tempcoord2 = [10,11,12];
distmat1 = calcdistmattemp(tempcoord1,tempcoord2);
distmat2 = calcdistmattemp(tempcoord2,tempcoord1);

disp(['tempcoord1: ',num2str(tempcoord1)]);
disp(['tempcoord2: ',num2str(tempcoord2)]);
disp('distmat(1to2): ')
disp(num2str(distmat1));
disp('distmat(2to1): ')
disp(num2str(distmat2));

tempcoord1 = zeros(1,0);
tempcoord2 = [10,11,12];
distmat1 = calcdistmattemp(tempcoord1,tempcoord2);
distmat2 = calcdistmattemp(tempcoord2,tempcoord1);

disp(['tempcoord1: ',num2str(tempcoord1)]);
disp(['tempcoord2: ',num2str(tempcoord2)]);
disp('distmat(1to2): ')
disp(num2str(distmat1));
disp('distmat(2to1): ')
disp(num2str(distmat2));
