function calclagbound_test()

load('testdata_2306_06.mat');
lonlat = [TestDataLon,TestDataLat];

distmat = calcdistmatdeg2km(lonlat,lonlat);

lagparam = 20;
[lagbound,numpairs] = calclagbound(distmat,lagparam);

if ~isequal(size(lagbound,2),lagparam)
    error('Number of lags is wrong');
end

idx1 = find(lagbound(1) < distmat & distmat<=lagbound(2));
idx2 = find(lagbound(2) < distmat & distmat<=lagbound(3));
idx3 = find(lagbound(3) < distmat & distmat<=lagbound(4));

disp(['numpairs(2) - calclagbound: ',num2str(numpairs(2))]);
disp(['numpairs(2) - distmat     : ',num2str(size(idx1,1)/2)]);
disp(['numpairs(3) - calclagbound: ',num2str(numpairs(3))]);
disp(['numpairs(3) - distmat     : ',num2str(size(idx2,1)/2)]);
disp(['numpairs(4) - calclagbound: ',num2str(numpairs(4))]);
disp(['numpairs(4) - distmat     : ',num2str(size(idx3,1)/2)]);

disp(['(Max Sep/2) - calclagbound: ',num2str(lagbound(end))]);
disp(['(Max Sep/2) - distmat     : ',num2str(max(max(distmat))/2)]);

disp('lag,numpairs')
disp(num2str([lagbound',numpairs']));
