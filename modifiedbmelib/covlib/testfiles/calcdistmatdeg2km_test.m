function calcdistmatdeg2km_test()

lonlat = [0,0;1,0;1,1;0,1];
distmat1 = calcdistmatdeg2km(lonlat,lonlat);
disp(['1 degree:',num2str(deg2km(1))]);
disp(['(lon,lat)=(0,0)-(1,0):',num2str(distmat1(1,2))]);
disp(['(lon,lat)=(0,0)-(1,1):',num2str(distmat1(1,3))]);
disp(['(lon,lat)=(0,0)-(0,1):',num2str(distmat1(1,4))]);
disp(['(lon,lat)=(0,1)-(1,1):',num2str(distmat1(4,3))]);

lonlat = [0,10;1,10;1,11;0,11];
distmat2 = calcdistmatdeg2km(lonlat,lonlat);
disp(['(lon,lat)=(0,10)-(1,10):',num2str(distmat2(1,2))]);
disp(['(lon,lat)=(0,10)-(1,11):',num2str(distmat2(1,3))]);
disp(['(lon,lat)=(0,10)-(0,11):',num2str(distmat2(1,4))]);
disp(['(lon,lat)=(0,11)-(1,11):',num2str(distmat2(4,3))]);

lonlat1 = [0,0;0,4];
lonlat2 = [0,1;0,2;0,3];
distmat3 = calcdistmatdeg2km(lonlat1,lonlat2);
disp(['(lon,lat)=(0,0)-(0,1):',num2str(distmat3(1,1))]);
disp(['(lon,lat)=(0,0)-(0,2):',num2str(distmat3(1,2))]);
disp(['(lon,lat)=(0,0)-(0,3):',num2str(distmat3(1,3))]);
disp(['(lon,lat)=(0,4)-(0,1):',num2str(distmat3(2,1))]);
disp(['(lon,lat)=(0,4)-(0,2):',num2str(distmat3(2,2))]);
disp(['(lon,lat)=(0,4)-(0,3):',num2str(distmat3(2,3))]);

latval = -90:90;
londist = zeros(size(latval)) * NaN;

for i = 1:size(latval,2)
    lonlat = [0,latval(i);1,latval(i)];
    distmat = calcdistmatdeg2km(lonlat,lonlat);
    londist(i) = distmat(1,2);
end

figure;
plot(latval,londist);
xlabel('Latitude')
ylabel('1degree(longitude) [km]');
xlim([min(latval),max(latval)]);
box on;

latval = -90:90;
latdist = zeros(size(latval)) * NaN;

for i = 1:size(latval,2)
    if i == size(latval,2)
        lonlat = [0,latval(i);0,latval(i-1)];
    else
        lonlat = [0,latval(i);0,latval(i+1)];
    end
    distmat = calcdistmatdeg2km(lonlat,lonlat);
    latdist(i) = distmat(1,2);
end

figure;
plot(latval,latdist);
xlabel('Latitude')
ylabel('1degree(latitude) [km]');
xlim([min(latval),max(latval)]);
box on;

lonval = -180:180;
londist1 = zeros(size(lonval)) * NaN;
londist2 = zeros(size(lonval)) * NaN;
londist3 = zeros(size(lonval)) * NaN;

for i = 1:size(lonval,2)
    if i == size(lonval,2)
        lonlat = [lonval(i),0;-179,0];
    else
        lonlat = [lonval(i),0;lonval(i+1),0];
    end
    distmat = calcdistmatdeg2km(lonlat,lonlat);
    londist1(i) = distmat(1,2);

    if i == size(lonval,2)
        lonlat = [lonval(i),30;lonval(i-1),30];
    else
        lonlat = [lonval(i),30;lonval(i+1),30];
    end
    distmat = calcdistmatdeg2km(lonlat,lonlat);
    londist2(i) = distmat(1,2);

    if i == size(lonval,2)
        lonlat = [lonval(i),60;lonval(i-1),60];
    else
        lonlat = [lonval(i),60;lonval(i+1),60];
    end
    distmat = calcdistmatdeg2km(lonlat,lonlat);
    londist3(i) = distmat(1,2);
end

figure;
hold on;
plot(lonval,londist1,'-b');
plot(lonval,londist2,'-r');
plot(lonval,londist3,'-g');
xlabel('Longitude')
ylabel('1degree(longitude) [km]');
xlim([min(latval),max(latval)]);
box on;
legend('0deg','30deg','60deg');
