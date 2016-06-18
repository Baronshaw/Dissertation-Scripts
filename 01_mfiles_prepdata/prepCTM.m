function [whiproj nad83 distCTM coordCTM yrmodaCTM dailyCTMg distCTMv ...
    yrmodaCTMv dailyCTMv] = prepCTM(filename,years)
% this function will prepare the CTM data for the soft data analysis by: 1)
% converting the CTM data to NAD83 and convert to the projection presented
% in the CTM data, 2) average hourly data across time, 3) have data in
% format: lat lon year mo da projlat projlon

if nargin < 1, filename = 'J3a_b313.36km.std.2001year_PM.conc'; end
if nargin < 2, years = 2001; end

% loading the CMAQ data
cd ../datafiles/Modeled_PM2p5
globalinfo = ncinfo(filename);
cd ../../01_mfiles_prepdata
for i = 1:length(globalinfo.Attributes)
    name{i} = globalinfo.Attributes(i).Name;
    value{i} = globalinfo.Attributes(i).Value;
end

% x and y coordinates of the projection given a centriod values
[vertdist,hortdist,vertgrid,hortgrid] = getVHdistgrid(name,value);

% defining the projection
whiproj.lat = [value{strcmp(name,'P_ALP')} value{strcmp(name,'P_BET')}];
whiproj.ellips = 'nad83';
whiproj.ORell = [value{strcmp(name,'XCENT')} value{strcmp(name,'YCENT')}];
whiproj.ORproj = [0 0];
if value{strcmp(name,'GDTYP')}==2, whiproj.type = 'lambertcc2'; end

% convert data from Lambert CC to NAD83
cd ../09_mfiles_projections
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj');
coordCTM = lambertcc2ell([hortgrid(:) vertgrid(:)],'whiproj'); 
cd ../01_mfiles_prepdata

% converting data from NAD83 back to the CTM projection of 2001
filename2001 = 'J3a_b313.36km.std.2001year_PM.conc';
cd ../datafiles/Modeled_PM2p5
globalinfo = ncinfo(filename2001);
cd ../../01_mfiles_prepdata
for i = 1:length(globalinfo.Attributes)
    name{i} = globalinfo.Attributes(i).Name;
    value{i} = globalinfo.Attributes(i).Value;
end
whiproj2001.lat = [value{strcmp(name,'P_ALP')} value{strcmp(name,'P_BET')}];
whiproj2001.ellips = 'nad83';
whiproj2001.ORell = [value{strcmp(name,'XCENT')} value{strcmp(name,'YCENT')}];
whiproj2001.ORproj = [0 0];
if value{strcmp(name,'GDTYP')}==2, whiproj2001.type = 'lambertcc2'; end

% convert data from NAD83 to Lambert CC 
cd ../09_mfiles_projections
load Projections.mat
save('Projections.mat','agk28','agk31','agk34','bev','bmn_gk','france_1',...
    'france_2','france_2_et','france_3','france_4','gk','lambert93',...
    'utm','whiproj','whiproj2001');
load Ellipsoids.mat
% from Wikipedia:
nad83.a = 6378137; nad83.b = 6356752.3141; nad83.f = 1/298.257222101; 
save('Ellipsoids.mat','airy1830','bessel1841','besseldhdn','clarke1880',...
    'grs80','hayford','wgs84','nad83');
distCTM = ell2lambertcc([coordCTM(:,1) coordCTM(:,2)],'whiproj2001'); 
cd ../01_mfiles_prepdata

% averaging all the CTM data from hourly to daily
cd ../datafiles/Modeled_PM2p5
if years == 2005
    % I must use the times from a different years, because there are
    % problems with 2005
    tflag = ncread('F1.2007yearly.PM2p5.D36.ncf','TFLAG');
else
    tflag = ncread(filename,'TFLAG');
end
if years == 2005,
    CTM = ncread(filename,'PM25_TOT');
elseif years == 2007
    CTM = ncread(filename,'PM2.5');
else
    CTM = ncread(filename,'PM25');
end
cd ../../01_mfiles_prepdata
CTM = squeeze(CTM);
CTM = double(CTM);
daysCTM = str2num(num2str(tflag(1,1,:) - floor(tflag(1,1,:)./1000).*1000));
yearCTM = str2num(num2str(floor(tflag(1,1,:)./1000)));
hoursCTM = str2num(num2str(floor(tflag(2,1,:)./10000)));

% average all values to daily values
[alldays uniidx] = unique(daysCTM);
yrmodaCTM = datevec(datenum(years,0,alldays));
[r c h] = size(CTM);
dailyCTMg = NaN*ones(c,r,length(alldays)); % r and c reverse on purpose

 % added 9/14/2014, correcting for time zone
fid = fopen('../for_timezone/CMAQ_36k_US_timezone_on_grid.txt');
Zone = textscan(fid,'%f%f%f%s%f', 'delimiter', '\t','headerLines',1);
fclose(fid);
Col = Zone{1};
Row = Zone{2};
ColRow = Zone{3};
TimeZn = Zone{4};
Zn = Zone{5};
Zn_shift = datenum(yearCTM',0,daysCTM',hoursCTM',0,0); 
Zn_human = unique(datevec(datenum(yearCTM',0,daysCTM')),'rows');
Zn_humans = Zn_human(:,1)*10^4 + Zn_human(:,2)*10^2 + Zn_human(:,3);
Zn_shifted = arrayfun(@(x) datevec(Zn_shift+(x/24)), Zn, 'UniformOutput', false);
Zn_shifteds = cellfun(@(x) x(:,1)*10^4 + x(:,2)*10^2 + x(:,3), Zn_shifted, 'UniformOutput', false);

for i = 1:length(Zn_shifteds)
    disp(i);
    temp = Zn_shifteds{i};
    idx = arrayfun(@(x) x==temp, Zn_humans, 'UniformOutput', false);
    dailyCTMg(Row(i),Col(i),:) = cell2mat( cellfun(@(x) mean(CTM(Col(i),Row(i),x)), idx, 'UniformOutput', false) );
end

[rC cC hC] = size(dailyCTMg);
dailyCTMg = reshape(dailyCTMg,rC*cC,hC);

% converting from space/time grid to space/time vector
dailyCTMv = dailyCTMg(:);
yrmodaCTMv = repmat(datenum(years,0,alldays),length(distCTM),1);
yrmodaCTMv = datevec(yrmodaCTMv(:));
distCTMv = repmat(distCTM,length(alldays),1);

end