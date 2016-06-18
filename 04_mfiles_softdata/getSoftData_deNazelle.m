function [meanGivMod,varGivMod,valMod,valObs,perctile_data,numpnts,...
    idxtME,idxcMS,modplots,dailyCTMg,mean_Mod,mean_Obs,var_Obs,cMSCTM] ...
    = getSoftData_deNazelle(dayWHI,n,deltaT,numbins,minpnts,modplots,negval, ...
        D,Obsg,Modg,cMS,tME,cMSCTM,tMECTM,distCTMv,yrmodaCTMv,dailyCTMv)
% this will create the mean and variance for a given modeled values within
% a time-centered interval

% parameters subject to change
if nargin < 1, dayWHI = [ 2001 06 15 ]; end
if nargin < 2, n = 3; end % number of closest stations
if nargin < 3, deltaT = 365; end % number of days in the interval
if nargin < 4, numbins = 10; end % number of bins for each plot
if nargin < 5, minpnts = 150; end % min number of point to each plot
if nargin < 6, modplots = 0:5:50; end % these are the modeled values you will see
if nargin < 7, negval = 0; end % 0 = there are no negative predicted values

% run program for certain days on interest
[rWHI cWHI] = size(dayWHI);
minT = datenum(dayWHI)-deltaT/2;
maxT = datenum(dayWHI)+deltaT/2;
yr = floor(tME'./10^4);
mo = floor(tME'./10^2)-yr*10^2;
da = tME'-yr*10^4-mo*10^2;
numtME = datenum([yr mo da]);

% first getting rid of all the days outside the window of interest
idxtME = numtME >= minT & numtME <= maxT;
tMEs = tME(idxtME);
tempObsg = Obsg(:,idxtME);
[r c] = size(tempObsg);
nancount = sum(isnan(tempObsg'))';
idxcMS = nancount < c;
cMSs = cMS(idxcMS);
Obsgs = Obsg(idxcMS,idxtME);
Modgs = Modg(idxcMS,idxtME);

Ds = D(:,idxcMS);

% getting the values for the n closest stations

[dummy sortidx] = sort(Ds,2);
whichindexes = sortidx(:,1:n);        

valObs = []; valMod = [];
for j = 1:n
    valObs(:,:,j) = Obsgs(whichindexes(:,j),:);
    valMod(:,:,j) = Modgs(whichindexes(:,j),:);
end
[r c h] = size(valObs);    

valObs = reshape(valObs,r,c*h);
valMod = reshape(valMod,r,c*h);

% making sure each site has enough points
numpnts = sum(~isnan(valObs'))';
idxpnts = numpnts < minpnts;    
[addr addc] = size(Obsgs);
addcol = n;
numstations = n*ones(length(valObs),1);

while sum(idxpnts) > 0 
    valObs(:,end+1:end+addc) = NaN;
    valMod(:,end+1:end+addc) = NaN;
    valObs(idxpnts,end-addc+1:end) = Obsgs(sortidx(idxpnts,addcol+1),:);
    valMod(idxpnts,end-addc+1:end) = Modgs(sortidx(idxpnts,addcol+1),:);
    %valObs(~idxpnts,end-addc+1:end) = NaN;
    %valMod(~idxpnts,end-addc+1:end) = NaN;
    addcol = addcol + 1;
    numpnts = sum(~isnan(valObs'))';
    idxpnts = numpnts < minpnts;
    numstations(idxpnts) = numstations(idxpnts)+1; % added in '_new'
end

% percentiles for each grid 
[rnew cnew] = size(valMod);
perctile_data = cell(rnew,1);
tempmax = max(valMod');
for j = 1:rnew
    if numbins > tempmax(j)
        perctile_data{j} = [0 tempmax(j)];
    else
        perctile_data{j} = 0:numbins:tempmax(j);
    end
end

% find the mean and variance of each bin
mean_Mod = cell(rnew,1);
var_Mod = cell(rnew,1);
mean_Obs = cell(rnew,1);
var_Obs = cell(rnew,1);
for j = 1:rnew
    for k = 1:size(perctile_data{j},2)-1
        idx_highlow = valMod(j,:) >= perctile_data{j}(k) & valMod(j,:) < perctile_data{j}(k+1);
        mean_Mod{j}(k) = nanmean(valMod(j,idx_highlow));
        mean_Obs{j}(k) = nanmean(valObs(j,idx_highlow));
        var_Obs{j}(k) = nanvar(valObs(j,idx_highlow),0);
    end
end 

% bring in modeled value and find mean and var given mod value
dayWHImod = dayWHI(:,1).*10^4 + dayWHI(:,2).*10^2 + dayWHI(:,3);
givenday = yrmodaCTMv == dayWHImod;
dailyCTMg = dailyCTMv(givenday);
[lia lib] = ismember(distCTMv(givenday,:),cMSCTM,'rows');
dailyCTMg = dailyCTMg(lib);

meanGivMod = NaN*ones(length(dailyCTMg),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMg),1+length(modplots));
for j = 1:length(meanGivMod)
    if size(mean_Mod{j},2)>1 & sum(isnan(mean_Mod{j}))==0
        meanGivMod(j,:) = interp1(mean_Mod{j},mean_Obs{j},...
            [dailyCTMg(j) modplots],'linear','extrap');
        varGivMod(j,:) = interp1(mean_Mod{j},var_Obs{j},...
            [dailyCTMg(j) modplots],'linear','extrap');
    else
        meanGivMod(j,:) = NaN;
        varGivMod(j,:) = NaN;
    end
end

% if values are negative, make them zero, but show S-curves
if negval == 0
    meanGivMod(meanGivMod<0) = 0;
end

end