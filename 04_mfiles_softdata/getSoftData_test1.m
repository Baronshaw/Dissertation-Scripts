function [meanGivMod,varGivMod,valMod,valObs,perctile_data,numpnts,...
    idxtME,idxcMS,modplots,dailyCTMg,mean_Mod,mean_Obs,var_Obs,cMSCTM] ...
    = getSoftData_test1(dayWHI,n,deltaT,numbins,minpnts,modplots,negval, ...
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
    addcol = addcol + 1;
    numpnts = sum(~isnan(valObs'))';
    idxpnts = numpnts < minpnts;
    numstations(idxpnts) = numstations(idxpnts)+1; 
end

% perciles for each grid 
[rnew cnew] = size(valMod);
perctile_data = prctile(valMod',linspace(0,100,numbins+1))';

% find the mean and variance of each bin
mean_Mod = NaN*ones(rnew,numbins);
var_Mod = NaN*ones(rnew,numbins);
mean_Obs = NaN*ones(rnew,numbins);
var_Obs = NaN*ones(rnew,numbins);
for j = 1:numbins
    perctile_low = repmat(perctile_data(:,j),[1 cnew]);
    perctile_high = repmat(perctile_data(:,j+1),[1 cnew]);
    valModtemp = valMod;
    valObstemp = valObs;

    idx_highlow = valModtemp >= perctile_low & valModtemp < perctile_high;
    tempmod = valModtemp(:);
    tempobs = valObstemp(:);
    temphighlow = idx_highlow(:);
    tempmod( temphighlow == 0 ) = NaN;
    tempobs( temphighlow == 0 ) = NaN;
    valModtemp = reshape(tempmod,rnew,cnew);
    valObstemp = reshape(tempobs,rnew,cnew);

    mean_Mod(:,j) = nanmean(valModtemp,2);
    mean_Obs(:,j) = nanmean(valObstemp,2);
    var_Obs(:,j) = nanvar(valObstemp,0,2);
end 

% bring in modeled value and find mean and var given mod value
dayWHImod = dayWHI(:,1).*10^4 + dayWHI(:,2).*10^2 + dayWHI(:,3);
givenday = yrmodaCTMv == dayWHImod;
dailyCTMg = dailyCTMv(givenday);
if length(dailyCTMg)~=length(cMSCTM)
    [lia lib] = ismember(cMSCTM,distCTMv(givenday,:),'rows');
    lib(lib==0)=[];
    dailyCTMg = cMSCTM(lib);
else
    [lia lib] = ismember(distCTMv(givenday,:),cMSCTM,'rows');
    dailyCTMg = dailyCTMg(lib);
end

meanGivMod = NaN*ones(length(dailyCTMg),1+length(modplots));
varGivMod = NaN*ones(length(dailyCTMg),1+length(modplots));
for j = 1:length(meanGivMod)
    idx = ~isnan(mean_Mod(j,:)); % sometimes there's not data between given percentiles
    meanGivMod(j,:) = interp1(mean_Mod(j,idx),mean_Obs(j,idx),...
        [dailyCTMg(j) modplots],'linear','extrap');
    varGivMod(j,:) = interp1(mean_Mod(j,idx),var_Obs(j,idx),...
        [dailyCTMg(j) modplots],'linear','extrap');
end

% if values are negative, make them zero, but show S-curves
if negval == 0
    meanGivMod(meanGivMod<0) = 0;
end

end