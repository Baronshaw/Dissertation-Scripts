function [meanGivMod,varGivMod,valMod,valObs,numpnts,...
    idxtME,idxcMS,modplots,dailyCTMg,cMSCTM,cMSCTMv, ...
    cMSObsx,cMSObsy,tMEObs] ...
    = getSoftData_LinReg(dayWHI,n,deltaT,numbins,minpnts,modplots,negval, ...
        D,Obsg,Modg,cMS,tME,cMSCTM,tMECTM,distCTMv,yrmodaCTMv,dailyCTMv)
% this will create the mean and variance for a given modeled values within
% a time-centered interval

% parameters subject to change
if nargin < 1, dayWHI = [ 2001 07 23 ]; end
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
cMSs = cMS(idxcMS,:);
Obsgs = Obsg(idxcMS,idxtME);
Modgs = Modg(idxcMS,idxtME);

Ds = D(:,idxcMS);

% getting the values for the n closest stations

[dummy sortidx] = sort(Ds,2);
whichindexes = sortidx(:,1:n);        

valObs = []; valMod = [];
cMSObs = []; tMEObs = [];
for j = 1:n
    valObs(:,:,j) = Obsgs(whichindexes(:,j),:);
    valMod(:,:,j) = Modgs(whichindexes(:,j),:);
    cMSObs(:,:,j) = cMSs(whichindexes(:,j),:); 
    tMEObs(:,j) = tMEs; 
end
[r c h] = size(valObs);    

valObs = reshape(valObs,r,c*h);
valMod = reshape(valMod,r,c*h);
cMSObsx = squeeze(cMSObs(:,1,:));
cMSObsy = squeeze(cMSObs(:,2,:));
tMEObs = reshape(tMEObs,1,c*h);

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
    cMSObsx(:,addcol+1) = NaN; cMSObsx(idxpnts,addcol+1) = cMSs(sortidx(idxpnts,addcol+1),1);
    cMSObsy(:,addcol+1) = NaN; cMSObsy(idxpnts,addcol+1) = cMSs(sortidx(idxpnts,addcol+1),2);
    tMEObs(end+1:end+addc) = tMEs; 
    
    addcol = addcol + 1;
    numpnts = sum(~isnan(valObs'))';
    idxpnts = numpnts < minpnts;
    numstations(idxpnts) = numstations(idxpnts)+1; % added in '_new'
end

% bring in modeled value and find mean and var given mod value
dayWHImod = dayWHI(:,1).*10^4 + dayWHI(:,2).*10^2 + dayWHI(:,3);
givenday = yrmodaCTMv == dayWHImod;
dailyCTMg = dailyCTMv(givenday);
if length(dailyCTMg)~=length(cMSCTM)    
    [lia lib] = ismember(cMSCTM,distCTMv(givenday,:),'rows');
    lib(lib==0)=[];
    cMSCTMv = distCTMv(givenday,:);
    dailyCTMg = dailyCTMg(lib);
else
    % [lia lib] = ismember(distCTMv(givenday,:),cMSCTM,'rows'); "corrected" 9/16/2014
    [lia lib] = ismember(cMSCTM,distCTMv(givenday,:),'rows');
    cMSCTMv = cMSCTM;
    dailyCTMg = dailyCTMg(lib);
end

meanGivMod = NaN*ones(length(dailyCTMg),1);
varGivMod = NaN*ones(length(dailyCTMg),1);

for i = 1:size(valObs,1)
    idx = ~isnan(valMod(i,:));
    fitresult = fit(valMod(i,idx)',valObs(i,idx)','poly1');
    [bounds meanGivMod(i,1)] = predint(fitresult,dailyCTMg(i),0.95,'observation','off');
    varGivMod(i,1) = ((bounds(2)-bounds(1))/4).^2;
end

end