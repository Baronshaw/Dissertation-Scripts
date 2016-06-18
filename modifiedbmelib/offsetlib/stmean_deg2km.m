function [ms,mss,mt,mts] = stmean_deg2km(Z,lonlat,tME,p)

% stmean_deg2km - Estimates space/time mean values from measurements
%                 (Nov 28,2012)
%
% [ms,mss,mt,mts] = stmean_deg2km(Z,lonlat,tME,p);
%
% INPUTS :
%  
% Z(nMS by nME): matrix of measurements for Z at the nMS monitoring
%                sites and nME measuring events. Z may have NaN values.
% lonlat(nMS by 2): Spatial coordinates
%                 1st column is longitude (x)
%                 2nd column is latitude (y)
% tME(1 by nME): vector with the time of the measuring events
% p(1 by 5): parameters to smooth the spatial and temporal average
%            p(1)=dNeib distance (radius) of spatial neighborhood
%            p(2)=ar    spatial range of exponential smoothing function
%            p(3)=tNeib time (radius) of temporal neighborhood
%            p(4)=at    temporal range of exponential smoothing function
%            p(5)=tloop is an optional input used for temporal smoothing
%            When tloop>0, the measuring events are looped in a
%            cycle of duration tloop.
%
% OUTPUT :
%
% ms(nMS by 1): vector of spatial average
% mss(nMS by 1): vector of smoothed spatial average
% mt(1 by nME): vector of temporal average
% mts(1 by nME): vector of smoothed temporal average
%

%
% Check the input arguments
%
nMS = size(Z,1);
nME = size(Z,2);

if size(lonlat,1) ~= nMS | size(lonlat,2) ~= 2
    error('cMS must be a nMS by 2 matrix'); 
end

if size(tME,1) ~= 1 | size(tME,2) ~= nME
    error('tME must be a 1 by nME vector'); 
end

%
% set parameters and variables
%
dNeib = p(1); % Radius use to select local neighborhood to average Z
ar = p(2);    % spatial range of exponential smoothing function
tNeib = p(3); % Radius use to select local neighborhood to average PM10
at = p(4);    % temporal range of exponential smoothing function

if length(p) > 4
    tloop = p(5);
else
    tloop = 0;
end

xMS = lonlat(:,1);
yMS = lonlat(:,2);

%
%  Calculate the spatial averages ms
%

for iMS = 1:nMS
    ms(iMS) = mean(Z(iMS,~isnan(Z(iMS,:))));
end

%
%  smooth the spatial average with an exponential filter to get mss
%
for iMS = 1:nMS
% The following line was originally used in stmean function
% to calculate Euclidian distance
%    d=sqrt((xMS-xMS(iMS)).^2+(yMS-yMS(iMS)).^2);

    distval = distance([yMS(iMS),xMS(iMS)],[yMS,xMS]);
    d = deg2km(distval);
    
    idxMSloc = find(d<dNeib);
    mss(iMS) = sum(ms(idxMSloc).*exp(-d(idxMSloc)/ar)');
    mss(iMS) = mss(iMS)/(sum(exp(-d(idxMSloc)/ar)));
end

%
%  Calculate the temporal averages mt
%
for iME = 1:nME
    mt(iME) = mean(Z(~isnan(Z(:,iME)),iME));
end

%
%  smooth the temporal average with an exponential filter to get mts
%
for iME = 1:nME
    t = abs(tME-tME(iME));
    if tloop>0
        t = min(t,tloop-t); 
    end
    idxMEloc = find(t<tNeib);
    mts(iME) = sum(mt(idxMEloc).*exp(-t(idxMEloc)/at));
    mts(iME) = mts(iME)/(sum(exp(-t(idxMEloc)/at)));
end

ms=ms';
mss=mss';
