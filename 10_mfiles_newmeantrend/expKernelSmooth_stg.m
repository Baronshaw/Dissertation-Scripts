function [mI]=expKernelSmooth_stg(Z,cMS,idMS,tME,smoothingParam,pI);

% expKernelSmooth_stg         - exponential kernel smoothing on space/time grid data (v.2011/01/04)
%
% Calculates smoothed values at interpolation points using an 
% exponential kernel on space/time grid data. The data consist of  
% measurements at fixed monitoring sites cMS and fixed event times tME
%
% SYNTAX :
%
% [mI]=expKernelSmooth_stg(Z,cMS,idMS,tME,smoothingParam,pI);
%
% INPUTS :
%  
%  Z     nMS by nME matrix of measurements for Z at the nMS monitoring
%                   sites and nME measuring events. Z may have NaN values.
%  cMS   nMS by 2   matrix of spatial x-y coordinates for the nMS monitoring
%                   sites
%  idMS  nMS by 1   vector with a unique id number for each monitoring site
%  tME   1 by nME   vector with the time of the measuring events
%  smoothingParam   1 by 5 parameters to smooth the spatial and temporal average
%                   p(1)=dNeib  distance (radius) of spatial neighborhood
%                   p(2)=ar     spatial range of exponential smoothing function
%                   p(3)=tNeib  time (radius) of temporal neighborhood
%                   p(4)=at     temporal range of exponential smoothing function
%                   p(5)=tloop  is an optional input used for temporal smoothing.
%                        When tloop>0, the measuring events are looped in a
%                    cycle of duration tloop.
% pI     nI by 3    matrix of space/time coordinates for the Interpolation points
%
% OUTPUT :
%
%  mI    nI by 1    vector of smoothed values
%
% EXAMPLE : 
%
% cMS=1000*rand(100,2);
% idMS=[1:100]';
% tME=[1:24];
% Z=1+randn(100,24);
% Z(Z<0)=NaN;
% smoothingParam=[200,100,12,4];
% [pd,zd]=valstg2stv(Z,cMS,tME);
% [mst]=expKernelSmooth_stg(Z,cMS,idMS,tME,smoothingParam,pd);
% [Mst]=valstv2stg(pd,mst,cMS,tME);
% for iMS=1:2, figure; plot(tME,Z(iMS,:),'o-'); hold on; plot(tME,Mst(iMS,:),'-r'); end
%
% NOTE : 
%
% See also expKernelSmooth_stv.m, expKernelSmoothTest.m, stmean_pI.m
% Use help stgridsyntax for help on space/time grid and vector format


nMS=size(cMS,1);
if size(cMS,2)~=2, error('cMS must have two columns'); end;
if size(idMS,1)~=nMS | size(idMS,2)~=1, error('idMS must be a nMS by 1 vector'); end;
if size(tME,1)~=1, error('tME must have one row'); end;
nME=size(tME,2);
if size(pI,2)~=3, error('pI must have 3 columns'); end;
nI=size(pI,1);

%
% set parameters and variables
%
dNeib=smoothingParam(1);      % Radius use to select local neighborhood to average Z
ar=smoothingParam(2);         % spatial range of exponential smoothing function
tNeib=smoothingParam(3);         % Radius use to select local neighborhood to average PM10
at=smoothingParam(4);           % temporal range of exponential smoothing function
if length(smoothingParam)>4
  tloop=smoothingParam(5);
else
  tloop=0;
end; 
xMS=cMS(:,1);
yMS=cMS(:,2);
xI=pI(:,1);
yI=pI(:,2);
tI=pI(:,3);

%
%  smooth with an exponential filter to get mst
%
for i=1:nI
  d2=(xMS-xI(i)).^2+(yMS-yI(i)).^2;
  t=abs(tME-tI(i));
  if tloop>0
    t=min(t,tloop-t); 
  end;
  idxMS=find(d2<dNeib^2);
  idxME=find(t<tNeib);
  [pd,zd]=valstg2stv(Z(idxMS,idxME),cMS(idxMS,:),tME(1,idxME));
  d2=(pd(:,1)-xI(i)).^2+(pd(:,2)-yI(i)).^2;
  t=abs(pd(:,3)-tI(i));
  if tloop>0
    t=min(t,tloop-t); 
  end;
  idx=find(~isnan(zd));
  weights=exp(-sqrt(d2(idx))/ar-t(idx)/at);
  mI(i,1)=sum(zd(idx).*weights(:))/sum(weights);
end;

