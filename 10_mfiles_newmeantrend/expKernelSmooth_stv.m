function [mI]=expKernelSmooth_stv(pd,zd,smoothingParam,pI);

% expKernelSmooth_stv         - exponential kernel smoothing on space/time vector data (v.2011/01/04)
%
% Calculates smoothed values at interpolation points using an 
% exponential kernel on space/time vector data. The data consist of  
% measurements at any arbitrary space/time points.
%
% SYNTAX :
%
% [mI]=expKernelSmooth_stv(pd,zd,smoothingParam,pI);
%
% INPUTS :
%  
%  pd    nd by 3    matrix of space/time coordinates for the data points
%  zd    nd by 1    vector of measurements at the data points
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
% [mI]=expKernelSmooth_stv(pd,zd,smoothingParam,pd);
% [MI]=valstv2stg(pd,mI,cMS,tME);
% for iMS=1:2, figure; plot(tME,Z(iMS,:),'o-'); hold on; plot(tME,MI(iMS,:),'-r'); end
%
% NOTE : 
%
% See also expKernelSmooth_stg.m, expKernelSmoothTest.m
% Use help stgridsyntax for help on space/time grid and vector format



%
% Check the input arguments
%
if size(pd,2)~=3, error('pd must have 3 columns '); end;
nd=size(pd,1);
if size(zd,2)~=1, error('zd must have 1 column'); end;
if size(zd,1)~=nd, error('zd must have the same number of rows as pd'); end;
if length(smoothingParam)<4, error('smoothingParam must be a vector of at least 4 elements'); end
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
xd=pd(:,1);
yd=pd(:,2);
td=pd(:,3);
xI=pI(:,1);
yI=pI(:,2);
tI=pI(:,3);

%
%  smooth with an exponential filter to get mst
%

idx = ~isnan(zd);
zd = zd(idx);
xd = xd(idx);
yd = yd(idx);

mI = NaN*ones(nI,1);
len = length(mI);
for i=1:nI
    %tic
  if mod(i,1000)==0, disp([i len]); end
  d2=(xd-xI(i)).^2+(yd-yI(i)).^2;
  t=abs(td-tI(i));
  if tloop>0
    t=min(t,tloop-t); 
  end;
  idx=d2<dNeib^2 & t<tNeib;
  weights=exp(-sqrt(d2(idx))/ar-t(idx)/at);
  mI(i,1)=sum(zd(idx).*weights(:))/sum(weights);
  %toc
end

end