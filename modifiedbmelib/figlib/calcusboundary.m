function [xmin,xmax,ymin,ymax] = calcusboundary(marginratio)

% calcusboundary - Calculate the limit of latitude and longitude of
%                  the contiguous United States (Nov 30,2012)
%
% [xmin,xmax,ymin,ymax] = calcusboundary(marginratio)
%
% INPUT:
%
% marginratio(scalar): Take the margin based on this number (optional)
%                      If omitted, no margin
%                      Margin length is given by the following eq.
%                      max([xmax - xmin, ymax - ymin] * marginratio
%
% OUTPUT :
%
% xmin, xmax, ymin, ymax: Boundary of longitude and latitude axis
%

if nargin < 1
    marginratio = 0;
end

s = shaperead('GISfiles/contig_dtl_st');
listx = [];
listy = [];
for j = 1:size(s,1)
    listx = [listx,s(j).X];
    listy = [listy,s(j).Y];
end

xmin = min(listx);
xmax = max(listx);
ymin = min(listy);
ymax = max(listy);

delx = (xmax - xmin) * marginratio;
dely = (ymax - ymin) * marginratio;
delval = max(delx, dely);

xmin = xmin - delval;
xmax = xmax + delval;
ymin = ymin - delval;
ymax = ymax + delval;
