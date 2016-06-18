function [xmin,xmax,ymin,ymax] = calcfigureboundary(sptlcoord,marginratio)

% calcfigureboundary - Calculate the limit of x and y axes based on
%                      spatial coordinate (Nov 28,2012)
%
% [xmin,xmax,ymin,ymax] = calcfigureboundary(sptlcoord,marginratio)
%
% INPUT:
%
% sptlcoord(n by 2): Spatial Coordinate
%                    1st column is x
%                    2nd column is y
% marginratio(scalar): Take the margin based on this number (optional)
%                      If omitted, no margin
%                      Margin length is given by the following eq.
%                      max([xmax - xmin, ymax - ymin] * marginratio
%
% OUTPUT :
%
% xmin, xmax, ymin, ymax: Boundary of x and y axis
%

if nargin < 2
    marginratio = 0;
end

xmin = min(sptlcoord(:,1));
xmax = max(sptlcoord(:,1));
ymin = min(sptlcoord(:,2));
ymax = max(sptlcoord(:,2));

delx = (xmax - xmin) * marginratio;
dely = (ymax - ymin) * marginratio;
delval = max(delx, dely);

xmin = xmin - delval;
xmax = xmax + delval;
ymin = ymin - delval;
ymax = ymax + delval;
