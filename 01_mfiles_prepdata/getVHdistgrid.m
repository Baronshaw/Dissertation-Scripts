function [vertdist,hortdist,vertgrid,hortgrid] = getVHdistgrid(name,value)
% this function will output the horitzontal and vertical distance between
% the origin and the centriods of each grid
 
% naming some variables needed for this function
idx = strcmp(name,'NCOLS'); ncols = value(idx);
idx = strcmp(name,'NROWS'); nrows = value(idx);
idx = strcmp(name,'XCENT'); xcent = value(idx);
idx = strcmp(name,'YCENT'); ycent = value(idx);
idx = strcmp(name,'XORIG'); xorig = value(idx);
idx = strcmp(name,'YORIG'); yorig = value(idx);
idx = strcmp(name,'XCELL'); xcell = value(idx);
idx = strcmp(name,'YCELL'); ycell = value(idx);

% converting into the right type
xorig = str2num(num2str(xorig{1}));
xcell = str2num(num2str(xcell{1}));
yorig = str2num(num2str(yorig{1}));
ycell = str2num(num2str(ycell{1}));
nrows = str2num(num2str(nrows{1}));
ncols = str2num(num2str(ncols{1}));
ycent = str2num(num2str(ycent{1}));
xcent = str2num(num2str(xcent{1}));

% find the grid row and column that's at the lower left hand corner of the
% origin
ncolorig = abs(xorig/xcell);
nroworig = abs(yorig/ycell);

% finding distances between origin and centriod of grids
nrowsmat = 1:nrows;  
vertdist = (nrowsmat - nroworig).*ycell; % vertical distance
vertdist = vertdist - (ycell/2); % centriod of the grid
ncolsmat = 1:ncols;
hortdist = (ncolsmat - ncolorig).*xcell; % horizontal distance
hortdist = hortdist - (xcell/2); % centriod of the gird
[hortgrid,vertgrid] = meshgrid(hortdist,vertdist);

end