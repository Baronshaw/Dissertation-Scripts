function plotusstate(boundarycolor)

% plotusstate - Plot US state boundary (Nov 30,2012)
%
% plotusstate(boundarycolor)
%
% INPUT:
%
% boundarycolor(string): Color of state boundary
%                        Default color is blue ('b')
%

if nargin < 1
    boundarycolor = 'b';
end

s = shaperead('GISfiles/contig_dtl_st');

for j = 1:size(s,1)
    plot(s(j).X,s(j).Y,'-','Color',boundarycolor);
    hold on;
end
