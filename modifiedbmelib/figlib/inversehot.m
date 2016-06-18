function cmap = inversehot(numcolor)

% inversehot - Inverse hot colormap (Nov 28,2012)
%
% cmap = inversehot(numcolor)
%
% INPUT:
%
% numcolor(scalar): Number of color scale
%
% OUTPUT :
%
% cmap(numcolor by 3): Colormap
%

if nargin < 1
    cmap = hot();
else
    cmap = hot(numcolor);
end

cmap = cmap(end:-1:1,:);
