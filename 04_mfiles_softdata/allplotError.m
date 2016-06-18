function [] = allplotError(soft,constant,gauss)
% this function will disp maps of error for every day in 2001

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

dayz = datevec(datenum(2001,1,1):datenum(2001,12,31));

for i = 1:length(dayz)
    plotError(soft,constant,gauss,[dayz(i,1) dayz(i,2) dayz(i,3)]);
    if length(findall(0,'type','figure')) >= 20, close all; end
end

end