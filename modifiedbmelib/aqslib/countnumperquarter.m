function [numperquarter,lengthquarter] = countnumperquarter(serialdate)

% countnumperquarter - Count the number of records in each quarter of
%                      the year (Nov 28,2012)
%
% [numperquarter,lengthquarter] = countnumperquarter(serialdate)
%
% INPUT:
%
% serialdate(n by 1): Serial date number
%
% OUTPUT :
%
% numperquarter(1 by 4): Number of records in each quarter of the year
% lengthquarter(1 by 4): Length (in days) of quarter of the year
%

if size(unique(year(serialdate)),1) > 1
    error('Input date must have the same single calender year');
end

monthval = month(serialdate);
numperquarter = [sum(monthval <= 3),...
                 sum(monthval > 3 & monthval <= 6),...
                 sum(monthval > 6 & monthval <= 9),...
                 sum(monthval > 9)];
yearval = unique(year(serialdate));
lengthquarter = [31 + (datenum(yearval,3,1) - datenum(yearval,2,1)) + 31,...
                 30 + 31 + 30, 31 + 31 + 30, 31 + 30 + 31];
