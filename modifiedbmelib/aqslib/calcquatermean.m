function [qmeanval,numperquarter,lengthquarter] = ...
             calcquatermean(serialdate,value)

% calcquatermean - Calculate the mean of the vlaues in each quarter of
%                  the year (Dec 10,2012)
%
% [qmeanval,numperquarter,lengthquarter] = ...
%              calcquatermean(serialdate,value)
%
% INPUT:
%
% serialdate(n by 1): Serial date number
% value(n by 1): Value
%
% OUTPUT :
%
% qmeanval(1 by 4): The mean of the value in each quarter of the year
% numperquarter(1 by 4): Number of records in each quarter of the year
% lengthquarter(1 by 4): Length (in days) of quarter of the year
%

if size(unique(year(serialdate)),1) > 1
    error('Input date must have the same single calender year');
end

monthval = month(serialdate);

logiQ1 = (monthval <= 3);
logiQ2 = (monthval > 3 & monthval <= 6);
logiQ3 = (monthval > 6 & monthval <= 9);
logiQ4 = (monthval > 9);

qmeanval = [mean(value(logiQ1)),mean(value(logiQ2)),...
            mean(value(logiQ3)),mean(value(logiQ4))];
numperquarter = [sum(logiQ1),sum(logiQ2),sum(logiQ3),sum(logiQ4)];
yearval = unique(year(serialdate));
lengthquarter = [31 + (datenum(yearval,3,1) - datenum(yearval,2,1)) + 31,...
                 30 + 31 + 30, 31 + 31 + 30, 31 + 30 + 31];
