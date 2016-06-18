function serialdate = calcnumdaysfromstr(numericdate)

% calcnumdaysfromstr - Calculate serial date number from numeric with
%                      'yyyymmdd' format (Nov 28,2012)
%
% serialdate = calcnumdaysfromstr(numericdate)
%
% INPUT:
%
% numericdate(n by 1): Numeric date with 'yyyymmdd' format
%
% OUTPUT :
%
% serialdate(n by 1): Serial date number
%

numdatestr = num2str(numericdate);
yearval = cellfun(@str2double,cellstr(numdatestr(:,1:4)));
monthval = cellfun(@str2double,cellstr(numdatestr(:,5:6)));
dayval = cellfun(@str2double,cellstr(numdatestr(:,7:8)));
serialdate = datenum(yearval,monthval,dayval);
