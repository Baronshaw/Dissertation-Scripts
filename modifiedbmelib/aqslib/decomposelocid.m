function [statecode,countycode,siteid] = decomposelocid(locid)

% decomposelocid - Decompose numeric location ID of AQS site into
%                  Numeric state code, county code and site ID 
%                  (Nov 28,2012)
%
% [statecode,countycode,siteid] = decomposelocid(locid)
%
% INPUT:
%
% locid(n by 1): Numeric location ID of AQS sites
%
% OUTPUT :
%
% statecode(n by 1): State code
% countycode(n by 1): County code
% siteid(n by 1): Site ID
%

strlocid = num2str(locid);
statecode = str2num(strlocid(:,1:2));
countycode = str2num(strlocid(:,3:5));
siteid = str2num(strlocid(:,6:9));
