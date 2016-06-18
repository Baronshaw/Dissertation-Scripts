function locid = composelocid(statecode,countycode,siteid)

% decomposelocid - Decompose numeric location ID of AQS site into
%                  Numeric state code, county code and site ID 
%                  (Nov 28,2012)
%
% locid = composelocid(statecode,countycode,siteid)
%
% INPUT:
%
% statecode(n by 1): State code
% countycode(n by 1): County code
% siteid(n by 1): Site ID
%
% OUTPUT :
%
% locid(n by 1): Numeric location ID of AQS sites
%

formatstatecode = arrayfun(@(x) sprintf('%02d',x), statecode,...
                           'UniformOutput', false);
formatcountycode = arrayfun(@(x) sprintf('%03d',x), countycode,...
                            'UniformOutput', false);
formatsiteid = arrayfun(@(x) sprintf('%04d',x), siteid,...
                        'UniformOutput', false);

formatlocid = strcat(formatstatecode,formatcountycode,formatsiteid);

locid = str2num(cell2mat(formatlocid));