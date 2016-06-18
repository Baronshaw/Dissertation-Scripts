function startup()

% startup - Add modifiedbmelib path (Dec 05,2012)
%
% startup
%

curpath = pwd;
addpath(curpath,'-end');
addpath([curpath,'/aqslib'],'-end');
addpath([curpath,'/figlib'],'-end');
addpath([curpath,'/covlib'],'-end');
addpath([curpath,'/offsetlib'],'-end');
addpath([curpath,'/GISfiles'],'-end');
addpath([curpath,'/softlib'],'-end');
% addpath([curpath,'\stcovlib'],'-end');

disp('modifiedbmelib is ready to use');
