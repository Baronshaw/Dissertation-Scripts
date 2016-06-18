function [] = generatelambda(dec2simu,yr2simu)
% this function will perform a simulation of RAMP and find lambda1 and
% lambda2 for given fixed modeled values from simulated lambda1s and
% lambda2s

% runall_Cluster_17.sh
% ##!/bin/bash
% for ((i=1;i<=10;i++))
% do
%   bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "generatelambda(${i})" -logfile "generatelambda${i}.out"
%   sleep 60
% done

if nargin < 1, dec2simu = 1; end
if nargin < 2, yr2simu = 2001; end

% loading BME function
cd ../BMELIB2.0b
startup
cd ../17_mfiles_validatelambda1

% taken from: 
% http://www.mathworks.com/matlabcentral/fileexchange/27613-random-field-simulation
% -3/ar1 = -1/2c0 => c0 = ar/6

% get some covariance information
load('../matfiles/covmod_r_long_joint exponential exponential_joint.mat');
xdir = f.alp*f.ar1+(1-f.alp)*f.ar2/6; xdir = xdir*(1+dec2simu/10);
ydir = f.alp*f.ar1+(1-f.alp)*f.ar2/6; ydir = xdir*(1+dec2simu/10);
tdir = f.alp*f.at1+(1-f.alp)*f.at2/6; tdir = tdir*(1+dec2simu/10);
sigmaz = f.Cr1;
corr.name = 'exp'; corr.sigma = sigmaz; 
corr.c0 = [xdir ydir tdir]; corr.c1 = [];
corr.A = []; corr.B = [];
    
%%% generate random field 

% get mean modeled space/time locations
load(sprintf('PM2p5_meanMod_yr%d.mat',yr2simu));
% use css
unidays = datenum(yr2simu,1:12,1);
[aidx bidx] = ismember(css(:,3),unidays);
cssSub = css(aidx,:);

% get observed space/time locations
load(sprintf('../matfiles/prepCTMandObs_%d.mat',yr2simu));
yr = floor(yrmodaObs./10000);
mo = floor((yrmodaObs - yr.*10000)./100);
da = yrmodaObs - yr.*10000 - mo.*100;
% use coordObs, yrmodaObs 
cObs = [coordObs datenum(yr,mo,da)];

% combine modeled and obs locations
cGen = [cssSub ; cObs];
vecGen = datevec(cGen(:,3));
unimo = unique(vecGen(:,2));

% loop through each month
FcGen = NaN*ones(length(cObs),1);
for i = 1:length(unimo)
    disp(i);      
    idx = vecGen(:,2) == unimo(i);

    tic
    cd('gp/file_exchange3')    
    C = correlation_fun(corr,cGen(idx,:));
    corr.C = C;   
    [Fday,dummy] = randomfield(corr,cGen(idx,:),'filter',0.95);
    cd('../..')
    toc     
    
    % putting the output of the RF in the correct final matrix
    [aidx bidx] = ismember(cGen,cGen(idx,:),'rows');
    FcGen(aidx) = Fday;
end

% save file
save(sprintf('RFlambda1_%d_decile_%0.2d.mat',yr2simu,dec2simu),'cGen','FcGen','-v7.3');

end