function [] = runall_Cluster_08b(ESTIDX)
% this function will run all the data in the 08_mfiles_QAQC folder

if nargin < 1, ESTIDX = 40; end % ESTIDX is from 3 to 160

% load BME
cd ../BMELIB2.0b
startup();
cd ../08_mfiles_QAQC

% calculating the BME mean and variance for WHIMS participants 
disp(datevec(datenum(1998,ESTIDX,1)));
WHIMS_est(ESTIDX);

end