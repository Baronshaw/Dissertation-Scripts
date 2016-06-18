function [] = runall_Cluster_05e_hard(FORIDX)
% this function will run all the data in the 05_mfiles_softdata folder

if nargin < 1, FORIDX = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../05_mfiles_crossvalidation

forcedoptions = 0:100000:1000000;
forceddist = forcedoptions(FORIDX);

% calculate the spatial LOOCV cross validation statistics
getXvalforcedisolation_LOOCV(0,0,1,forceddist);
%getXvalforcedisolation_LOOCV(1,0,1,forceddist);

end