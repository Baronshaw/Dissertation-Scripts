function [] = runall_test1(FOLDIDX)
% this function will run all the data in the 05_mfiles_softdata folder

if nargin < 1, FOLDIDX = 1; end

% load BME
cd ../BMELIB2.0b
startup();
cd ../05_mfiles_crossvalidation

% calculate the spatial 10 fold cross validation 
% getXvalidation_10fold_test1(1,0,1,FOLDIDX);
% getXvalidation_LOOCV_test1(0,0,1,FOLDIDX);
% getXvalforcedisolation_10fold_test1(0,0,1,500000,FOLDIDX);
getXvalforcedisolation_LOOCV_test1(1,0,1,500000,FOLDIDX);

end