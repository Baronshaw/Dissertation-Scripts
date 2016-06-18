function [] = runall_Cluster_05f(yearz)
% this function will run all the data in the 05_mfiles_softdata folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../05_mfiles_crossvalidation

% calculate the spatial 10 fold cross validation 
 SeeEffectOfLambda1(1,0,1,yearz);

end