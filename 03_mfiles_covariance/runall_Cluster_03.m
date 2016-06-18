function [] = runall_Cluster_03()
% this function will run all the data in the 03_mfiles_covariance folder

% load BME
cd ../BMELIB2.0b
startup();
cd ../03_mfiles_covariance

% get experimental covariance model
getExpCovariance();

% explore covariance models
getCovModel();

% plot covariance models
plotCovModel();

% explore joint covariance models
getCovModelJoint();

% plot join covariance models
plotCovModelJoint();

% calculate covariance model fittings
displayCovModel();

% dominance plots
plotDominance();

% calculate and plot constant mean trend
CovModelConstant();











end