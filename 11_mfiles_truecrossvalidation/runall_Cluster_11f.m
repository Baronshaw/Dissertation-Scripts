function [] = runall_Cluster_11f()
% this function will display results of the true 10foldCV

% load BME
cd ../BMELIB2.0b
startup();
cd ../11_mfiles_truecrossvalidation

% mean trend maps
for FOLDIDX = 1:10
    displayMeanTrend(FOLDIDX);
    close all;
end

% mean trend maps of soft data
for FOLDIDX = 1:10
    displayMeanTrendSoft(FOLDIDX);
    close all;
end

% covariance plots/tables
for FOLDIDX = 1:10
    displayCovariancePlots(FOLDIDX);
    close all;
end
displayCovarianceVariables();

end