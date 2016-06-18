function [] = runall()
% this function will run all the scripts in the 13_mfiles_modelperformance
% folder

% loading BME function
cd ../BMELIB2.0b
startup
cd ../13_mfiles_modelperformance

% get network and region of data 
get_info();

% calculating measures of traditional model performance of CMAQ 
stat_MP_cmaq();
stat_MP_cmaqII();

% create measures of traditional model perofrmance on a grid
stat_MP_cmaq_grid();
stat_MP_cmaq_grid_2(); % used to see more of a temporal trend in performance
stat_MP_cmaq_gridbin(); 
patch_MP_cmaq_grid(); % run on cluster
par_stat_MP_cmaq_gridbin(); % run on cluster
stat_MP_cmaq_extendbin(); % run on cluster

% create tables of model performance
table_MP_cmaq();

% create maps of model performance
map_MP_cmaq();

% create time series of model performance
ts_MP_cmaq();

% create maps of model performance on a grid
map_MP_cmaq_grid();
map_MP_cmaq_grid_2(); % used to see more of a temporal trend in performance
map_MP_cmaq_gridlamper();
map_MP_cmaq_gridPatch(); % work in this next

% creating a global S-curve (at least for 2001)
globalSCurve();
plotGlobalSCurve();

% comparing measures with approximations
lambdacompare();
plotSCurve([2001 1 1],20010101);
plotSCurve_2([2001 1 1],20010101); % make SCurves on a grid

% compare model performance with BME
for i = 0:100000:900000
    disp(i);
    stat_MP_bme(i);
end
map_MP_bme();

% gathering all results together
antiPerform();

% applying exclusion criteria
exclusionCriteria();
exclusionCriteria_continuous_singular_cutoffs();
exclusionCriteria_continuous_singular_intervals();
exclusionCriteria_continuous_double_cutoffs();

% maps to show Will
map_MP_will();
findMax_will();
SCurve_will();
findMS_will();

% anticipating performance of RAMP using CMAQ characteristics
for i = 0:100000:900000
    disp(i);
    %antiPerform(i);
    antiPerform_2(i);
end

end