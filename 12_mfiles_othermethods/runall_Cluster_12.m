function [] = runall_Cluster_12()
% this function will run all the functions in this folder

% parameters estimation for static downscaler
staticDownScaler();

% static cross-validation
staticDSXval();

% put DS results in a table
staticDSstats();

% space only RAMP cross-validation
staticRAMPXval();

% put RAMP results in a table
staticRAMPstats();

% static kriging
staticKrigingXval();

% put static kriging results in a table
staticKrigingstats();

% visualize covariance
visualize_covmodel_DS();

% make maps of static DS (estimates and parameters)
for i = 1:3
    staticDS_estimate(i);
    for j = 1:2    
        staticDS_maps(i,j);
        close all;
    end
end

% make maps of static Kriging
staticKriging_estimate();
for i = 1:2    
    staticKriging_maps(i);
    close all;
end

% make maps of static RAMP
staticRAMP_estimate();
for i = 1:2
    staticRAMP_maps(i);  
    close all;
end

% parameter estimation for space/time downscaler
spacetimeDownScaler('ind','ind');
spacetimeDownScaler('dyn','ind');

% cross-validation of the space/time downscaler
spacetimeDSXval('ind','ind');
spacetimeDSXval('dyn','ind');

% put space/time downscaler results in a table
spacetimeDSstats('ind','ind');
spacetimeDSstats('dyn','ind');

% visualize covariance
visualize_covmodel_stDS('ind','ind');
visualize_covmodel_stDS('dyn','ind');

% space/time DS (estimates and parameters) (through cluster)
for i = 1:3
    spacetimeDS_estimate(i,'ind','ind');
    spacetimeDS_estimate(i,'dyn','ind');
end

% maps of space/time DS (estimates and parameters)
for i = 1:3
    for j = 1:2
        spacetimeDS_maps(i,j,'ind','ind');
        close all;
        spacetimeDS_maps(i,j,'dyn','ind');
        close all;
    end
end

% make maps of space/time Kriging
spacetimeKriging_estimate();
for i = 1:2    
    spacetimeKriging_maps(i);
    close all;
end

% make maps of space/time RAMP
spacetimeRAMP_estimate();
for i = 1:3
    spacetimeRAMP_maps(i,1); 
    close all;
end
spacetimeRAMP_maps(1,2); 
close all;

end