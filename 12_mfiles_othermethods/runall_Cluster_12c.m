function [] = runall_Cluster_12c(FORIDX)
% loop through each distance in static DS

if nargin < 1, FORIDX = 1; end

forcedoptions = 0:100000:1000000;
forceddist = forcedoptions(FORIDX);

spacetimeDSXvalforcedisolation_LOOCV('dyn','ind',forceddist);

end