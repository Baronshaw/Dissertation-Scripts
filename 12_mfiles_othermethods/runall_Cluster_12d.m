function [] = runall_Cluster_12d(FORIDX)
% loop through each distance in static DS

if nargin < 1, FORIDX = 1; end

forcedoptions = 0:100000:1000000;
forceddist = forcedoptions(FORIDX);

CAMP_LOOCV(forceddist);

end