function [] = todeliver()
% this function will deliver what I need to Yasu to perform the QA analysis

% boundaries for the US
lax = [-2377710.9671949,2431053.02030718,-1674346.09967391,1368406.15170322];

% the 500 random locations
rand('seed',0);
simulocs = 500;
ax = [-1969760 1852823 -922563 663551];
simux = rand(simulocs,1)*(ax(2)-ax(1)) + ax(1); % x coordinates
simuy = rand(simulocs,1)*(ax(4)-ax(3)) + ax(3); % y coordinates

% the 20 locations from the simulated data
rand('seed',0);
simuQA = randsample(simulocs,20);
simuxQA = simux(simuQA);
simuyQA = simuy(simuQA);

% the 40 locations from the participant data
load('../matfiles_estimation/partdata.mat');
rand('seed',0);
partQA = randsample(length(partid),40);
partxQA = partx(partQA);
partyQA = party(partQA);

save('../matfiles_QAQC/forQA.mat','lax','ax','simulocs','simuQA','simux','simuy','simuxQA', ...
    'simuyQA','partQA','partxQA','partyQA');

end