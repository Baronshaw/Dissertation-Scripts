function [] = analysis()
% this function will call all the functions in the 
% 17_mfiles_validatelambda1 folder

% loading BME function
cd ../BMELIB2.0b
startup
cd ../17_mfiles_validatelambda1

%%% Compare RAMP with CAMP

% do an implimentation of CAMP
for i = 2001:2002
    runCAMP();
end

% do an implimentation of a regionalized CAMP
runsubCAMP();

% do an implimentation of constant error correction
runConstant();

% gathering all info about space/time locations of obs including: day,
% season, cmaq, lambda1, lambda2, monitoring network, etc.
performanceInfo();

% calculating statistics for all these different space/time characteristics
% and stratifications
stratCMAQRAMP(); 
stratCMAQRAMP_testing();

% put results of stratCMAQRAMP into tables and figures
tableCMAQRAMP();
tableCMAQRAMP_testing();

%%% simulate RAMP

% getting all the mean modeled data needed for the simulation
gatherMeanMod();

% randomly generate "true" lambda1 
for i = 1:10
    for j = 2001:2002
        generatelambda(i,j);
    end
end

% select lambda1 and lambda2 for fixed modeled values from the data files I
% already have
% submitted to cluster - 1 hour
for i = 2001:2002
    selectlambda(i);
end

% plotting the selected lambda1s for fixed modeled values
for i = 2001:2002
    plotselectlambda(i);
end

% plot random fields
for i = 1:10
    dispRF_2001(i);
end

% simulate obs
for i = 2001:2002
    simulateObs(i);
end

% plot simulated obs
for i = 2001:2002
    plotsimuObs(i);
end

% recalculate lambda1/lambda2 across percentiles for each CMAQ grid for 2001
% submitted to cluster - 8 hours
recalculatelambdaGivPer();

recalculatelambdaGivPer_CAMP();

recalculatelambdaGivPer_Constant();

% recalculate lambda1/lambda2 across fixed modeled values for each CMAQ grid for 2001
% submitted to cluster - 1 hour
recalculatelambdaGivMod();

% submitted to cluster - 1 hour
recalculatelambdaGivMod_CAMP();

% submitted to cluster - 1 hour
recalculatelambdaGivMod_Constant();

% compare simulated lambda1 with calculated lambda1
compareSimulambda();
compareSimulambda_CAMP();
compareSimulambda_Constant();

% create table showing R2 for Constant, CAMP, RAMP
tableSimu();

end