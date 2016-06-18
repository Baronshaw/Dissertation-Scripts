function decomposecovmodel_test()

covmodel = {'gaussianC/sphericalC','sphericalC/exponentialC'};
covparam = {[1,2,4],[2,4,8]};

[sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam);

disp('-----------set 1-----------');
disp(['Sptl Cov Model 1:',sptlcovmodel{1}]);
disp(['Sptl Cov Param 1:',num2str(sptlcovparam{1})]);
disp(['Sptl Cov Model 2:',sptlcovmodel{2}]);
disp(['Sptl Cov Param 2:',num2str(sptlcovparam{2})]);

disp(['Temp Cov Model 1:',tempcovmodel{1}]);
disp(['Temp Cov Param 1:',num2str(tempcovparam{1})]);
disp(['Temp Cov Model 2:',tempcovmodel{2}]);
disp(['Temp Cov Param 2:',num2str(tempcovparam{2})]);

covmodel = {'nuggetC/sphericalC','sphericalC/exponentialC'};
covparam = {[1,4],[2,4,8]};

[sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam);
 
disp('-----------set 2-----------');
disp(['Sptl Cov Model 1:',sptlcovmodel{1}]);
disp(['Sptl Cov Param 1:',num2str(sptlcovparam{1})]);
disp(['Sptl Cov Model 2:',sptlcovmodel{2}]);
disp(['Sptl Cov Param 2:',num2str(sptlcovparam{2})]);

disp(['Temp Cov Model 1:',tempcovmodel{1}]);
disp(['Temp Cov Param 1:',num2str(tempcovparam{1})]);
disp(['Temp Cov Model 2:',tempcovmodel{2}]);
disp(['Temp Cov Param 2:',num2str(tempcovparam{2})]);

covmodel = {'nuggetC/nuggetC','sphericalC/exponentialC'};
covparam = {1,[2,4,8]};

[sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam);
 
disp('-----------set 3-----------');
disp(['Sptl Cov Model 1:',sptlcovmodel{1}]);
disp(['Sptl Cov Param 1:',num2str(sptlcovparam{1})]);
disp(['Sptl Cov Model 2:',sptlcovmodel{2}]);
disp(['Sptl Cov Param 2:',num2str(sptlcovparam{2})]);

disp(['Temp Cov Model 1:',tempcovmodel{1}]);
disp(['Temp Cov Param 1:',num2str(tempcovparam{1})]);
disp(['Temp Cov Model 2:',tempcovmodel{2}]);
disp(['Temp Cov Param 2:',num2str(tempcovparam{2})]);

covmodel = {'gaussianC/nuggetC','sphericalC/exponentialC'};
covparam = {[1,2],[2,4,8]};

[sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam);
 
disp('-----------set 4-----------');
disp(['Sptl Cov Model 1:',sptlcovmodel{1}]);
disp(['Sptl Cov Param 1:',num2str(sptlcovparam{1})]);
disp(['Sptl Cov Model 2:',sptlcovmodel{2}]);
disp(['Sptl Cov Param 2:',num2str(sptlcovparam{2})]);

disp(['Temp Cov Model 1:',tempcovmodel{1}]);
disp(['Temp Cov Param 1:',num2str(tempcovparam{1})]);
disp(['Temp Cov Model 2:',tempcovmodel{2}]);
disp(['Temp Cov Param 2:',num2str(tempcovparam{2})]);

covmodel = 'sphericalC/exponentialC';
covparam = [2,4,8];

[sptlcovmodel,sptlcovparam,tempcovmodel,tempcovparam] = ...
             decomposecovmodel(covmodel,covparam);
 
disp('-----------set 5-----------');
disp(['Sptl Cov Model 1:',sptlcovmodel{1}]);
disp(['Sptl Cov Param 1:',num2str(sptlcovparam{1})]);

disp(['Temp Cov Model 1:',tempcovmodel{1}]);
disp(['Temp Cov Param 1:',num2str(tempcovparam{1})]);
