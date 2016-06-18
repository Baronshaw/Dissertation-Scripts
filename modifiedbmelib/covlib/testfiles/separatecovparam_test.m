function separatecovparam_test()

covmodel = {'gaussianC','sphericalC'};
covparam = {[1,2],[2,4]};

[sill,covrange] = separatecovparam(covmodel,covparam);

disp('-----------set 1-----------');
disp(['Sill    :',num2str(sill)]);
disp(['Covrange:',num2str(covrange)]);

covmodel = {'nuggetC','sphericalC'};
covparam = {[1],[2,4]};

[sill,covrange] = separatecovparam(covmodel,covparam);

disp('-----------set 2-----------');
disp(['Sill    :',num2str(sill)]);
disp(['Covrange:',num2str(covrange)]);

covmodel = {'sphericalC','nuggetC'};
covparam = {[2,4],[1]};

[sill,covrange] = separatecovparam(covmodel,covparam);

disp('-----------set 3-----------');
disp(['Sill    :',num2str(sill)]);
disp(['Covrange:',num2str(covrange)]);
