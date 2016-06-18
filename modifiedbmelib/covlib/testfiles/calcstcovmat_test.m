function calcstcovmat_test()

c1 = [0 0 0;0 1 0;0 2 0;0 0 1;0 1 1;0 2 1;0 0 2;0 1 2;0 2 2];
sptldistmat = calcdistmateuclid(c1(:,1:2),c1(:,1:2));
tempdistmat = calcdistmattemp(c1(:,3)',c1(:,3)');

covmodel = {'exponentialC/exponentialC','exponentialC/gaussianC'};
covparam = {[2,1,2],[3,7,8]};
covmat1 = coord2K_test(c1,c1,covmodel,covparam);
covmat2 = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam);
disp(['Set1: ',num2str(isequal(covmat1,covmat2))]);

covmodel = {'exponentialC/nuggetC','exponentialC/gaussianC'};
covparam = {[2,1],[3,7,8]};
covmat1 = coord2K_test(c1,c1,covmodel,covparam);
covmat2 = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam);
disp(['Set2: ',num2str(isequal(covmat1,covmat2))]);

covmodel = {'nuggetC/exponentialC','exponentialC/gaussianC'};
covparam = {[2,2],[3,7,8]};
covmat1 = coord2K_test(c1,c1,covmodel,covparam);
covmat2 = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam);
disp(['Set3: ',num2str(isequal(covmat1,covmat2))]);

covmodel = {'nuggetC/nuggetC','exponentialC/gaussianC'};
covparam = {2,[3,7,8]};
covmat1 = coord2K_test(c1,c1,covmodel,covparam);
covmat2 = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam);
disp(['Set4: ',num2str(isequal(covmat1,covmat2))]);

c1 = [0 0 0;0 1 0;0 2 0;0 0 1;0 1 1;0 2 1;0 0 2;0 1 2;0 2 2];
c2 = [5 0 0;5 1 0;5 2 0];
sptldistmat = calcdistmateuclid(c1(:,1:2),c2(:,1:2));
tempdistmat = calcdistmattemp(c1(:,3)',c2(:,3)');

covmodel = {'exponentialC/exponentialC','exponentialC/gaussianC'};
covparam = {[2,1,2],[3,7,8]};
covmat1 = coord2K_test(c1,c2,covmodel,covparam);
covmat2 = calcstcovmat(sptldistmat,tempdistmat,covmodel,covparam);
disp(['Set5: ',num2str(isequal(covmat1,covmat2))]);
