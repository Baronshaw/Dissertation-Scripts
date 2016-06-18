function lscovfitsttwoexp_test()

sptllag = 0:15;
sptlcovval = exponentialC(sptllag,[0.7,4]) + ...
             exponentialC(sptllag,[0.3,8]);
sptlweight  = ones(size(sptllag)) + randn(size(sptllag)) * 0.01;

templag = 0:10:200;
tempcovval = exponentialC(templag,[0.7,40]) + ...
             exponentialC(templag,[0.3,80]);
tempweight  = ones(size(templag)) + randn(size(templag)) * 0.01;

covparam = lscovfitsttwoexp(sptllag,sptlcovval,sptlweight,...
                            templag,tempcovval,tempweight);

disp(['two exponential : ',num2str(covparam)]);
