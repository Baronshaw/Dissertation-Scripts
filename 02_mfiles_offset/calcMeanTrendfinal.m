function [] = calcMeanTrendfinal()
% this function will calculate the final mean trends for ALL the
% observational data

% getting mean trend parameters
ars = [20000 50000 300000 1000000];
ats = [10 20 50 200];

% loading all the data
for i = 1999:2010
    load(sprintf('../matfiles/prepObs_%d.mat',i)); 
    ch{i-1998,1} = coordObs;
    yrs = floor(yrmodaObs./10000);
    mos = floor( (yrmodaObs-10000*yrs) ./ 100 );
    das = yrmodaObs - 10000*yrs - 100*mos;
    cht{i-1998,1} = datenum(yrs,mos,das);
    zh{i-1998,1} = Obs;
end
ch = [cell2mat(ch) cell2mat(cht)];
zh = cell2mat(zh);
pd = ch;
zd = zh;
pI = pd;

for i = 1:length(ars)
    cd ../10_mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv_parallel(pd,zd,[900000 ars(i) 2*ats(i) ats(i)],pI);
    toc
    cd ../02_mfiles_offset 
    % saving results
    save(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',[900000 ars(i) 2*ats(i) ats(i)]), ...
        'pd','zd','pI','mI');
end

end