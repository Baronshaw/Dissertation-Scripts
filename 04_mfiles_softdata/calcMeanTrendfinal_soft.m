function [] = calcMeanTrendfinal_soft(MEANIDX)
% this function will calculate the final mean trends for ALL the
% soft data

if nargin < 1, MEANIDX = 4; end

% getting mean trend parameters
ars = 300000;
ats = 50;

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

% finding all locations where mean trend needs to be calculated
% loading all the modeled data
CTMyears = [2001 2002 2005 2006 2007];
CTMyears = CTMyears(MEANIDX);
for i = 1:length(CTMyears);
    
    if CTMyears(i) == 2006
        load(sprintf('../matfiles/prepCTM_%d.mat',CTMyears(i)));  
        pI{1,1} = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];
        load(sprintf('../matfiles/prepCTM_%dsum.mat',CTMyears(i)));  
        pI{2,1} = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];
        pI = cell2mat(pI);
    else
        load(sprintf('../matfiles/prepCTM_%d.mat',CTMyears(i)));  
        pI = [distCTMv datenum(yrmodaCTMv(:,1),yrmodaCTMv(:,2),yrmodaCTMv(:,3))];
    end

    cd ../10_mfiles_newmeantrend
    tic
    [mI]=expKernelSmooth_stv_parallel(pd,zd,[2*900000 ars 2*ats ats],pI);  % '2*' added later
    toc
    cd ../04_mfiles_softdata  
    % saving results
    save(sprintf('../matfiles/meanTrend_%d_%d_%d_%d_soft_yr%d.mat',[900000 ars 2*ats ats],CTMyears(i)), ...
        'pd','zd','pI','mI');

end

end