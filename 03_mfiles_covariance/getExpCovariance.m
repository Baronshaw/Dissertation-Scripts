function [] = getExpCovariance()
% this function will calculate the experimental covariance for all of the
% observational data

% getting mean trend parameters
ars = [20000 50000 300000 1000000];
ats = [10 20 50 200];
meanNow = {'short','intermediate','long','very_long'};

for i = 1:length(ars)
    
    load(sprintf('../matfiles/meanTrend_%d_%d_%d_%d.mat',900000,ars(i),2*ats(i),ats(i)));
    
    ch = pd;
    zh = zd;
    zh_mrmvd = zh - mI;
    [Z_mrmvd,cMS,tME,nanratio]=valstv2stg(ch,zh_mrmvd);

    % in order to come up with rLag, I need to know what the spatial
    % distances are and to equal spacing in a log space
    DMS = sqrt(bsxfun(@plus,dot(cMS,cMS,2),dot(cMS,cMS,2)')-2*(cMS*cMS'));
    DME = abs(bsxfun(@minus,tME,tME'))';

    % spatial covariance, tic/toc ~= 20 seconds
    rLag = [0 prctile(unique(DMS(:)),[0.25 0.5 0.75 1 1.5 2:10 12.5 15:5:50])];
    rTol = [0 (rLag(2:end)-rLag(1:end-1))/2];

    % temporal covariance, tic/toc ~= 15 seconds
    tLag = [0:10 15:5:40 50:25:150];
    tTol = [0 repmat(0.5,1,10) repmat(2.5,1,6) repmat(25,1,5)];

    % trying it a new way that's mine to compare to stcov
    tic
    [Crtest nprtest]=stcov(Z_mrmvd,cMS,tME,Z_mrmvd,cMS,tME,rLag,rTol,0,0);
    toc
    tic
    [Cttest npttest]=stcov(Z_mrmvd,cMS,tME,Z_mrmvd,cMS,tME,0,0,tLag,tTol);
    toc

    % diplaying results
    figure; hold on;
    plot(rLag,Crtest,'bo');
    figure; hold on;
    plot(tLag,Cttest,'bo');

    % saving results of experimental covariance
    save(sprintf('../matfiles/expcov_%s.mat',meanNow{i}), ...
        'Z_mrmvd','cMS','tME','rLag','rTol','tLag','tTol','Crtest','nprtest',...
        'Cttest','npttest');

end

end