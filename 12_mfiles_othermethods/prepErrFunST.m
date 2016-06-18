function [rLags,Crtests,problems] = prepErrFunST(chsub,chtsub,zhosub,zhmsub,lenp)
% prep objective function

rLags = cell(lenp,lenp);
Crtests = cell(lenp,lenp);
zhmp = prctile(zhmsub,0:lenp:100);

for i = 1:lenp
    for j = i:lenp
        
        idxp = zhmsub >= zhmp(i) & zhmsub < zhmp(i+1);
        [Z1,cMS1,tME1,nanratio]=valstv2stg([chsub(idxp,:) chtsub(idxp)],zhosub(idxp));
        idxp = zhmsub >= zhmp(j) & zhmsub < zhmp(j+1);
        [Z2,cMS2,tME2,nanratio]=valstv2stg([chsub(idxp,:) chtsub(idxp)],zhosub(idxp));

        % spatial distances 
        DMS = sqrt(bsxfun(@plus,dot(cMS1,cMS1,2),dot(cMS2,cMS2,2)')-2*(cMS1*cMS2'));
        rLags{i,j} = [0 prctile(unique(DMS(:)),[1 2:10 12.5 15:5:50])];
        rTol = [0 (rLags{i,j}(2:end)-rLags{i,j}(1:end-1))/2];

        % experimental covariance
        if size(Z1,1) > 1 & size(Z2,1) > 1
            [Crtests{i,j} nprtest]=stcov(Z1,cMS1,tME1,Z2,cMS2,tME2,rLags{i,j},rTol,0,0);
            idxnan = isnan(Crtests{i,j}); Crtests{i,j}(idxnan) = []; rLags{i,j}(idxnan) = [];
        end

    end
end

for i = 1:lenp
    problems(i) = (zhmp(i+1)-zhmp(i))./2 + zhmp(i); 
end

end
