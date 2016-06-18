function [slag,scovval,snumpair,sweight,...
          tlag,tcovval,tnumpair,tweight] = ...
             calcstcovmom(sdistmat,tdistmat,val,slagbound,tlagbound)

% calcstcovmom - Calculate space/time covariance using method of moments
%                (Nov 28,2012)
%
% [slag,scovval,snumpair,sweight,tlag,tcovval,tnumpair,tweight] = ...
%     calcstcovmom(sdistmat,tdistmat,val,slagbound,tlagbound)
%
% INPUT:
%
% sdistmat(n by n): Diatance matrix
% tdistmat(m by m): Temporal distance matrix
% val(n by m): Value
% slagbound(1 by x): Bound of lag
% tlagbound(1 by y): Bound of lag
%
% OUTPUT :
%
% slag(1 by x): Spatial lag
% scovval(1 by x): Experimental covariance
% snumpair(1 by x): Number of pairs used to calculate experimental covariance
% sweight(1 by x): Weight for general least square covarinace model fit
%
% tlag(1 by y): Temporal lag
% tcovval(1 by y): Experimental covariance
% tnumpair(1 by y): Number of pairs used to calculate experimental covariance
% tweight(1 by y): Weight for general least square covarinace model fit
%

tempval = val(:);
meanval = mean(tempval(~isnan(tempval)));
varval = var(tempval(~isnan(tempval)),1);
numval = size(tempval(~isnan(tempval)),1);
zerocovvals = (tempval(~isnan(tempval)) - meanval).^2;
zerocovvar = var(zerocovvals);

tlag = zeros(1,size(tlagbound,2)-1);
tcovval = zeros(1,size(tlagbound,2)-1);
tnumpair = zeros(1,size(tlagbound,2)-1);
tweight = zeros(1,size(tlagbound,2)-1);

for i = 1:size(tlagbound,2) - 1
    [idxrow,idxcol] = find(tlagbound(i) < tdistmat & tdistmat <= tlagbound(i+1));
    tdistvals = tdistmat(tlagbound(i) < tdistmat & tdistmat <= tlagbound(i+1));

    if isempty(idxrow)
        tlag(i) = NaN;
        tcovval(i) = NaN;
        tnumpair(i) = NaN;
        tweight(i) = NaN;
    else
        idxuniquepair = find(idxrow>=idxcol);
        idxrow = idxrow(idxuniquepair);
        idxcol = idxcol(idxuniquepair);
        tdistvals = tdistvals(idxuniquepair);

        valhead = val(:,idxrow);
        valtail = val(:,idxcol);

        valhead = valhead(:);
        valtail = valtail(:);
        tdistlist = kron(tdistvals,ones(size(val,1),1));

        idxvalid = find(~isnan(valhead) & ~isnan(valtail));

        if isempty(idxvalid)
            tlag(i) = NaN;
            tcovval(i) = NaN;
            tnumpair(i) = NaN;
            tweight(i) = NaN;
        else
            valhead_valid = valhead(idxvalid);
            valtail_valid = valtail(idxvalid);

            tdistpair = tdistlist(idxvalid);
            tcovpair = (valhead_valid - meanval).*(valtail_valid - meanval);

            tlag(i) = mean(tdistpair);
            tcovval(i) = mean(tcovpair);
            tnumpair(i) = size(idxvalid,1);
            tweight(i) = 1/var(tcovpair);
        end
    end
end

tlag = [0,tlag];
tcovval = [varval,tcovval];
tnumpair = [numval,tnumpair];
tweight = [1/zerocovvar,tweight];

slag = zeros(1,size(slagbound,2)-1);
scovval = zeros(1,size(slagbound,2)-1);
snumpair = zeros(1,size(slagbound,2)-1);
sweight = zeros(1,size(slagbound,2)-1);

for i = 1:size(slagbound,2) - 1 
    [idxrow,idxcol] = find(slagbound(i) < sdistmat & sdistmat <= slagbound(i+1));
    sdistvals = sdistmat(slagbound(i) < sdistmat & sdistmat <= slagbound(i+1));

    if isempty(idxrow)
        slag(i) = NaN;
        scovval(i) = NaN;
        snumpair(i) = NaN;
        sweight(i) = NaN;
    else
        idxuniquepair = find(idxrow>=idxcol);
        idxrow = idxrow(idxuniquepair);
        idxcol = idxcol(idxuniquepair);
        sdistvals = sdistvals(idxuniquepair);

        valhead = val(idxrow,:);
        valtail = val(idxcol,:);
        sdistlist = kron(sdistvals',ones(1,size(val,2)));

        thead = valhead';
        ttail = valtail';
        sdistlist = sdistlist';

        valhead = thead(:);
        valtail = ttail(:);

        idxvalid = find(~isnan(valhead) & ~isnan(valtail));

        if isempty(idxvalid)
            slag(i) = NaN;
            scovval(i) = NaN;
            sweight(i) = NaN;
        else
            valhead_valid = valhead(idxvalid);
            valtail_valid = valtail(idxvalid);

            sdistpair = sdistlist(idxvalid);
            scovpair = (valhead_valid - meanval).*(valtail_valid - meanval);

            slag(i) = mean(sdistpair);
            scovval(i) = mean(scovpair);
            snumpair(i) = size(idxvalid,1);
            sweight(i) = 1/var(scovpair);
        end
    end
end

slag = [0,slag];
scovval = [varval,scovval];
snumpair = [numval,snumpair];
sweight = [1/zerocovvar,sweight];
