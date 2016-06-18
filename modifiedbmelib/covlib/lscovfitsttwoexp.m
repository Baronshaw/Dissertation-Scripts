function [covparam,aicval] = lscovfitsttwoexp(slag,scovval,sweight,...
                                              tlag,tcovval,tweight)

% lscovfitsttwoexp - Calculate covariance parameter for space/time two
%                    exponential model using weighted least square fit
%                    (Nov 28,2012)
%
% [covparam,aicval] = lscovfitsttwoexp(slag,scovval,sweight,...
%                                      tlag,tcovval,tweight)
%
% INPUT:
%
% slag(1 by n): Spatial lag
% scovval(1 by n): Sample spatial covariance
% sweigh(1 by n): Weight for spatial component
% tlag(1 by m): Temporal lag
% tcovval(1 by m): Sample temporal covariance
% tweigh(1 by m): Weight for temporal component
%
% OUTPUT :
%
% covparam(1 by 6): Covariance parameter
%                   covparam(1) = Sill for 1st exponenital model
%                   covparam(2) = Spatial range  for 1st exponenital model
%                   covparam(3) = Temporal range  for 1st exponenital model
%                   covparam(4) = Sill for 2nd exponenital model
%                   covparam(5) = Spatial range  for 2nd exponenital model
%                   covparam(6) = Temporal range  for 2nd exponenital model
% aicval(scalar): AIC
%

sslopeval = ((scovval(3) - scovval(2)) / (slag(3) - slag(2)));
sintcep = scovval(2) - sslopeval * slag(2);

tslopeval = ((tcovval(3) - tcovval(2)) / (tlag(3) - tlag(2)));
tintcep = tcovval(2) - tslopeval * tlag(2);

init_sill1 = min([scovval(1) - sintcep,tcovval(1) - tintcep]);
if init_sill1 < 0
    init_sill1 = 0;
end
init_sill2 = scovval(1) - init_sill1;
init_srange1 = max(slag)/2;
init_srange2 = max(slag)/2;
init_trange1 = max(tlag)/2;
init_trange2 = max(tlag)/2;

initval = [init_sill1,init_srange1,init_trange1,...
           init_sill2,init_srange2,init_trange2];
lb = [realmin,realmin,realmin,realmin,realmin,realmin];
ub = [scovval(1),slag(end)*10000,tlag(end)*10000,...
      scovval(1),slag(end)*10000,tlag(end)*10000];
numparam = 6;

fmops = optimset('Algorithm','interior-point','DerivativeCheck','off',...
                 'Diagnostic','off','Display','off','FunValCheck','on',...
                 'MaxFunEvals',10000,'MaxIter',10000,...
                 'TolX',1e-16,'TolFun',1e-16);

lsfunc = @objfunsttwoexp;

[covparam,fval,exitflag,output] = ...
    fmincon(@(x) lsfunc(x,slag,scovval,sweight,tlag,tcovval,tweight),...
            initval,[],[],[],[],lb,ub,[],fmops);

% if exitflag < 1
%     disp('--- msg ---')
%     disp(output.message);
% end

numsample = size(slag,2) + size(tlag,2);
aicval = numsample * log(fval/numsample) + 2 * numparam;
