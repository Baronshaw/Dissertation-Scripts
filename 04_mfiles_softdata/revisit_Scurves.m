function [] = revisit_Scurves(soft,constant,gauss)
% created 4/15/2014 this function will explore S-curves for stations with
% poor X-valdiation statistics

if nargin < 1, soft = 1; end % soft data or not
if nargin < 2, constant = 0; end % constant offset or not
if nargin < 3, gauss = 1; end % gaussian soft data or not

if soft == 0, softstr = '_nosoft'; else softstr = '_soft'; end
if constant == 1, constr = '_constant'; else constr = '_long'; end
if gauss == 1, gaussstr = '_gauss'; else gaussstr = '_nongauss'; end

% load cross-valdiation results
load(sprintf('../05_mfiles_crossvalidation/Xval_10fold%s%s%s_results_test1.mat', ...
    softstr,constr,gaussstr));

% variables to work with
unisID = unique(ckall(:,1:2),'rows');
for i = 1:length(unisID)
    idx = ckall(:,1) == unisID(i,1) & ckall(:,2) == unisID(i,2);
    obs(i) = mean(zhall(idx));
end

% plot histogram of cross-valdiation results
figure; hold on;
hist(r2_sID,100);

% poor preformers
lb = prctile(r2_sID,1);
len = sum(r2_sID <= lb);
unisIDlb = unisID(r2_sID <= lb,:);

% obs v pred for the bad stations
for i = 1:len
    figure; hold on;
    idx = unisIDlb(i,1) == ckall(:,1) & unisIDlb(i,2) == ckall(:,2);
    plot(zkall(idx),zhall(idx),'b.');
    plot([0 max(zkall(idx))],[0 max(zkall(idx))],'k--'); % 1-1 line
    title(sprintf('mod v pred for %d values',sum(idx)));
    xlabel('pred');
    ylabel('obs');
end

% display the S curves near this station
toshowsID = unisIDlb(2,:); 
toshowidx = toshowsID(1) == ckall(:,1) & toshowsID(2) == ckall(:,2);
toshowckall = ckall(toshowidx,:);
dayz = datevec(ckall(toshowidx,3));
dayzdisp = dayz(:,1).*10000 + dayz(:,2).*100 + dayz(:,3);
deltaT = 365; n = 3; minpnts = 150; numbins = 10; negval = 0;

% three S curves (nsmax = 3) for the really bad values?

% S curves
for i = 1:length(dayz)
    
    % loading soft data: soft data NOT from test1
    load(sprintf('../matfiles/PM2p5_%d_%d_%d_%d_%d_neg%d.mat', ...
        dayzdisp(i),deltaT,n,minpnts,numbins,negval));

    figure; hold on;

    distz = sqrt( (CTMlocs(:,1)-toshowsID(1)).^2 + (CTMlocs(:,2)-toshowsID(2)).^2 );
    idxMod = distz == min(distz);
    plot(valMod(idxMod,:),valObs(idxMod,:),'b.');
    plot([get(gca,'XLim')],[get(gca,'YLim')],'k--');
    xlabel('Modeled Values of PM_{2.5} (\mug/m^3)');
    ylabel('Observed Values of PM_{2.5} (\mug/m^3)'); 
    title(sprintf('PM_{2.5} (\\mug/m^3) on %d, \\DeltaT=%d days, n=%d\nx=%0.0f km, y=%0.0f km',...        
        dayzdisp(i),deltaT,n,CTMlocs(idxMod,1)./1000,CTMlocs(idxMod,2)./1000));

    % add bins
    for j = 1:numbins+1
        plot([perctile_data(idxMod,j) perctile_data(idxMod,j)],get(gca,'YLim'),'r-');
    end

    % calculate means in each bin
    plot(mean_Mod(idxMod,:),mean_Obs(idxMod,:),'ro','MarkerFaceColor','r');

    % linear interpolation in each bin
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    interpx = linspace(xlim(1),xlim(2));
    interpy = interp1(mean_Mod(idxMod,:),mean_Obs(idxMod,:),interpx,'linear','extrap'); 
    plot(interpx,interpy,'r-');
    
    figure; hold on;
    idx = toshowsID(1) == ckall(:,1) & toshowsID(2) == ckall(:,2);
    plot(zkall(idx),zhall(idx),'b.');
    plot([0 max(zkall(idx))],[0 max(zkall(idx))],'k--'); % 1-1 line
    idx = toshowsID(1) == ckall(:,1) & toshowsID(2) == ckall(:,2) & datenum(dayz(i,:)) == ckall(:,3);
    plot(zkall(idx),zhall(idx),'ro');
    title(sprintf('mod v pred for %d values',sum(idx)));
    xlabel('pred');
    ylabel('obs');

end

end